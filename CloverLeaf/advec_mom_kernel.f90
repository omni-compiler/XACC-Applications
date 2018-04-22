!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran momentum advection kernel
!>  @author Wayne Gaudin
!>  @details Performs a second order advective remap on the vertex momentum
!>  using van-Leer limiting and directional splitting.
!>  Note that although pre_vol is only set and not used in the update, please
!>  leave it in the method.

MODULE advec_mom_kernel_mod

CONTAINS

  SUBROUTINE advec_mom_kernel(x_min,x_max,y_min,y_max,   &
                              vel1,                      &
                              mass_flux_x,               &
                              vol_flux_x,                &
                              mass_flux_y,               &
                              vol_flux_y,                &
                              volume,                    &
                              density1,                  &
                              node_flux,                 &
                              node_mass_post,            &
                              node_mass_pre,             &
                              mom_flux,                  &
                              pre_vol,                   &
                              post_vol,                  &
                              celldx,                    &
                              celldy,                    &
                              which_vel,                 &
                              sweep_number,              &
                              direction                  )

    IMPLICIT NONE
    use xmp_defs_module

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: which_vel,sweep_number,direction

    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: vel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: node_flux
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: node_mass_post
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: node_mass_pre
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: mom_flux
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy

!$xmp nodes p(*,*)
!$xmp template t(x_min-2:x_max+3,y_min-2:y_max+3)
!$xmp distribute t(gblock(blocksize_x), gblock(blocksize_y)) onto p
!$xmp align (i,j) with t(i,j) :: vol_flux_x, mass_flux_x, vol_flux_y, mass_flux_y, volume, density1, &
!$xmp                            vel1, node_flux, node_mass_post, node_mass_pre, mom_flux, pre_vol, post_vol
!$xmp align (i) with t(i,*) :: celldx
!$xmp align (j) with t(*,j) :: celldy
!$xmp shadow (2:2, 2:2) :: volume, density1
!$xmp shadow (2:3, 2:2) :: vol_flux_x, mass_flux_x
!$xmp shadow (2:2, 2:3) :: vol_flux_y, mass_flux_y
!$xmp shadow (2:3, 2:3) :: vel1, node_flux, node_mass_post, node_mass_pre, mom_flux, pre_vol, post_vol
!$xmp shadow (2:2) :: celldx, celldy
!$xmp save_desc :: p, t, vol_flux_x, mass_flux_x, vol_flux_y, mass_flux_y, volume, density1, celldx, celldy, vel1, node_flux, node_mass_post, node_mass_pre, mom_flux, pre_vol, post_vol
    INTEGER :: j,k,mom_sweep
    INTEGER :: upwind,donor,downwind,dif
    REAL(KIND=8) :: sigma,wind,width
    REAL(KIND=8) :: vdiffuw,vdiffdw,auw,adw,limiter
    REAL(KIND=8) :: advec_vel_s

    mom_sweep=direction+2*(sweep_number-1)

!$ACC DATA          &
!$ACC PCOPY(vel1)    &
!$ACC PCOPY(volume,mass_flux_x,mass_flux_y,vol_flux_x,vol_flux_y,density1,celldx,celldy) &
!$ACC PCOPY(mom_flux,node_flux,node_mass_post,node_mass_pre,post_vol,pre_vol)


    IF(mom_sweep.EQ.1)THEN ! x 1
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
        DO j=x_min-2,x_max+2
          post_vol(j,k)= volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
          pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k)
        ENDDO
      ENDDO
    ELSEIF(mom_sweep.EQ.2)THEN ! y 1
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
        DO j=x_min-2,x_max+2
          post_vol(j,k)= volume(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k)
          pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
        ENDDO
      ENDDO
    ELSEIF(mom_sweep.EQ.3)THEN ! x 2
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
        DO j=x_min-2,x_max+2
          post_vol(j,k)=volume(j,k)
          pre_vol(j,k)=post_vol(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
        ENDDO
      ENDDO
    ELSEIF(mom_sweep.EQ.4)THEN ! y 2
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
        DO j=x_min-2,x_max+2
          post_vol(j,k)=volume(j,k)
          pre_vol(j,k)=post_vol(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k)
        ENDDO
      ENDDO
    ENDIF

    IF(direction.EQ.1)THEN
      IF(which_vel.EQ.1)THEN
!$xmp loop (j,k) on t(j,k) expand(2,0:1)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min,y_max+1
!$ACC LOOP INDEPENDENT
          DO j=x_min-2,x_max+2
            ! Find staggered mesh mass fluxes, nodal masses and volumes.
            node_flux(j,k)=0.25_8*(mass_flux_x(j,k-1  )+mass_flux_x(j  ,k)  &
              +mass_flux_x(j+1,k-1)+mass_flux_x(j+1,k))
          ENDDO
        ENDDO

!$xmp loop (j,k) on t(j,k) expand(1:2,0:1)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min,y_max+1
!$ACC LOOP INDEPENDENT
          DO j=x_min-1,x_max+2
            ! Staggered cell mass post advection
            node_mass_post(j,k)=0.25_8*(density1(j  ,k-1)*post_vol(j  ,k-1)                   &
              +density1(j  ,k  )*post_vol(j  ,k  )                   &
              +density1(j-1,k-1)*post_vol(j-1,k-1)                   &
              +density1(j-1,k  )*post_vol(j-1,k  ))
            node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j-1,k)+node_flux(j,k)
          ENDDO
        ENDDO
      ENDIF

!$xmp loop (j,k) on t(j,k) expand(1,0:1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max+1
!$ACC LOOP INDEPENDENT PRIVATE(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,advec_vel_s)
        DO j=x_min-1,x_max+1
          IF(node_flux(j,k).LT.0.0)THEN
            upwind=j+2
            donor=j+1
            downwind=j
            dif=donor
          ELSE
            upwind=j-1
            donor=j
            downwind=j+1
            dif=upwind
          ENDIF
          sigma=ABS(node_flux(j,k))/(node_mass_pre(donor,k))
          width=celldx(j)
          vdiffuw=vel1(donor,k)-vel1(upwind,k)
          vdiffdw=vel1(downwind,k)-vel1(donor,k)
          limiter=0.0
          IF(vdiffuw*vdiffdw.GT.0.0)THEN
            auw=ABS(vdiffuw)
            adw=ABS(vdiffdw)
            wind=1.0_8
            IF(vdiffdw.LE.0.0) wind=-1.0_8
            limiter=wind*MIN(width*((2.0_8-sigma)*adw/width+(1.0_8+sigma)*auw/celldx(dif))/6.0_8,auw,adw)
          ENDIF
          advec_vel_s=vel1(donor,k)+(1.0-sigma)*limiter
          mom_flux(j,k)=advec_vel_s*node_flux(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0:1,0:1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max+1
!$ACC LOOP INDEPENDENT
        DO j=x_min,x_max+1
          vel1 (j,k)=(vel1 (j,k)*node_mass_pre(j,k)+mom_flux(j-1,k)-mom_flux(j,k))/node_mass_post(j,k)
        ENDDO
      ENDDO

    ELSEIF(direction.EQ.2)THEN
      IF(which_vel.EQ.1)THEN
!$xmp loop (j,k) on t(j,k) expand(0:1,2:2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min,x_max+1
            ! Find staggered mesh mass fluxes and nodal masses and volumes.
            node_flux(j,k)=0.25_8*(mass_flux_y(j-1,k  )+mass_flux_y(j  ,k  ) &
              +mass_flux_y(j-1,k+1)+mass_flux_y(j  ,k+1))
          ENDDO
        ENDDO

!$xmp loop (j,k) on t(j,k) expand(0:1,1:2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-1,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min,x_max+1
            node_mass_post(j,k)=0.25_8*(density1(j  ,k-1)*post_vol(j  ,k-1)                     &
              +density1(j  ,k  )*post_vol(j  ,k  )                     &
              +density1(j-1,k-1)*post_vol(j-1,k-1)                     &
              +density1(j-1,k  )*post_vol(j-1,k  ))
            node_mass_pre(j,k)=node_mass_post(j,k)-node_flux(j,k-1)+node_flux(j,k)
          ENDDO
        ENDDO
      ENDIF

!$xmp loop (j,k) on t(j,k) expand(0:1,1:1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-1,y_max+1
!$ACC LOOP INDEPENDENT PRIVATE(upwind,donor,downwind,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,advec_vel_s)
        DO j=x_min,x_max+1
          IF(node_flux(j,k).LT.0.0)THEN
            upwind=k+2
            donor=k+1
            downwind=k
            dif=donor
          ELSE
            upwind=k-1
            donor=k
            downwind=k+1
            dif=upwind
          ENDIF

          sigma=ABS(node_flux(j,k))/(node_mass_pre(j,donor))
          width=celldy(k)
          vdiffuw=vel1(j,donor)-vel1(j,upwind)
          vdiffdw=vel1(j,downwind)-vel1(j,donor)
          limiter=0.0
          IF(vdiffuw*vdiffdw.GT.0.0)THEN
            auw=ABS(vdiffuw)
            adw=ABS(vdiffdw)
            wind=1.0_8
            IF(vdiffdw.LE.0.0) wind=-1.0_8
            limiter=wind*MIN(width*((2.0_8-sigma)*adw/width+(1.0_8+sigma)*auw/celldy(dif))/6.0_8,auw,adw)
          ENDIF
          advec_vel_s=vel1(j,donor)+(1.0_8-sigma)*limiter
          mom_flux(j,k)=advec_vel_s*node_flux(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0:1,0:1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max+1
!$ACC LOOP INDEPENDENT
        DO j=x_min,x_max+1
          vel1 (j,k)=(vel1(j,k)*node_mass_pre(j,k)+mom_flux(j,k-1)-mom_flux(j,k))/node_mass_post(j,k)
        ENDDO
      ENDDO
    ENDIF



!$ACC END DATA

  END SUBROUTINE advec_mom_kernel

END MODULE advec_mom_kernel_mod
