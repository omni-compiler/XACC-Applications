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

!>  @brief Fortran cell advection kernel.
!>  @author Wayne Gaudin
!>  @details Performs a second order advective remap using van-Leer limiting
!>  with directional splitting.

MODULE advec_cell_kernel_module

CONTAINS

  SUBROUTINE advec_cell_kernel(x_min,       &
                               x_max,       &
                               y_min,       &
                               y_max,       &
                               dir,         &
                               sweep_number,&
                               vertexdx,    &
                               vertexdy,    &
                               volume,      &
                               density1,    &
                               energy1,     &
                               mass_flux_x, &
                               vol_flux_x,  &
                               mass_flux_y, &
                               vol_flux_y,  &
                               pre_vol,     &
                               post_vol,    &
                               pre_mass,    &
                               post_mass,   &
                               advec_vol,   &
                               post_ener,   &
                               ener_flux    )

    IMPLICIT NONE
    use xmp_defs_module

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: sweep_number,dir
    INTEGER :: g_xdir=1,g_ydir=2

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: advec_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_ener
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

    REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy

    INTEGER :: j,k,upwind,donor,downwind,dif

    REAL(KIND=8) :: wind,sigma,sigmat,sigmav,sigmam,sigma3,sigma4
    REAL(KIND=8) :: diffuw,diffdw,limiter
    REAL(KIND=8) :: one_by_six=1.0_8/6.0_8
    REAL(KIND=8) :: pre_mass_s,post_mass_s,post_ener_s,advec_vol_s

!$xmp nodes p(*,*)
!$xmp template t(x_min-2:x_max+3,y_min-2:y_max+3)
!$xmp distribute t(gblock(blocksize_x), gblock(blocksize_y)) onto p
!$xmp align (i,j) with t(i,j) :: volume, density1, energy1, vol_flux_x, mass_flux_x, vol_flux_y, mass_flux_y, &
!$xmp                            pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux
!$xmp align (i) with t(i,*) :: vertexdx
!$xmp align (j) with t(*,j) :: vertexdy
!$xmp shadow (2:2, 2:2) :: volume, density1, energy1
!$xmp shadow (2:3, 2:2) :: vol_flux_x, mass_flux_x
!$xmp shadow (2:2, 2:3) :: vol_flux_y, mass_flux_y
!$xmp shadow (2:3, 2:3) :: pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux
!$xmp shadow (2:3) :: vertexdx, vertexdy
!$xmp save_desc :: p, t, volume, density1, energy1, vol_flux_x, vol_flux_y, mass_flux_x, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux
    INTEGER :: l_x_max, l_y_max
    l_x_max = g_l_max_x
    l_y_max = g_l_max_y
!$ACC DATA &
!$ACC PCOPY(density1,energy1) &
!$ACC PCOPY(vol_flux_x,vol_flux_y,volume,mass_flux_x,mass_flux_y,vertexdx,vertexdy) &
!$ACC PCOPY(pre_vol,post_vol,ener_flux)

    IF(dir.EQ.g_xdir) THEN
      IF(sweep_number.EQ.1)THEN
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min-2,x_max+2
            pre_vol(j,k)=volume(j,k)+(vol_flux_x(j+1,k  )-vol_flux_x(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
            post_vol(j,k)=pre_vol(j,k)-(vol_flux_x(j+1,k  )-vol_flux_x(j,k))
          ENDDO
        ENDDO
      ELSE
!$xmp loop (j,k) on t(j,k) expand(2,2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min-2,x_max+2
            pre_vol(j,k)=volume(j,k)+vol_flux_x(j+1,k)-vol_flux_x(j,k)
            post_vol(j,k)=volume(j,k)
          ENDDO
        ENDDO
      ENDIF

!$xmp loop (j,k) on t(j,k) expand(0:2, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max
!$ACC LOOP INDEPENDENT PRIVATE(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam, diffuw,diffdw,limiter,wind)
        DO j=x_min,x_max+2
          IF(vol_flux_x(j,k).GT.0.0)THEN
            upwind   =j-2
            donor    =j-1
            downwind =j
            dif      =donor
          ELSE
            upwind   =MIN(j+1,l_x_max+2)
            donor    =j
            downwind =j-1
            dif      =upwind
          ENDIF

          sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)
          sigma3=(1.0_8+sigmat)*(vertexdx(j)/vertexdx(dif))
          sigma4=2.0_8-sigmat

          sigma=sigmat
          sigmav=sigmat

          diffuw=density1(donor,k)-density1(upwind,k)
          diffdw=density1(downwind,k)-density1(donor,k)
          wind=1.0_8
          IF(diffdw.LE.0.0) wind=-1.0_8
          IF(diffuw*diffdw.GT.0.0)THEN
            limiter=(1.0_8-sigmav)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
          ELSE
            limiter=0.0
          ENDIF
          mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter)

          sigmam=ABS(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k))
          diffuw=energy1(donor,k)-energy1(upwind,k)
          diffdw=energy1(downwind,k)-energy1(donor,k)
          wind=1.0_8
          IF(diffdw.LE.0.0) wind=-1.0_8
          IF(diffuw*diffdw.GT.0.0)THEN
            limiter=(1.0_8-sigmam)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
          ELSE
            limiter=0.0
          ENDIF

          ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max
!$ACC LOOP INDEPENDENT PRIVATE(pre_mass_s,post_mass_s,post_ener_s,advec_vol_s)
        DO j=x_min,x_max
          pre_mass_s=density1(j,k)*pre_vol(j,k)
          post_mass_s=pre_mass_s+mass_flux_x(j,k)-mass_flux_x(j+1,k)
          post_ener_s=(energy1(j,k)*pre_mass_s+ener_flux(j,k)-ener_flux(j+1,k))/post_mass_s
          advec_vol_s=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k)
          density1(j,k)=post_mass_s/advec_vol_s
          energy1(j,k)=post_ener_s
        ENDDO
      ENDDO
    ELSEIF(dir.EQ.g_ydir) THEN
      IF(sweep_number.EQ.1)THEN
!$xmp loop (j,k) on t(j,k) expand(2, 2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min-2,x_max+2
            pre_vol(j,k)=volume(j,k)+(vol_flux_y(j  ,k+1)-vol_flux_y(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k))
            post_vol(j,k)=pre_vol(j,k)-(vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
          ENDDO
        ENDDO
      ELSE
!$xmp loop (j,k) on t(j,k) expand(2, 2)
!$ACC KERNELS LOOP INDEPENDENT
        DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
          DO j=x_min-2,x_max+2
            pre_vol(j,k)=volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
            post_vol(j,k)=volume(j,k)
          ENDDO
        ENDDO
      ENDIF

!$xmp loop (j,k) on t(j,k) expand(0, 0:2)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max+2
!$ACC LOOP INDEPENDENT PRIVATE(upwind,donor,downwind,dif,sigmat,sigma3,sigma4,sigmav,sigma,sigmam,diffuw,diffdw,limiter,wind)
        DO j=x_min,x_max
          IF(vol_flux_y(j,k).GT.0.0)THEN
            upwind   =k-2
            donor    =k-1
            downwind =k
            dif      =donor
          ELSE
            upwind   =MIN(k+1,l_y_max+2)
            donor    =k
            downwind =k-1
            dif      =upwind
          ENDIF

          sigmat=ABS(vol_flux_y(j,k))/pre_vol(j,donor)
          sigma3=(1.0_8+sigmat)*(vertexdy(k)/vertexdy(dif))
          sigma4=2.0_8-sigmat

          sigma=sigmat
          sigmav=sigmat

          diffuw=density1(j,donor)-density1(j,upwind)
          diffdw=density1(j,downwind)-density1(j,donor)
          wind=1.0_8
          IF(diffdw.LE.0.0) wind=-1.0_8
          IF(diffuw*diffdw.GT.0.0)THEN
            limiter=(1.0_8-sigmav)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
          ELSE
            limiter=0.0
          ENDIF
          mass_flux_y(j,k)=vol_flux_y(j,k)*(density1(j,donor)+limiter)

          sigmam=ABS(mass_flux_y(j,k))/(density1(j,donor)*pre_vol(j,donor))
          diffuw=energy1(j,donor)-energy1(j,upwind)
          diffdw=energy1(j,downwind)-energy1(j,donor)
          wind=1.0_8
          IF(diffdw.LE.0.0) wind=-1.0_8
          IF(diffuw*diffdw.GT.0.0)THEN
            limiter=(1.0_8-sigmam)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
          ELSE
            limiter=0.0
          ENDIF

          ener_flux(j,k)=mass_flux_y(j,k)*(energy1(j,donor)+limiter)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min,y_max
!$ACC LOOP INDEPENDENT PRIVATE(pre_mass_s,post_mass_s,post_ener_s,advec_vol_s)
        DO j=x_min,x_max
          pre_mass_s=density1(j,k)*pre_vol(j,k)
          post_mass_s=pre_mass_s+mass_flux_y(j,k)-mass_flux_y(j,k+1)
          post_ener_s=(energy1(j,k)*pre_mass_s+ener_flux(j,k)-ener_flux(j,k+1))/post_mass_s
          advec_vol_s=pre_vol(j,k)+vol_flux_y(j,k)-vol_flux_y(j,k+1)
          density1(j,k)=post_mass_s/advec_vol_s
          energy1(j,k)=post_ener_s
        ENDDO
      ENDDO
    ENDIF

!$ACC END DATA
  END SUBROUTINE advec_cell_kernel

END MODULE advec_cell_kernel_module

