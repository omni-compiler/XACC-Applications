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

!>  @brief Fortran kernel to update the external halo cells in a chunk.
!>  @author Wayne Gaudin
!>  @details Updates halo cells for the required fields at the required depth
!>  for any halo cells that lie on an external boundary. The location and type
!>  of data governs how this is carried out. External boundaries are always
!>  reflective.

MODULE update_halo_kernel_module

  USE data_module
  use xmp_defs_module

CONTAINS

  SUBROUTINE update_halo_kernel(x_min,x_max,y_min,y_max, &
    density0,                                            &
    energy0,                                             &
    pressure,                                            &
    viscosity,                                           &
    soundspeed,                                          &
    density1,                                            &
    energy1,                                             &
    xvel0,                                               &
    yvel0,                                               &
    xvel1,                                               &
    yvel1,                                               &
    vol_flux_x,                                          &
    vol_flux_y,                                          &
    mass_flux_x,                                         &
    mass_flux_y,                                         &
    fields,                                              &
    depth                                                )
    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density0,energy0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure,viscosity,soundspeed
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1,energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x,mass_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y,mass_flux_y
    INTEGER :: fields(:),depth

!$xmp nodes p(*,*)
!$xmp template t(x_min-2:x_max+3,y_min-2:y_max+3)
!$xmp distribute t(gblock(blocksize_x), gblock(blocksize_y)) onto p
!$xmp align (i,j) with t(i,j) :: density0, density1, energy0, energy1, pressure, viscosity, soundspeed,&
!$xmp                            xvel0, xvel1, yvel0, yvel1, vol_flux_x, mass_flux_x, vol_flux_y, mass_flux_y
!$xmp shadow (2:2, 2:2) :: density0, density1, energy0, energy1, pressure, viscosity, soundspeed
!$xmp shadow (2:3, 2:3) :: xvel0, xvel1, yvel0, yvel1
!$xmp shadow (2:3, 2:2) :: vol_flux_x, mass_flux_x
!$xmp shadow (2:2, 2:3) :: vol_flux_y, mass_flux_y
!$xmp save_desc :: p,t,density0,density1,energy0,energy1,pressure,viscosity,soundspeed,&
!$xmp              xvel0,xvel1,yvel0,yvel1,vol_flux_x,mass_flux_x,vol_flux_y,mass_flux_y

    ! These need to be kept consistent with the data module to avoid use statement
    !INTEGER,      PARAMETER :: CHUNK_LEFT   =1    &
    !                          ,CHUNK_RIGHT  =2    &
    !                          ,CHUNK_BOTTOM =3    &
    !                          ,CHUNK_TOP    =4    &
    !                          ,EXTERNAL_FACE=-1

    !INTEGER,      PARAMETER :: FIELD_DENSITY0   = 1         &
    !                          ,FIELD_DENSITY1   = 2         &
    !                          ,FIELD_ENERGY0    = 3         &
    !                          ,FIELD_ENERGY1    = 4         &
    !                          ,FIELD_PRESSURE   = 5         &
    !                          ,FIELD_VISCOSITY  = 6         &
    !                          ,FIELD_SOUNDSPEED = 7         &
    !                          ,FIELD_XVEL0      = 8         &
    !                          ,FIELD_XVEL1      = 9         &
    !                          ,FIELD_YVEL0      =10         &
    !                          ,FIELD_YVEL1      =11         &
    !                          ,FIELD_VOL_FLUX_X =12         &
    !                          ,FIELD_VOL_FLUX_Y =13         &
    !                          ,FIELD_MASS_FLUX_X=14         &
    !                          ,FIELD_MASS_FLUX_Y=15

    ! Update values in external halo cells based on depth and fields requested
    ! Even though half of these loops look the wrong way around, it should be noted
    ! that depth is either 1 or 2 so that it is more efficient to always thread
    ! loop along the mesh edge.
    IF(fields(FIELD_DENSITY0).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        density0,                                    &
        depth                                        )
    ENDIF

    IF(fields(FIELD_DENSITY1).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        density1,                                    &
        depth                                        )
    ENDIF

    IF(fields(FIELD_ENERGY0).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        energy0,                                     &
        depth                                        )
    ENDIF

    IF(fields(FIELD_ENERGY1).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        energy1,                                     &
        depth                                        )
    ENDIF

    IF(fields(FIELD_PRESSURE).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        pressure,                                    &
        depth                                        )
    ENDIF

    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        viscosity,                                   &
        depth                                        )
    ENDIF

    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
      CALL update_halo_cell(x_min,x_max,y_min,y_max, &
        soundspeed,                                  &
        depth                                        )
    ENDIF

    IF(fields(FIELD_XVEL0).EQ.1) THEN
      CALL update_halo_vertex(x_min,x_max,y_min,y_max, &
        xvel0,                                         &
        -1,1,                                          &
        depth                                          )
    ENDIF

    IF(fields(FIELD_XVEL1).EQ.1) THEN
      CALL update_halo_vertex(x_min,x_max,y_min,y_max, &
        xvel1,                                         &
        -1,1,                                          &
        depth                                          )
    ENDIF

    IF(fields(FIELD_YVEL0).EQ.1) THEN
      CALL update_halo_vertex(x_min,x_max,y_min,y_max, &
        yvel0,                                         &
        1,-1,                                          &
        depth                                          )
    ENDIF

    IF(fields(FIELD_YVEL1).EQ.1) THEN
      CALL update_halo_vertex(x_min,x_max,y_min,y_max, &
        yvel1,                                         &
        1,-1,                                          &
        depth                                          )
    ENDIF

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
      CALL update_halo_x_face(x_min,x_max,y_min,y_max, &
        vol_flux_x,                                    &
        depth                                          )
    ENDIF

    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
      CALL update_halo_x_face(x_min,x_max,y_min,y_max, &
        mass_flux_x,                                   &
        depth                                          )
    ENDIF

    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
      CALL update_halo_y_face(x_min,x_max,y_min,y_max, &
        vol_flux_y,                                    &
        depth                                          )
    ENDIF

    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
      CALL update_halo_y_face(x_min,x_max,y_min,y_max, &
        mass_flux_y,                                   &
        depth                                          )
    ENDIF

  CONTAINS

    SUBROUTINE update_halo_cell(x_min,x_max,y_min,y_max, &
      cell,                                              &
      depth                                              )
      IMPLICIT NONE

      INTEGER :: x_min,x_max,y_min,y_max
      REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: cell
      INTEGER :: depth
      INTEGER :: j,k
!$xmp align (i,j) with t(i,j) :: cell
!$xmp shadow (2:2, 2:2) :: cell
!$xmp save_desc :: cell

!$ACC DATA PRESENT(cell)
!$xmp loop (j,k) on t(j,k) expand(depth, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+depth
!$ACC LOOP INDEPENDENT
        DO k=1,depth
          cell(j,1-k)=cell(j,0+k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(depth, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+depth
!$ACC LOOP INDEPENDENT
        DO k=y_max-depth+1,y_max
          cell(j,2*y_max+1-k)=cell(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+depth
!$ACC LOOP INDEPENDENT
        DO j=1,depth
          cell(1-j,k)=cell(0+j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+depth
!$ACC LOOP INDEPENDENT
        DO j=x_max-depth+1,x_max
          cell(2*x_max+1-j,k)=cell(j,k)
        ENDDO
      ENDDO
!$ACC END DATA

    END SUBROUTINE update_halo_cell

    SUBROUTINE update_halo_vertex(x_min,x_max,y_min,y_max, &
      vertex,                                              &
      x_inv,y_inv,                                         &
      depth                                                )
      IMPLICIT NONE

      INTEGER :: x_min,x_max,y_min,y_max
      INTEGER :: x_inv,y_inv
      REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: vertex
      INTEGER :: depth
      INTEGER :: j,k
!$xmp align (i,j) with t(i,j) :: vertex
!$xmp shadow (2:3, 2:3) :: vertex
!$xmp save_desc :: vertex

!$ACC DATA PRESENT(vertex)
!$xmp loop (j,k) on t(j,k) expand(depth:depth+1, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+1+depth
!$ACC LOOP INDEPENDENT
        DO k=1,depth
          vertex(j,1-k)=y_inv*vertex(j,1+k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(depth:depth+1, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+1+depth
!$ACC LOOP INDEPENDENT
        DO k=y_max-depth+1,y_max
          vertex(j,2*(y_max+1)-k)=y_inv*vertex(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth+1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+1+depth
!$ACC LOOP INDEPENDENT
        DO j=1,depth
          vertex(1-j,k)=x_inv*vertex(1+j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth+1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+1+depth
        !$ACC LOOP INDEPENDENT
        DO j=x_max-depth+1,x_max
          vertex(2*(x_max+1)-j,k)=x_inv*vertex(j,k)
        ENDDO
      ENDDO
!$ACC END DATA

    END SUBROUTINE update_halo_vertex


    SUBROUTINE update_halo_x_face(x_min,x_max,y_min,y_max, &
      x_face,                                              &
      depth                                                )
      IMPLICIT NONE

      INTEGER :: x_min,x_max,y_min,y_max
      REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: x_face
      INTEGER :: depth
      INTEGER :: j,k
!$xmp align (i,j) with t(i,j) :: x_face
!$xmp shadow (2:3, 2:2) :: x_face
!$xmp save_desc :: x_face

!$ACC DATA PRESENT(x_face)
!$xmp loop (j,k) on t(j,k) expand(depth:depth+1, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+1+depth
!$ACC LOOP INDEPENDENT
        DO k=1,depth
          x_face(j,1-k)=x_face(j,1+k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(depth:depth+1, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+1+depth
!$ACC LOOP INDEPENDENT
        DO k=y_max-depth,y_max-1
          x_face(j,2*(y_max+0)-k)=x_face(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+depth
!$ACC LOOP INDEPENDENT
        DO j=1,depth
          x_face(1-j,k)=-x_face(1+j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+depth
!$ACC LOOP INDEPENDENT
        DO j=x_max-depth+1, x_max
          x_face(2*(x_max+1)-j,k)=-x_face(j,k)
        ENDDO
      ENDDO
!$ACC END DATA

    END SUBROUTINE update_halo_x_face


    SUBROUTINE update_halo_y_face(x_min,x_max,y_min,y_max, &
      y_face,                                              &
      depth                                                )
      IMPLICIT NONE

      INTEGER :: x_min,x_max,y_min,y_max
      REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: y_face
      INTEGER :: depth
      INTEGER :: j,k
!$xmp align (i,j) with t(i,j) :: y_face
!$xmp shadow (2:2, 2:3) :: y_face
!$xmp save_desc :: y_face

!$ACC DATA PRESENT(y_face)
!$xmp loop (j,k) on t(j,k) expand(depth:depth, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+depth
!$ACC LOOP INDEPENDENT
        DO k=1,depth
          y_face(j,1-k)=-y_face(j,1+k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(depth:depth, 0)
!$ACC KERNELS LOOP INDEPENDENT
      DO j=x_min-depth,x_max+depth
!$ACC LOOP INDEPENDENT
        DO k=y_max-depth+1,y_max
          y_face(j,2*(y_max+1)-k)=-y_face(j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth+1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+1+depth
!$ACC LOOP INDEPENDENT
        DO j=1,depth
          y_face(1-j,k)=y_face(1+j,k)
        ENDDO
      ENDDO

!$xmp loop (j,k) on t(j,k) expand(0, depth:depth+1)
!$ACC KERNELS LOOP INDEPENDENT
      DO k=y_min-depth,y_max+1+depth
!$ACC LOOP INDEPENDENT
        DO j=x_max-depth, x_max-1
          y_face(2*(x_max+0)-j,k)=y_face(j,k)
        ENDDO
      ENDDO
!$ACC END DATA

    END SUBROUTINE update_halo_y_face

  END SUBROUTINE update_halo_kernel

END  MODULE update_halo_kernel_module
