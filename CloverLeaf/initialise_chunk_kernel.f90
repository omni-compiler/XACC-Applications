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

!>  @brief Fortran chunk initialisation kernel.
!>  @author Wayne Gaudin
!>  @details Calculates mesh geometry for the mesh chunk based on the mesh size.

MODULE initialise_chunk_kernel_module

CONTAINS

  SUBROUTINE initialise_chunk_kernel(x_min,x_max,y_min,y_max, &
    xmin,ymin,dx,dy,                                          &
    vertexx,                                                  &
    vertexdx,                                                 &
    vertexy,                                                  &
    vertexdy,                                                 &
    cellx,                                                    &
    celldx,                                                   &
    celly,                                                    &
    celldy,                                                   &
    volume,                                                   &
    xarea,                                                    &
    yarea )

    IMPLICIT NONE
    use xmp_defs_module

    INTEGER      :: x_min,x_max,y_min,y_max
    REAL(KIND=8) :: xmin,ymin,dx,dy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
    REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
    REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: xarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: yarea
!$xmp nodes p(*,*)
!$xmp template t(x_min-2:x_max+3,y_min-2:y_max+3)
!$xmp distribute t(gblock(blocksize_x), gblock(blocksize_y)) onto p
!$xmp align (i) with t(i,*) :: vertexx, vertexdx, cellx, celldx
!$xmp align (j) with t(*,j) :: vertexy, vertexdy, celly, celldy
!$xmp align (i,j) with t(i,j) :: volume, xarea, yarea
!$xmp shadow (2:3) :: vertexx, vertexdx, vertexy, vertexdy
!$xmp shadow (2:2) :: cellx, celldx, celly, celldy
!$xmp shadow (2:2, 2:2) :: volume
!$xmp shadow (2:3, 2:2) :: xarea
!$xmp shadow (2:2, 2:3) :: yarea

    INTEGER      :: j,k
!$ACC DATA

!$xmp loop (j) on t(j,*) expand(2:3,0)
!$ACC KERNELS LOOP INDEPENDENT
    DO j=x_min-2,x_max+3
      vertexx(j)=xmin+dx*float(j-x_min)
    ENDDO

!$xmp loop (j) on t(j,*) expand(2:3,0)
!$ACC KERNELS LOOP INDEPENDENT
    DO j=x_min-2,x_max+3
      vertexdx(j)=dx
    ENDDO

!$xmp loop (k) on t(*,k) expand(0,2:3)
!$ACC KERNELS LOOP INDEPENDENT
    DO k=y_min-2,y_max+3
      vertexy(k)=ymin+dy*float(k-y_min)
    ENDDO

!$xmp loop (k) on t(*,k) expand(0,2:3)
!$ACC KERNELS LOOP INDEPENDENT
    DO k=y_min-2,y_max+3
      vertexdy(k)=dy
    ENDDO

!$xmp loop (j) on t(j,*) expand(2:2,0)
!$ACC KERNELS LOOP INDEPENDENT
    DO j=x_min-2,x_max+2
      cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
    ENDDO

!$xmp loop (j) on t(j,*) expand(2:2,0)
!$ACC KERNELS LOOP INDEPENDENT
    DO j=x_min-2,x_max+2
      celldx(j)=dx
    ENDDO

!$xmp loop (k) on t(*,k) expand(0,2:2)
!$ACC KERNELS LOOP INDEPENDENT
    DO k=y_min-2,y_max+2
      celly(k)=0.5*(vertexy(k)+vertexy(k+1))
    ENDDO

!$xmp loop (k) on t(*,k) expand(0,2:2)
!$ACC KERNELS LOOP INDEPENDENT
    DO k=y_min-2,y_max+2
      celldy(k)=dy
    ENDDO

!$xmp loop (j,k) on t(j,k) expand(2:2,2:2)
!$ACC KERNELS LOOP INDEPENDENT PRIVATE(j)
    DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
      DO j=x_min-2,x_max+2
        volume(j,k)=dx*dy
      ENDDO
    ENDDO

!$xmp loop (j,k) on t(j,k) expand(2:2,2:2)
!$ACC KERNELS LOOP INDEPENDENT PRIVATE(j)
    DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
      DO j=x_min-2,x_max+2
        xarea(j,k)=celldy(k)
      ENDDO
    ENDDO

!$xmp loop (j,k) on t(j,k) expand(2:2,2:2)
!$ACC KERNELS LOOP INDEPENDENT PRIVATE(j)
    DO k=y_min-2,y_max+2
!$ACC LOOP INDEPENDENT
      DO j=x_min-2,x_max+2
        yarea(j,k)=celldx(j)
      ENDDO
    ENDDO
!$ACC END DATA

  END SUBROUTINE initialise_chunk_kernel

END MODULE initialise_chunk_kernel_module
