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

!>  @brief  Allocates the data for each mesh chunk
!>  @author Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

SUBROUTINE build_field()

  USE clover_module

  IMPLICIT NONE

  INTEGER :: j,k

  call init_xmp_defs(1, grid%x_cells, 1, grid%y_cells)
  !$xmp template_fix (gblock(blocksize_x), gblock(blocksize_y)) tmplt_dyn(1-2:grid%x_cells+3, 1-2:grid%y_cells+3)

    ALLOCATE(density0  (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(density1  (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(energy0   (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(energy1   (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(pressure  (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(viscosity0(chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(soundspeed(chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))

    ALLOCATE(xvel0(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(xvel1(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(yvel0(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(yvel1(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))

    ALLOCATE(vol_flux_x (chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(mass_flux_x(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(vol_flux_y (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(mass_flux_y(chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+3))

    ALLOCATE(work_array1(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array2(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array3(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array4(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array5(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array6(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(work_array7(chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+3))

    ALLOCATE(cellx   (chunk%x_min-2:chunk%x_max+2))
    ALLOCATE(celly   (chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(vertexx (chunk%x_min-2:chunk%x_max+3))
    ALLOCATE(vertexy (chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(celldx  (chunk%x_min-2:chunk%x_max+2))
    ALLOCATE(celldy  (chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(vertexdx(chunk%x_min-2:chunk%x_max+3))
    ALLOCATE(vertexdy(chunk%y_min-2:chunk%y_max+3))
    ALLOCATE(volume  (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(xarea   (chunk%x_min-2:chunk%x_max+3, &
      chunk%y_min-2:chunk%y_max+2))
    ALLOCATE(yarea   (chunk%x_min-2:chunk%x_max+2, &
      chunk%y_min-2:chunk%y_max+3))

    ! Zeroing isn't strictly neccessary but it ensures physical pages
    ! are allocated. This prevents first touch overheads in the main code
    ! cycle which can skew timings in the first step

    !$xmp loop (j,k) on tmplt_dyn(j,k) expand(2:3,2:3)
    DO k=chunk%y_min-2,chunk%y_max+3
      DO j=chunk%x_min-2,chunk%x_max+3
        work_array1(j,k)=0.0
        work_array2(j,k)=0.0
        work_array3(j,k)=0.0
        work_array4(j,k)=0.0
        work_array5(j,k)=0.0
        work_array6(j,k)=0.0
        work_array7(j,k)=0.0

        xvel0(j,k)=0.0
        xvel1(j,k)=0.0
        yvel0(j,k)=0.0
        yvel1(j,k)=0.0
      ENDDO
    ENDDO

    !$xmp loop (j,k) on tmplt_dyn(j,k) expand(2:2,2:2)
    DO k=chunk%y_min-2,chunk%y_max+2
      DO j=chunk%x_min-2,chunk%x_max+2
        density0(j,k)=0.0
        density1(j,k)=0.0
        energy0(j,k)=0.0
        energy1(j,k)=0.0
        pressure(j,k)=0.0
        viscosity0(j,k)=0.0
        soundspeed(j,k)=0.0
        volume(j,k)=0.0
      ENDDO
    ENDDO

    !$xmp loop (j,k) on tmplt_dyn(j,k) expand(2:3,2:2)
    DO k=chunk%y_min-2,chunk%y_max+2
      DO j=chunk%x_min-2,chunk%x_max+3
        vol_flux_x(j,k)=0.0
        mass_flux_x(j,k)=0.0
        xarea(j,k)=0.0
      ENDDO
    ENDDO

    !$xmp loop (j,k) on tmplt_dyn(j,k) expand(2:2,2:3)
    DO k=chunk%y_min-2,chunk%y_max+3
      DO j=chunk%x_min-2,chunk%x_max+2
        vol_flux_y(j,k)=0.0
        mass_flux_y(j,k)=0.0
        yarea(j,k)=0.0
      ENDDO
    ENDDO

    !$xmp loop (j) on tmplt_dyn(j,*) expand(2:2,0)
    DO j=chunk%x_min-2,chunk%x_max+2
      cellx(j)=0.0
      celldx(j)=0.0
    ENDDO

    !$xmp loop (k) on tmplt_dyn(*,k) expand(0,2:2)
    DO k=chunk%y_min-2,chunk%y_max+2
      celly(k)=0.0
      celldy(k)=0.0
    ENDDO

    !$xmp loop (j) on tmplt_dyn(j,*) expand(2:3,0)
    DO j=chunk%x_min-2,chunk%x_max+3
      vertexx(j)=0.0
      vertexdx(j)=0.0
    ENDDO

    !$xmp loop (k) on tmplt_dyn(*,k) expand(0,2:3)
    DO k=chunk%y_min-2,chunk%y_max+3
      vertexy(k)=0.0
      vertexdy(k)=0.0
    ENDDO

END SUBROUTINE build_field
