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

!>  @brief Cell centred advection driver.
!>  @author Wayne Gaudin
!>  @details Invokes the user selected advection kernel.

MODULE  advec_cell_driver_module

CONTAINS

  SUBROUTINE advec_cell_driver(sweep_number,dir)

    USE clover_module
    USE advec_cell_kernel_module

    IMPLICIT NONE

    INTEGER :: sweep_number,dir

    CALL advec_cell_kernel(chunk%x_min, &
      chunk%x_max,                      &
      chunk%y_min,                      &
      chunk%y_max,                      &
      dir,                              &
      sweep_number,                     &
      vertexdx,                         &
      vertexdy,                         &
      volume,                           &
      density1,                         &
      energy1,                          &
      mass_flux_x,                      &
      vol_flux_x,                       &
      mass_flux_y,                      &
      vol_flux_y,                       &
      work_array1,                      &
      work_array2,                      &
      work_array3,                      &
      work_array4,                      &
      work_array5,                      &
      work_array6,                      &
      work_array7                       )

  END SUBROUTINE advec_cell_driver

END MODULE  advec_cell_driver_module
