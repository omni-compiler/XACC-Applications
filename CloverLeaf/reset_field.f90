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

!>  @brief Reset field driver
!>  @author Wayne Gaudin
!>  @details Invokes the user specified field reset kernel.

MODULE reset_field_module

CONTAINS

  SUBROUTINE reset_field()

    USE clover_module
    USE reset_field_kernel_module

    IMPLICIT NONE


    REAL(KIND=8) :: kernel_time,timer

    IF(profiler_on) kernel_time=timer()

    CALL reset_field_kernel(chunk%x_min, &
      chunk%x_max,                       &
      chunk%y_min,                       &
      chunk%y_max,                       &
      density0,                          &
      density1,                          &
      energy0,                           &
      energy1,                           &
      xvel0,                             &
      xvel1,                             &
      yvel0,                             &
      yvel1 )

    IF(profiler_on) profiler%reset=profiler%reset+(timer()-kernel_time)

  END SUBROUTINE reset_field

END MODULE reset_field_module
