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

!>  @brief Driver for the timestep kernels
!>  @author Wayne Gaudin
!>  @details Invokes the user specified timestep kernel.

MODULE calc_dt_module

CONTAINS

  SUBROUTINE calc_dt(local_dt,local_control,xl_pos,yl_pos,jldt,kldt)

    USE clover_module
    USE calc_dt_kernel_module

    IMPLICIT NONE

    REAL(KIND=8)     :: local_dt
    CHARACTER(LEN=8) :: local_control
    REAL(KIND=8)     :: xl_pos,yl_pos
    INTEGER          :: jldt,kldt

    INTEGER          :: l_control
    INTEGER          :: small

    local_dt=g_big

    small = 0

    CALL calc_dt_kernel(chunk%x_min, &
      chunk%x_max,                   &
      chunk%y_min,                   &
      chunk%y_max,                   &
      g_small,                       &
      g_big,                         &
      dtmin,                         &
      dtc_safe,                      &
      dtu_safe,                      &
      dtv_safe,                      &
      dtdiv_safe,                    &
      xarea,                         &
      yarea,                         &
      cellx,                         &
      celly,                         &
      celldx,                        &
      celldy,                        &
      volume,                        &
      density0,                      &
      energy0,                       &
      pressure,                      &
      viscosity0,                    &
      soundspeed,                    &
      xvel0,                         &
      yvel0,                         &
      work_array1,                   &
      local_dt,                      &
      l_control,                     &
      xl_pos,                        &
      yl_pos,                        &
      jldt,                          &
      kldt,                          &
      small                          )

    IF(l_control.EQ.1) local_control='sound'
    IF(l_control.EQ.2) local_control='xvel'
    IF(l_control.EQ.3) local_control='yvel'
    IF(l_control.EQ.4) local_control='div'

  END SUBROUTINE calc_dt

END MODULE calc_dt_module
