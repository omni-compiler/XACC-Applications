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

!>  @brief Driver for the halo updates
!>  @author Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_halo_module

CONTAINS

  SUBROUTINE update_halo(fields,depth)

    USE clover_module
    USE update_halo_kernel_module

    IMPLICIT NONE

    INTEGER :: fields(NUM_FIELDS),depth
    REAL(KIND=8) :: kernel_time,timer

    !TODO: fix the chunk comms phase

    IF(profiler_on) kernel_time=timer()

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
       !$xmp reflect(density0) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
       !$xmp reflect(density1) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
       !$xmp reflect(energy0) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
       !$xmp reflect(energy1) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
       !$xmp reflect(pressure) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
       !$xmp reflect(viscosity0) width(depth,depth) acc
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
       !$xmp reflect(soundspeed) width(depth,depth) acc
    ENDIF

    IF(fields(FIELD_XVEL0).EQ.1) THEN
       !$xmp reflect(xvel0) width(depth:depth+1, depth:depth+1) acc
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
       !$xmp reflect(xvel1) width(depth:depth+1, depth:depth+1) acc
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
       !$xmp reflect(yvel0) width(depth:depth+1, depth:depth+1) acc
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
       !$xmp reflect(yvel1) width(depth:depth+1, depth:depth+1) acc
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
       !$xmp reflect(vol_flux_x) width(depth:depth+1, depth:depth) acc
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
       !$xmp reflect(mass_flux_x) width(depth:depth+1, depth:depth) acc
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
       !$xmp reflect(vol_flux_y) width(depth:depth, depth:depth+1) acc
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
       !$xmp reflect(mass_flux_y) width(depth:depth, depth:depth+1) acc
    ENDIF

    IF(profiler_on) profiler%mpi_halo_exchange=profiler%mpi_halo_exchange+(timer()-kernel_time)
 
    IF(profiler_on) kernel_time=timer()

    CALL update_halo_kernel(chunk%x_min, &
      chunk%x_max,                       &
      chunk%y_min,                       &
      chunk%y_max,                       &
      density0,                          &
      energy0,                           &
      pressure,                          &
      viscosity0,                        &
      soundspeed,                        &
      density1,                          &
      energy1,                           &
      xvel0,                             &
      yvel0,                             &
      xvel1,                             &
      yvel1,                             &
      vol_flux_x,                        &
      vol_flux_y,                        &
      mass_flux_x,                       &
      mass_flux_y,                       &
      fields,                            &
      depth                              )

    IF(profiler_on) profiler%self_halo_exchange=profiler%self_halo_exchange+(timer()-kernel_time)

  END SUBROUTINE update_halo

END MODULE update_halo_module
