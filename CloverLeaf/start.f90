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

!>  @brief Main set up routine
!>  @author Wayne Gaudin
!>  @details Invokes the mesh decomposer and sets up chunk connectivity. It then
!>  allocates the communication buffers and call the chunk initialisation and
!>  generation routines. It calls the equation of state to calculate initial
!>  pressure before priming the halo cells and writing an initial field summary.

SUBROUTINE start

  USE clover_module
  USE parse_module
  USE update_halo_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: c

  INTEGER :: x_cells,y_cells
  INTEGER:: right,left,top,bottom

  INTEGER :: fields(NUM_FIELDS) !, chunk_task_responsible_for

  LOGICAL :: profiler_off

  !$xmp task on p(1,1) nocomm
  WRITE(g_out,*) 'Setting up initial geometry'
  WRITE(g_out,*)
  !$xmp end task

  time  = 0.0
  step  = 0
  dtold = dtinit
  dt    = dtinit

  !$xmp barrier

  CALL clover_get_num_chunks(number_of_chunks)


  !create the chunks

  chunk%task = xmp_node_num() - 1

  !chunk_task_responsible_for = parallel%task+1

  chunk%left   = 1
  chunk%bottom = 1
  chunk%right  = grid%x_cells
  chunk%top    = grid%y_cells
  chunk%x_min = 1
  chunk%y_min = 1
  chunk%x_max = grid%x_cells
  chunk%y_max = grid%y_cells

  CALL build_field()

  !$xmp barrier

  !$xmp task on p(1,1) nocomm
  WRITE(g_out,*) 'Generating chunks'
  !$xmp end task

  !$ACC ENTER DATA &
  !$ACC COPYIN(density0)   &
  !$ACC COPYIN(density1)   &
  !$ACC COPYIN(energy0)    &
  !$ACC COPYIN(energy1)    &
  !$ACC COPYIN(pressure)   &
  !$ACC COPYIN(soundspeed) &
  !$ACC COPYIN(viscosity0) &
  !$ACC COPYIN(xvel0)      &
  !$ACC COPYIN(yvel0)      &
  !$ACC COPYIN(xvel1)      &
  !$ACC COPYIN(yvel1)      &
  !$ACC COPYIN(vol_flux_x) &
  !$ACC COPYIN(vol_flux_y) &
  !$ACC COPYIN(mass_flux_x)&
  !$ACC COPYIN(mass_flux_y)&
  !$ACC COPYIN(volume)     &
  !$ACC COPYIN(work_array1)&
  !$ACC COPYIN(work_array2)&
  !$ACC COPYIN(work_array3)&
  !$ACC COPYIN(work_array4)&
  !$ACC COPYIN(work_array5)&
  !$ACC COPYIN(work_array6)&
  !$ACC COPYIN(work_array7)&
  !$ACC COPYIN(cellx)      &
  !$ACC COPYIN(celly)      &
  !$ACC COPYIN(celldx)     &
  !$ACC COPYIN(celldy)     &
  !$ACC COPYIN(vertexx)    &
  !$ACC COPYIN(vertexdx)   &
  !$ACC COPYIN(vertexy)    &
  !$ACC COPYIN(vertexdy)   &
  !$ACC COPYIN(xarea)      &
  !$ACC COPYIN(yarea)

  CALL initialise_chunk()
  CALL generate_chunk()

  advect_x=.TRUE.

  !$xmp barrier

  ! Do no profile the start up costs otherwise the total times will not add up
  ! at the end
  profiler_off=profiler_on
  profiler_on=.FALSE.


  CALL ideal_gas(.FALSE.)

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY0)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_PRESSURE)=1
  fields(FIELD_VISCOSITY)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1

  CALL update_halo(fields,2)

  !$xmp task on p(1,1) nocomm
  WRITE(g_out,*)
  WRITE(g_out,*) 'Problem initialised and generated'
  !$xmp end task

  CALL field_summary()

  IF(visit_frequency.NE.0) CALL visit()

  !$xmp barrier

  profiler_on=profiler_off

END SUBROUTINE start
