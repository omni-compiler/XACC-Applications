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

!>  @brief Controls the main hydro cycle.
!>  @author Wayne Gaudin
!>  @details Controls the top level cycle, invoking all the drivers and checks
!>  for outputs and completion.

SUBROUTINE hydro

  USE clover_module
  USE timestep_module
  USE viscosity_module
  USE PdV_module
  USE accelerate_module
  USE flux_calc_module
  USE advection_module
  USE reset_field_module

  IMPLICIT NONE

  INTEGER         :: loc(1)
  REAL(KIND=8)    :: timer,timerstart,wall_clock,step_clock

  REAL(KIND=8)    :: grind_time,cells,rstep
  REAL(KIND=8)    :: step_time,step_grind
  REAL(KIND=8)    :: first_step,second_step
  REAL(KIND=8)    :: kernel_total,totals(num_nodes)

  ! this is warmup and not necessary
  !$xmp reflect(density0, energy0, pressure, viscosity0) width(1,1) acc
  !$xmp reflect(density1, energy1) width(2,2) acc
  !$xmp reflect(xvel0, yvel0) width(1:2,1:2) acc
  !$xmp reflect(xvel1, yvel1) width(2:3,2:3) acc
  !$xmp reflect(vol_flux_x, mass_flux_x) width(2:3,2:2) acc
  !$xmp reflect(vol_flux_y, mass_flux_y) width(2:2,2:3) acc
!  call xacc_reflect_reset_time()
  timerstart = timer()

  !FIXME it can be removed
  !$ACC DATA &
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
  DO

    step_time = timer()

    step = step + 1

    CALL timestep()

    CALL PdV(.TRUE.)

    CALL accelerate()

    CALL PdV(.FALSE.)

    CALL flux_calc()

    CALL advection()

    CALL reset_field()

    advect_x = .NOT. advect_x

    time = time + dt

    IF(summary_frequency.NE.0) THEN
      IF(MOD(step, summary_frequency).EQ.0) CALL field_summary()
    ENDIF
    IF(visit_frequency.NE.0) THEN
      IF(MOD(step, visit_frequency).EQ.0) CALL visit()
    ENDIF

    ! Sometimes there can be a significant start up cost that appears in the first step.
    ! Sometimes it is due to the number of MPI tasks, or OpenCL kernel compilation.
    ! On the short test runs, this can skew the results, so should be taken into account
    !  in recorded run times.
    IF(step.EQ.1) first_step=(timer() - step_time)
    IF(step.EQ.2) second_step=(timer() - step_time)

    IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

      complete=.TRUE.
      CALL field_summary()
      IF(visit_frequency.NE.0) CALL visit()

      wall_clock=timer() - timerstart
      !$xmp task on p(1,1) nocomm
      WRITE(g_out,*)
      WRITE(g_out,*) 'Calculation complete'
      WRITE(g_out,*) 'Clover is finishing'
      WRITE(g_out,*) 'Wall clock ', wall_clock
      WRITE(g_out,*) 'First step overhead', first_step-second_step
      WRITE(    0,*) 'Wall clock ', wall_clock
      WRITE(    0,*) 'First step overhead', first_step-second_step
      !$xmp end task

      IF ( profiler_on ) THEN
        ! First we need to find the maximum kernel time for each task. This
        ! seems to work better than finding the maximum time for each kernel and
        ! adding it up, which always gives over 100%. I think this is because it
        ! does not take into account compute overlaps before syncronisations
        ! caused by halo exhanges.
        kernel_total=profiler%timestep+profiler%ideal_gas+profiler%viscosity+profiler%PdV          &
          +profiler%revert+profiler%acceleration+profiler%flux+profiler%cell_advection   &
          +profiler%mom_advection+profiler%reset+profiler%summary+profiler%visit         &
          +profiler%self_halo_exchange+profiler%mpi_halo_exchange
        CALL clover_allgather(kernel_total,totals)
        ! So then what I do is use the individual kernel times for the
        ! maximum kernel time task for the profile print
        loc=1 !MAXLOC(totals)
        kernel_total=totals(loc(1))
        CALL clover_allgather(profiler%timestep,totals)
        profiler%timestep=totals(loc(1))
        CALL clover_allgather(profiler%ideal_gas,totals)
        profiler%ideal_gas=totals(loc(1))
        CALL clover_allgather(profiler%viscosity,totals)
        profiler%viscosity=totals(loc(1))
        CALL clover_allgather(profiler%PdV,totals)
        profiler%PdV=totals(loc(1))
        CALL clover_allgather(profiler%revert,totals)
        profiler%revert=totals(loc(1))
        CALL clover_allgather(profiler%acceleration,totals)
        profiler%acceleration=totals(loc(1))
        CALL clover_allgather(profiler%flux,totals)
        profiler%flux=totals(loc(1))
        CALL clover_allgather(profiler%cell_advection,totals)
        profiler%cell_advection=totals(loc(1))
        CALL clover_allgather(profiler%mom_advection,totals)
        profiler%mom_advection=totals(loc(1))
        CALL clover_allgather(profiler%reset,totals)
        profiler%reset=totals(loc(1))
        CALL clover_allgather(profiler%self_halo_exchange,totals)
        profiler%self_halo_exchange=totals(loc(1))
        CALL clover_allgather(profiler%mpi_halo_exchange,totals)
        profiler%mpi_halo_exchange=totals(loc(1))
!        call xacc_reflect_get_time_cuda(profiler%mpi_halo_exchange_packunpack)
!        call xacc_reflect_get_time_mpi(profiler%mpi_halo_exchange_sendrecv)
        CALL clover_allgather(profiler%mpi_halo_exchange_packunpack,totals)
        profiler%mpi_halo_exchange_packunpack=totals(loc(1))
        CALL clover_allgather(profiler%mpi_halo_exchange_sendrecv,totals)
        profiler%mpi_halo_exchange_sendrecv=totals(loc(1))
        CALL clover_allgather(profiler%summary,totals)
        profiler%summary=totals(loc(1))
        CALL clover_allgather(profiler%visit,totals)
        profiler%visit=totals(loc(1))

        !$xmp task on p(1,1) nocomm
        WRITE(g_out,*)
        WRITE(g_out,'(a58,2f16.4)')"Profiler Output                 Time            Percentage"
        WRITE(g_out,'(a23,2f16.4)')"Timestep              :",profiler%timestep,&
          100.0*(profiler%timestep/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Ideal Gas             :",profiler%ideal_gas,&
          100.0*(profiler%ideal_gas/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Viscosity             :",profiler%viscosity,&
          100.0*(profiler%viscosity/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"PdV                   :",profiler%PdV,&
          100.0*(profiler%PdV/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Revert                :",profiler%revert,&
          100.0*(profiler%revert/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Acceleration          :",profiler%acceleration,&
          100.0*(profiler%acceleration/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Fluxes                :",profiler%flux,&
          100.0*(profiler%flux/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Cell Advection        :",profiler%cell_advection,&
          100.0*(profiler%cell_advection/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Momentum Advection    :",profiler%mom_advection,&
          100.0*(profiler%mom_advection/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Reset                 :",profiler%reset,&
          100.0*(profiler%reset/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Summary               :",profiler%summary,&
          100.0*(profiler%summary/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Visit                 :",profiler%visit,&
          100.0*(profiler%visit/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Tile Halo Exchange    :",0.0,&
          100.0*(0.0/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Self Halo Exchange    :",profiler%self_halo_exchange,&
          100.0*(profiler%self_halo_exchange/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"MPI Halo Exchange     :",profiler%mpi_halo_exchange,&
          100.0*(profiler%mpi_halo_exchange/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"__pack/unpack         :",profiler%mpi_halo_exchange_packunpack,&
          100.0*(profiler%mpi_halo_exchange_packunpack/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"__send/recv           :",profiler%mpi_halo_exchange_sendrecv,&
          100.0*(profiler%mpi_halo_exchange_sendrecv/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"Total                 :",kernel_total,&
          100.0*(kernel_total/wall_clock)
        WRITE(g_out,'(a23,2f16.4)')"The Rest              :",wall_clock-kernel_total,&
          100.0*(wall_clock-kernel_total)/wall_clock
        !$xmp end task
      ENDIF

      CALL clover_finalize

      EXIT

    END IF

    !$xmp task on p(1,1) nocomm
    wall_clock=timer()-timerstart
    step_clock=timer()-step_time
    WRITE(g_out,*)"Wall clock ",wall_clock
    WRITE(0    ,*)"Wall clock ",wall_clock
    cells = grid%x_cells * grid%y_cells
    rstep = step
    grind_time   = wall_clock/(rstep * cells)
    step_grind   = step_clock/cells
    WRITE(0    ,*)"Average time per cell ",grind_time
    WRITE(g_out,*)"Average time per cell ",grind_time
    WRITE(0    ,*)"Step time per cell    ",step_grind
    WRITE(g_out,*)"Step time per cell    ",step_grind
    !$xmp end task

  END DO
  !FIXME it can be removed
  !$ACC END DATA

END SUBROUTINE hydro
