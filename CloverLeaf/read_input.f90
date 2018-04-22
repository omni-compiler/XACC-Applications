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

!>  @brief Reads the user input
!>  @author Wayne Gaudin
!>  @details Reads and parses the user input from the processed file and sets
!>  the variables used in the generation phase. Default values are also set
!>  here.

SUBROUTINE read_input()

  USE clover_module
  USE parse_module
  USE report_module

  IMPLICIT NONE

  INTEGER            :: state,stat,state_max,n

  REAL(KIND=8) :: dx,dy

  CHARACTER(LEN=500) :: word

  test_problem=0

  state_max=0

  grid%xmin=  0.0_8
  grid%ymin=  0.0_8
  grid%xmax=100.0_8
  grid%ymax=100.0_8

  grid%x_cells=10
  grid%y_cells=10

  end_time=10.0_8
  end_step=g_ibig
  complete=.FALSE.

  visit_frequency=0
  summary_frequency=10


  dtinit=0.1_8
  dtmax=1.0_8
  dtmin=0.0000001_8
  dtrise=1.5_8
  dtc_safe=0.7_8
  dtu_safe=0.5_8
  dtv_safe=0.5_8
  dtdiv_safe=0.7_8

  use_fortran_kernels=.TRUE.
  use_C_kernels=.FALSE.
  profiler_on=.FALSE.
  profiler%timestep=0.0
  profiler%acceleration=0.0
  profiler%PdV=0.0
  profiler%cell_advection=0.0
  profiler%mom_advection=0.0
  profiler%viscosity=0.0
  profiler%ideal_gas=0.0
  profiler%visit=0.0
  profiler%summary=0.0
  profiler%reset=0.0
  profiler%revert=0.0
  profiler%flux=0.0
  profiler%self_halo_exchange=0.0
  profiler%mpi_halo_exchange=0.0
  profiler%mpi_halo_exchange_packunpack=0.0
  profiler%mpi_halo_exchange_sendrecv=0.0

  !$xmp task on p(1,1) nocomm
  WRITE(g_out,*) 'Reading input file'
  WRITE(g_out,*)
  !$xmp end task

  stat=parse_init(g_in,'*clover')

  DO
    stat=parse_getline(dummy)
    IF (stat.ne.0) exit
    DO
      word=parse_getword(.FALSE.)
      IF(word.EQ.'')EXIT
      IF (word.EQ.'state') THEN
        state_max=MAX(state_max,parse_getival(parse_getword(.TRUE.)))
        EXIT
      ENDIF
    ENDDO
  ENDDO

  number_of_states=state_max

  IF(number_of_states.LT.1) CALL report_error('read_input','No states defined.')

  stat=parse_init(g_in,'*clover')

  ALLOCATE(states(number_of_states))
  states(:)%defined=.FALSE.
  states(:)%energy=0.0
  states(:)%density=0.0
  states(:)%xvel=0.0
  states(:)%yvel=0.0

  DO
    stat=parse_getline(dummy)

    IF(stat.NE.0)EXIT

    DO
      word=parse_getword(.FALSE.)

      IF(word.EQ.'')EXIT
      SELECT CASE(word)
      CASE('initial_timestep')
        dtinit=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'initial_timestep ',dtinit
        !$xmp end task
      CASE('max_timestep')
        dtmax=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'max_timestep',dtinit
        !$xmp end task
      CASE('timestep_rise')
        dtrise=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'timestep_rise',dtrise
        !$xmp end task
      CASE('end_time')
        end_time=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'end_time',end_time
        !$xmp end task
      CASE('end_step')
        end_step=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'end_step',end_step
        !$xmp end task
      CASE('xmin')
        grid%xmin=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'xmin',grid%xmin
        !$xmp end task
      CASE('xmax')
        grid%xmax=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'xmax',grid%xmax
        !$xmp end task
      CASE('ymin')
        grid%ymin=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'ymin',grid%ymin
        !$xmp end task
      CASE('ymax')
        grid%ymax=parse_getrval(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,e12.4)")'ymax',grid%ymax
        !$xmp end task
      CASE('x_cells')
        grid%x_cells=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'x_cells',grid%x_cells
        !$xmp end task
      CASE('y_cells')
        grid%y_cells=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'y_cells',grid%y_cells
        !$xmp end task
      CASE('visit_frequency')
        visit_frequency=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'visit_frequency',visit_frequency
        !$xmp end task
      CASE('summary_frequency')
        summary_frequency=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'summary_frequency',summary_frequency
        !$xmp end task
      CASE('use_fortran_kernels')
        use_fortran_kernels=.TRUE.
        use_C_kernels=.FALSE.
      CASE('use_c_kernels')
        use_fortran_kernels=.FALSE.
        use_C_kernels=.TRUE.
      CASE('profiler_on')
        profiler_on=.TRUE.
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25)")'Profiler on'
        !$xmp end task
      CASE('test_problem')
        test_problem=parse_getival(parse_getword(.TRUE.))
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,"(1x,a25,i12)")'test_problem',test_problem
        !$xmp end task
      CASE('state')

        state=parse_getival(parse_getword(.TRUE.))

        !$xmp task on p(1,1) nocomm
        WRITE(g_out,*)'Reading specification for state ',state
        !$xmp end task
        IF (states(state)%defined) CALL report_error('read_input','State defined twice.')
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,*)
        !$xmp end task

        states(state)%defined=.TRUE.
        DO
          word=parse_getword(.FALSE.)
          IF(word.EQ.'') EXIT

          SELECT CASE(word)

          CASE('xvel')
            states(state)%xvel=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'xvel ',states(state)%xvel
            !$xmp end task
          CASE('yvel')
            states(state)%yvel=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'yvel ',states(state)%yvel
            !$xmp end task
          CASE('xmin')
            states(state)%xmin=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state xmin ',states(state)%xmin
            !$xmp end task
          CASE('ymin')
            states(state)%ymin=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state ymin ',states(state)%ymin
            !$xmp end task
          CASE('xmax')
            states(state)%xmax=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state xmax ',states(state)%xmax
            !$xmp end task
          CASE('ymax')
            states(state)%ymax=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state ymax ',states(state)%ymax
            !$xmp end task
          CASE('radius')
            states(state)%radius=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state radius ',states(state)%radius
            !$xmp end task
          CASE('density')
            states(state)%density=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state density ',states(state)%density
            !$xmp end task
          CASE('energy')
            states(state)%energy=parse_getrval(parse_getword(.TRUE.))
            !$xmp task on p(1,1) nocomm
            WRITE(g_out,"(1x,a25,e12.4)")'state energy ',states(state)%energy
            !$xmp end task
          CASE('geometry')
            word=TRIM(parse_getword(.TRUE.))
            SELECT CASE(word)
            CASE("rectangle")
              states(state)%geometry=g_rect
              !$xmp task on p(1,1) nocomm
              WRITE(g_out,"(1x,a26)")'state geometry rectangular'
              !$xmp end task
            CASE("circle")
              states(state)%geometry=g_circ
              !$xmp task on p(1,1) nocomm
              WRITE(g_out,"(1x,a25)")'state geometry circular'
              !$xmp end task
            CASE("point")
              states(state)%geometry=g_point
              !$xmp task on p(1,1) nocomm
              WRITE(g_out,"(1x,a25)")'state geometry point'
              !$xmp end task
            END SELECT
          END SELECT
        ENDDO
        !$xmp task on p(1,1) nocomm
        WRITE(g_out,*)
        !$xmp end task
      END SELECT
    ENDDO
  ENDDO

  !$xmp task on p(1,1) nocomm
  WRITE(g_out,*)
  IF(use_fortran_kernels) THEN
    WRITE(g_out,"(1x,a25)")'Using Fortran Kernels'
  ELSEIF(use_c_kernels) THEN
    WRITE(g_out,"(1x,a25)")'Using C Kernels'
  ENDIF
  WRITE(g_out,*)
  WRITE(g_out,*) 'Input read finished.'
  WRITE(g_out,*)
  !$xmp end task

  ! If a state boundary falls exactly on a cell boundary then round off can
  ! cause the state to be put one cell further that expected. This is compiler
  ! /system dependent. To avoid this, a state boundary is reduced/increased by a 100th
  ! of a cell width so it lies well with in the intended cell.
  ! Because a cell is either full or empty of a specified state, this small
  ! modification to the state extents does not change the answers.
  dx=(grid%xmax-grid%xmin)/float(grid%x_cells)
  dy=(grid%ymax-grid%ymin)/float(grid%y_cells)
  DO n=2,number_of_states
    states(n)%xmin=states(n)%xmin+(dx/100.0_8)
    states(n)%ymin=states(n)%ymin+(dy/100.0_8)
    states(n)%xmax=states(n)%xmax-(dx/100.0_8)
    states(n)%ymax=states(n)%ymax-(dy/100.0_8)
  ENDDO

END SUBROUTINE read_input
