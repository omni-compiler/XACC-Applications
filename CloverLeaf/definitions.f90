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

!>  @brief Holds the high level Fortran data types
!>  @author Wayne Gaudin
!>  @details The high level data types used to store the mesh and field data
!>  are defined here.
!>
!>  Also the global variables used for defining the input and controlling the
!>  scheme are defined here.

MODULE definitions_module

  USE data_module
  USE xmp_defs_module

!  IMPLICIT NONE

  TYPE state_type
    LOGICAL            :: defined

    REAL(KIND=8)       :: density          &
      ,energy           &
      ,xvel             &
      ,yvel

    INTEGER            :: geometry

    REAL(KIND=8)       :: xmin               &
      ,xmax               &
      ,ymin               &
      ,ymax               &
      ,radius
  END TYPE state_type

  TYPE(state_type), ALLOCATABLE             :: states(:)
  INTEGER                                   :: number_of_states

  TYPE grid_type
    REAL(KIND=8)       :: xmin            &
      ,ymin            &
      ,xmax            &
      ,ymax

    INTEGER            :: x_cells              &
      ,y_cells
  END TYPE grid_type

  INTEGER      :: step

  LOGICAL      :: advect_x



  INTEGER      :: error_condition

  INTEGER      :: test_problem
  LOGICAL      :: complete

  LOGICAL      :: use_fortran_kernels
  LOGICAL      :: use_C_kernels

  LOGICAL      :: profiler_on ! Internal code profiler to make comparisons across systems easier

  TYPE profiler_type
    REAL(KIND=8)       :: timestep           &
      ,acceleration       &
      ,PdV                &
      ,cell_advection     &
      ,mom_advection      &
      ,viscosity          &
      ,ideal_gas          &
      ,visit              &
      ,summary            &
      ,reset              &
      ,revert             &
      ,flux               &
      ,self_halo_exchange &
      ,mpi_halo_exchange_packunpack &
      ,mpi_halo_exchange_sendrecv &
      ,mpi_halo_exchange

  END TYPE profiler_type
  TYPE(profiler_type)  :: profiler

  REAL(KIND=8) :: end_time

  INTEGER      :: end_step

  REAL(KIND=8) :: dtold          &
    ,dt             &
    ,time           &
    ,dtinit         &
    ,dtmin          &
    ,dtmax          &
    ,dtrise         &
    ,dtu_safe       &
    ,dtv_safe       &
    ,dtc_safe       &
    ,dtdiv_safe     &
    ,dtc            &
    ,dtu            &
    ,dtv            &
    ,dtdiv

  INTEGER      :: visit_frequency   &
    ,summary_frequency

  INTEGER         :: jdt,kdt


  TYPE chunk_type

    INTEGER         :: task   !mpi task


    ! Idealy, create an array to hold the buffers for each field so a commuincation only needs
    !  one send and one receive per face, rather than per field.
    ! If chunks are overloaded, i.e. more chunks than tasks, might need to pack for a task to task comm
    !  rather than a chunk to chunk comm. See how performance is at high core counts before deciding


    INTEGER         :: x_min  &
      ,y_min  &
      ,x_max  &
      ,y_max

    INTEGER         :: left &
      ,right  &
      ,bottom &
      ,top
  END TYPE chunk_type



  TYPE(chunk_type)       :: chunk

  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: density0,density1
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: energy0,energy1
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: pressure
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: viscosity0
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: soundspeed
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: xvel0,xvel1
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: yvel0,yvel1
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_x,mass_flux_x
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_y,mass_flux_y
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array1 !node_flux, stepbymass, volume_change, pre_vol
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array2 !node_mass_post, post_vol
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array3 !node_mass_pre,pre_mass
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array4 !advec_vel, post_mass
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array5 !mom_flux, advec_vol
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array6 !pre_vol, post_ener
  REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array7 !post_vol, ener_flux

  REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: cellx    &
    ,celly    &
    ,vertexx  &
    ,vertexy  &
    ,celldx   &
    ,celldy   &
    ,vertexdx &
    ,vertexdy

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: volume  &
    ,xarea   &
    ,yarea

  INTEGER                :: number_of_chunks

  TYPE(grid_type)        :: grid

!$xmp nodes p(*,*)
!$xmp template tmplt_dyn(:,:)
!$xmp distribute tmplt_dyn(gblock(*), gblock(*)) onto p
!$xmp align (i,j) with tmplt_dyn(i,j) :: density0, density1, energy0, energy1, pressure, viscosity0, soundspeed, volume, &
!$xmp                                    vol_flux_x, mass_flux_x, xarea, vol_flux_y, mass_flux_y, yarea, xvel0, xvel1, yvel0,yvel1,&
!$xmp                                    work_array1, work_array2, work_array3, work_array4, work_array5, work_array6, work_array7

!$xmp shadow (2:2, 2:2) :: density0, density1, energy0, energy1, pressure, viscosity0, soundspeed, volume
!$xmp shadow (2:3, 2:2) :: vol_flux_x, mass_flux_x, xarea
!$xmp shadow (2:2, 2:3) :: vol_flux_y, mass_flux_y, yarea
!$xmp shadow (2:3, 2:3) :: xvel0, xvel1, yvel0, yvel1, &
!$xmp                      work_array1, work_array2, work_array3, work_array4, work_array5, work_array6, work_array7

!$xmp align (i) with tmplt_dyn(i,*) :: cellx, celldx, vertexx, vertexdx
!$xmp shadow (2:2) :: cellx, celldx
!$xmp shadow (2:3) :: vertexx, vertexdx

!$xmp align (j) with tmplt_dyn(*,j) :: celly, celldy, vertexy, vertexdy
!$xmp shadow (2:2) :: celly, celldy
!$xmp shadow (2:3) :: vertexy, vertexdy

contains
  subroutine init_xmp_defs(xmin, xmax, ymin, ymax)
    INTEGER :: xmin, xmax, ymin, ymax
    INTEGER :: xmp_node_xsize, xmp_node_xidx, xmp_node_ysize, xmp_node_yidx
    external xmp_nodes_index, xmp_nodes_size

    num_nodes = xmp_num_nodes()

    call xmp_nodes_index(xmp_desc_of(p), 1, xmp_node_xidx)
    call xmp_nodes_index(xmp_desc_of(p), 2, xmp_node_yidx)
    call xmp_nodes_size(xmp_desc_of(p), 1, xmp_node_xsize)
    call xmp_nodes_size(xmp_desc_of(p), 2, xmp_node_ysize)

    g_l_min_x = xmax / xmp_node_xsize * (xmp_node_xidx - 1) + 1
    g_l_min_y = ymax / xmp_node_ysize * (xmp_node_yidx - 1) + 1
    g_l_max_x = xmax / xmp_node_xsize * (xmp_node_xidx)
    g_l_max_y = ymax / xmp_node_ysize * (xmp_node_yidx)

    allocate(blocksize_x(xmp_node_xsize))
    allocate(blocksize_y(xmp_node_ysize))
    blocksize_x = xmax / xmp_node_xsize
    blocksize_y = ymax / xmp_node_ysize
    blocksize_x(1) = blocksize_x(1) + 2
    blocksize_x(xmp_node_xsize) = blocksize_x(xmp_node_xsize) + 3
    blocksize_y(1) = blocksize_y(1) + 2
    blocksize_y(xmp_node_ysize) = blocksize_y(xmp_node_ysize) + 3

    ! IF (xmp_node_xidx == 1 .and. xmp_node_yidx == 1) THEN
    ! write(*,*) "#nodesX", xmp_node_xsize, "#nodesY", xmp_node_ysize
    ! write(*,*) "xsize/node", g_x_max / xmp_node_xsize, "ysize/node", g_y_max / xmp_node_ysize
    ! write(*,*) "blocksize_x", blocksize_x
    ! write(*,*) "blocksize_y", blocksize_y
    ! ENDIF
  end subroutine init_xmp_defs

END MODULE definitions_module
