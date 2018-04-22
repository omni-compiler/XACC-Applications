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

!>  @brief Driver for chunk initialisation.
!>  @author Wayne Gaudin
!>  @details Invokes the user specified chunk initialisation kernel.

SUBROUTINE initialise_chunk()

  USE clover_module
  USE initialise_chunk_kernel_module

  IMPLICIT NONE

  REAL(KIND=8) :: xmin,ymin,dx,dy

  dx=(grid%xmax-grid%xmin)/float(grid%x_cells)
  dy=(grid%ymax-grid%ymin)/float(grid%y_cells)

  xmin=grid%xmin
  ymin=grid%ymin

  CALL initialise_chunk_kernel(chunk%x_min, &
    chunk%x_max,                            &
    chunk%y_min,                            &
    chunk%y_max,                            &
    xmin,ymin,dx,dy,                        &
    vertexx,                                &
    vertexdx,                               &
    vertexy,                                &
    vertexdy,                               &
    cellx,                                  &
    celldx,                                 &
    celly,                                  &
    celldy,                                 &
    volume,                                 &
    xarea,                                  &
    yarea )

END SUBROUTINE initialise_chunk
