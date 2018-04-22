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

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE clover_abort

    INTEGER :: ierr,err

  END SUBROUTINE clover_abort

  SUBROUTINE clover_finalize

    CLOSE(g_out)
    CALL FLUSH(0)
    CALL FLUSH(6)
    CALL FLUSH(g_out)

  END SUBROUTINE clover_finalize

  SUBROUTINE clover_get_num_chunks(count)

    IMPLICIT NONE

    INTEGER :: count

    ! Should be changed so there can be more than one chunk per mpi task

    count=1

  END SUBROUTINE clover_get_num_chunks

  SUBROUTINE clover_allgather(value,values)

    IMPLICIT NONE

    REAL(KIND=8) :: value
    REAL(KIND=8) :: values(1)

    values(:) = 0.0
    values(xmp_node_num())=value ! Just to ensure it will work in serial
    !$xmp reduction(+:values)

  END SUBROUTINE clover_allgather

  SUBROUTINE clover_check_error(error)

    IMPLICIT NONE

    INTEGER :: error

    !$xmp reduction(max:error)
  END SUBROUTINE clover_check_error

END MODULE clover_module
