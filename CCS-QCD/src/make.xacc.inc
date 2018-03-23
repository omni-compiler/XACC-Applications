#$Id: make.inc,v 1.1 2009/12/02 10:24:37 ishikawa Exp $

REVISION = -D_REVISION_='"$(shell cat .svnrevision)"'

#
# specify Fortran90 MPI compiler FC, compiler options FFLAGS and linker options LDFLAGS
#

FC        = xmpf90
CC        = mpicc

FFLAGS   += -cpp -O2 -Mpreprocess -m64 
CFLAGS    = $(FFLAGS)

ACCFLAGS  = -D_OPENACC -Minfo=accel -xacc --Wn-acc --Wl-acc 
FFLAGS   += $(ACCFLAGS)

LDFLAGS   = $(FFLAGS)

INCLUDE   = -I./
