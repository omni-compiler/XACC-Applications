#Crown Copyright 2012 AWE.
#
# This file is part of CloverLeaf.
#
# CloverLeaf is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version.
#
# CloverLeaf is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# CloverLeaf. If not, see http://www.gnu.org/licenses/.

#  @brief Makefile for CloverLeaf
#  @author Wayne Gaudin, Andy Herdman
#  @details Agnostic, platform independent makefile for the Clover Leaf benchmark code.

# It is not meant to be clever in anyway, just a simple build out of the box script.
# Just make sure mpif90 is in your path. It uses mpif90 even for all builds because this abstracts the base
#  name of the compiler. If you are on a system that doesn't use mpif90, just replace mpif90 with the compiler name
#  of choice. The only mpi dependencies in this non-MPI version are mpi_wtime in timer.f90.

# There is no single way of turning OpenMP compilation on with all compilers.
# The known compilers have been added as a variable. By default the make
#  will use no options, which will work on Cray for example, but not on other
#  compilers.
# To select a OpenMP compiler option, do this in the shell before typing make:-
#
#  export COMPILER=INTEL       # to select the Intel flags
#  export COMPILER=SUN         # to select the Sun flags
#  export COMPILER=GNU         # to select the Gnu flags
#  export COMPILER=CRAY        # to select the Cray flags
#  export COMPILER=PGI         # to select the PGI flags
#  export COMPILER=PATHSCALE   # to select the Pathscale flags
#  export COMPILER=XL          # to select the IBM Xlf flags

# or this works as well:-
#
# make COMPILER=INTEL
# make COMPILER=SUN
# make COMPILER=GNU
# make COMPILER=CRAY
# make COMPILER=PGI
# make COMPILER=PATHSCALE
# make COMPILER=XL
#

# Don't forget to set the number of threads you want to use, like so
# export OMP_NUM_THREADS=4

# usage: make                     # Will make the binary
#        make clean               # Will clean up the directory
#        make DEBUG=1             # Will select debug options. If a compiler is selected, it will use compiler specific debug options
#        make IEEE=1              # Will select debug options as long as a compiler is selected as well
# e.g. make COMPILER=INTEL MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc DEBUG=1 IEEE=1 # will compile with the intel compiler with intel debug and ieee flags included

ifndef COMPILER
	COMPILER=PGI
#  MESSAGE=select a compiler to compile in OpenMP, e.g. make COMPILER=INTEL
endif

OMNI_HOME=$(dir $(shell which xmpcc))/..

OMP_INTEL     = -openmp
OMP_SUN       = -xopenmp=parallel -vpara
OMP_GNU       = -fopenmp
OMP_CRAY      =
OMP_PGI       = -mp=nonuma
OMP_PATHSCALE = -mp
OMP_XL        = -qsmp=omp -qthreaded
OMP=$(OMP_$(COMPILER))

FLAGS_INTEL     = -O3 -no-prec-div
FLAGS_SUN       = -fast -xipo=2 -Xlistv4
FLAGS_GNU       = -O3 -march=native -funroll-loops
FLAGS_CRAY      = -em -ra -h acc_model=fast_addr:no_deep_copy:auto_async_all
FLAGS_PGI       = -fastsse -Mipa=fast -Mlist -Minfo=acc -ta=nvidia,cc35,nordc,pinned #"-acc" is not required
FLAGS_PATHSCALE = -O3
FLAGS_XL        = -O5 -qipa=partition=large -g -qfullpath -Q -qsigtrap -qextname=flush:ideal_gas_kernel_c:viscosity_kernel_c:pdv_kernel_c:revert_kernel_c:accelerate_kernel_c:flux_calc_kernel_c:advec_cell_kernel_c:advec_mom_kernel_c:reset_field_kernel_c:timer_c:unpack_top_bottom_buffers_c:pack_top_bottom_buffers_c:unpack_left_right_buffers_c:pack_left_right_buffers_c:field_summary_kernel_c:update_halo_kernel_c:generate_chunk_kernel_c:initialise_chunk_kernel_c:calc_dt_kernel_c:clover_unpack_message_bottom_c:clover_pack_message_bottom_c:clover_unpack_message_top_c:clover_pack_message_top_c:clover_unpack_message_right_c:clover_pack_message_right_c:clover_unpack_message_left_c:clover_pack_message_left_c -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036FLAGS_          = -O3
CFLAGS_INTEL     = -O3 -no-prec-div -restrict -fno-alias
CFLAGS_SUN       = -fast -xipo=2
CFLAGS_GNU       = -O3 -march=native -funroll-loops
CFLAGS_CRAY      = -em -h list=a
CFLAGS_PGI       = -fastsse -Mipa=fast -Mlist
CFLAGS_PATHSCALE = -O3
CFLAGS_XL       = -O5 -qipa=partition=large -g -qfullpath -Q -qlistopt -qattr=full -qlist -qreport -qxref=full -qsource -qsuppress=1506-224:1500-036 -qsrcmsg
CFLAGS_          = -O3

ifdef DEBUG
  FLAGS_INTEL     = -O0 -g -debug all -check all -traceback -check noarg_temp_created
  FLAGS_SUN       = -g -xopenmp=noopt -stackvar -u -fpover=yes -C -ftrap=common
  FLAGS_GNU       = -O0 -g -O -Wall -Wextra -fbounds-check
  FLAGS_CRAY      = -O0 -g -em -eD
  FLAGS_PGI       = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk -Mchkptr -Minfo=acc -ta=nvidia,cc35
  FLAGS_PATHSCALE = -O0 -g
  FLAGS_XL       = -O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qinit=f90ptr -qsigtrap -qextname=flush:ideal_gas_kernel_c:viscosity_kernel_c:pdv_kernel_c:revert_kernel_c:accelerate_kernel_c:flux_calc_kernel_c:advec_cell_kernel_c:advec_mom_kernel_c:reset_field_kernel_c:timer_c:unpack_top_bottom_buffers_c:pack_top_bottom_buffers_c:unpack_left_right_buffers_c:pack_left_right_buffers_c:field_summary_kernel_c:update_halo_kernel_c:generate_chunk_kernel_c:initialise_chunk_kernel_c:calc_dt_kernel_c
  FLAGS_          = -O0 -g
  CFLAGS_INTEL    = -O0 -g -debug all -traceback
  CFLAGS_SUN      = -g -O0 -xopenmp=noopt -stackvar -u -fpover=yes -C -ftrap=common
  CFLAGS_GNU       = -O0 -g -O -Wall -Wextra -fbounds-check
  CFLAGS_CRAY     = -O0 -g -em -eD
  CFLAGS_PGI      = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk
  CFLAGS_PATHSCALE= -O0 -g
  CFLAGS_XL      = -O0 -g -qfullpath -qcheck -qflttrap=ov:zero:invalid:en -qsource -qinitauto=FF -qmaxmem=-1 -qsrcmsg
endif

ifdef IEEE
  I3E_INTEL     = -fp-model strict -fp-model source -prec-div -prec-sqrt
  I3E_SUN       = -fsimple=0 -fns=no
  I3E_GNU       = -ffloat-store
  I3E_CRAY      = -hflex_mp=intolerant
  I3E_PGI       = -Kieee
  I3E_PATHSCALE = -mieee-fp
  I3E_XL       = -qfloat=nomaf
  I3E=$(I3E_$(COMPILER))
endif

FLAGS=$(FLAGS_$(COMPILER)) $(OMP) $(I3E) $(OPTIONS) -xacc -lstdc++
CFLAGS=$(CFLAGS_$(COMPILER)) $(OMP) $(I3E) $(C_OPTIONS) -c
MPI_FLAGS=$(FLAGS_$(COMPILER)) $(OMP) $(I3E) $(OPTIONS) -I$(OMNI_HOME)/include
MPI_CFLAGS=$(CFLAGS_$(COMPILER)) $(OMP) $(I3E) $(C_OPTIONS) -c
XMP_COMPILER=xmpf90
C_XMP_COMPILER=xmpcc
MPI_COMPILER=mpif90
C_MPI_COMPILER=mpicc

OBJS=	data.o				\
	definitions.o			\
	clover.o			\
	report.o			\
	timer.o				\
	parse.o				\
	read_input.o			\
	initialise_chunk_kernel.o	\
	initialise_chunk.o		\
	build_field.o			\
	update_halo_kernel.o		\
	update_halo.o			\
	ideal_gas_kernel.o		\
	ideal_gas.o			\
	start.o				\
	generate_chunk_kernel.o		\
	generate_chunk.o		\
	initialise.o			\
	field_summary_kernel.o		\
	field_summary.o			\
	viscosity_kernel.o		\
	viscosity.o			\
	calc_dt_kernel.o		\
	calc_dt.o			\
	timestep.o			\
	accelerate_kernel.o		\
	accelerate.o			\
	revert_kernel.o			\
	revert.o			\
	PdV_kernel.o			\
	PdV.o				\
	flux_calc_kernel.o		\
	flux_calc.o			\
	advec_cell_kernel.o		\
	advec_cell_driver.o		\
	advec_mom_kernel.o		\
	advec_mom_driver.o		\
	advection.o			\
	reset_field_kernel.o		\
	reset_field.o			\
	hydro.o				\
	visit.o				\
	clover_leaf.o			\
	xmp_defs.o \
	timer_c.o

clover_leaf: ${OBJS} Makefile
	$(XMP_COMPILER) $(FLAGS) ${OBJS} \
	-o clover_leaf; echo $(MESSAGE)

.SUFFIXES: .f90 .c

.f90.o:
	rm -f `basename $< .f90`.s
	$(XMP_COMPILER) -c $(FLAGS) $<
	rm -f `basename $< .f90`.s

# read_input.o: read_input.f90
# 	rm -f `basename $< .f90`.s
# 	$(MPI_COMPILER) -c $(MPI_FLAGS) $<
# 	rm -f `basename $< .f90`.s

%.o: %.c
	$(C_XMP_COMPILER) $(CFLAGS) $<

%_module.mod: %.f90 %.o
	@:

%_mod.mod: %.f90 %.o
	@:

pdv_module.mod: PdV.f90 PdV.o
	@:

pdv_kernel_module.mod: PdV_kernel.f90 PdV_kernel.o
	@:

clean:
	rm -f *.o *.mod *genmod* *cuda* *hmd* *.cu *.oo *.hmf *.lst *.cub *.ptx *.cl clover_leaf *.xmod *.s


PdV.o: clover_module.mod report_module.mod pdv_kernel_module.mod revert_module.mod update_halo_module.mod ideal_gas_module.mod xmp_defs_module.mod
PdV_kernel.o: xmp_defs_module.mod
accelerate.o: clover_module.mod accelerate_kernel_module.mod xmp_defs_module.mod
accelerate_kernel.o: xmp_defs_module.mod
advec_cell_driver.o: clover_module.mod advec_cell_kernel_module.mod xmp_defs_module.mod
advec_cell_kernel.o: xmp_defs_module.mod
advec_mom_driver.o: clover_module.mod advec_mom_kernel_mod.mod xmp_defs_module.mod
advec_mom_kernel.o: xmp_defs_module.mod
advection.o: clover_module.mod advec_cell_driver_module.mod advec_mom_driver_module.mod update_halo_module.mod
build_field.o: clover_module.mod xmp_defs_module.mod
calc_dt.o: clover_module.mod calc_dt_kernel_module.mod xmp_defs_module.mod
calc_dt_kernel.o: xmp_defs_module.mod
clover.o: data_module.mod definitions_module.mod xmp_defs_module.mod
clover_leaf.o: clover_module.mod
data.o:
definitions.o: data_module.mod xmp_defs_module.mod
field_summary.o: clover_module.mod ideal_gas_module.mod field_summary_kernel_module.mod xmp_defs_module.mod
field_summary_kernel.o: xmp_defs_module.mod
flux_calc.o: clover_module.mod flux_calc_kernel_module.mod xmp_defs_module.mod
flux_calc_kernel.o: xmp_defs_module.mod
generate_chunk.o: clover_module.mod generate_chunk_kernel_module.mod xmp_defs_module.mod
generate_chunk_kernel.o: xmp_defs_module.mod
hydro.o: clover_module.mod timestep_module.mod viscosity_module.mod pdv_module.mod accelerate_module.mod flux_calc_module.mod advection_module.mod reset_field_module.mod
ideal_gas.o: clover_module.mod ideal_gas_kernel_module.mod xmp_defs_module.mod
ideal_gas_kernel.o: xmp_defs_module.mod
initialise.o: clover_module.mod parse_module.mod report_module.mod
initialise_chunk.o: clover_module.mod initialise_chunk_kernel_module.mod xmp_defs_module.mod
initialise_chunk_kernel.o: xmp_defs_module.mod
parse.o: data_module.mod report_module.mod clover_module.mod
read_input.o: clover_module.mod parse_module.mod report_module.mod
report.o: data_module.mod clover_module.mod
reset_field.o: clover_module.mod reset_field_kernel_module.mod xmp_defs_module.mod
reset_field_kernel.o: xmp_defs_module.mod
revert.o: clover_module.mod revert_kernel_module.mod xmp_defs_module.mod
revert_kernel.o: xmp_defs_module.mod
start.o: clover_module.mod parse_module.mod update_halo_module.mod ideal_gas_module.mod
timer.o:
timestep.o: clover_module.mod report_module.mod update_halo_module.mod viscosity_module.mod calc_dt_module.mod ideal_gas_module.mod definitions_module.mod
update_halo.o: clover_module.mod update_halo_kernel_module.mod xmp_defs_module.mod #update_tile_halo_module.mod
update_halo_kernel.o: data_module.mod xmp_defs_module.mod
# update_tile_halo.o: clover_module.mod update_tile_halo_kernel_module.mod
# update_tile_halo_kernel.o: data_module.mod
viscosity.o: clover_module.mod viscosity_kernel_module.mod xmp_defs_module.mod
viscosity_kernel.o: xmp_defs_module.mod
visit.o: clover_module.mod update_halo_module.mod viscosity_module.mod ideal_gas_module.mod xmp_defs_module.mod
xmp_defs.o:
