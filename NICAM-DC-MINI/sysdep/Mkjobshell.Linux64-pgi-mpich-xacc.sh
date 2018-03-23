#! /bin/bash -x

GLEV=${1}
RLEV=${2}
NMPI=${3}
ZL=${4}
VGRID=${5}
TOPDIR=${6}
BINNAME=${7}
RUNCONF=${8}

# System specific
MPIEXEC="mpirun -np ${NMPI}"

GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${NMPI} -ge 10000 ]; then
	NP=`printf %05d ${NMPI}`
elif [ ${NMPI} -ge 1000 ]; then
	NP=`printf %04d ${NMPI}`
elif [ ${NMPI} -ge 100 ]; then
	NP=`printf %03d ${NMPI}`
else
	NP=`printf %02d ${NMPI}`
fi

dir2d=gl${GL}rl${RL}pe${NP}
dir3d=gl${GL}rl${RL}z${ZL}pe${NP}
res2d=GL${GL}RL${RL}
res3d=GL${GL}RL${RL}z${ZL}

MNGINFO=rl${RL}-prc${NP}.info

outdir=${dir3d}
cd ${outdir}

#=- HA-PACS/TCA ---------------------------------------------------------------=
cat << EOFHAPACS1 > run_hapacs.sh
#PBS -S /bin/bash
#PBS -N xacc-nicamdc
#PBS -A XMPTCA
#PBS -q tcaq
#PBS -l select=3:ncpus=4:mpiprocs=4
#PBS -l walltime=00:01:00

. /opt/Modules/default/init/bash

module purge
module load pgi mvapich2/2.2_pgi_medium_cuda-8.0.44
module list

cd \$PBS_O_WORKDIR
export FORT_FMT_RECL=400
#export OMP_NUM_THREADS=1
export PGI_ACC_TIME=1

ln -s ${TOPDIR}/bin/${BINNAME} .
ln -s ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -s ${TOPDIR}/data/grid/vgrid/${VGRID} .

EOFHAPACS1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -s ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run_hapacs.sh
done

cat << EOFHAPACS2 >> run_hapacs.sh

# run
### OPT="MV2_ENABLE_AFFINITY=0 MV2_SHOW_CPU_BINDING=1 MV2_USE_CUDA=1 MV2_USE_GPUDIRECT=1 MV2_NUM_PORTS=1"
OPT="MV2_ENABLE_AFFINITY=0 MV2_SHOW_CPU_BINDING=1 MV2_USE_CUDA=1 MV2_USE_GPUDIRECT=1"
mpirun_rsh -np ${NMPI} -hostfile \$PBS_NODEFILE \$OPT ./numa.sh ./${BINNAME}

exit \$?

EOFHAPACS2
#=-----------------------------------------------------------------------------=
cat << EOFHAPACS3 > numa.sh
LOCAL_RANK=\$MV2_COMM_WORLD_LOCAL_RANK
SOCKET=\$(expr \$LOCAL_RANK / 2)
numactl --cpunodebind=\$SOCKET --localalloc \$@
EOFHAPACS3
#=-----------------------------------------------------------------------------=
chmod +x numa.sh
cat << EOF1 > run.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & mpich2 -----
#
################################################################################
export FORT_FMT_RECL=400


ln -sv ${TOPDIR}/bin/${BINNAME} .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/grid/vgrid/${VGRID} .
EOF1

for f in $( ls ${TOPDIR}/data/grid/boundary/${dir2d} )
do
   echo "ln -sv ${TOPDIR}/data/grid/boundary/${dir2d}/${f} ." >> run.sh
done

cat << EOF2 >> run.sh

# run
${MPIEXEC} ./${BINNAME} || exit

################################################################################
EOF2

exit
cat << EOFICO2LL1 > ico2ll.sh
#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & intel C&fortran & mpich2 -----
#
################################################################################
export FORT_FMT_RECL=400


ln -sv ${TOPDIR}/bin/fio_ico2ll_mpi .
ln -sv ${TOPDIR}/data/mnginfo/${MNGINFO} .
ln -sv ${TOPDIR}/data/zaxis .
EOFICO2LL1

for f in $( ls ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/ )
do
   echo "ln -sv ${TOPDIR}/data/grid/llmap/gl${GL}/rl${RL}/${f} ." >> ico2ll.sh
done

cat << EOFICO2LL2 >> ico2ll.sh

# run
${MPIEXEC} ./fio_ico2ll_mpi \
history \
glevel=${GLEV} \
rlevel=${RLEV} \
mnginfo="./${MNGINFO}" \
layerfile_dir="./zaxis" \
llmap_base="./llmap" \
-lon_swap \
-comm_smallchunk

################################################################################
EOFICO2LL2
