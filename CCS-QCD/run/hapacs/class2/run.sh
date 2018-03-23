#PBS -S /bin/bash
#PBS -N xacc-ccs-qcd
#PBS -A XMPTCA
#PBS -q tcaq-q1
#PBS -l select=16:ncpus=4:mpiprocs=4
#PBS -l walltime=00:01:00

. /opt/Modules/default/init/bash

module purge
module load pgi/17.10
module load mvapich2/2.2_pgi_cuda-8.0.44
module load omnicompiler-20180226/pgi17.10
module list

cd $PBS_O_WORKDIR

export PGI_ACC_TIME=1
OPT="MV2_ENABLE_AFFINITY=0 MV2_SHOW_CPU_BINDING=1 MV2_USE_CUDA=1 MV2_USE_GPUDIRECT=1"
LD=../../../src/ccs_qcd_solver_bench_class2
time -p mpirun_rsh -np 64 -hostfile $PBS_NODEFILE $OPT ./numa.sh ${LD}

exit $?
