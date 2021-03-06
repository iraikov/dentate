#!/bin/bash
#
#SBATCH -J dentate_Test_GC_slice_300um
#SBATCH -o ./results/dentate_Test_GC_slice_300um.%j.o
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=56
#SBATCH -p normal
#SBATCH -t 2:00:00
#SBATCH --mail-user=pmoolcha@stanford.edu
#SBATCH --mail-type=END
#

#module load python3
#module load phdf5
module load intel19

set -x

export MODEL_HOME=/scratch1/04119/pmoolcha/HDM
export DG_HOME=$MODEL_HOME/dentate
export LD_PRELOAD=$MKLROOT/lib/intel64_lin/libmkl_core.so:$MKLROOT/lib/intel64_lin/libmkl_sequential.so

#export I_MPI_EXTRA_FILESYSTEM=enable
#export I_MPI_ADJUST_ALLGATHER=4
#export I_MPI_ADJUST_ALLGATHERV=4
#export I_MPI_ADJUST_ALLTOALL=4
#export I_MPI_ADJUST_ALLTOALLV=2

export RAIKOVSCRATCH=/scratch1/03320/iraikov

results_path=results/Test_GC_slice_300um_$SLURM_JOB_ID
export results_path

mkdir -p $results_path

#git ls-files | tar -zcf ${results_path}/dentate.tgz --files-from=/dev/stdin
#git --git-dir=../dgc/.git ls-files | grep Mateos-Aparicio2014 | tar -C ../dgc -zcf ${results_path}/dgc.tgz --files-from=/dev/stdin

#ibrun -n 1 gdb -batch -x ../../debug/gdb_backtrace --args python3 ./scripts/main.py  \
ibrun python3 ./scripts/main.py  \
    --arena-id=A --trajectory-id=Diag \
    --config-file=Test_Slice_300um_IN_Izh.yaml \
    --config-prefix=./config \
    --template-paths=../dgc/Mateos-Aparicio2014:templates \
    --dataset-prefix "$RAIKOVSCRATCH/striped/dentate" \
    --results-path=$results_path \
    --io-size=8 \
    --tstop=9500 \
    --v-init=-75 \
    --max-walltime-hours=1.9 \
    --spike-input-path "$RAIKOVSCRATCH/striped/dentate/Full_Scale_Control/DG_input_spike_trains_20200910_compressed.h5" \
    --spike-input-namespace='Input Spikes A Diag' \
    --spike-input-attr='Spike Train' \
    --microcircuit-inputs \
    --checkpoint-interval 0. \
    --recording-fraction 0.01 \
    --use-coreneuron \
    --verbose

