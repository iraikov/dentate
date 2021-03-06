#!/bin/bash
#
#SBATCH -J dentate_Test_GC_slice_600um
#SBATCH -o ./results/dentate_Test_GC_slice_600um.%j.o
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=56
#SBATCH -p normal
#SBATCH -t 2:00:00
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#

module load python3
module load phdf5

set -x

export NEURONROOT=$HOME/bin/nrnpython3
export PYTHONPATH=$HOME/model:$NEURONROOT/lib/python:$SCRATCH/site-packages:$PYTHONPATH
export PATH=$NEURONROOT/x86_64/bin:$PATH
export MODEL_HOME=$HOME/model
export DG_HOME=$MODEL_HOME/dentate
export LD_PRELOAD=$MKLROOT/lib/intel64_lin/libmkl_core.so:$MKLROOT/lib/intel64_lin/libmkl_sequential.so

export I_MPI_EXTRA_FILESYSTEM=enable
export I_MPI_ADJUST_ALLGATHER=4
export I_MPI_ADJUST_ALLGATHERV=4
export I_MPI_ADJUST_ALLTOALL=4
export I_MPI_ADJUST_ALLTOALLV=2

results_path=$SCRATCH/striped/dentate/results/Test_GC_slice_600um_$SLURM_JOB_ID
export results_path

mkdir -p $results_path

#git ls-files | tar -zcf ${results_path}/dentate.tgz --files-from=/dev/stdin
#git --git-dir=../dgc/.git ls-files | grep Mateos-Aparicio2014 | tar -C ../dgc -zcf ${results_path}/dgc.tgz --files-from=/dev/stdin

ibrun env PYTHONPATH=$PYTHONPATH python3 ./scripts/main.py  \
    --arena-id=A --trajectory-id=Diag \
    --config-file=Test_Slice_600um_IN_Izh.yaml \
    --config-prefix=./config \
    --template-paths=../dgc/Mateos-Aparicio2014:templates \
    --dataset-prefix="$SCRATCH/striped/dentate" \
    --results-path=$results_path \
    --io-size=16 \
    --tstop=9500 \
    --v-init=-75 \
    --max-walltime-hours=1.9 \
    --spike-input-path "$SCRATCH/striped/dentate/Full_Scale_Control/DG_input_spike_trains_20200910_compressed.h5" \
    --spike-input-namespace='Input Spikes A Diag' \
    --spike-input-attr='Spike Train' \
    --microcircuit-inputs \
    --checkpoint-interval 0. \
    --use-coreneuron \
    --verbose

