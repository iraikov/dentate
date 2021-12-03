#!/bin/bash
#
#SBATCH -J generate_distance_normal_weights_DG_GC
#SBATCH -o ./results/generate_normal_weights_DG_GC.%j.o
#SBATCH -N 40
#SBATCH --ntasks-per-node=56            # # of mpi tasks per node
#SBATCH -p normal      # Queue (partition) name
#SBATCH -t 1:00:00
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

module load gcc/9.1.0
module load python3
module load phdf5/1.10.4

set -x

export NEURONROOT=$SCRATCH/bin/nrnpython3_gcc9
export PYTHONPATH=$HOME/model:$NEURONROOT/lib/python:$SCRATCH/site-packages/gcc9:$PYTHONPATH
export DATA_PREFIX=$SCRATCH/striped2/dentate

export I_MPI_HYDRA_TOPOLIB=ipl
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=off
export I_MPI_HYDRA_BRANCH_COUNT=0
export UCX_TLS="knem,dc_x"

export I_MPI_ADJUST_ALLTOALL=4
export I_MPI_ADJUST_ALLTOALLV=2

cd $SLURM_SUBMIT_DIR

mpirun -rr python3 ./scripts/generate_normal_weights_as_cell_attr.py \
    -d GC -s MC -s ConMC \
    --config=Full_Scale_Basis.yaml \
    --config-prefix=./config \
    --weights-path=$DATA_PREFIX/Full_Scale_Control/DG_GC_syn_weights_LN_20210920.h5 \
    --connections-path=$DATA_PREFIX/Full_Scale_Control/DG_GC_connections_20210920_compressed.h5 \
    --io-size=40 --cache-size=40  --value-chunk-size=10000 --chunk-size=10000 --write-size=0 -v

