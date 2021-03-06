#!/bin/bash

#SBATCH -J dentate_Full_Scale_GC_Exc_Sat_LNN_dbg           # Job name
#SBATCH -o ./results/dentate.o%j       # Name of stdout output file
#SBATCH -e ./results/dentate.e%j       # Name of stderr error file
#SBATCH -p normal      # Queue (partition) name
#SBATCH -N 256            # Total # of nodes 
#SBATCH -n 14336           # Total # of mpi tasks
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=all    # Send email at begin and end of job

module unload python2
module unload intel
module load intel/18.0.5
module load python3
module load phdf5

set -x

export NEURONROOT=$HOME/bin/nrnpython3
export PYTHONPATH=$HOME/model:$NEURONROOT/lib/python:$PYTHONPATH
export PATH=$NEURONROOT/x86_64/bin:$PATH

results_path=$SCRATCH/dentate/results/Full_Scale_GC_Exc_Sat_LNN_$SLURM_JOB_ID
export results_path

cd $SLURM_SUBMIT_DIR

mkdir -p $results_path

git ls-files | tar -zcf ${results_path}/dentate.tgz --files-from=/dev/stdin
git --git-dir=../dgc/.git ls-files | grep Mateos-Aparicio2014 | tar -C ../dgc -zcf ${results_path}/dgc.tgz --files-from=/dev/stdin

export I_MPI_EXTRA_FILESYSTEM=enable
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre
export I_MPI_ADJUST_ALLGATHER=4
export I_MPI_ADJUST_ALLGATHERV=4

ibrun python3 ./scripts/main.py  \
    --config-file=Full_Scale_GC_Exc_Sat_LNN.yaml  \
    --arena-id=A --trajectory-id=Diag \
    --template-paths=../dgc/Mateos-Aparicio2014:templates \
    --dataset-prefix="$SCRATCH/dentate" \
    --results-path=$results_path \
    --io-size=256 \
    --tstop=150 \
    --v-init=-75 \
    --results-write-time=600 \
    --stimulus-onset=50.0 \
    --max-walltime-hours=1.0 \
    --vrecord-fraction=0.001 \
    --verbose

