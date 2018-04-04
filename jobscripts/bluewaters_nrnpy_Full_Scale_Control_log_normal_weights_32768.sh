#!/bin/bash

### set the number of nodes and the number of PEs per node
#PBS -l nodes=2048:ppn=16:xe
### which queue to use
#PBS -q high
### set the wallclock time
#PBS -l walltime=18:00:00
### set the job name
#PBS -N dentate_Full_Scale_Control_log_normal_weights_32768
### set the job stdout and stderr
#PBS -e ./results/dentate.$PBS_JOBID.err
#PBS -o ./results/dentate.$PBS_JOBID.out
### set email notification
##PBS -m bea
### Set umask so users in my group can read job stdout and stderr files
#PBS -W umask=0027


module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel
module load bwpy 
module load bwpy-mpi
module load atp

set -x

export ATP_ENABLED=1 
export LD_LIBRARY_PATH=/sw/bw/bwpy/0.3.0/python-mpi/usr/lib:/sw/bw/bwpy/0.3.0/python-single/usr/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$HOME/bin/nrn/lib/python:/projects/sciteam/baef/site-packages:$PYTHONPATH
export PATH=$HOME/bin/nrn/x86_64/bin:$PATH
export DARSHAN_LOGPATH=$PBS_O_WORKDIR/darshan-logs

echo python is `which python`
results_path=./results/Full_Scale_Control_$PBS_JOBID
export results_path

cd $PBS_O_WORKDIR

mkdir -p $results_path

git ls-files | tar -zcf ${results_path}/dentate.tgz --files-from=/dev/stdin
git --git-dir=../dgc/.git ls-files | grep Mateos-Aparicio2014 | tar -C ../dgc -zcf ${results_path}/dgc.tgz --files-from=/dev/stdin

## Necessary for correct loading of Darshan with LD_PRELOAD mechanism
export PMI_NO_FORK=1
export PMI_NO_PREINITIALIZE=1
export LD_PRELOAD=/opt/cray/hdf5-parallel/1.8.16/GNU/4.9/lib/libhdf5_parallel_gnu_49.so.10

aprun -n 32768 \
    python main.py  \
    --config-file=config/Full_Scale_Control_log_normal_weights.yaml  \
    --template-paths=../dgc/Mateos-Aparicio2014 \
    --dataset-prefix="/projects/sciteam/baef" \
    --results-path=$results_path \
    --io-size=256 \
    --tstop=12050 \
    --v-init=-75 \
    --max-walltime-hours=5.8 \
    --node-rank-file=parts.32768 \
    --vrecord-fraction=0.001 \
    --verbose

