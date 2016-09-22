#!/bin/bash

### set the number of nodes and the number of PEs per node
#PBS -l nodes=1024:ppn=2:xe
### which queue to use
#PBS -q high
### set the wallclock time
#PBS -l walltime=6:00:00
### set the job name
#PBS -N dentate_Full_Scale_ReducedGCtoHC
### set the job stdout and stderr
#PBS -e ./results/$PBS_JOBID.err
#PBS -o ./results/$PBS_JOBID.out
### set email notification
##PBS -m bea
### Set umask so users in my group can read job stdout and stderr files
#PBS -W umask=0027
### Get darsan profile data
#PBS -lgres=darshan

module swap PrgEnv-cray PrgEnv-gnu

set -x

cd $PBS_O_WORKDIR

results_path=./results/Full_Scale_ReducedGCtoHC_$PBS_JOBID
export results_path

mkdir -p $results_path

git ls-files | tar -zcf ${results_path}/dentate.tgz --files-from=/dev/stdin
git --git-dir=../dgc/.git ls-files | grep Mateos-Aparicio2014 | tar -C ../dgc -zcf ${results_path}/dgc.tgz --files-from=/dev/stdin

runhoc="./jobscripts/bluewaters_Full_Scale_ReducedGCtoHC_run_${PBS_JOBID}.hoc"

sed -e "s/JOB_ID/$PBS_JOBID/g" ./jobscripts/bluewaters_Full_Scale_ReducedGCtoHC_run.hoc > $runhoc

##export LD_PRELOAD=/sw/xe/darshan/2.3.0/darshan-2.3.0_cle52/lib/libdarshan.so 
export DARSHAN_LOGPATH=$PWD/darshan-logs

aprun -n 2048 ./mechanisms/x86_64/special -mpi -nobanner -nogui $runhoc


