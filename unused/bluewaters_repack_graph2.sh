#!/bin/bash
### set the number of nodes and the number of PEs per node
#PBS -l nodes=1:ppn=8:xe
### which queue to use
#PBS -q high
### set the wallclock time
#PBS -l walltime=12:00:00
### set the job name
#PBS -N dentate_repack_graph
### set the job stdout and stderr
#PBS -e ./results/dentate_repack_graph.$PBS_JOBID.err
#PBS -o ./results/dentate_repack_graph.$PBS_JOBID.out
### set email notification
##PBS -m bea
### Set umask so users in my group can read job stdout and stderr files
#PBS -W umask=0027

module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel

set -x

cd $PBS_O_WORKDIR

export prefix=/projects/sciteam/baqc/Full_Scale_Control
export input=$prefix/DG_GC_connections_20180721.h5
export output=$prefix/DG_GC_connections_20180721_compressed.h5

aprun h5repack -v -f SHUF -f GZIP=9 -i $input -o $output
