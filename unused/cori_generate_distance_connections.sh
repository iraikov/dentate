#!/bin/bash
#
#SBATCH -J generate_distance_connections
#SBATCH -o ./results/generate_distance_connections.%j.o
#SBATCH -N 64
#PBS -N generate_distance_connections
#SBATCH -p regular
#SBATCH -t 6:00:00
#SBATCH -L SCRATCH   # Job requires $SCRATCH file system
#SBATCH -C haswell   # Use Haswell nodes
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#

module swap PrgEnv-intel PrgEnv-gnu
module unload darshan
module load cray-hdf5-parallel/1.8.16
module load python/2.7-anaconda

set -x

export PYTHONPATH=$HOME/model/dentate:/opt/python/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=$HOME/bin/nrn/lib/python:$PYTHONPATH
export PYTHONPATH=$HOME/.local/cori/2.7-anaconda/lib/python2.7/site-packages:$PYTHONPATH

srun -n 2048 -c 1 python ./scripts/generate_distance_connections.py \
       --config=./config/Full_Scale_Control.yaml \
       --forest-path=$SCRATCH/dentate/Full_Scale_Control/MC_forest_syns_20171013.h5 \
       --connectivity-path=$SCRATCH/dentate/Full_Scale_Control/DG_MC_connections_20171026.h5 \
       --connectivity-namespace=Connections \
       --coords-path=$SCRATCH/dentate/Full_Scale_Control/dentate_Full_Scale_Control_coords_20171005.h5 \
       --coords-namespace=Coordinates \
       --io-size=320 --cache-size=1 --value-chunk-size=50000 --chunk-size=10000 --quick
