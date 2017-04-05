#!/bin/bash
#
#SBATCH -J place_cell_spikes
#SBATCH -o ./results/ca3_place_cell_spikes.%j.o
#SBATCH -t 12:00:00
#SBATCH --mem 8192
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

module load matlab

BATCH_SIZE=400
BATCH_INDEX=$SLURM_ARRAY_TASK_ID
OUTPUT_PATH=$SCRATCH/ca3_placecells
INPUT_DATA_PATH=$SCRATCH/ca3_placecells/linear_ca3_placecell_data.bin

export BATCH_SIZE
export BATCH_INDEX
export OUTPUT_PATH
export INPUT_DATA_PATH

cd $HOME/model/dentate/scripts/SpatialCells
./run_gen_ca3_placecell_spikes.sh /share/sw/licensed/MATLAB-R2016b

