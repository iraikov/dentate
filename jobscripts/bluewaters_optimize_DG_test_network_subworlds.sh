### set the number of nodes and the number of PEs per node
#PBS -l nodes=72:ppn=16:xe
### which queue to use
#PBS -q normal
### set the wallclock time
#PBS -l walltime=5:00:00
### set the job name
#PBS -N optimize_DG_test_network_subworlds
### set the job stdout and stderr
#PBS -e ./results/optimize_DG_test_network_subworlds.$PBS_JOBID.err
#PBS -o ./results/optimize_DG_test_network_subworlds.$PBS_JOBID.out
### set email notification
##PBS -m bea
### Set umask so users in my group can read job stdout and stderr files
#PBS -W umask=0027
#PBS -A bayj


module load bwpy/2.0.1
module load craype-hugepages2M

set -x

export SCRATCH=/projects/sciteam/bayj
export NEURONROOT=$SCRATCH/nrnintel
export PYTHONPATH=$HOME/model:$NEURONROOT/lib/python:$SCRATCH/site-packages:$PYTHONPATH
export PATH=$NEURONROOT/x86_64/bin:$PATH
export MODEL_HOME=$HOME/model
export DG_HOME=$MODEL_HOME/dentate

results_path=$DG_HOME/results/optimize_DG_test_network_subworlds.$PBS_JOBID
export results_path

cd $PBS_O_WORKDIR

mkdir -p $results_path

cd tests

aprun -n 1152 -b -- bwpy-environ -- \
    python2.7 -m nested.optimize  \
     --config-file-path=$DG_HOME/config/DG_test_network_subworlds_config.yaml \
     --output-dir=$results_path \
     --pop_size=4 \
     --max_iter=4 \
     --path_length=1 \
     --framework=pc \
     --disp \
     --procs-per-worker=288 \
     --no-cleanup \
     --verbose \
     --template_paths=$MODEL_HOME/dgc/Mateos-Aparicio2014:$DG_HOME/templates \
     --dataset_prefix="$SCRATCH" \
     --config_prefix=$DG_HOME/config \
     --results_path=$results_path \
     --cell_selection_path=$DG_HOME/datasets/DG_slice.yaml \
     --spike_input_path=$DG_HOME/results/Full_Scale_GC_Exc_Sat_LNN_9870802.bw/dentatenet_Full_Scale_GC_Exc_Sat_LNN_results.h5 \
     --spike_input_namespace='Spike Events' \
     --max-walltime-hours=3.75 \
     --io-size=64 \
     -v
