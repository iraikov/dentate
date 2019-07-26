#!/bin/bash

dataset_prefix=/media/igr/d865f900-7fcd-45c7-a7a7-bd2a7391bc40/Data/DG/Full_Scale_Control

mpirun.mpich -np 8 python3 ./scripts/generate_DG_input_selectivity_features.py \
             --config=Full_Scale_Basis.yaml \
             --config-prefix=./config \
             --coords-path=${dataset_prefix}/DG_coords_20190717_compressed.h5 \
             --output-path=${dataset_prefix}/DG_input_features_20190724.h5 \
             -p LPP -p MPP -p CA3c -v --plot --save-fig DG_input_selectivity_features

