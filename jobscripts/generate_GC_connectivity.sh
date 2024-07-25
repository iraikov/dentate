#!/bin/bash
#             --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_gid_601886.yaml \
#             --forest-path=./datasets/Single/data_601886_20240713.h5 \
#             --connectivity-path=./datasets/Single/data_601886_20240713.h5 \
#             --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_center_pf.yaml \
#             --forest-path=./datasets/Slice/dentatenet_Slice_SLN_center_pf_20240715.h5 \
#             --connectivity-path=./datasets/Slice/dentatenet_Slice_SLN_center_pf_20240715.h5 \

mpirun.mpich -np 1 python3 ./scripts/generate_distance_connections.py \
             --config-prefix=./config \
             --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_gid_720884.yaml \
             --forest-path=./datasets/Single/data_720884_20240718.h5 \
             --connectivity-path=./datasets/Single/data_720884_20240718.h5 \
             --connectivity-namespace=Connections \
             --coords-path=./datasets/Full_Scale_Control/DG_coords_20190717_compressed.h5 \
             --coords-namespace=Coordinates \
             --io-size=1 -v
