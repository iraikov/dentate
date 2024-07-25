#        --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_gid_601886.yaml \
#        --forest-path=./datasets/Single/tree_601886.h5 \
#        --output-path=./datasets/Single/data_601886_20240713.h5 \
#        --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_center_pf.yaml \
#        --forest-path=./datasets/Slice/dentatenet_Slice_SLN_center_pf_20240715.h5 \
#        --output-path=./datasets/Slice/dentatenet_Slice_SLN_center_pf_20240715.h5 \

mpirun.mpich -n 1 python3 ./scripts/distribute_synapse_locs.py \
        --config-prefix=./config \
        --config=Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfd_SLN_CLS_IN_PR_gid_720884.yaml \
        --forest-path=./datasets/Single/tree_720884.h5 \
        --output-path=./datasets/Single/data_720884_20240718.h5 \
        --template-path=templates \
        --populations=GC \
        --distribution=poisson \
        --io-size=1 -v


