export DATA_PREFIX=/media/igr/d865f900-7fcd-45c7-a7a7-bd2a7391bc40/Data/DG/

mpirun -n 1 python3  network_clamp.py go -c Network_Clamp_GC_Exc_Sat_SLN_IN_Izh_extent.yaml \
       --template-paths templates --dt 0.001 \
       -p MC -g 1015624 -t 9500 \
       --dataset-prefix $DATA_PREFIX \
       --config-prefix config \
       --input-features-path "$DATA_PREFIX/Full_Scale_Control/DG_input_features_20200910_compressed.h5" \
       --input-features-namespaces 'Place Selectivity' \
       --input-features-namespaces 'Grid Selectivity' \
       --input-features-namespaces 'Constant Selectivity' \
       --arena-id A --trajectory-id Diag \
       --results-path results/netclamp \

       --params-path config/20201105_Izhi_compiled.yaml



