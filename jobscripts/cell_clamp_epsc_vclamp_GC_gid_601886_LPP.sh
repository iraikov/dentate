#!/bin/bash

export DATA_PREFIX=/media/igr/d865f900-7fcd-45c7-a7a7-bd2a7391bc40/Data/DG
export DATA_PREFIX=./datasets

#        --section-index 11 \
#        --section-index 10 \
#        --section-index 12 \
#        --section-index 13 \
#        --syn-weight 1 \
    
mpirun -n 1 python3 ./cell_clamp.py -v \
        -g 601886 \
        --population GC \
        --config-prefix $HOME/src/model/dentate/config \
        --config Network_Clamp_GC_Exc_Sat_SynExp3NMDA2SGfdb_SLN_CLS_IN_PR_gid_601886.yaml \
        --dataset-prefix $DATA_PREFIX \
        --template-paths=templates \
        --results-path=results/cell_clamp \
        --load-weights \
        --presyn-name LPP \
        --syn-mech-name AMPA \
        --syn-mech-name NMDA \
        --syn-count 1 \
        --erev 0 \
        --v-init -75 \
        --stim-count 3 \
        --stim-interval 5 \
        --vclamp-hold -90 \
        --vclamp-hold -80 \
        --vclamp-hold -70 \
        --vclamp-hold -60 \
        --vclamp-hold -50 \
        --vclamp-hold -40 \
        --vclamp-hold -30 \
        --vclamp-hold -20 \
        --vclamp-hold -10 \
        --vclamp-hold 0 \
        -d ina -d ik -d i_Kir21 -d i_KA_Aradi -d i_sKDR_Aradi -d i_fKDR_Aradi \
        -m psc_vclamp
