#!/bin/bash
#
#SBATCH -J netclamp_single_cell 
#SBATCH -o ./results/netclamp_single_cell.%j.o
#SBATCH --nodes=25
#SBATCH --ntasks-per-node=56
#SBATCH -p normal
#SBATCH -t 4:00:00
#SBATCH --mail-user=pmoolcha@stanford.edu
#SBATCH --mail-type=END
#


export LD_PRELOAD=$MKLROOT/lib/intel64_lin/libmkl_core.so:$MKLROOT/lib/intel64_lin/libmkl_sequential.so
export FI_MLX_ENABLE_SPAWN=yes
ml load intel19

fil=(
" HC 1030000 15879716 "
" HC 1030000 45419272 "
" HC 1030000 28682721 "
" HC 1030000 53736785 "
" HC 1030000 63599789 "
" HC 1032250 39786206 "
" HC 1032250 29585629 "
" HC 1032250 80350627 "
" HC 1032250 13571683 "
" HC 1032250 204284 "
" HC 1034500 17324444 "
" HC 1034500 19663111 "
" HC 1034500 48663281 "
" HC 1034500 85653587 "
" HC 1034500 61645166 "
" HC 1036750 96962730 "
" HC 1036750 96494881 "
" HC 1036750 9389064 "
" HC 1036750 72215563 "
" HC 1036750 86292756 "
" HC 1038999 69758687 "
" HC 1038999 19285510 "
" HC 1038999 57153876 "
" HC 1038999 45197076 "
" HC 1038999 58867820 "
" BC 1039000 93454042 "
" BC 1039000 74865768 "
" BC 1039000 1503844 "
" BC 1039000 52357252 "
" BC 1039000 28135771 "
" BC 1039950 27431042 "
" BC 1039950 90485672 "
" BC 1039950 49779904 "
" BC 1039950 97486157 "
" BC 1039950 67539344 "
" BC 1040900 32603857 "
" BC 1040900 98230850 "
" BC 1040900 97616150 "
" BC 1040900 63255735 "
" BC 1040900 4643442 "
" BC 1041850 87742037 "
" BC 1041850 63875616 "
" BC 1041850 21074287 "
" BC 1041850 67275914 "
" BC 1041850 19191189 "
" BC 1042799 183453 "
" BC 1042799 36739208 "
" BC 1042799 62574877 "
" BC 1042799 57586813 "
" BC 1042799 67308587 "
" AAC 1042800 4893658 "
" AAC 1042800 27137089 "
" AAC 1042800 36010476 "
" AAC 1042800 53499406 "
" AAC 1042800 49937004 "
" AAC 1042913 291547 "
" AAC 1042913 8379552 "
" AAC 1042913 80515073 "
" AAC 1042913 38840365 "
" AAC 1042913 9111622 "
" AAC 1043025 62387369 "
" AAC 1043025 52829774 "
" AAC 1043025 59206615 "
" AAC 1043025 82956063 "
" AAC 1043025 7818268 "
" AAC 1043138 19281943 "
" AAC 1043138 40133402 "
" AAC 1043138 82470709 "
" AAC 1043138 85264434 "
" AAC 1043138 70337332 "
" AAC 1043249 54652217 "
" AAC 1043249 43400137 "
" AAC 1043249 95905199 "
" AAC 1043249 66598438 "
" AAC 1043249 26662642 "
" HCC 1043250 33236209 "
" HCC 1043250 92055940 "
" HCC 1043250 71407528 "
" HCC 1043250 17609813 "
" HCC 1043250 12260638 "
" HCC 1043600 42402504 "
" HCC 1043600 89433777 "
" HCC 1043600 60991105 "
" HCC 1043600 64290895 "
" HCC 1043600 17293770 "
" HCC 1043950 99434948 "
" HCC 1043950 57660249 "
" HCC 1043950 54593731 "
" HCC 1043950 72125941 "
" HCC 1043950 41566230 "
" HCC 1044300 97569363 "
" HCC 1044300 66834161 "
" HCC 1044300 94613541 "
" HCC 1044300 63592626 "
" HCC 1044300 92910319 "
" HCC 1044649 84121621 "
" HCC 1044649 94560988 "
" HCC 1044649 46366417 "
" HCC 1044649 24805208 "
" HCC 1044649 59396015 "
)
NGFC=(
" NGFC 1044650 12740157 "
" NGFC 1044650 97895890 "
" NGFC 1044650 93872787 "
" NGFC 1044650 95844113 "
" NGFC 1044650 96772370 "
" NGFC 1045900 67428613 "
" NGFC 1045900 95436908 "
" NGFC 1045900 6112188 "
" NGFC 1045900 71039025 "
" NGFC 1045900 89814943 "
" NGFC 1047150 59071557 "
" NGFC 1047150 77901687 "
" NGFC 1047150 27400566 "
" NGFC 1047150 50965365 "
" NGFC 1047150 48744644 "
" NGFC 1048400 80347840 "
" NGFC 1048400 38650070 "
" NGFC 1048400 62046131 "
" NGFC 1048400 35472841 "
" NGFC 1048400 35297351 "
" NGFC 1049649 77179804 "
" NGFC 1049649 26628153 "
" NGFC 1049649 99082330 "
" NGFC 1049649 89481705 "
" NGFC 1049649 10249569 "
)
lif=(
" IS 1049650 4259860 "
" IS 1049650 11745958 "
" IS 1049650 75940072 "
" IS 1049650 49627038 "
" IS 1049650 84013649 "
" IS 1050400 63796673 "
" IS 1050400 69320701 "
" IS 1050400 7843435 "
" IS 1050400 10084233 "
" IS 1050400 93591428 "
" IS 1051150 22725943 "
" IS 1051150 21032749 "
" IS 1051150 1339500 "
" IS 1051150 83916441 "
" IS 1051150 49587749 "
" IS 1051900 82185961 "
" IS 1051900 27654574 "
" IS 1051900 23672271 "
" IS 1051900 70119958 "
" IS 1051900 51871840 "
" IS 1052649 45707385 "
" IS 1052649 37549278 "
" IS 1052649 18680556 "
" IS 1052649 60814941 "
" IS 1052649 82004212 "
" MOPP 1052650 31571230 "
" MOPP 1052650 45373570 "
" MOPP 1052650 85763600 "
" MOPP 1052650 68839073 "
" MOPP 1052650 29079471 "
" MOPP 1053650 35281577 "
" MOPP 1053650 82093235 "
" MOPP 1053650 78038978 "
" MOPP 1053650 39888091 "
" MOPP 1053650 59550066 "
" MOPP 1054650 60567645 "
" MOPP 1054650 94967765 "
" MOPP 1054650 21247851 "
" MOPP 1054650 26628185 "
" MOPP 1054650 3611780 "
" MOPP 1055650 34097792 "
" MOPP 1055650 44866707 "
" MOPP 1055650 61810606 "
" MOPP 1055650 83145544 "
" MOPP 1055650 79924848 "
" MOPP 1056649 17666981 "
" MOPP 1056649 88486608 "
" MOPP 1056649 92808036 "
" MOPP 1056649 73504121 "
" MOPP 1056649 68347478 "
)

pil=(
" HC 1030000 15879716 "
" HC 1030000 45419272 "
" HC 1030000 53736785 "
" HC 1030000 63599789 "
" HC 1032250 39786206 "
" BC 1042799 67308587 "
" AAC 1043249 26662642 "
" HCC 1043600 42402504 "
" MOPP 1053650 78038978 "
)

N_cores=35

IFS='
'
counter=0
for f in ${NGFC[@]}
do

set -- "$f" 
IFS=" " ; declare -a tempvar=($*) 


#ibrun -n 56 -o  0 task_affinity ./mycode.exe input1 &   # 56 tasks; offset by  0 entries in hostfile.
#ibrun -n 56 -o 56 task_affinity ./mycode.exe input2 &   # 56 tasks; offset by 56 entries in hostfile.
#wait                                                    # Required; else script will exit immediately.


#ibrun -n $N_cores -o $counter python3 network_clamp.py go -c 20201022_Network_Clamp_GC_Exc_Sat_SLN_IN_Izh.yaml \
#    --template-paths templates \
#    -p ${tempvar[0]:1:-1} -g ${tempvar[1]} -t 9500 --dt 0.001 \
#    --dataset-prefix /scratch1/03320/iraikov/striped/dentate \
#    --config-prefix config \
#    --input-features-path /scratch1/03320/iraikov/striped/dentate/Full_Scale_Control/DG_input_features_20200910_compressed.h5 \
#    --input-features-namespaces 'Place Selectivity' \
#    --input-features-namespaces 'Grid Selectivity' \
#    --input-features-namespaces 'Constant Selectivity' \
#    --arena-id A --trajectory-id Diag \
#    --results-path results/netclamp \
#    --opt-seed ${tempvar[3]} \
#    --params-path ${tempvar[2]:1:-1} &


#pop=${tempvar[0]:1:-1}
pop=${tempvar[0]}
gid=${tempvar[1]}
seed=${tempvar[2]}

#ibrun -n 8 python3  network_clamp.py optimize -c 20201022_Network_Clamp_GC_Exc_Sat_SLN_IN_Izh.yaml \
ibrun -n $N_cores -o $((counter * 56))  python3  network_clamp.py optimize -c 20201022_Network_Clamp_GC_Exc_Sat_SLN_IN_Izh.yaml \
    --template-paths templates \
    -p $pop -g $gid -t 9500 --dt 0.001 \
    --dataset-prefix /scratch1/03320/iraikov/striped/dentate \
    --config-prefix config \
    --input-features-path /scratch1/03320/iraikov/striped/dentate/Full_Scale_Control/DG_input_features_20200910_compressed.h5 \
    --input-features-namespaces 'Place Selectivity' \
    --input-features-namespaces 'Grid Selectivity' \
    --input-features-namespaces 'Constant Selectivity' \
    --arena-id A --trajectory-id Diag \
    --results-path results/netclamp \
    --param-config-name "Weight exc inh microcircuit" \
    --opt-seed $seed \
    --opt-iter 4000 rate & 


#    --results-file network_clamp.optimize.$pop\_$gid\_$(date +'%Y%m%d_%H%M%S')\_$seed.h5 \

counter=$((counter + 1))

done
wait
