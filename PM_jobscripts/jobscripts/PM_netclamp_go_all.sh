now=$(date +"%Y%m%d_%H%M%S")
Izhiseries='_Izhi_20201022_'
#Optseries='20201028_233742'
expfx=$now$Optseries$Izhiseries'netclamp_go_'
ext='.log'
logdir=results/netclamp/logs/
logdirexpfx=$logdir$expfx

cp jobscripts/PM_netclamp_go_all.sh $logdirexpfx'all'.sh 


#ACyaml='network_clamp.optimize.AAC_1042800_20201028_233742.yaml'
#BCyaml='network_clamp.optimize.BC_1039000_20201028_233742.yaml'
#HCyaml='network_clamp.optimize.HCC_1043250_20201028_233742.yaml'
#HPyaml='network_clamp.optimize.HC_1030000_20201028_233742.yaml'
#ISyaml='network_clamp.optimize.IS_1049650_20201028_233742.yaml'
#MOyaml='network_clamp.optimize.MOPP_1052650_20201028_233742.yaml'
#NGyaml='network_clamp.optimize.NGFC_1044650_20201028_233742.yaml'


HCyaml='network_clamp.optimize.HCC_1043250_20201029_213934.yaml'  
ISyaml='network_clamp.optimize.IS_1049650_20201029_213934.yaml'  
NGyaml='network_clamp.optimize.NGFC_1044650_20201029_213934.yaml'

#ibrun -n 1 ./jobscripts/PM_netclamp_go_AAC.sh    $ACyaml 2>&1 | tee $logdirexpfx'AAC'$ext  &
#ibrun -n 1 ./jobscripts/PM_netclamp_go_BC.sh     $BCyaml 2>&1 | tee $logdirexpfx'BC'$ext  &
ibrun -n 1 ./jobscripts/PM_netclamp_go_HICAP.sh  $HCyaml 2>&1 | tee $logdirexpfx'HICAP'$ext  &
#ibrun -n 1 ./jobscripts/PM_netclamp_go_HIPP.sh   $HPyaml 2>&1 | tee $logdirexpfx'HIPP'$ext  &
ibrun -n 1 ./jobscripts/PM_netclamp_go_IS.sh     $ISyaml 2>&1 | tee $logdirexpfx'IS'$ext  &
#ibrun -n 1 ./jobscripts/PM_netclamp_go_MOPP.sh   $MOyaml 2>&1 | tee $logdirexpfx'MOPP'$ext  &
ibrun -n 1 ./jobscripts/PM_netclamp_go_NGFC.sh   $NGyaml 2>&1 | tee $logdirexpfx'NGFC'$ext  

#ibrun -n 8 ./jobscripts/netclamp_opt_MC.sh*
