#! /bin/sh


cubmethod="cubfits" # change method
sdlog="0.5 1 2 4" # change sdlog in command line, have to figure out why it is not working with variable
#
#genome="../data/ecoli_K12_MG1655_genome_filtered.fasta"
#genome="../S.cervisiae.REU13/yeast.sim.fasta"
genome="../S.cervisiae.REU13/section1.fasta"
#
#empphi="../data/ecoli_X_obs.csv"
#empphi="../S.cervisiae.REU13/yeast.sim.Xobs.csv"
empphi="../S.cervisiae.REU13/section1.csv"
#
#prefix="ecoli_fits_" # change name
prefix="yeast_NSEfits2_" # change name
#prefix="without_xobs_" # change name
suffix="singlechain"
folder="../results/conv_test/" # change folder
pinit="wxobs_pinit.csv" # file with initial values for multiple parameters such as Ax, s_phi^{(0)}, s_\varepsilon
#pinit="woxobs_pinit.csv" # file with initial values for multiple parameters such as Ax, s_phi^{(0)}, s_\varepsilon
#pinit="woxobs_pinit.csv" # file with initial values for multiple parameters such as Ax, s_phi^{(0)}, s_\varepsilon
                ##NOTE: With and without Xobs fits use different files.

if [ ! -d "$folder" ]; then
  mkdir -p $folder
fi

foutname="$prefix$suffix"
logfile="$folder$prefix$suffix.log"

if [ "$cubmethod" = "cubfits" ]; then
	echo "running with x obs run\n"
       nohup Rscript run_nsef.r -c $cubmethod -s "0.5 1 2 4" -f $genome -p $empphi -o $folder -n $foutname -i $pinit >> $logfile &
#	nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -p $empphi -o $folder -n $foutname -i $pinit >> $logfile &
fi
if [ "$cubmethod" = "cubappr" ]; then
	if [ "$sdlog" = "0" ]; then #double check that sdlog value is properly defined
		echo "sdog is 0, execution stopped! \n" 
		exit 1 
	fi

	echo "running without x obs run"
	nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -o $folder -n $foutname -i $pinit >> $logfile &
#	nohup Rscript run_nsef.r -c $cubmethod -s "0.5 1 2 4" -f $genome -o $folder -n $foutname -i $pinit >> $logfile
fi

