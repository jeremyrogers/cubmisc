#! /bin/sh


cubmethod="cubappr" # change method
sdlog="0.5 1 2 4" # change sdlog in command line, have to figure out why it is not working with variable
genome="../data/ecoli_K12_MG1655_genome_filtered.fasta"
empphi="../data/ecoli_X_obs.csv"
prefix="without_xobs_" # change name
suffix="multichain_lapply_4cores"
folder="../results/conv_test/" # change folder
pinit="woxobs_pinit.csv" # make sure correct file is loaded (with X obs, w/o X obs)

if [ ! -d "$folder" ]; then
  mkdir $folder
fi

foutname="$prefix$suffix"
logfile="$folder$prefix$suffix.log"

if [ "$cubmethod" = "cubfits" ]; then
	echo "running with x obs run\n"
	nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -p $empphi -o $folder -n $foutname -i $pinit >> $logfile &
fi
if [ "$cubmethod" = "cubappr" ]; then
	if [ "$sdlog" = "0" ]; then 
		echo "sdog is 0, execution stopped! \n" 
		exit 1 
	fi

	echo "running without x obs run"
	nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -o $folder -n $foutname -i $pinit >> $logfile &
fi

