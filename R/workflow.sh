#! /bin/sh


cubmethod="cubappr" # change method
sdlog="0.5 1 2 4" # change sdlog 
genome="../data/ecoli_K12_MG1655_genome_filtered.fasta"
empphi="../data/ecoli_X_obs.csv"
prefix="without_xobs_" # change name
suffix="multichain_lapply_4cores"
folder="../results/benchmark/" # change folder

if [ ! -d "$folder" ]; then
  mkdir $folder
fi

fnout="$folder$prefix$suffix"
logfile="$folder$prefix$suffix.log"

if [ "$cubmethod" = "cubfits" ]; then
	echo "running with x obs run\n"
	nohup Rscript run_roc.r -c $cubmethod -s $sdlog -f $genome -p $empphi -o $fnout > $logfile &
fi
if [ "$cubmethod" = "cubappr" ]; then
	if [ $sdlog = 0 ]; then 
		echo "sdog is 0, execution stopped! \n" 
		exit 1 
	fi

	echo "running without x obs run\n"
	nohup Rscript run_roc.r -c $cubmethod -s $sdlog -f $genome -o $fnout > $logfile &
fi

