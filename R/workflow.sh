#! /bin/sh


cubmethod="cubfits" # change method
sdlog="0.5 1 2 4" # change sdlog in command line, have to figure out why it is not working with variable

#genome="../K12_genome/Genome/ecoli_K12_MG1655_genome_filtered.fasta"
#genome="../brewersYeast/s288c.genome.fasta"
#genome="../S.cervisiae.REU13/section1.fasta"
genome="../jeremyYeast/Genome/jeremyYeast.fasta"
#genome="../loganYeast/CubfitsValues1.fasta"
#genome="../prestonYeast/S.cerevisiae.S288c.fasta"

#empphi="../K12_genome/Genome/ecoli_X_obs.csv"
#empphi="../brewersYeast/yassour2009.phi.tsv"
#empphi="../S.cervisiae.REU13/section1.csv"
empphi="../jeremyYeast/jyeast.phi.tsv"
#empphi="../loganYeast/pyeast.phi.tsv"
#empphi="../prestonYeast/pyeast.phi.tsv"

prefix="jeremyyeast" # change name
#prefix="yeast_NSEfits_" # change name
#prefix="S.cervisiae.REU13" # change name
suffix="withPhi"

folder="../results/jeremyyeast/scaledphi3/" # change folder

pinit="wxobs_pinit.csv" # file with initial values for multiple parameters such as Ax, s_phi^{(0)}, s_\varepsilon
#pinit="woxobs_pinit.csv" # file with initial values for multiple parameters such as Ax, s_phi^{(0)}, s_\varepsilon
                ##NOTE: With and without Xobs fits use different files.

if [ ! -d "$folder" ]; then
  mkdir -p $folder
fi

foutname="$prefix$suffix"
logfile="$folder$prefix$suffix.log"

if [ "$cubmethod" = "cubfits" ]; then
	echo "running with x obs run\n"
      #nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -p $empphi -o $folder -n $foutname -i $pinit >> $logfile &
      nohup Rscript run_nsef.r -c $cubmethod -s "0.5 1 2 4" -f $genome -p $empphi -o $folder -n $foutname -i $pinit >> $logfile &
fi
if [ "$cubmethod" = "cubappr" ]; then
	if [ "$sdlog" = "0" ]; then #double check that sdlog value is properly defined
		echo "sdlog is 0, execution stopped! \n" 
		exit 1 
	fi

	echo "running without x obs run\n"
	#nohup Rscript run_roc.r -c $cubmethod -s "0.5 1 2 4" -f $genome -o $folder -n $foutname -i $pinit >> $logfile &
	nohup Rscript run_nsef.r -c $cubmethod -s "0.5 1 2 4" -f $genome -o $folder -n $foutname -i $pinit >> $logfile &

fi

