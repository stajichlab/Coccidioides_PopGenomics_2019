#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --time 2:00:00 -p short --mem 64G --out logs/mosdepth.log
#SBATCH -J modepth
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
module unload python/2.7.5

for cfg in C*_config.txt
do
	source $cfg
	mkdir -p coverage/mosdepth/$GENOMENAME
	export PATH="/bigdata/stajichlab/jstajich/miniconda3/bin:$PATH"
	
	WINDOW=10000
	D=${ASMALNFOLDER}
	echo "$D $GENOMENAME"
	parallel --jobs $CPU mosdepth -f ${REFGENOME} -T 1,10,50,100,200 -n --by ${WINDOW}  -t 2 "{= s:${D}\/:coverage/mosdepth/${GENOMENAME}/:; s:\.cram:.${WINDOW}bp: =}" {} ::: ${D}/*.cram
	bash scripts/mosdepth_prep_ggplot.sh $GENOMENAME
	mkdir -p plots
	Rscript Rscripts/plot_mosdepth_CNV.R $GENOMENAME
done
