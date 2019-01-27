#!/bin/bash
#SBATCH -p batch --time 1-0:0:0
module load aspera
module load sratoolkit

for acc in $(cat acc.txt);
do
 prefetch -a "ascp|$ASPERAKEY" $acc
 if [ ! -f ${acc}_1.fastq.gz ]; then
  fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files  --gzip $acc
 fi
done


