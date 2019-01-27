#!/bin/bash

NAME=$1
for WINDOW in 10000;
do
	for file in $(ls coverage/mosdepth/$NAME/*.${WINDOW}bp.regions.bed.gz)
 	do
     		b=$(basename $file .${WINDOW}bp.regions.bed.gz)
     		GROUP=$(echo $b | perl -p -e 'my $in = $_; chomp($in); open(my $fh => "cat *samples.csv |") || die $!; while(my $n = <$fh>) { chomp($n); my @row = split(",",$n); $lookup{$row[0]} = [@row];} 
my $new = $lookup{$in}->[1]; 
$_ = $new."\n"')
 # echo "$GROUP $b"
 		mean=$(zcat $file | awk '{total += $4} END { print total/NR}') 
 		zcat $file | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4/'${mean}','\"$GROUP\"','\"$b\"'}' 
	done  > coverage/$NAME.mosdepth.${WINDOW}bp.gg.tab
 	pigz -f coverage/$NAME.mosdepth.${WINDOW}bp.gg.tab
done


