#!/usr/bin/env bash

FULLNAME=( $(ls | grep flagstat) )

i=0 
for n in "${FULLNAME[@]}"; do
 	y=$(awk 'NR==1{print $1}' $n)
	arrayNumber[$i]=$y
	if [ $i -eq 0 ];then 
 		min=$y
 	fi
	(( i++ ))
 	((y < min)) && min=$y
done
j=0
for i in "${arrayNumber[@]}"; do
    id="$( cut -d '.' -f1 <<< "${FULLNAME[j]}")"   
    sub=$(echo "scale=2; $min/$i" | bc -l)
    	if [ $(echo $sub'<'1 | bc -l) -eq 1 ]; then
        	java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar DownsampleSam  I=${id}.merged.bam O=${id}.subset.bam P=$sub
	else
		cp  ${id}.merged.bam ${id}.subset.bam
	fi
	
	(( j++ ))
done
 

