#!/usr/bin/env bash

FULLNAME=( $(ls | grep flagstat) ) # more than one file 
bamname=( $(ls | grep merged.bam) ) # only one file
bamid=("$( cut -d '.' -f1 <<< $bamname)") #common part

### first extract stat number 

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

### check if current bam should be down sampled and if so do it

j=0
for i in "${arrayNumber[@]}"; do
    id="$( cut -d '.' -f1 <<< "${FULLNAME[j]}")" # this is id from flagstat
    if [ "$id" = "$bamid" ]; then      
	sub=$(echo "scale=2; $min/$i" | bc -l) # check if sequencing depth is largere than min
    	if [ $(echo $sub'<'1 | bc -l) -eq 1 ]; then
        	java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar DownsampleSam  I=${id}.merged.bam O=${id}.subset.bam P=$sub
	else
		cp  ${id}.merged.bam ${id}.subset.bam
	fi
    fi
	
	(( j++ ))
done
 

