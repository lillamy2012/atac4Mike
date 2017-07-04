
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt
workpath=$(grep workpath $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')

names=( $(cut -d$'\t' -f1 $PARAM_FILE ) )
inputs=( $(cut -d$'\t' -f6 $PARAM_FILE ) )

##########################
### tables for individule files:
## statistics, counts 
###########################

## statistics

filetypes=(  uniq_filtered.bam )
table=( "IS" "numbers" "depth.table" )

for j in "${filetypes[@]}"
do
    for k in "${table[@]}"
    do
        echo R_$j.$k.tab
        if [ -f $workpath/$outname/R_$j.$k.tab ];then
            echo "rm"
            rm $workpath/$outname/R_$j.$k.tab
        fi
        for i in "${names[@]}"
        do
            echo $i
            printf "%s\n" $i.$j.$k.tab >> $workpath/$outname/R_$j.$k.tab
        done
	#rm $workpath/$outname/R_$j.$k.tab
    done
done

###

echo "next"

for i in "${inputs[@]}";do
    echo $i
done

## peaks
p_filetypes=( uniq_filtered )
p_table=( "MACS2_paired"  "MACS2_noInput" )

isNotInput () {
    local e
    for e in "${@:2}"; do
        if [[ "$e" == "$1" ]];then
            return 1
        fi
    done
    return 0
}


for j in "${p_filetypes[@]}"
do
    for k in "${p_table[@]}"
    do
        echo R_$j.$k.tab
        if [ -f $workpath/$outname/R_$j.$k.tab ];then
            echo "rm"
            rm $workpath/$outname/R_$j.$k.tab
        fi
        for i in "${names[@]}"
        do
            echo $i
            isNotInput "$i" "${inputs[@]}"
            if [ $? -eq 0 ]; then
                printf "%s\n" ${i}.${j}_${k}_peaks.narrowPeak >> $workpath/$outname/R_$j.$k.tab
            fi
        done
	rm $workpath/$outname/R_$j.$k.tab
    done
done

## deseq2 files
## merged macs output

de_table=( "MACS2_paired" "MACS2_noInput")

for k in "${de_table[@]}"
   do
      if [ -f $workpath/$outname/R_${k}_deseq.tab ];then
            echo "rm"
            rm $workpath/$outname/R_${k}_deseq.tab
        fi
        for i in "${names[@]}"
        do
            echo $i
	type=$(grep ^$i $PARAM_FILE | awk '{print $4}')
            printf "%s\t" ${i}_${k}.gff.counts.tab $type $(grep ^$i $PARAM_FILE | awk '{print $3}') >> $workpath/$outname/R_${k}_deseq.tab
            printf "\n" >>  $workpath/$outname/R_${k}_deseq.tab
        done
done





##########################
### tables for merged files:
## run macs, R peak, 
###########################
## merged bam files
FILEA=$workpath/files.tab
types=($(awk '{print $4}' $FILEA | uniq))
if [ -f $workpath/files_macs.tab ];then
	echo "rm"
	rm $workpath/files_macs.tab
fi
for t in "${types[@]}"
do
    printf "%s\n" $workpath/$outname/$t/${t}.uniq_filtered.bam.subset.bam >> $workpath/files_macs.tab
done




p_table=( "MACS2_paired"  "MACS2_noInput" )

for k in "${p_table[@]}"
do
	echo R_$k.tab
        #	if [ -f $workpath/$outname/R_$j.$k.tab ];then
         #   echo "rm"
         #   rm $workpath/$outname/R_$j.$k.tab
        #fi
        for t in "${types[@]}"
        do
            echo $t
                printf "%s\n" $workpath/$outname/macs_results/${t}.uniq_filtered.bam.subset.bam_${k}_peaks.narrowPeak >> $workpath/$outname/R_$k.tab
        done
done

