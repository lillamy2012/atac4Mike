# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load Picard/1.141
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
ALIGN=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $2}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
bampath=$(grep bampath $settingsfile | awk '{print $2}')
index=$(grep index $settingsfile | awk '{print $2}')
alignpath=$(grep alignpath $settingsfile | awk '{print $2}')
quality=$(grep Q_threshold $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')
tmp=$(grep tmp_path $settingsfile | awk '{print $2}')


echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "ALIGN: $ALIGN"

# === START ===
## This script starts with a aligned and sorted bam file that has been downsampled if needed


## mark duplicates for alignment statistics
export TMPDIR=$WORKDIR/Mike/picard_tmp
# make sure the directory exists
mkdir -p $TMPDIR

### make tmp folder
mkdir -p ${workpath}/${tmp}
mkdir -p ${workpath}/${outname}/results

###### MARKED DUPLICATES
input=${alignpath}/${ALIGN}.bam
output=${workpath}/${tmp}/${NAME}.marked_duplicates.bam
output_keep=${workpath}/${outname}/results/${NAME}.marked_duplicates.bam

if [ ! -f ${output} ];then
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${input}  O= ${output} M=${NAME}.dup_metrics.txt AS=true REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR
## alignment statistics with marked duplicates
    samtools flagstat ${output} > ${output_keep}.flagstat.txt
fi

samtools depth -a ${output} > ${output_keep}.depth.txt
samtools stats -c 0,2000,1 ${output} > ${output_keep}.stats.txt

####### QUALITY FILTER

## quality filter
input=${output}
output=${workpath}/${tmp}/${NAME}.marked_duplicates.filtered.bam
output_keep=${workpath}/${outname}/results/${NAME}.marked_duplicates.filtered.bam

if [ ! -f ${output} ];then
    echo "filter"
    echo $quality
    samtools view -bSq $quality ${input}> ${output}
    ## alignment statistics filtered
    samtools flagstat ${output} > ${output_keep}.flagstat.txt
fi

samtools stats -c 0,2000,1 ${output} > ${output_keep}.stats.txt
samtools depth -a ${output} > ${output_keep}.depth.txt

###### REMOVE DUPLICATES

input=${output}
output=${workpath}/${tmp}/${NAME}.uniq_filtered.bam
output_keep=${workpath}/${outname}/results/${NAME}.uniq_filtered.bam

if [ ! -f ${output} ];then
    echo "rm"
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${input}  O= ${output} M=${NAME}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
## alignment statistics dup rm
    samtools flagstat ${output} > ${output_keep}.flagstat.txt
fi
samtools depth -a ${output} > ${output_keep}.depth.txt
samtools stats -c 0,2000,1 ${output} > ${output_keep}.stats.txt



