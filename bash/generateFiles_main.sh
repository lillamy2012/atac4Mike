
# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load cutadapt/1.9-ictce-5.3.0-Python-2.7.6
module load Bowtie2/2.1.0-goolf-1.4.10
module load Picard/1.141
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
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

#subn=$(grep subn $settingsfile | awk '{print $2}')

echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "ALIGN: $ALIGN"

# === START ===
## This script starts with a aligned and sorted bam file


## mark duplicates for alignment statistics
export TMPDIR=$WORKDIR/Mike/picard_tmp
# make sure the directory exists
mkdir -p $TMPDIR

input=${alignpath}/${ALIGN}.bam
output=${workpath}/${outname}/${NAME}.marked_duplicates.bam

if [ ! -f ${output} ];then
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${input}  O= ${output} M=${NAME}.dup_metrics.txt AS=true REMOVE_DUPLICATES=false TMP_DIR=$TMPDIR
    ## alignment statistics with marked duplicates
    samtools flagstat ${output} > ${output}.flagstat.txt
    samtools stats ${output} > ${output}.stats.txt
fi
samtools depth -a ${output} > ${output}.depth.txt
samtools stats -c 0,2000,1 ${output} > ${output}.stats.txt

## quality filter 
input=${output}
output=${workpath}/${outname}/${NAME}.marked_duplicates.filtered.bam
if [ ! -f ${output} ];then
    echo "filter"
echo $quality
    samtools view -bSq $quality ${input}> ${output}
## alignment statistics filtered
    samtools flagstat ${output} > ${output}.flagstat.txt
    samtools stats ${output}> ${output}.stats.txt
fi
samtools stats -c 0,2000,1 ${output} > ${output}.stats.txt
samtools depth -a ${output} > ${output}.depth.txt

## remove duplicates

input=${output}
output=${workpath}/${outname}/${NAME}.uniq_filtered.bam

if [ ! -f ${output} ];then
    echo "rm"
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${input}  O= ${output} M=${NAME}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
## alignment statistics dup rm
    samtools flagstat ${output} > ${output}.flagstat.txt
    samtools stats -c 0,2000,1 ${output} > ${output}.stats.txt

fi
samtools depth -a ${output} > ${output}.depth.txt

input=${output}
output1=${workpath}/${outname}/${NAME}.filter.duprm.0mer.bam
output2=${workpath}/${outname}/${NAME}.filter.duprm.mer.bam

## split
if [ ! -f $output1 ] || [ ! -f $output2 ];then
    echo "split"
#split bam, use reads with insert size <100
    samtools view -f 2 -h ${input}  | perl -ane '(m/^@/ || abs($F[8])<=100) && print' | samtools view -Sb - > ${output1}

# split bam, use reads with insert size >180
    samtools view -f 2 -h ${input}  | perl -ane '(m/^@/ || abs($F[8])>=180) && print' | samtools view -Sb - > ${output2}
fi

## index files
for i in ${workpath}/${outname}/${NAME}.filter.duprm.mer.bam ${workpath}/${outname}/${NAME}.filter.duprm.0mer.bam ${workpath}/${outname}/${NAME}.uniq_filtered.bam ${workpath}/${outname}/${NAME}.marked_duplicates.filtered.bam; do
    if [ ! -f ${i}.bai ];then
    echo "indeex"
        samtools index ${i} ${i}.bai
    fi
done



