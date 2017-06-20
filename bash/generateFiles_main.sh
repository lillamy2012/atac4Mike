
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
type=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
bampath=$(grep bampath $settingsfile | awk '{print $2}')
index=$(grep index $settingsfile | awk '{print $2}')
alignpath=$(grep alignpath $settingsfile | awk '{print $2}')
quality=$(grep Q_threshold $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')
tmp=$(grep tmp_file $settingsfile | awk '{print $2}')


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


## quality filter and remove duplicates

mkdir -p ${workpath}/${outname}/${type}

input=${alignpath}/${ALIGN}.bam
output=${workpath}/${outname}/${type}/${NAME}.uniq_filtered.bam
tmpfile=${workpath}/${tmp}/${NAME}.sam_filtered.bam

if [ ! -f ${output} ];then
    samtools view -bSq $quality ${input}> ${tmpfile}
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${tmpfile}  O= ${output} M=${NAME}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
fi



