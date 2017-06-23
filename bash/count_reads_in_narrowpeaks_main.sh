

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load  HTSeq/0.6.1-goolf-1.4.10-Python-2.7.3
module load pysam/0.8.2.1-goolf-1.4.10-Python-2.7.3
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
FULL=$PBS_O_WORKDIR/files_macs.tab
settingsfile=$PBS_O_WORKDIR/setting.txt
outname=$(grep outname $settingsfile | awk '{print $2}')

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
ALIGN=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $2}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
bampath=$(grep bampath $settingsfile | awk '{print $2}')
index=$(grep index $settingsfile | awk '{print $2}')
alignpath=$(grep alignpath $settingsfile | awk '{print $2}')
#subn=$(grep subn $settingsfile | awk '{print $2}')
type=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')
FULLPATH=${workpath}/${outname}/${type}/${type}.uniq_filtered.bam

echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "ALIGN: $ALIGN"
echo $(pwd)
echo $FULLPATH

filetype=uniq_filtered
peaktype=( "MACS2_paired" "MACS2_noInput" )
path="${FULLPATH%/*}"
filename="${FULLPATH##*/}"


for j in "${peaktype[@]}"
do
    echo ${path}/$NAME.$filetype_${j}_peaks.narrowPeak_counts.tab
    gff=${FULLPATH}_${j}_peaks.narrowPeak.gff
    echo $gff

# === START ===

    if [ ! -f ${FULLPATH}_${j}_peaks.narrowPeak_counts.tab ];
    then
        echo "htseq"
        htseq-count -f bam -s no ${path}/$NAME.$filetype.bam ${gff} > ${path}/$NAME.$filetype_${j}_peaks.narrowPeak_counts.tab
    fi
done

