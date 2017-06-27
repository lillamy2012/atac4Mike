

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
type=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')
#iPEAKPATH=${workpath}/${outname}/macs_results

echo "NAME: $NAME"
echo "workpath: $workpath"
echo "ALIGN: $ALIGN"

filetype=uniq_filtered
peaktype=( "MACS2_paired" "MACS2_noInput" )

#filename="${FULLPATH##*/}"

mkdir -p ${workpath}/${outname}/counts

for j in "${peaktype[@]}"
do
    gff=${workpath}/${outname}/${type}.${filetype}.bam_${j}_peaks.narrowPeak.gff
    echo $gff

# === START ===

    if [ ! -f $workpath/${outname}/counts/${NAME}_${j}_peaks.narrowPeak_counts.tab ];
    then
        echo "htseq"
        htseq-count -f bam -s no ${workpath}/${outname}/${type}/$NAME.$filetype.bam ${gff} > ${workpath}/${outname}/counts/${NAME}_${j}_peaks.narrowPeak_counts.tab
    fi
done

