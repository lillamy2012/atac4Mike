

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
mytype=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')
type=( $(awk '{print $4}' $PARAM_FILE | sort | uniq) )

echo "NAME: $NAME"
echo "workpath: $workpath"
echo "ALIGN: $ALIGN"

filetype=uniq_filtered
peaktype=( "MACS2_paired" "MACS2_noInput" )

#filename="${FULLPATH##*/}"

mkdir -p ${workpath}/${outname}/counts

for j in "${peaktype[@]}"
do
for t in "${type[@]}"
    do
    gffname=${t}.${filetype}.bam.subset.bam_${j}_peaks.narrowPeak.gff
    gff=${workpath}/${outname}/${gffname}
    echo $gff

# === START ===

    if [ ! -f $workpath/${outname}/counts/${NAME}_${gffname}.counts.tab ];
    then
        echo "htseq"
        htseq-count -f bam -s no ${workpath}/${outname}/${mytype}/$NAME.$filetype.bam ${gff} > $workpath/${outname}/counts/${NAME}_${gffname}.counts.tab 
    fi
done
done
