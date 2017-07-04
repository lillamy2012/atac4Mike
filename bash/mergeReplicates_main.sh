# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt

type=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')
rep=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $5}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')
tmp=$(grep tmp_file $settingsfile | awk '{print $2}')

if [ $rep -eq 1 ]; then

    output=${workpath}/${outname}/${type}/${type}.uniq_filtered.bam
    if [ ! -f ${output} ];then
        cd ${workpath}/${outname}/${type}/
        samtools merge $output *.bam
	samtools flagstat ${output} > ${output}.flagstat.txt 
   fi
fi
