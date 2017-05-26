

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
input=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $6}')
#subn=$(grep subn $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')


echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "ALIGN: $ALIGN"
echo "input: $input"

# === START ===
## This script starts with a aligned and sorted bam file

export TMPDIR=$WORKDIR/Mike/picard_tmp
# make sure the directory exists
mkdir -p $TMPDIR


#if [ "$input" != "NA" ]; then

# call peaks MACS2 - on all 3 bam files (or only combined)

    for i in uniq_filtered.bam; do

#${workpath}/${NAME}.marked_duplicates.filtered.bam; do
        full=${workpath}/$outname/${NAME}.$i
        filename=$(basename "${full}")
        filename="${filename%.*}"
 #       input_file=${workpath}/$outname/${input}.$i
        if [ ! -f ${workpath}/$outname/${filename}_MACS2peaks.done ];then
            echo "macs"
            echo ${filename}
            macs2 callpeak -t $full -g 1.2e8 -n ${filename}_MACS2_noInput --outdir ${workpath}/$outname --nomodel --shift -100 --extsize 200 -B -q 0.01 
            macs2 callpeak -t $full -g 1.2e8 -n ${filename}_MACS2_paired  --outdir ${workpath}/$outname  -B -q 0.01 -f BAMPE
#macs2 callpeak -t $full -c $input_file -g 1.2e8 -n ${filename}_MACS2_withInput --outdir ${workpath}/$outname --nomodel --shift -100 --extsize 200 -B -q 0.01 -f BAMPE
            touch ${workpath}/$outname/${filename}_MACS2peaks.done
        fi
    done
#fi

## next we should stage out the files (MACS peaks and final bam files), run script to retrieve statistsic

## stage out selected files - ${workpath}/${NAME}.sorted.duprm.mer.bam ${workpath}/${NAME}.sorted.duprm.0mer.bam ${workpath}/${NAME}.marked_duplicates.bam  ``


