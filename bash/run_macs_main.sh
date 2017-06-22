

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load cutadapt/1.9-ictce-5.3.0-Python-2.7.6
module load Bowtie2/2.1.0-goolf-1.4.10
module load Picard/1.141
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files_macs.tab
settingsfile=$PBS_O_WORKDIR/setting.txt

FULLNAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')


echo "NAME: $FULLNAME"
echo "workpath: $workpath"

path="${FULLNAME%/*}"
filename="${FULLNAME##*/}"


# === START ===
## This script starts with a aligned and sorted bam file

export TMPDIR=$WORKDIR/Mike/picard_tmp
# make sure the directory exists
mkdir -p $TMPDIR




if [ ! -f ${FULLNAME}_MACS2.done ]; then 
	 macs2 callpeak -t $FULLNAME -g 1.2e8 -n ${filename}_MACS2_noInput --outdir ${path} --nomodel --shift -100 --extsize 200 -B -q 0.01
         macs2 callpeak -t $full -g 1.2e8 -n ${filename}_MACS2_paired  --outdir ${path}  -B -q 0.01 -f BAMPE  

	 touch ${FULLNAME}_MACS2.done
fi



