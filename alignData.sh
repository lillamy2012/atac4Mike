#! /bin/bash

#PBS -N generate_atac_files
#PBS -P berger_common
#PBS -q workq
#PBS -l walltime=40:00:00
#PBS -l select=1:ncpus=6:mem=20gb
#PBS -J 2-4


# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load cutadapt/1.9-ictce-5.3.0-Python-2.7.6
module load Bowtie2/2.1.0-goolf-1.4.10
module load Picard/1.141
module load FASTX-Toolkit/0.0.13.2-goolf-1.4.10
module load MACS/2.1.0.20150420.1-goolf-1.4.10-Python-2.7.5
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files2.tab
settingsfile=$PBS_O_WORKDIR/setting.txt

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
bampath=$(grep bampath $settingsfile | awk '{print $2}')
index=$(grep index $settingsfile | awk '{print $2}')
alignpath=$(grep alignpath $settingsfile | awk '{print $2}')
l_trim=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $6}')
echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "l_trim: $l_trim"
# === START ===



export TMPDIR=$PBS_O_WORKDIR/${NAME}_tmp
# make sure the directory exists
mkdir -p $TMPDIR
echo $TMPDIR

#output folder
myout=$PBS_O_WORKDIR/${NAME}
mkdir -p $myout
echo $myout

input=${bampath}/${NAME}.bam
output1=${myout}/${NAME}_R1_.fastq

#trim to same length if different sequencing lengths 
if [ ! -f ${output1} ];then
    echo "fastq"
# fastx_trimmer [-l N]       = Last base to keep. Default is entire read.
    samtools view -f 0x40 -b ${input} | samtools bam2fq - | fastx_trimmer -l $l_trim -Q33 -o ${output1}
fi

output2=${myout}/${NAME}_R2_.fastq
if [ ! -f ${output2} ];then
# fastx_trimmer [-l N]       = Last base to keep. Default is entire read.
   samtools view -f 0x80 -b ${input} | samtools bam2fq -| fastx_trimmer -l $l_trim -Q33 -o ${output2}
fi


## cutadapt, if less than 5 bp left remove read

input1=${output1}
input2=${output2}
output1=${myout}/${NAME}_R1_cut1.fastq
output2=${myout}/${NAME}_R2_cut1.fastq
output11=${myout}/${NAME}_R1_cut3.fastq
output22=${myout}/${NAME}_R2_cut3.fastq



if [ ! -f ${output1} ] || [ ! -f  ${output2} ];then
   echo "cutadapt"
    cutadapt --minimum-length 5 --overlap 1 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o ${output1} -p ${output2} ${input1} ${input2}
cutadapt --minimum-length 5 --overlap 3 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -o ${output11} -p ${output22} ${input1} ${input2}
fi

## (align with trim to avoid problem with short residual adapter seq)
input1=${output1}
input2=${output2}
input11=${output11}
input22=${output22}
output=${myout}/${NAME}.aligned_cut1.sam
output2=${myout}/${NAME}.aligned_cut3.sam

if [ ! -f ${output} ];then
    echo "bowtie2"
   bowtie2 --end-to-end -x ${index} -1 ${input1} -2 ${input2} -S ${output}
   bowtie2 --end-to-end --trim3 3 -x ${index} -1 ${input11} -2 ${input22} -S ${output2}
fi

#fi

## sam to sorted bam
input=$output
input2=$output2
output=${myout}/${NAME}.aligned_sorted_cut1.bam
output2=${myout}/${NAME}.aligned_sorted_cut3.bam
if [ ! -f ${output} ];then
	echo "sort"
	samtools view -b ${input} | samtools sort -T $TMPDIR - > ${output}
	samtools view -b ${input2} | samtools sort -T $TMPDIR - > ${output2}
fi

