#! /bin/bash

#PBS -N subsamp
#PBS -P berger_common
#PBS -q workq
#PBS -l walltime=01:00:00
#PBS -J 1-7


# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
# === end ENVIRONMENT SETUP ===

PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files2.tab
settingsfile=$PBS_O_WORKDIR/setting.txt

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
bampath=$(grep bampath $settingsfile | awk '{print $2}')
index=$(grep index $settingsfile | awk '{print $2}')
alignpath=$(grep alignpath $settingsfile | awk '{print $2}')

# === START ===

myout=$PBS_O_WORKDIR/bam

input=${myout}/${NAME}.aligned_sorted_cut1.bam

output=${myout}/${NAME}.aligned_sorted_cut1_small.bam

samtools view -s 0.01 -b $input > $output
