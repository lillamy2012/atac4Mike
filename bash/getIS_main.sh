
# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
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
out=$(grep outname $settingsfile | awk '{print $2}')
tmp=$(grep tmp_path $settingsfile | awk '{print $2}')

echo "NAME: $NAME"
echo "workpath: $workpath"
echo "bampath: $bampath"
echo "index: $index"
echo "ALIGN: $ALIGN"
echo $(pwd)

filetypes=(  uniq_filtered.bam )

# === START ===

#if [ ! -f ${workpath}/$NAME.mult_stat.tab ]; then
#printf "#%s\n" "${filetypes[@]}"  > ${workpath}/$NAME.mult_stat.tab

for f in "${filetypes[@]}"
do
## multi (proper only)
echo $NAME.$f
samtools view -f 66 ${workpath}/$tmp/$NAME.$f | awk '{print $9}' | python  $PBS_O_WORKDIR/python/my_counter.py > ${workpath}/$out/$NAME.$f.IS.tab
done

