
# === begin ENVIRONMENT SETUP ===
# Load the required modules
# === end ENVIRONMENT SETUP ===


PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt
workpath=$(grep workpath $settingsfile | awk '{print $2}')
NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
out=$(grep outname $settingsfile | awk '{print $2}')

filetypes=( uniq_filtered.bam )

for f in "${filetypes[@]}"
do
awk '{print $3}' ${workpath}/$out/$NAME.$f.depth.txt | python $PBS_O_WORKDIR/python/my_counter.py > ${workpath}/$out/$NAME.$f.depth.table.tab
done








