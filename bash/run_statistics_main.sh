


# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
# === end ENVIRONMENT SETUP ===

#BS_O_WORKDIR="../"
echo $PBS_O_WORKDIR
PARAM_FILE=$PBS_O_WORKDIR/files.tab
nr=$(wc -l $PARAM_FILE | awk '{print $1}')
echo $nr
settingsfile=$PBS_O_WORKDIR/setting.txt
workpath=$(grep workpath $settingsfile | awk '{print $2}')
out=$(grep outname $settingsfile | awk '{print $2}')

bashpath="bash"
echo $workpath

for id in $(seq 1 $nr)
do
    NAME=$(sed "${id}q;d" $PARAM_FILE | awk '{print $1}')
    echo $NAME
    filetypes=( uniq_filtered.bam )

    ## extract/parse the statistics files (stats for barplot,...)
    for f in "${filetypes[@]}"
    do
        echo $f
        bash ${workpath}/$bashpath/collectStats.sh ${workpath}/$out/results/$NAME.$f.stats.txt > ${workpath}/$out/$NAME.$f.numbers.tab
    done
done


## stats = numbers.tab

###############
## all files are generated and the uniq mult mapping + insert sizes are run

## extract/parse the statistics files (stats for barplot,...)
## generare files with filenames to be read into R
