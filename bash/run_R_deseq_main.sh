#

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load R/3.2.3-foss-2016a
module load Python/2.7.11-intel-2016a
# === end ENVIRONMENT SETUP ===

settingsfile=$PBS_O_WORKDIR/setting.txt
workpath=$(grep workpath $settingsfile | awk '{print $2}')
out=$(grep outname $settingsfile | awk '{print $2}')
cd $workpath/$out


##
files=( $(ls | grep "^R_" | grep "deseq") )
for f in "${files[@]}"
do
    echo $f
    Rscript --vanilla $PBS_O_WORKDIR/R/deseq2.R $f
done


