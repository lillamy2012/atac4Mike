

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load R/3.2.3-foss-2016a
module load Python/2.7.11-intel-2016a
# === end ENVIRONMENT SETUP ===

settingsfile=$PBS_O_WORKDIR/setting.txt
workpath=$(grep workpath $settingsfile | awk '{print $2}')
out=$(grep outname $settingsfile | awk '{print $2}')
mkdir -p $workpath/$out/plots
cd $workpath/$out


## barplots.R
files=( $(ls | grep "R_" | grep "numbers") )
for f in "${files[@]}"
do
    Rscript --vanilla $PBS_O_WORKDIR/R/barplot.R $f
done

#insertplots.R

files=( $(ls | grep "R_" | grep "IS") )
for f in "${files[@]}"
do
    Rscript --vanilla $PBS_O_WORKDIR/R/insertplots.R $f
done

#analyseCov.R

files=( $(ls | grep "R_" | grep "depth.table") )
for f in "${files[@]}"
do
    Rscript --vanilla $PBS_O_WORKDIR/R/analyseCov.R $f
done

mv *.pdf plots
mv *.png plots
