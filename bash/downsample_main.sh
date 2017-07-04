

# === begin ENVIRONMENT SETUP ===
# Load the required modules
module load SAMtools/1.3-goolf-1.4.10
module load picard/2.3.0
# === end ENVIRONMENT SETUP ===


export TMPDIR=$WORKDIR/Mike/picard_tmp
# make sure the directory exists
mkdir -p $TMPDIR
settingsfile=$PBS_O_WORKDIR/setting.txt
PARAM_FILE=$PBS_O_WORKDIR/files_macs.tab #$PBS_O_WORKDIR

FULLNAME=( $(awk '{print $1}' $PARAM_FILE) )
workpath=$(grep workpath $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')

i=0
max=0
for n in "${FULLNAME[@]}"; do 
 	y=$(awk 'NR==1{print $1}' $n.flagstat.txt)
	arrayNumber[$i]=$y
	(( i++ ))
 	((y > max)) && max=$y
done
j=0
for i in "${arrayNumber[@]}"; do
	sub=$(echo "scale=2; $i/$max" | bc -l)
	if [ $(echo $sub'<'1 | bc -l) -eq 1 ]; then 
		java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar DownsampleSam  I=${FULLNAME[j]} O=${FULLNAME[j]}.subset.bam P=$sub
	else
		ln -s ${FULLNAME[j]} ${FULLNAME[j]}.subset.bam
	fi	
	(( j++ ))
done
