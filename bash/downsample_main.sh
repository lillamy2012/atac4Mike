

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

FULLNAME=( $(awk '{print $1}' $PARAM_FILE | sed s/".subset.bam"//) )
echo $FULLNAME
workpath=$(grep workpath $settingsfile | awk '{print $2}')
outname=$(grep outname $settingsfile | awk '{print $2}')

i=0
for n in "${FULLNAME[@]}"; do 
 	y=$(awk 'NR==1{print $1}' $n.flagstat.txt)
	arrayNumber[$i]=$y
	if [ $i -eq 0 ];then 
 		min=$y
 	fi
	(( i++ ))
 	((y < min)) && min=$y
done
j=0
for i in "${arrayNumber[@]}"; do
	if [ ! -f ${FULLNAME[j]}.subset.bam ]; then 
		sub=$(echo "scale=2; $min/$i" | bc -l)
		if [ $(echo $sub'<'1 | bc -l) -eq 1 ]; then 
			java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar DownsampleSam  I=${FULLNAME[j]} O=${FULLNAME[j]}.subset.bam P=$sub
		else
			ln -s ${FULLNAME[j]} ${FULLNAME[j]}.subset.bam
		fi
	fi	
	(( j++ ))
done

