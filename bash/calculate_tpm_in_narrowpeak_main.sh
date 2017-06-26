


PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt
outname=$(grep outname $settingsfile | awk '{print $2}')
type=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $4}')

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
ALIGN=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $2}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')
TOPATH=${workpath}/${outname}/${type}

filetype="uniq_filtered"
peaktype=( "MACS2_paired" "MACS2_noInput" )

for j in "${peaktype[@]}"
do
    gff=${TOPATH}/${type}.${filetype}.bam_${j}_peaks.narrowPeak.gff
echo $gff
    counts=${TOPATH}/$NAME.${j}_peaks.narrowPeak_counts.tab 
echo $counts
    out=${TOPATH}/$NAME.${j}_length.tab
echo $out
    joined=${TOPATH}/$NAME.${j}_joined.tab
echo $joined
    tpm=${TOPATH}/$NAME.${j}_tpm.tab
echo $tmp

    awk '{ $11 = $5 - $4 } { print $10 "\t" $11} ' $gff | tr -d '""' | tr -d ";" > $out
    join -j1 <(grep -v __ $counts | sort -t$'\t' -n --key=1.5 ) <(sort -t$'\t' -n --key=1.5 $out) > $joined
    awk '{$4=$2/($3/1000)} 1' $joined > $joined.tmp
    mv $joined.tmp $joined
    total=$(awk 'BEGIN{ total=0 } { total=total+$4 } END{ printf total }' $joined)
	echo $total
    awk -v total=$total '{ print $1 "\t" 1000000*$4/total }' $joined > $tpm

done
