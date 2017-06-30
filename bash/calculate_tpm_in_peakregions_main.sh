


PARAM_IDX=$PBS_ARRAY_INDEX
PARAM_FILE=$PBS_O_WORKDIR/files.tab
settingsfile=$PBS_O_WORKDIR/setting.txt
outname=$(grep outname $settingsfile | awk '{print $2}')

NAME=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $1}')
ALIGN=$(sed "${PARAM_IDX}q;d" $PARAM_FILE | awk '{print $2}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')

filetype="uniq_filtered"
peaktype=( "MACS2_paired" "MACS2_noInput" )

for j in "${peaktype[@]}"
do
    gff=${workpath}/$outname/$j.gff
echo $gff
    counts=${workpath}/$outname/${NAME}_${j}.gff.counts.tab 
echo $counts
    out=${workpath}/$outname/$NAME.${j}_length.tab
echo $out
    joined=${workpath}/$outname/$NAME.${j}_joined.tab
echo $joined
    tpm=${workpath}/$outname/${NAME}_${j}.gff.tpm.tab
echo $tmp
    fpkm=${workpath}/$outname/${NAME}_${j}.gff.fpkm.tab
echo $fpkm

    awk '{ $11 = $5 - $4 } { print $10 "\t" $11} ' $gff | tr -d '""' | tr -d ";" > $out
    join -j1 <(grep -v __ $counts | sort -t$'\t' -n --key=1.5 ) <(sort -t$'\t' -n --key=1.5 $out) > $joined
## fpkm     
    tot_read=$(awk 'BEGIN{ tot_read=0} {tot_read=tot_read+$2 } END{printf tot_read/1000000}' $counts)
    awk -v tot_read=$tot_read '{ print $1 "\t" ($2/tot_read)/($3/1000) }' $joined > $fpkm 
## tpm 
    awk '{$4=$2/($3/1000)} 1' $joined > $joined.tmp
    mv $joined.tmp $joined
    total=$(awk 'BEGIN{ total=0 } { total=total+$4 } END{ printf total }' $joined)
	echo $total
    awk -v total=$total '{ print $1 "\t" 1000000*$4/total }' $joined > $tpm
   
rm $out
rm $joined
done
