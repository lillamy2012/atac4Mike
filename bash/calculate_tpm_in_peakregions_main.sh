

awk '{ $11 = $5 - $4 } { print $10 "\t" $11} ' $gff | tr -d '""' | tr -d ";" > $out
join -j1 <(sort -k 1,1 $counts) <(sort -k1,1 $out) > $joined
awk '{$4=$2/($3/1000)} 1' $joined > $joined.tmp
mv $joined.tmp $joined
total=$(awk 'BEGIN{ total=0 } { total=total+$4 } END{ printf total }' $joined)
awk -v total=$total '{ print $1 "\t" 1000000*$4/total }' $joined > $tpm