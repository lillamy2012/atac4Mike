


gff=$(find ./ -name "*.gff")
counts=$(find ./ -name "*counts.tab")
out=length.tab
joined=joined.tab


awk '{ $11 = $5 - $4 } { print $10 "\t" $11} ' $gff | tr -d '""' | tr -d ";" > $out
    join -j1 <(grep -v __ $counts | sort -t$'\t' -n --key=1.5 ) <(sort -t$'\t' -n --key=1.5 $out) > $joined
## fpkm     
    tot_read=$(awk 'BEGIN{ tot_read=0} {tot_read=tot_read+$2 } END{printf tot_read/1000000}' $counts)
    awk -v tot_read=$tot_read '{ print $1 "\t" ($2/tot_read)/($3/1000) }' $joined




