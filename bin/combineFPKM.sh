#!/usr/bin/env bash


gffs=( $(find . -name "*.gff") )

for g in "${gffs[@]}"; do
    gf=$(basename $g)
    outfile=${gf}_summary.tab
    awk '{print $10"\t"$4"\t"$5}' $g | sed 's/"//g' | sed 's/;//g'> $outfile
    printf "%s\t" "peaks" > header_${gf}.tab    
    echo $g
    if [ "$gf" = "master.gff" ];then
        echo "yes" 
        fpkm=( $(find . -name "*master_fpkm.tab") ) 
    else
    	fpkm=( $(find . -name "*${gf}_fpkm.tab") )
    fi 
    for f in "${fpkm[@]}"; do
         echo $f 
         join $f $outfile > test.tmp
	 mv test.tmp $outfile
     done
     for (( idx=${#fpkm[@]}-1 ; idx>=0 ; idx-- )) ; do
   	printf "%s\t" $(basename "${fpkm[idx]}") >> header_${gf}.tab
     done
     printf "%s\t" "start" "end"  >>  header_${gf}.tab
     printf "\n" >> header_${gf}.tab
     cat header_${gf}.tab $outfile  > tmpfile2 
     mv tmpfile2 $outfile
     rm header_${gf}.tab 

done

deseq="deseq_results.csv"

## handle comments 
grep '#' $deseq > comments
grep -v '#' $deseq > data


fpkm="master.gff_summary.tab"
out="deseq_total_results.csv"

join -j1 -t ";" <(awk 'FNR > 1 {print}' data | sort -n --key=1.5 ) <(perl -wnlp -e 's/\s+/;/g;' $fpkm | awk 'FNR > 1 {print}' | sort -n --key=1.5) > merge_tmp


## header
awk 'BEGIN {FS=","}; NR==1 {print}' data > h1
awk 'NR==1 {print}' $fpkm > h2

cat comments > $out
printf "Peak"\; >> $out
printf $(cat h1) >> $out
printf \; >> $out
printf $(cat h2 | perl -wnlp -e 's/\s+/;/g;' |  cut -d';' -f 2-)>> $out
printf "\n" >> $out
cat merge_tmp >> $out
cp $out out.tmp
nrC=$(awk -F";" '{print NF}' out.tmp | sort -nu | tail -n 1)
nn=$((nrC-3))
cut -d';' -f1-$nn out.tmp > $out

# rm h1 h2 merge_tmp out.tmp data comments







