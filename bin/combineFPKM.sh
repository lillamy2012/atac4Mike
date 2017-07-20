#!/usr/bin/env bash


gffs=( $(find . -name "*.gff") )

for g in "${gffs[@]}"; do
    gf=$(basename $g)
    #echo $g
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


