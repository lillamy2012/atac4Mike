

settingsfile=$PBS_O_WORKDIR/setting.txt
#settingsfile="../setting.txt"
outname=$(grep outname $settingsfile | awk '{print $2}')
workpath=$(grep workpath $settingsfile | awk '{print $2}')

gffPath=$workpath/$outname/gff
tpmPath=$workpath/$outname/tpm
countPath=$workpath/$outname/counts


gfs=( $(find $gffPath -type f -name *.gff ) )

for g in "${gfs[@]}"; do 
	echo $g
	gf=$(basename $g)
	echo $gf
	outfile=$workpath/$outname/outfiles/${gf}_summary.tab
#echo $outfile
	awk '{print $10"\t"$4"\t"$5"\t"$5-$4}' $g | sed 's/"//g' | sed 's/;//g'> $outfile
	printf "%s\t" "peaks" > header_${gf}.tab
	
	tps=( $(find $tpmPath -type f -name *$gf.fpkm.tab) ) 

	for t in "${tps[@]}"; do 
		join $t $outfile > test.tmp
		mv test.tmp $outfile
		echo $t 
	done
	for (( idx=${#tps[@]}-1 ; idx>=0 ; idx-- )) ; do
   		printf "%s\t" $(basename "${tps[idx]}") >> header_${gf}.tab
	done
	printf "%s\t" "start" "end" "width" >>  header_${gf}.tab
	printf "\n" >> header_${gf}.tab
	cat header_${gf}.tab $outfile  > tmpfile2 
	mv tmpfile2 $outfile
	rm header_${gf}.tab 
done

cd $workpath/$outname/outfiles/

anno="_anno.csv"
deseq="_deseq_type_sperm_veg.csv"
fpkm=".gff_summary.tab"

for type in "MACS2_paired" "MACS2_noInput"; do

    out=$type.total.csv

### deseq and fpkm

    join -j1 -t ";" <(awk 'FNR > 1 {print}' $type$deseq | sort -n --key=1.5 ) <(perl -wnlp -e 's/\s+/;/g;' $type$fpkm | awk 'FNR > 1 {print}' | sort -n --key=1.5) > $type.merge_tmp

## header

    awk 'BEGIN {FS=","}; NR==1 {print}' $type$deseq > $type.h1
    awk 'NR==1 {print}' $type$fpkm > $type.h2

    printf "Peak"\; > $out
    printf $(cat $type.h1) >> $out
    printf \; >> $out
    printf $(cat $type.h2 | perl -wnlp -e 's/\s+/;/g;' |  cut -d';' -f 2-)>> $out
    printf "\n" >> $out
    cat $type.merge_tmp >> $out
    cp $out $type.tmp
    cut -d';' -f1-23 $type.tmp > $out

    rm $type.h1 $type.h2 $type.merge_tmp $type.tmp
done



 
