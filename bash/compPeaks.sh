

settingsfile="../setting.txt"
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
