

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
	#index=

	tps=( $(find $tpmPath -type f -name *_$gf.tpm.tab) ) 

	for t in "${tps[@]}"; do 
		echo $t
	done
done 
