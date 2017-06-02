mkdir -p bam


nr=$(wc -l files2.tab | awk '{print $1}')
for i in `seq 1 $nr`; do 
cp $(sed "${i}q;d" files2.tab | awk '{print $1}')/$(sed "${i}q;d" files.tab | awk '{print $2}').bam bam
done
