#! /bin/bash


. functions.sh

r="yes" # run main work
s="yes" # run statistics

#step 1
#create output folder, scripts folder  and set path to bash files

out=$(grep outname setting.txt | awk '{print $2}')
mkdir -p $out
bashpath="bash/"
mkdir -p $out/scripts

#step 2
# generatate PBS header for files

for ff in cleanUpDir compPeaks calculate_tpm_in_peakregions calculate_tpm_in_narrowpeak generateFiles generateTables generateStatistics mergeReplicates run_macs run_statistics run_samtools_mult run_samtools_uniq getIS tableCov run_R_plots run_R_peaks count_reads_in_peakregions run_R_deseq count_reads_in_narrowpeaks; 
do
    createHead -f files.tab -s setting.txt -o $(pwd)/$out -w $ff > $out/scripts/${ff}.sh
    cat ${bashpath}${ff}_main.sh >> $out/scripts/${ff}.sh
done

#step 3
# setup submission
# main work

if [ "${r}" = "yes" ]; then
	echo "running" 
	a1=$(qsub -h $out/scripts/generateFiles.sh) # create the uniq filtered bam files
	b1=$(qsub -W depend=afterok:${a1} $out/scripts/mergeReplicates.sh) # merge bam files for replicates, depends on a1
	c1=$(qsub -W depend=afterok:${b1} $out/scripts/generateTables.sh) # generate file tables, depends on the merging of replicates (only!)
	d1=$(qsub -W depend=afterok:${c1} $out/scripts/run_macs.sh) # run macs on merged bam files, depends on a1, b1 and c1!
	e1=$(qsub -W depend=afterok:${d1} $out/scripts/run_R_peaks.sh) # merge peaks , generate gff (both merged and for each narrow)
#counts and tpm/fpkm 
	f1=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_narrowpeaks.sh) # narrow peaks 
	g1=$(qsub -W depend=afterok:${f1} $out/scripts/calculate_tpm_in_narrowpeak.sh)
	f2=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_peakregions.sh) # in the merged peaks
	g2=$(qsub -W depend=afterok:${f2} $out/scripts/calculate_tpm_in_peakregions.sh)
#deseq on merged peaks
	g3=$(qsub -W depend=afterok:${f2} $out/scripts/run_R_deseq.sh)

#start!
	qrls ${a1}
fi

# stats

if [ "${s}" = "yes" ]; then 
	echo "stats"
# statistics:
	st1=$(qsub $out/scripts/generateStatistics.sh)
	i=0
	for ff in  run_statistics run_samtools_mult run_samtools_uniq getIS tableCov;
 	do
  		st2[i]=$(qsub -W depend=afterok:${st1} $out/scripts/$ff.sh )
 		((i+=1))
	done

	st3=$(qsub -W depend=afterok:${st2[0]}:${st2[1]}:${st2[2]}:${st2[3]}:${st2[4]} $out/scripts/run_R_plots.sh)
fi

if [ "${r}" = "yes" ]; then 
	if [ "${s}" = "yes" ]; then 
		clen=$(qsub -W depend=afterok:${st3}:${g3}:${f1}:${f2} $out/scripts/cleanUpDir.sh)
	else
		clen=$(qsub -W depend=afterok:${g3}:${f1}:${f2} $out/scripts/cleanUpDir.sh)

	fi
	merge=$(qsub -W depend=afterok:${clen} $out/scripts/compPeaks.sh)
fi









## add summary - R mark down report, change script files, output path etc
