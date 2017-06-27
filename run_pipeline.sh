#! /bin/bash


. functions.sh

#step 1
#create output folder, scripts folder  and set path to bash files

out=$(grep outname setting.txt | awk '{print $2}')
mkdir -p $out
bashpath="bash/"
mkdir -p $out/scripts

#step 2
# generatate PBS header for files

for ff in calculate_tpm_in_peakregions calculate_tpm_in_narrowpeak generateFiles generateTables generateStatistics mergeReplicates run_macs run_statistics run_samtools_mult run_samtools_uniq getIS tableCov run_R_plots run_R_peaks count_reads_in_peakregions run_R_deseq count_reads_in_narrowpeaks; 
do
    createHead -f files.tab -s setting.txt -o $(pwd)/$out -w $ff > $out/scripts/${ff}.sh
    cat ${bashpath}${ff}_main.sh >> $out/scripts/${ff}.sh
done

r="yes"
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
#counts and tpm 
f1=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_narrowpeaks.sh) # narrow peaks FIXME count in all peaks!
g1=$(qsub -W depend=afterok:${f1} $out/scripts/calculate_tpm_in_narrowpeak.sh)
f2=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_peakregions.sh) # in the merged peaks
g2=$(qsub -W depend=afterok:${f2} $out/scripts/calculate_tpm_in_peakregions.sh)

#deseq on merged peaks
g3=$(qsub -W depend=afterok:${f2} $out/scripts/run_R_deseq.sh)

#combine tables (add tpm)
#clean up script

#start!
qrls ${a1}

fi

# statistics:
st1=$(qsub $out/scripts/generateStatistics.sh)
i=0
for ff in  run_statistics run_samtools_mult run_samtools_uniq getIS tableCov;
 do
  st2[i]=$(qsub -W depend=afterok:${st1} $out/scripts/$ff.sh )
 ((i+=1))
done

st3=$(qsub -W depend=afterok:${st2[0]}:${st2[1]}:${st2[2]}:${st2[3]}:${st2[4]} $out/scripts/run_R_plots.sh)




## level c ## rplots
#ff=run_R_plots
#cl=$(qsub -W depend=afterok:${b[0]}:${b[1]}:${b[2]}:${b[3]}:${b[4]}:${b[5]} $out/$ff.sh)

## level d ## could be level c, merging MACS2 peaks - add so also gff of individual peaks
#ff=run_R_peaks
#d=$(qsub -W depend=afterok:$cl $out/$ff.sh)


## level e ## count - add also ind
#ff=count_reads_in_peakregions
#e=$(qsub -W depend=afterany:$d $out/$ff.sh)

## level f ## run deseq and calculate tpm (should tpm be added to ouput)
#ff=run_R_deseq
#fl=$(qsub -W depend=afterok:$e $out/$ff.sh)
#ff2=run_tpm
#ffl=$(qsub -W depend=afterok:$e $out/$ff2.sh)


## level g
#summarize


## release first level

#qrls ${a[0]}
#qrls ${a[1]}

#qrls ${b[0]}
#qrls ${b[1]}
#qrls ${b[2]}
#qrls ${b[3]}
#qrls ${b[4]}
#qrls ${b[5]}







## add summary - R mark down report, change script files, output path etc
