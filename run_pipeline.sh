#! /bin/bash


. functions.sh
## read in 3 things
# 1, settings file, 2, files file 3, output location, overwrite=FALSE, steps to run

# variables needed #jobs, location for output/error, #conditions, colors...,

#step 0 read in arg


#step 1

#create output folder

out=$(grep outname setting.txt | awk '{print $2}')
mkdir -p $out

bashpath="bash/"

#step 2

# generatate PBS header for files

mkdir -p $out/scripts

for ff in calculate_tpm_in_peakregions calculate_tpm_in_narrowpeak generateFiles generateTables generateStatistics mergeReplicates run_macs run_statistics run_samtools_mult run_samtools_uniq getIS tableCov run_R_plots run_R_peaks count_reads_in_peakregions run_R_deseq count_reads_in_narrowpeaks; do
    createHead -f files.tab -s setting.txt -o $(pwd)/$out -w $ff > $out/scripts/${ff}.sh
    cat ${bashpath}${ff}_main.sh >> $out/scripts/${ff}.sh
done


#step 3

# setup submission

## level a: ## generate all alignment files and tables

#ind1=$(qsub $out/scripts/generateStatistics.sh)
#a1=$(qsub -h $out/scripts/generateFiles.sh)
#b1=$(qsub -W depend=afterok:${a1} $out/scripts/mergeReplicates.sh)
#c1=$(qsub -W depend=afterok:${b1} $out/scripts/generateTables.sh)
qsub $out/scripts/generateTables.sh 
#d1=$(qsub -W depend=afterok:${c1} $out/scripts/run_macs.sh)
#e1=$(qsub -W depend=afterok:${d1} $out/scripts/run_R_peaks.sh)
#f1=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_narrowpeaks.sh)
#f2=$(qsub -W depend=afterok:${e1} $out/scripts/count_reads_in_peakregions.sh)
#g1=$(qsub -W depend=afterok:${f1} $out/scripts/calculate_tpm_in_narrowpeak.sh)
#g2=$(qsub -W depend=afterok:${f2} $out/scripts/calculate_tpm_in_peakregions.sh) 

#i=0
#for ff in generateFiles generateStatistics; do
 #  	a[i]=$(qsub -h $out/scripts/{ff}.sh)
#	((i+=1))
#done

#qrls ${a1}


####ae=$(qsub $out/generateStatistic.sh) ## only need to check in end that it's finished

## level b: ## run statistics on alignment files, run macs2

#$(qsub -W depend=afterok:${a1} $out/mergeReplicates.sh)
#i=0
#for ff in  run_statistics run_macs run_samtools_mult run_samtools_uniq getIS tableCov; do
#    b[i]=$(qsub -W depend=afterok:${a[1]}:${a[0]} $out/$ff.sh )
#b[i]=$(qsub -h $ff.sh )
#	((i+=1))
#done


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
