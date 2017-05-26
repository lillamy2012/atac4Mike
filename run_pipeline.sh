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

for ff in generateFiles generateTables run_macs run_statistics run_samtools_mult run_samtools_uniq getIS tableCov run_R_plots run_R_peaks count_reads_in_peakregions run_R_deseq; do
    createHead -f files.tab -s setting.txt -o $(pwd)/$out -w $ff > $out/${ff}.sh
    cat ${bashpath}${ff}_main.sh >> $out/${ff}.sh
done


#step 3

# setup submission

## level a:

i=0
for ff in generateFiles generateTables; do
   	a[i]=$(qsub -h $out/${ff}.sh)
	((i+=1))
done

## level b

i=0
for ff in  run_statistics run_macs run_samtools_mult run_samtools_uniq getIS tableCov; do
    b[i]=$(qsub -W depend=afterok:${a[1]}:${a[0]} $out/$ff.sh )
#b[i]=$(qsub -h $ff.sh )
	((i+=1))
done


## level c
ff=run_R_plots
cl=$(qsub -W depend=afterok:${b[0]}:${b[1]}:${b[2]}:${b[3]}:${b[4]}:${b[5]} $out/$ff.sh)

## level d
ff=run_R_peaks
d=$(qsub -W depend=afterok:$cl $out/$ff.sh)


## level e
ff=count_reads_in_peakregions
e=$(qsub -W depend=afterany:$d $out/$ff.sh)

## level f
ff=run_R_deseq
fl=$(qsub -W depend=afterok:$e $out/$ff.sh)

## level g
#summarize


## release first level

qrls ${a[0]}
qrls ${a[1]}

#qrls ${b[0]}
#qrls ${b[1]}
#qrls ${b[2]}
#qrls ${b[3]}
#qrls ${b[4]}
#qrls ${b[5]}







## add summary - R mark down report, change script files, output path etc
