#! /bin/bash

## create head - generateFiles
createHead() { local OPTIND;

    while getopts ":f:s:o:w:" opt; do
        case $opt in
        f)
        local file=$OPTARG
        ;;
        s)
        local setting=$OPTARG
        ;;
        o)
        local out=$OPTARG
        ;;
        w)
        local which=$OPTARG
        esac
    done

    nrBJ=$(wc -l $file | awk '{print $1}') # simply the number of rows in files.tab
    nrPJ=$(awk '{print $3}' files.tab | grep sample | wc -l)
    nrMJ=$(awk '{print $4}' files.tab | sort | uniq | wc -l )


    case $which in
    generateFiles) # generateFiles
        printf %"s\n" "#!/bin/bash" "#PBS -N generate_files" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=05:00:00" "#PBS -e $out/generate_files_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=62gb" "#PBS -J 1-$nrBJ"
    ;;
    generateStatistics) # generateStatistics
        printf %"s\n" "#!/bin/bash" "#PBS -N generate_statistics" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=05:00:00" "#PBS -e $out/generate_statistics_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=62gb" "#PBS -J 1-$nrBJ"
    ;;
    mergeReplicates) #
        printf %"s\n" "#!/bin/bash" "#PBS -N merge_reps" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=02:00:00" "#PBS -e $out/merge_reps_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;
    generateTables) # generateTables
        printf %"s\n" "#!/bin/bash" "#PBS -N generate_tables" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:02:00" "#PBS -e $out/generate_tables_error_out" "#PBS -j eo"
    ;;
    run_macs) ## fix J to not include input
        printf %"s\n" "#!/bin/bash" "#PBS -N run_macs2" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=10:00:00" "#PBS -e $out/run_macs_error_out" "#PBS -j eo" "#PBS -J 1-$nrMJ" "#PBS -l select=1:ncpus=16:mem=62gb"
    ;;
    run_statistics)
        printf %"s\n" "#!/bin/bash" "#PBS -N run_stats" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:10:00" "#PBS -e $out/run_statistics_error_out" "#PBS -j eo"
    ;;
    run_samtools_uniq)
        printf %"s\n" "#!/bin/bash" "#PBS -N run_samt_un" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=03:00:00" "#PBS -e $out/run_samtools_uniq_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=20gb" "#PBS -J 1-$nrBJ"
    ;;
    run_samtools_mult)
        printf %"s\n" "#!/bin/bash" "#PBS -N run_samt_mult" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=03:00:00" "#PBS -e $out/run_samtools_mult_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=20gb" "#PBS -J 1-$nrBJ"
    ;;
    getIS)
        printf %"s\n" "#!/bin/bash" "#PBS -N getIS" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=01:00:00" "#PBS -e $out/getIS_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=20gb" "#PBS -J 1-$nrBJ"
    ;;
    tableCov)
        printf %"s\n" "#!/bin/bash" "#PBS -N tableCov" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:30:00" "#PBS -e $out/tableCov_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;
    run_R_plots)
         printf %"s\n" "#!/bin/bash" "#PBS -N run_R_plots" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=01:30:00" "#PBS -e $out/run_R_plots_error_out" "#PBS -j eo" "#PBS -l select=1:ncpus=16:mem=20gb"
    ;;
    run_R_peaks)
        printf %"s\n" "#!/bin/bash" "#PBS -N run_R_peaks" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:30:00" "#PBS -e $out/run_R_peaks_error_out" "#PBS -j eo"
    ;;
    count_reads_in_peakregions)
        printf %"s\n" "#!/bin/bash" "#PBS -N count_reads_peaks" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=05:00:00" "#PBS -e $out/count_reads_peaks_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;
    count_reads_in_narrowpeaks)
        printf %"s\n" "#!/bin/bash" "#PBS -N count_reads_narrowpeaks" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=05:00:00" "#PBS -e $out/count_reads_narrowpeaks_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;
    run_R_deseq)
        printf %"s\n" "#!/bin/bash" "#PBS -N run_R_deseq" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:10:00" "#PBS -e $out/run_R_deseq_error_out" "#PBS -j eo"
;;
   calculate_tpm_in_peakregions)
        printf %"s\n" "#!/bin/bash" "#PBS -N calc_tpm_peaks" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:10:00" "#PBS -e $out/calc_tpm_peaks_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;
   calculate_tpm_in_narrowpeak)
        printf %"s\n" "#!/bin/bash" "#PBS -N calc_tpm_narrowpeaks" "#PBS -P berger_common" "#PBS -q workq" "#PBS -l walltime=00:10:00" "#PBS -e $out/calc_tpm_narrowpeaks_error_out" "#PBS -j eo" "#PBS -J 1-$nrBJ"
    ;;



esac


}



