#!/usr/bin/env nextflow

//params.design     =
//params.macs_call  = 
params.bams       = "../bams/*small.bam" 
params.quality    = 10  

bamset = Channel
        .fromPath(params.bams)
        .map { file -> tuple(file.baseName, file) }


process generateFiles {
tag "bam : {name}"  
  
    input:
    set name, file(bam) from bamset

    output:
    set name, file("${name}_uniq_filtered.bam") into uniq_filtered

    script:
    """
    samtools view -bSq ${params.quality} ${bam} > ${name}_tmp.bam
    java -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${name}_tmp.bam  O= ${name}_uniq_filtered.bam} M=${NAME}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true TMP_DIR=$TMPDIR
    """
}

process mergeReplicates {
}

process subsampleMerged {
}

process callMACS2 {
}

process makeMasterPeaks {
}

process count_reads_in_master {
}

process count_reads_in_narrow {
}

process calculate_fpkm_in_master {

}
process calculate_fpkm_in_narrow {
}

process deseq2 {
}

