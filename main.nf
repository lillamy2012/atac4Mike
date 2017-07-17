#!/usr/bin/env nextflow

params.design     = 'exp.tab'
params.macs_call  = '-B -q 0.01 -f BAMPE'
params.genomesize = '1.2e8'
params.bams       = "bam/*small.bam" 
params.quality    = 10  

bamset = Channel
        .fromPath(params.bams)
        .map { file -> tuple(file.baseName, file) }

design = Channel
        .fromPath(params.design)
        .splitCsv()
        .map { row -> [ id:row[1],cond:row[2] ] }


process generateFiles {
tag "bam : $name"  
  
    input:
    set name, file(bam) from bamset

    output:
    set name, file("${name}_uniq_filtered.bam") into uniq_filtered

    script:
    """
    export TMPDIR=\$(pwd)
    samtools view -bSq ${params.quality} ${bam} > ${name}_tmp.bam
    java -Djava.io.tmpdir=\$TMPDIR -jar \$EBROOTPICARD/picard.jar MarkDuplicates I=${name}_tmp.bam  O= ${name}_uniq_filtered.bam  M=${name}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true 
   TMP_DIR=\$TMPDIR 
    """
}

uniq_filtered.into { for_rep; uniq_filtered }

comb = for_rep.map { it -> [ id:it[0], file:it[1] ] }.phase(design) {it -> it.id}
     .map { it -> tuple( it.cond[1] , it.file[0])  }
     .groupTuple(by: 0)


process mergeReplicates {
    
    input:
    set id, file(ff) from comb

    output:
    set id, file("${id}.merged.bam") into merged_bam_count
    set id, file("${id}.merged.bam") into  merged_bam_sub

    script:
    """
    samtools merge ${id}.merged.bam ${ff.collect { it }.join(' ')}
    """
}

process mergeStats {
 
    input:
    set type, file(merge) from merged_bam_count

    output:
    set type, file("${type}.flagstat.txt") into stats

    script:
    """
    samtools flagstat ${merge} > ${type}.flagstat.txt
    """
}


process subsampleMerged {

    input:
    file(stat) from stats.collect()
    set type, file(bam) from merged_bam_sub
  
    output:
    set type, file("${type}.subset.bam") into subbams

    script:
    """
    $baseDir/bin/downsample.sh 
    """
}

/*process callMACS2 {

 input: 
 set type, file(sbam) from subbams

 output:
 file("${type}.peaks.narrowPeak") into narrowPeaks

script:
"""
macs2 callpeak -t ${sbam} -g ${param.genomesize} -n ${type} ${param.macs_call}
"""
}*/



/*process makeMasterPeaks {
}*/

/*process count_reads_in_master {
}*/

/*process count_reads_in_narrow {
}*/

/*process calculate_fpkm_in_master {

}*/

/*process calculate_fpkm_in_narrow {
}*/

/*process deseq2 {
}*/

