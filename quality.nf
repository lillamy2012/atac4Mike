#!/usr/bin/env nextflow

params.bams          = "bam/*.bam"
params.quality       = 10

bamset = Channel
        .fromPath(params.bams)
        .map { file -> tuple(file.baseName, file) }

process markDuplicates {
 
   input: 
   set name, file(bam) from bamset

   output:
   set name, file("${name}_mark_bam") into markedBam

   script:
   """
    export TMPDIR=\$(pwd)
    java -Djava.io.tmpdir=\$TMPDIR -jar \$EBROOTPICARD/picard.jar MarkDuplicates I=${bam}  O=${name}_mark_bam M=filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=false
    """
}

markedBam=markedBam.map{ it -> tuple(it[0],"marked",it[1])}

process qualityFilter {

  input: 
  set name, type, file(mark_bam) from markedBam

  output:
  set name, file("${name}_filter_mark_bam") into filter_markedBam

  script:
  """
  samtools view -bSq ${params.quality} ${mark_bam} > ${name}_filter_mark_bam
  """
}

filter_markedBam=filter_markedBam.map{ it -> tuple(it[0],"filter_mark",it[1])}


process removeDuplicates {

   input:
   set name,type, file(filter_mark_bam) from filter_markedBam

   output:
   set name, file("${name}_uniq_filtered.bam") into uniq_filteredBam
   
   script:
   """
    export TMPDIR=\$(pwd)
    java -Djava.io.tmpdir=\$TMPDIR -jar \$EBROOTPICARD/picard.jar MarkDuplicates I=${filter_mark_bam}  O=${name}_uniq_filtered.bam M=filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true
   """
}

uniq_filteredBam=uniq_filteredBam.map{ it -> tuple(it[0],"uniq_filtered",it[1])}

/*
/*process flagstat {

    input: 
    set name, file(bam) from allbams

    output:
    set name, file("${name}.
}

process stats {
}

process depth {
}
*/
