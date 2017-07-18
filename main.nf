#!/usr/bin/env nextflow

params.design        = 'exp.tab'
params.macs_call     = '-B -q 0.01 -f BAMPE'
params.genomesize    = '1.2e8'
params.bams          = "bam/*small.bam" 
params.quality       = 10  
params.output        = "results/"
params.anno_distance = 900


// set up start channels, from bam and design file
bamset = Channel
        .fromPath(params.bams)
        .map { file -> tuple(file.baseName, file) }

design = Channel
        .fromPath(params.design)
        .splitCsv()
        .map { row -> [ id:row[1],cond:row[2] ] }

designFile = file(params.design)

// bamset will be used again later
bamset.into { bamset; bamset_count }


// generate uniq filtered files
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


// set up channel with design info for each bam - needed for merging of replicates
comb = uniq_filtered.map { it -> [ id:it[0], file:it[1] ] }.phase(design) {it -> it.id}
     .map { it -> tuple( it.cond[1] , it.file[0])  }
     .groupTuple(by: 0)


// merged filtered, uniq bam replicates
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

// run statistics on merged files (needed for down sampling) 
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

// sub sample to smallest condition
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

// call macs2
process callMACS2 {

    input: 
    set type, file(sbam) from subbams

    output:
    file("${type}_peaks.narrowPeak") into narrowPeaks_to_merge
    set type, file("${type}_peaks.narrowPeak") into narrowPeaks_to_gff

    script:
    """
    export TMPDIR=\$(pwd)
    macs2 callpeak -t ${sbam} -g ${params.genomesize} -n ${type} ${params.macs_call}
    """
}

// combine macs peaks to master peaks set 
process makeMasterPeaks {
    
   input:
   file(narrowPeaks) from narrowPeaks_to_merge.collect() 

   output:
   file("master.gff") into master
   file("master_anno.csv") into anno

   script:
   """
   $baseDir/bin/masterPeaks.R
   """

}

// gff file from each narrow peak
process makeFileGff {

   input:
   set type, file(narrowPeaks) from narrowPeaks_to_gff

   output:
   file("${type}_peaks.narrowPeak.gff") into fileGFF

   script:
   """
   $baseDir/bin/peaks2gff.R
   """
}  

// count reads in master peaks and calculate FPKM
process count_reads_in_master {
    
   publishDir 'data/counts', mode: 'copy'  
 
   input:
   file("master.gff") from master
   set name, file(bams) from bamset_count

   output:
   file ("${bams}_master_counts.tab") into master_counts
   file ("${bams}_master_fpkm.tab") into master_fpkm

   """
   htseq-count -f bam -s no ${bams} master.gff > ${bams}_master_counts.tab
   $baseDir/bin/calcFPKM.sh > ${bams}_master_fpkm.tab
   """
}

//bam files again 
bamset_count_n=Channel.fromPath(params.bams)

// all combination of gff file and bam file
xch=bamset_count_n.combine(fileGFF)

// count reads in narrowpeaks and calc FPKM
process count_reads_in_narrow {
tag "gff: $gff"

   input:
   set file(bams), file(gff) from xch

   output:
   file("${bams}_${gff}_counts.tab") into narrow_counts
   file("${bams}_${gff}_fpkm.tab") into narrow_fpkm  

   """
   htseq-count -f bam -s no ${bams} ${gff} > ${bams}_${gff}_counts.tab
   $baseDir/bin/calcFPKM.sh > ${bams}_${gff}_fpkm.tab
   """
}


process deseq2 {

   input: 
   file("master_anno.csv") from anno
   file(counts) from master_counts.collect()
   file(designFile)

   output:
   file("dds.Rdata")   

   script:
   """
   $baseDir/bin/deseq2.R ${designFile}
   """
}

/*process combineResults {
}*/

