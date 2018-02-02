#!/usr/bin/env nextflow

/*********************
* Parameters
*********************/

params.design		= 'exp_1.tab'
params.macs_call	= '-B -q 1e-5 -f BAMPE'
params.genome		= 'Athaliana'
params.bams		= "bam_1/*.bam" 
params.quality		= 10  
params.output		= "results_1/"
params.anno_distance	= 900
// params.txdb		= "TxDb.Athaliana.BioMart.plantsmart28" 
// "mm10ref_seq_txdb.sqlite" "TxDb.Mmusculus.UCSC.mm10.knownGene"
params.peak_merge_dist	= 50
params.deseq_p		= 0.01
params.deseq_fc		= 2
params.bw_binsize	= 10

//int val = new BigDecimal(stringValue).intValue();
/**********************
* set genome size, txdb
***********************/

if(params.genome == 'Athaliana'){
	params.genomesize = 1.2e8
	params.txdb = "TxDb.Athaliana.BioMart.plantsmart28"
}
	int tonorm = new BigDecimal(params.genomesize).intValue()


/*************************************************
*************************************************/

log.info "ATAC-SEQ PIPE  N F  ~  version 0.1"
log.info "====================================="
log.info "====================================="
log.info "bam files		: ${params.bams}"
log.info "design		: ${params.design}"
log.info "quality threshold	: ${params.quality}"
log.info "output		: ${params.output}"
log.info "****************************************"
log.info "macs call		: ${params.macs_call}"
log.info "genome		: ${params.genome}"
log.info "genome size		: ${params.genomesize}"
log.info "****************************************"
log.info "peak merge distance	: ${params.peak_merge_dist}"
log.info "annotation distance	: ${params.anno_distance}"
log.info "tx db			: ${params.txdb}"
log.info "****************************************"
log.info "deseq2 p-value	: ${params.deseq_p}"
log.info "deseq2 fc		: ${params.deseq_fc}"
log.info "****************************************"
log.info "bigwig bin size	: ${params.bw_binsize}"


/*************************************************
*************************************************/


// set up start channels, from bam and design file
bamset = Channel
        .fromPath(params.bams)
        .map { file -> tuple(file.baseName, file) }

design = Channel
        .fromPath(params.design)
        .splitCsv()
        .map { row -> [ id:row[1],cond:row[2] ] }

// for deseq script
designFile = file(params.design)


// generate uniq filtered files, needed for replicate merge and htseq counts, For the narrowPeaks we don't want the name in the output

process generateFiles {
tag "bam : $name"  
publishDir "${params.output}/bam", mode: 'copy'
  
    input:
    set name, file(bam) from bamset

    output:
    set name, file("${name}_uniq_filtered.bam") into uniq_filtered
    set name, file("${name}_uniq_filtered.bam") into final_bamset_count
    file("${name}_uniq_filtered.bam") into final_bamset_count_n

    script:
    """
    export TMPDIR=\$(pwd)
    samtools view -bSq ${params.quality} ${bam} > ${name}_tmp.bam
    java -Djava.io.tmpdir=\$TMPDIR -jar \$EBROOTPICARD/picard.jar MarkDuplicates I=${name}_tmp.bam  O= ${name}_uniq_filtered.bam  M=${name}_filtered.dup_metrics.txt AS=true REMOVE_DUPLICATES=true 
   TMP_DIR=\$TMPDIR 
    """
}

uniq_filtered.into { uniq_filtered ; bw_bam }

// bw of individual bam files

process bwIndBam {
publishDir "${params.output}/ds_bam", mode: 'copy'
tag "name: $name"

    input:
    set name, file(bam) from bw_bam

    output:
    set name, file("${name}.bw")

    script:
    """
    export TMPDIR=\$(pwd)
    samtools index ${bam}
    bamCoverage -b ${bam} -o ${name}.bw --normalizeTo1x ${tonorm} --binSize=${params.bw_binsize} 
    """	
}


// set up channel with design info for each bam - needed for merging of replicates

comb = uniq_filtered.map { it -> [ id:it[0], file:it[1] ] }.phase(design) {it -> it.id}
     .map { it -> tuple( it.cond[1] , it.file[0])  }
     .groupTuple(by:0 )



// merged filtered, uniq bam replicates. Needed for both flagstat and subsampling

process mergeReplicates {
tag "id: $id"

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
tag "type: $type" 
    input:
    set type, file(merge) from merged_bam_count

    output:
    set type, file("${type}.flagstat.txt") into stats

    script:
    """
    samtools flagstat ${merge} > ${type}.flagstat.txt
    """
}

// sub sample to smallest condition. This will be used to call peaks

process subsampleMerged {
tag "type: $type"
publishDir "${params.output}/ds_bam", mode: 'copy'
    
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

subbams.into { subbams; subbams_bw }

process sortBam {
tag "bam: $bam"

   input:
   set type, file(bam) from subbams_bw

   output:
   set type, file("${type}.subset_sort.bam") into sortbams_bw

   script:
   """
   samtools sort -o ${type}.subset_sort.bam ${bam}
   """

}

process bam2bw {
tag "bam: $bam"
publishDir "${params.output}/ds_bam", mode: 'copy'

   input:
   set type, file(bam) from sortbams_bw

   output:
   set type, file("${type}.subset.bw") 

   script:
   """
   export TMPDIR=\$(pwd)
   samtools index ${bam}
   bamCoverage -b ${bam} -o ${type}.subset.bw --normalizeTo1x ${tonorm} --binSize=${params.bw_binsize}
   """
}

// call macs2, narrowPeaks are needed for masterpeak and counting in individual peaks

process callMACS2 {
tag "type: $type"
publishDir "${params.output}/macs2", mode: 'copy'

    input: 
    set type, file(sbam) from subbams

    output:
    file("${type}_peaks.narrowPeak") into narrowPeaks_to_merge
    set type, file("${type}_peaks.narrowPeak") into narrowPeaks_to_gff
    file("${type}_peaks.narrowPeak") into narrowPeaks_to_plot

    script:
    """
    export TMPDIR=\$(pwd)
    macs2 callpeak -t ${sbam} -g ${params.genomesize} -n ${type} ${params.macs_call}
    """
}

// combine macs peaks to master peaks set 

process makeMasterPeaks {
publishDir "${params.output}/gff", mode: 'copy'   

   input:
   file(narrowPeaks) from narrowPeaks_to_merge.collect() 

   output:
   file("master.gff") into master
   file("master.gff") into master2 
   file("master_anno.csv") into anno
   file("master_anno.csv") into anno2

   script:
   """
   $baseDir/bin/masterPeaks.R ${params.anno_distance} ${params.txdb} ${params.peak_merge_dist}
   """

}
process plotPeaks {
publishDir "${params.output}/macs2", mode: 'copy'

  input:
  file(narrowPeaks) from narrowPeaks_to_plot.collect()
  file("master_anno.csv") from anno2

  output:
  file("peakOV.fig1.pdf")
  file("peakOV.fig2.pdf")
  file("peakOV.fig3.pdf")
  file("peakOV.fig1.png")
  file("peakOV.fig2.png")
  file("peakOV.fig3.png")


  script:
  """
  $baseDir/bin/peakplots.R 
  """
}
// gff file from each narrow peak

process makeFileGff {
tag "peak: $type"
publishDir "${params.output}/gff", mode: 'copy' 
   
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
tag "name: $name"
publishDir "${params.output}/counts", mode: 'copy'  
 
   input:
   file("master.gff") from master
   set name, file(bams) from final_bamset_count

   output:
   file ("${bams}_master_counts.tab") into master_counts
   file ("${bams}_master_fpkm.tab") into master_fpkm

   """
   htseq-count -f bam -s no ${bams} master.gff > ${bams}_master_counts.tab
   $baseDir/bin/calcFPKM.sh > ${bams}_master_fpkm.tab
   """
}


fileGFF.into {fileGFF; GFFS}

// all combination of gff file and bam file
xch=final_bamset_count_n.combine(fileGFF)

// count reads in narrowpeaks and calc FPKM
process count_reads_in_narrow {
publishDir "${params.output}/counts", mode: 'copy'
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
publishDir "${params.output}/deseq", mode: 'copy'

   input: 
   file("master_anno.csv") from anno
   file(counts) from master_counts.collect()
   file(designFile)

   output:
   file("dds.Rdata")   
   file("master_peaks_deseq_fig*")
   file("deseq_results.csv") into deseq

   script:
   """
   $baseDir/bin/deseq2.R ${designFile} ${params.deseq_p} ${params.deseq_fc}
   """
}

process combineResults {
publishDir "${params.output}/fpkm", mode: 'copy'

  input: 
  file(narrow) from GFFS.collect()
  file(fpkm) from narrow_fpkm.collect() 
  file("master.gff") from master2
  file(m_fpkm) from master_fpkm.collect()
  file("deseq_results.csv") from deseq
 
  output:
  file("*_summary.tab")
  file("deseq_total_results.csv")   

  script:
  """
  $baseDir/bin/combineFPKM.sh 
  """


}


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}



