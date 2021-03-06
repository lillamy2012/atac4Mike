#!/usr/bin/env nextflow

// inparameters
params.bams          = "/lustre/scratch/users/elin.axelsson/Mike-25-09-2017/demult/*.bam"
params.output        = "bam/"
params.stats         = "align_logs/"
params.min_length    = 5
params.overlap       = 1
params.read_length   = 50
params.A             = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
params.a             = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

params.index         = "/lustre/scratch/projects/berger_common/backup_berger_common/Bowtie2Index/genome"

//TAIR10/Bowtie2Index/Arabidopsis_thaliana.TAIR10.31.dna_rm.genome"

//params.index         = "/lustre/scratch/projects/berger_common/mm10/mm10_index_4Bowtie/mm10"


// index
index=file(params.index)

// start channel
bams = Channel 
       .fromPath(params.bams)
       .map { file -> tuple(file.baseName, file) }

bams.into { bam_read1; bam_read2 }

// bam2fastq 1st and 2nd read in parallel , reads trimmed if read_length not 0

process bam_to_fastq1 {
tag "name: $name"

   input:
   set name, file(bam) from bam_read1

   output:
   set name, file("${name}_R1_.fastq") into fq1

   script:
   if( params.read_length == 0 )
       """
       samtools view -f 0x40 -b ${bam} | samtools bam2fq - > "${name}_R1_.fastq" 
       """

   else if( params.read_length > 0 )
       """
       samtools view -f 0x40 -b ${bam} | samtools bam2fq - | fastx_trimmer -l $params.read_length -Q33 -o "${name}_R1_.fastq"
       """

   else 
       error "Invalid read_length argument"
}

process bam_to_fastq2 {
tag "name: $name"

   input:
   set name, file(bam) from bam_read2

   output:
   set name, file("${name}_R2_.fastq") into fq2

   script:
   if( params.read_length == 0 )
       """
       samtools view -f 0x80 -b ${bam} | samtools bam2fq - > "${name}_R2_.fastq"
       """

   else if( params.read_length > 0 )
       """
       samtools view -f 0x80 -b ${bam} | samtools bam2fq - | fastx_trimmer -l ${params.read_length} -Q33 -o "${name}_R2_.fastq"
       """
   
   else 
       error "Invalid read_length argument"
}

// downstream pipeline needs the first and second read at once. 

fqs = fq1.cross(fq2).map{ it -> tuple( it[0][0], it[0][1],it[1][1] )}


process cutadapt {
tag "name: $name"

   input:
   set name, file(fq1), file(fq2) from fqs
 
   output:
   set name, file("${name}_cutadapt_R1_fastq"), file("${name}_cutadapt_R2_fastq") into trimmed
   file("${name}.cutadapt.log") into calogs 

   script:
   """
   cutadapt --minimum-length ${params.min_length} --overlap ${params.overlap} -a ${params.a}  -A ${params.A} -o ${name}_cutadapt_R1_fastq -p ${name}_cutadapt_R2_fastq ${fq1} ${fq2} > ${name}.cutadapt.log 
   """
}

process bowtie2 {
tag "name: $name"

   input:
   set name, file(cut1), file(cut2) from trimmed

   output:
   set name, file("${name}.aligned_cut.sam") into sams  
   file("${name}.bowtie2.log") into btlogs

   script:
   """
   bowtie2 --end-to-end -X 2000 -x ${index} -1 ${cut1} -2 ${cut2} -S ${name}.aligned_cut.sam 2> ${name}.bowtie2.log 

   """
}

process sort {
tag "name: $name"
publishDir params.output, mode: 'copy'

   input:
   set name, file(sam) from sams

   output:
   set name, file("${name}.aligned_cut_sorted.bam")

   script:
   """
   export TMPDIR=./${name}_tmp
# make sure the directory exists
   mkdir -p \$TMPDIR
   samtools view -b ${sam} | samtools sort -T \$TMPDIR - > ${name}.aligned_cut_sorted.bam
   """
}

logs=calogs.mix(btlogs)


process savelog {
tag "name: $name"
publishDir params.stats, mode: 'copy'

 input: 
 file(mylog) from logs

 output:
 file "${mylog.baseName}.logfile" 

  script:
  """
  fname=\$(basename ${mylog})
  fbname=\${fname%.*}
  cp ${mylog} \${fbname}.logfile
  """
}
 

println "Project : $workflow.projectDir/$params.stats"
