#!/usr/bin/env nextflow


params.bams          = "/lustre/scratch/users/elin.axelsson/demult/test.bam"
params.output        = "bam/"
params.min_length    = 5
params.overlap       = 1
params.read_length   = 0
params.A             = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"
params.a             = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
params.index         = "/lustre/scratch/projects/berger_common/TAIR10/Bowtie2Index/Arabidopsis_thaliana.TAIR10.31.dna_rm.genome"

index=file(params.index)


bams = Channel 
       .fromPath(params.bams)
       .map { file -> tuple(file.baseName, file) }

bams.into { bam_read1; bam_read2 }

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

fqs = fq1.cross(fq2).map{ it -> tuple( it[0][0], it[0][1],it[1][1] )}

process cutadapt {
tag "name: $name"

   input:
   set name, file(fq1), file(fq2) from fqs
 
   output:
   set name, file("${name}_cutadapt_R1_fastq"), file("${name}_cutadapt_R2_fastq") into trimmed

   script:
   """
   cutadapt --minimum-length ${params.min_length} --overlap ${params.overlap} -a ${params.a}  -A ${params.A} -o ${name}_cutadapt_R1_fastq -p ${name}_cutadapt_R2_fastq ${fq1} ${fq2}  
   """
}


process bowtie2 {
tag "name: $name"

   input:
   set name, file(cut1), file(cut2) from trimmed

   output:
   set name, file("${name}.aligned_cut.sam") into sams 

   script:
   """
   bowtie2 --end-to-end -x ${index} -1 ${cut1} -2 ${cut2} -S ${name}.aligned_cut.sam
   """
}

process sort {
tag "name: $name"

   input:
   set name, file(sam) from sams

   output:
   set name, file("${name}.aligned_cut_sorted.bam") into aligned

   script:
   """
   export TMPDIR=./${name}_tmp
# make sure the directory exists
   mkdir -p \$TMPDIR
   samtools view -b ${sam} | samtools sort -T \$TMPDIR - > ${name}.aligned_cut_sorted.bam
   """
}

