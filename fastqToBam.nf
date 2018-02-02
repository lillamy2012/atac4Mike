#!/usr/bin/env nextflow

params.in = 'fastq/*end{1,2}.fq'
params.output = 'results'



Channel 
	.fromFilePairs( params.in, size: -1)
	.ifEmpty { error "Cannot find any reads matching: ${params.in}" }
	.flatten()
	.collate(3)
	.set { read_files }

//read_files.subscribe{println it}


process fastqSAM {
	module='Picard/1.141'
	memory='62 GB'

	input:
	set name, file(read1), file(read2)  from read_files 

	output:
	set name, file("${name}.sam") into sam 

	script:
	"""
	export TMPDIR=\$(pwd)	
        java -Djava.io.tmpdir=\$TMPDIR -jar \$EBROOTPICARD/picard.jar FastqToSam F1=${read1} F2=${read2}  O= ${name}.sam  SM=for_alignment	
        """	

}

process sam2bam {
	module='SAMtools/1.3-goolf-1.4.10'
	publishDir "${params.output}/bam", mode: 'copy'	

	input:
	set name, file(sam) from sam

	output:
	file("${name}.bam")

	script:
	"""
	samtools view -b -o ${name}.bam  ${sam}
        """
  }


