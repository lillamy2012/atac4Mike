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
markedBam.into { markedBam; markedBam2 }

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
filter_markedBam.into { filter_markedBam; filter_markedBam2 }

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

allbams=uniq_filteredBam.mix(filter_markedBam2).mix(markedBam2)
allbams.into { allbams; allbams2; allbams3; allbams4 }

process flagstat {

    input: 
    set name, type, file(bam) from allbams

    output:
    set name, file("${name}.${type}.flagstat.txt") into flagstats
    set name, file("${name}.${type}.stats.txt") into stats
    set name, type, file("${name}.${type}.depth.txt") into depths
    set name, type, file("${name}.${type}.numbers.txt") into numbers

    script:
    """
    samtools flagstat ${bam} > ${name}.${type}.flagstat.txt
    samtools stats -c 0,2000,1 ${bam} > ${name}.${type}.stats.txt
    for stat in "raw total sequences:" "reads mapped:" "reads mapped and paired:" "reads properly paired:" "reads duplicated:" "average length:" "maximum length:" "average quality:" "insert size average:" "insert size standard deviation:" "inward oriented pairs:" "outward oriented pairs:" "pairs with other orientation:" "pairs on different chromosomes:" ; do
     grep "^SN" ${name}.${type}.stats.txt | grep "\$stat" | awk -F \$'\t' '{ print \$2,"\t",\$3}' >> ${name}.${type}.numbers.txt
     done   
 samtools depth -a ${bam} > ${name}.${type}.depth.txt
    """

}

//allbams2.into { allbams2; allbams3 }

process samtools_mult {

    input:
    set name, type, file(bam) from allbams2
  
    output:
    file( "${name}.${type}.mult_stat.tab" )
    file( "${name}.${type}.uniq_stat.tab" )

"""
samtools view -f 66  ${bam} |grep "XS:i"|grep -E "@|AS:i" | wc -l >> ${name}.${type}.mult_stat.tab
samtools view -f 66  ${bam} |grep -v "XS:i"|grep -E "@|AS:i" | wc -l >> ${name}.${type}.uniq_stat.tab
"""
}

process getIS {

  input: 
  set name, type, file(bam) from allbams3

  output:
  set type, file("${name}.${type}.IS.tab") into is

  script:
  """
samtools view -f 66 ${bam} | awk '{print \$9}' | python  $baseDir/bin/my_counter.py > ${name}.${type}.IS.tab
  """
}

process coverTable {

  input:
  set name, type, file(depth) from depths

  output:
  file("${name}.${type}.depth.table.tab")

  script:
  """
  awk '{print \$3}' ${depth} | python $baseDir/bin/my_counter.py > ${name}.${type}.depth.table.tab
  """

}


//orderedNum=numbers.groupTuple(by:1).map{ it -> [ it[1] , it[2]]}.collect().subscribe{ println it }
orderedNum = numbers.map{it -> tuple( it[1] , it[2])}.groupTuple(by:0)
orderedIs =  is.groupTuple(by:0)

process R_barplots {
tag "type: $type"

   input: 
   set type, file(num) from orderedNum

   output:
   file("${type}_fig1.png")
   file("${type}_fig2.png")
   file("${type}_fig3.png")
   file("${type}_fig1.pdf")
   file("${type}_fig2.pdf")
   file("${type}_fig3.pdf")

   script:
   """
   Rscript --vanilla $baseDir/bin/barplot.R ${type}  ${num.collect { it }.join(' ')}  
   """
   
}

process R_insertplots {
tag "type: $type"

   input:
   set type, file(is) from orderedIs

   output:
   file("${type}_is_fig1.pdf")
   file("${type}_is_fig1.png")
   file("${type}_is_fig2.pdf")
   file("${type}_is_fig2.png")

   script:
   """
   Rscript --vanilla $baseDir/bin/insertplot.R ${type} ${is.collect { it }.join(' ')}
   """

}





