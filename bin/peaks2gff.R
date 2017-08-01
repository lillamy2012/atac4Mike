#!/usr/bin/env Rscript

library(ChIPseeker)
library(GenomicRanges)


writeGFF = function(GRanges,file){
df = data.frame(seqnames(GRanges),rep("rtracklayer",length(GRanges)),"exon",start(GRanges),end(GRanges),rep(".",length(GRanges)),rep("+",length(GRanges)),rep(".",length(GRanges)),paste0(paste0("gene_id \"peak",1:length(GRanges)),"\";"))
    write.table(df,sep="\t",file=paste0(file,".gff"),row.names = FALSE,col.names = FALSE,quote = FALSE)
}


#########################

files = grep("narrowPeak",dir(),value=TRUE)


peaks=list()
for (i in files){
    if(length(readLines(i))>0){
        peaks[[i]] = readPeakFile(i,header=F)
        peaks[[i]] = peaks[[i]][which(!seqnames(peaks[[i]])%in%c("mitochondria","chloroplast","Mt","Pt"))]
        writeGFF(peaks[[i]],file=i) 
   }
}
