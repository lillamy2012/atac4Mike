#!/usr/bin/env Rscript

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb = TxDb.Athaliana.BioMart.plantsmart28


#####################
## functions
######################

## merge peaks 
mergeClosePeaks=function(grange, length=50){
  extended = grange
  start(extended)=start(extended)-length/2
  end(extended)=end(extended)+length/2
  extended_r = reduce(extended)
  start(extended_r)=start(extended_r)+length/2
  end(extended_r)=end(extended_r)-length/2
  return(extended_r)
}

## create overlap mat 
createOverlapMat = function(peaks,tot){
  nsamp=length(peaks)
  ovs = list()
  allRegions= matrix(0,ncol=nsamp,nrow=length(tot))
  colnames(allRegions)=names(peaks)
  rownames(allRegions)=paste(seqnames(tot),start(tot),end(tot),sep="_")
  for (i in 1:nsamp){
    ovs[[i]] = findOverlaps(peaks[[i]],tot)
    allRegions[unique(subjectHits(ovs[[i]])),i] = 1 
  }
  return(allRegions)
}
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
    }
}

tot = reduce(Reduce(c, peaks)) # combine all peaks called
tot_merge = mergeClosePeaks(tot)
print(head(tot_merge))
overlaps = createOverlapMat(peaks,tot_merge)

## annotate

tot_anno = annotatePeak(tot_merge, tssRegion=c(-900, 900),TxDb=txdb)

tot_df_anno = cbind(rownames(overlaps),as.data.frame(tot_anno),overlaps)


## create gff
df = data.frame(seqnames(tot_merge),rep("rtracklayer",length(tot_merge)),"exon",start(tot_merge),end(tot_merge),rep(".",length(tot_merge)),rep("+",length(tot_merge)),rep(".",length(tot_merge)),paste0(paste0("gene_id \"peak",1:length(tot_merge)),"\";"))

write.table(df,sep="\t",file="master.gff",row.names = FALSE,col.names = FALSE,quote = FALSE)

write.table(tot_df_anno, file ="master_anno.csv",row.names = FALSE,quote=FALSE,sep=";")


