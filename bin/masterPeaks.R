#!/usr/bin/env Rscript

library(ChIPseeker)
library(GenomicRanges)



args = commandArgs(trailingOnly=TRUE)

### get which txdb to use:
tx=args[2]

library(tx,character.only = TRUE)
txdb=get(tx)

#library(TxDb.Athaliana.BioMart.plantsmart28)
#txdb = TxDb.Athaliana.BioMart.plantsmart28

bp_distance = as.numeric(args[1])


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
  rownames(allRegions)=paste(seqnames(tot),start(tot),end(tot),sep="__")
  for (i in 1:nsamp){
    ovs[[i]] = findOverlaps(peaks[[i]],tot)
    allRegions[unique(subjectHits(ovs[[i]])),i] = 1 
  }
  tot_base = as.data.frame(tot)
  allRegions = data.frame(cbind(tot_base,allRegions))
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
print("anno")
tot_anno = annotatePeak(tot_merge, tssRegion=c(-bp_distance, bp_distance),TxDb=txdb)
tmp_anno = as.data.frame(tot_anno)
rownames(tmp_anno) = paste(tmp_anno$seqnames, tmp_anno$start, tmp_anno$end,sep="__")

toRm = c("seqnames","start","end","width","strand")
toRm = intersect(colnames(overlaps),colnames(tmp_anno))

tmp_anno = tmp_anno[,!colnames(tmp_anno)%in%toRm]
tot_df_anno = merge(tmp_anno,overlaps,by="row.names",all=TRUE,sort=FALSE)

#if (sum(tot_df_anno$seqnames.x==tot_df_anno$seqnames.y,na.rm=T)!=sum(!is.na(tot_df_anno$seqnames.x))){
#stop("error with seqnames")
#}

#tot_df_anno$seqnames.x=tot_df_anno$seqnames.y
rownames(tot_df_anno) = tot_df_anno$Row.names
toRemove=c("Row.names")
tot_df_anno=tot_df_anno[,!colnames(tot_df_anno)%in%toRemove]
#colnames(tot_df_anno)=sub("seqnames.x","seqnames",colnames(tot_df_anno))

## create gff
df = data.frame(seqnames(tot_merge),rep("rtracklayer",length(tot_merge)),"exon",start(tot_merge),end(tot_merge),rep(".",length(tot_merge)),rep("+",length(tot_merge)),rep(".",length(tot_merge)),paste0(paste0("gene_id \"peak",1:length(tot_merge)),"\";"))

write.table(df,sep="\t",file="master.gff",row.names = FALSE,col.names = FALSE,quote = FALSE)

write.table(tot_df_anno, file ="master_anno.csv",row.names = FALSE,quote=FALSE,sep=";")


