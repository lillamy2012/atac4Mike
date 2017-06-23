args = commandArgs(trailingOnly=TRUE)
#source("../R/functions.R")

if (length(args)==0){
  stop("Provide tab filewith files to read")
} else if (length(args)>1) {
  print ("more than one input file, only first will be processed")
} else
  input=args[1]


### in argument with files to use
#input = "R_marked_duplicates.filtered.MACS2_noInput.tab"
files=read.table(input)
files  = files[,1]
#files = checkInput(input=files,file=colors)



###
name = basename(input)
name = sub("R_","",name)
lname = strsplit(name,"\\.")[[1]]
exclude = c("bam","tab","txt")
new = paste(lname[!lname%in%exclude],collapse = "_")


library(geneplotter)
library(ChIPseeker)
library(GenomicRanges)
library(VennDiagram)
library(wesanderson)
library(rtracklayer)
library(Cairo)
library(ChIPseeker)
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
df = data.frame(seqnames(GRanges),rep("rtracklayer",length(GRanges)),"exon",start(GRanges),end(GRanges),rep(".",length(GRanges)),rep("+",length(GRanges)),rep(".",length(GRanges)),paste0(paste0("gene_id \"peak",1:length(GRanges)),"\";")
    write.table(df,sep="\t",file=paste0(file,".gff"),row.names = FALSE,col.names = FALSE,quote = FALSE)
}

## summarize overlaps want to be up to 3 groups
#summarizeOverlaps = function(allRegions,groups){
#  all_groups = apply(allRegions,1,paste,collapse=".")
#  nG = length(unique(groups))
#  print(nG)
#  gs = list()
#  for ( i in 1:nG){
#      print(which(groups==i)) ## this is why it breaks,
#    gs[[i]] = apply(allRegions[, which(groups==i),drop=FALSE],1,paste,collapse=".")
#  }
#  groups=gs
#  return(list(all_groups,groups))
#}

######################################


### define colors
#colors = wes_palette(n=4, name="Zissou")[c(2:1,3:4)]
#cc =  wes_palette(n=5, name="Royal2")[c(3:5)]

### read in peaks, overlap, barplot 

## read in function ## write also gff for counting
peaks=list()
for (i in files){
    if(length(readLines(i))>0){
        peaks[[i]] = readPeakFile(i,header=F)
        peaks[[i]] = peaks[[i]][which(!seqnames(peaks[[i]])%in%c("mitochondria","chloroplast"))]
        writeGFF(peaks[[i]],file=i)
    }
}

tot = reduce(Reduce(c, peaks)) # combine all peaks called
tot_merge = mergeClosePeaks(tot)
print(head(tot_merge))
overlaps = createOverlapMat(peaks,tot_merge)


## annotate

tot_anno = annotatePeak(tot_merge, tssRegion=c(-900, 900),TxDb=txdb)

print(head(as.data.frame(tot_anno)))


tot_df_anno = cbind(rownames(overlaps),as.data.frame(tot_anno),overlaps)


## create gff
df = data.frame(seqnames(tot_merge),rep("rtracklayer",length(tot_merge)),"exon",start(tot_merge),end(tot_merge),rep(".",length(tot_merge)),rep("+",length(tot_merge)),rep(".",length(tot_merge)),paste0(paste0("gene_id \"peak",1:length(tot_merge)),"\";"))

write.table(df,sep="\t",file=paste0(new,".gff"),row.names = FALSE,col.names = FALSE,quote = FALSE)
#write.table(as.data.frame(tot_anno), file = paste0(new,"_anno",".txt"),row.names = FALSE,quote=FALSE)
write.table(tot_df_anno, file = paste0(new,"_anno",".csv"),row.names = FALSE,quote=FALSE,sep=";")

