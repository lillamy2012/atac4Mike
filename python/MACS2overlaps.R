args = commandArgs(trailingOnly=TRUE)
source("../R/functions.R")

if (length(args)==0){
  stop("Provide tab filewith files to read")
} else if (length(args)>1) {
  print ("more than one input file, only first will be processed")
} else
  input=args[1]


### in argument with files to use
#input = "R_marked_duplicates.filtered.MACS2_noInput.tab"
files=read.table(input)
fileF=read.table("../files.tab",comment.char="")
colors=defineColors(fileF)
c_order=unique(colors$V4[-nrow(colors)])
print(c_order)
colors = droplevels(subset(colors,colors$V3=="sample"))
files  = files[,1]
files = checkInput(input=files,file=colors)
fileF = colors[,-7]
colors = colors$col



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
mergeClosePeaks=function(grange, length=100){
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

## summarize overlaps want to be up to 3 groups
summarizeOverlaps = function(allRegions,groups){
  all_groups = apply(allRegions,1,paste,collapse=".")
  nG = length(unique(groups))
  print(nG)
  gs = list()
  for ( i in 1:nG){
      print(which(groups==i)) ## this is why it breaks, 
    gs[[i]] = apply(allRegions[, which(groups==i),drop=FALSE],1,paste,collapse=".")
  }
  groups=gs
  return(list(all_groups,groups))
}

## approach 1 
plotApp1 = function(sumList,cols){
  all = table(sumList[[1]])
  common = all["1.1.1.1"]
  un1 = all["1.1.0.0"]+ all["1.1.1.0"]+  all["1.1.0.1"]
  un2 = all["0.0.1.1"]+ all["0.1.1.1"]+  all["1.0.1.1"]
  vp = draw.pairwise.venn(un1+common,un2+common,common,fill=cols)  #c("red", "green"))
  grid.draw(vp)
}

## approach 2
plotApp2 = function(sumList,cols){
  in1 = grep("1",sumList[[2]][[1]])
  in2 = grep("1",sumList[[2]][[2]])
  common = length(intersect(in1,in2))
  un1 = length(setdiff(in1,in2))
  un2 = length(setdiff(in2,in1))
  vp = draw.pairwise.venn(un1+common,un2+common,common,fill=cols)
  grid.draw(vp)
}

######################################


### define colors
#colors = wes_palette(n=4, name="Zissou")[c(2:1,3:4)]
cc =  wes_palette(n=5, name="Royal2")[c(3:5)]

### read in peaks, overlap, barplot 

## read in function 
peaks=list()
for (i in files){
    if(length(readLines(i))>0){
        peaks[[i]] = readPeakFile(i,header=F)
        peaks[[i]] = peaks[[i]][which(!seqnames(peaks[[i]])%in%c("mitochondria","chloroplast"))]
    }
}

tot = reduce(Reduce(c, peaks)) # combine all peaks called
tot_merge = mergeClosePeaks(tot)
print(head(tot_merge))
overlaps = createOverlapMat(peaks,tot_merge)


gr=findRepCondOrder(files,fileF)
s_over = summarizeOverlaps(overlaps,gr$group)  

print(head(overlaps))
## annotate

tot_anno = annotatePeak(tot_merge, tssRegion=c(-900, 900),TxDb=txdb)

print(head(as.data.frame(tot_anno)))


tot_df_anno = cbind(rownames(overlaps),as.data.frame(tot_anno),overlaps)


## create gff
df = data.frame(seqnames(tot_merge),rep("rtracklayer",length(tot_merge)),"exon",start(tot_merge),end(tot_merge),rep(".",length(tot_merge)),rep("+",length(tot_merge)),rep(".",length(tot_merge)),paste0(paste0("gene_id \"peak",1:length(tot_merge)),"\";"))

write.table(df,sep="\t",file=paste0(new,".gff"),row.names = FALSE,col.names = FALSE,quote = FALSE)
#write.table(as.data.frame(tot_anno), file = paste0(new,"_anno",".txt"),row.names = FALSE,quote=FALSE)
write.table(tot_df_anno, file = paste0(new,"_anno",".csv"),row.names = FALSE,quote=FALSE,sep=";")

for(i in 1:max(gr$group)){
  pdf(paste(paste0(paste(new,"fig1",sep="_"),i),"pdf",sep="."))
  plotReplicateVenn(s_over,v_colors=colors,gr=gr$group,group=i,cat=paste0("rep",gr$rep))
  dev.off()
  CairoPNG(paste(paste0(paste(new,"fig1",sep="_"),i),"png",sep="."))
  plotReplicateVenn(s_over,v_colors=colors,gr=gr$group,group=i,cat=paste0("rep",gr$rep))
  dev.off()

  }

## approcach 1 : common if in all, unique if both rep (1 or zero other)
pdf(paste(new,"fig3.pdf",sep="_"))
grid.newpage()
nrR = table(gr$group)
try(plotAppVenn(defineA(overlaps,gr$group,nrR),cc,c_order))
dev.off()
CairoPNG(paste(new,"fig3.png",sep="_"))
grid.newpage()
nrR = table(gr$group)
try(plotAppVenn(defineA(overlaps,gr$group,nrR),cc,c_order))
dev.off()
## approcach 2 : common if in both , unique if only one

pdf(paste(new,"fig4.pdf",sep="_"))
grid.newpage() 
try(plotAppVenn(defineA(overlaps,gr$group,rep(1,length(nrR))),cc,c_order))
dev.off()
CairoPNG(paste(new,"fig4.png",sep="_"))
grid.newpage()
try(plotAppVenn(defineA(overlaps,gr$group,rep(1,length(nrR))),cc,c_order))
dev.off()


rle_w=lapply(peaks, function(x) width(x))
pdf(paste(new,"fig5.pdf",sep="_"))
#boxplot(lapply(rle_w,log10),lwd=2,col=colors,ylab="log10(length)",xlab="",main="length")
multidensity(lapply(rle_w,log10),lwd=3,col=colors,xlab="log10(length)",main="")
dev.off()
CairoPNG(paste(new,"fig5.png",sep="_"))
multidensity(lapply(rle_w,log10),lwd=3,col=colors,xlab="log10(length)",main="")
dev.off()


#export.gff(tot_merge,"peakregions.gff")
