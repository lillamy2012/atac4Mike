#!/usr/bin/env Rscript

library(Cairo)
library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library("gsubfn")

## functions
############
runVisPlot=function(dds,prename){
  dds <- estimateSizeFactors(dds)
  rld <- rlog(dds, blind=FALSE)
  pdf(paste(prename,"deseq_fig1.pdf",sep="_"))
  plot(plotPCA(rld,c("type")))
  dev.off()
  
  CairoPNG(paste(prename,"deseq_fig1.png",sep="_"))
  plot(plotPCA(rld,c("type")))
  dev.off()
  
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
  pdf(paste(prename,"deseq_fig2.pdf",sep="_"),onefile=FALSE)
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=colData)
  dev.off()
  
  CairoPNG(paste(prename,"deseq_fig2.png",sep="_"))
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
  cluster_cols=FALSE, annotation_col=colData)
  dev.off()
  
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <-  colnames(dds)
#print(rld$type)
#print(colnames(dds))
  #colnames(sampleDistMatrix) <- colnames(dds)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf(paste(prename,"deseq_fig3.pdf",sep="_"),onefile=FALSE)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  CairoPNG(paste(prename,"deseq_fig3.png",sep="_"))
  pheatmap(sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors)
  dev.off()
  
}

plotScatter = function(dds,contrast,th,fc){
  #print(log10(fc))
print(class(fc))
  res = results(dds,contrast)
  ncounts = counts(dds,normalized=T)
  groups= colData(dds)[contrast[1]]
  print(groups)
  xcounts = rowMeans(ncounts[,which(groups[,1]==contrast[2]),drop=F])
  ycounts = rowMeans(ncounts[,which(groups[,1]==contrast[3]),drop=F])
  c_max = max(max(xcounts),max(ycounts))
  c_max=(round(log10(c_max))+1)
  plot(log10(xcounts),log10(ycounts),pch=".",col="grey",axes=FALSE,ylab=contrast[3],xlab=contrast[2])
  axis(1,at=seq(0,c_max,1),labels = 10^seq(0,c_max,1))
  axis(2,at=seq(0,c_max,1),labels = 10^seq(0,c_max,1))
  abline(0,1,col="black",lwd=3)
  sig=which(res$padj<th)
  points(log10(xcounts)[sig],log10(ycounts)[sig],col="red",pch=".",cex=2)
  abline(-log10(fc),1,col="black",lty=2,lwd=2)
  abline(log10(fc),1,col="black",lty=2,lwd=2)
}

panel.smoothScatter = function(x,y){
  smoothScatter(x,y,add=TRUE)
  abline(0,1,col="red")
} 

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

plotRep = function(dds){
  dds <- estimateSizeFactors(dds)
  ncounts = counts(dds,normalized=T)
  pairs(log10(ncounts+1),lower.panel = panel.smoothScatter,upper.panel = panel.cor)
  
}

resTable = function(res,dds,th){
  resOrder = res[order(res$padj<th,abs(res$log2FoldChange),decreasing = T),]
  ncounts = counts(dds,normalized=T)
  groups= gsub(" ","",unlist(strsplit(strsplit(strsplit(res@elementMetadata[5,2],":")[[1]][2],c("type"))[[1]][2],"vs")))
  print(groups)
  xcounts = rowMeans(ncounts[,which(colData(dds)$type==groups[1]),drop=F])
  ycounts = rowMeans(ncounts[,which(colData(dds)$type==groups[2]),drop=F])
  tot = cbind(as.data.frame(resOrder),xcounts[rownames(resOrder)],ycounts[rownames(resOrder)])
  colnames(tot)[7:8]=paste(groups,"norm.counts",sep="_")
  tot
  }

############################
args = commandArgs(trailingOnly=TRUE)
print(args)

if (length(args)<3){
  stop("Provide tab file with files to read, p-value and fold change")
} else if (length(args)>3) {
  print ("too many args,  only 3 first will be processed")
} else
  design=args[1]
  th = args[2]
  fc = as.numeric(args[3])

setup=read.table(design,sep=",",comment.char="")

print(setup[,2])
files = paste(setup[,2],"uniq_filtered.bam_master_counts.tab",sep="_")
f_names = setup[,1]

print(f_names)
annoF = "master_anno.csv"

anno = read.table(annoF,sep=";",header=T,comment.char="#",quote = "") 

n=0
s=-1
x=TRUE
while(x==TRUE){
print (x)
  n = n+1
  s = s+1
   loi = scan(annoF,nlines=n,skip=s, what = character(),sep=";")
   x = grepl("^#",loi)[1]
print(x)
 }

print(s)	
comment = read.table(annoF,nrow=s,sep=";",comment.char="")


print(comment)
print(head(anno))
rownames(anno)=paste0("peak",1:nrow(anno))
peaks = grep("narrowPeak",colnames(anno),value=T)
print("ok")
ofInt = c(c("seqnames", "start", "end", "width", "annotation", "geneId" , "transcriptId" , "distanceToTSS"),c(peaks))
print(ofInt)
anno = anno[,ofInt]

tmp=list()
for (i in 1:length(files)){
  print(files[i])
  tmp[[files[i]]]=read.table(files[i],row.names = 1,header=F,sep="\t")
  colnames(tmp[[files[i]]])=f_names[i]
}

countTab=do.call("cbind",tmp)
countTab = countTab[-grep("__",rownames(countTab)),]

colData = data.frame(type=as.factor(setup[,3]))
rownames(colData)=colnames(countTab)
contr=c("type",levels(colData$type)[1],levels(colData$type)[2])

print("dds")
dds <- DESeqDataSetFromMatrix(countData = countTab,
                              colData = colData,
                              design = ~ type)

save(dds,file="dds.Rdata")

runVisPlot(dds,"master_peaks")

pdf("master_peaks_deseq_fig4.pdf")
plotRep(dds)
dev.off()

CairoPNG("master_peaks_deseq_fig4.png")
plotRep(dds)
dev.off()


dds = DESeq(dds)
res = results(dds,contrast=contr)

pdf("master_peaks_deseq_fig5.pdf")
plotMA(res,main="ma plot")
dev.off()

CairoPNG("master_peaks_deseq_fig5.png")
plotMA(res,main="ma plot")
dev.off()

pdf("master_peaks_deseq_fig6.pdf")
plotScatter(dds,contrast = contr,th=th,fc=fc)
dev.off()
  
CairoPNG("master_peaks_deseq_fig6.png")
plotScatter(dds,contrast = contr, th=th, fc=fc)
dev.off()

resOrder=resTable(res,dds,th=th)
resOrder = cbind(resOrder,anno[rownames(resOrder),])
write.table(comment,file="deseq_results.csv",sep=";",quote = FALSE,row.names=FALSE,col.names=FALSE)
write.table(resOrder,file="deseq_results.csv",sep=";",quote = FALSE,append=TRUE)














