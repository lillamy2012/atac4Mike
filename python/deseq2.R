args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
  stop("Provide tab filewith files to read")
} else if (length(args)>1) {
  print ("more than one input file, only first will be processed")
} else
  input=args[1]





## first get some info from R_file: type, path, 
## read in 
## assign conditions
## plots on all data (including input): PCA, pheatmap, heatmap, replicate plots
## remove input 
## plots, PCA, pheatmap, heatmap
## 'anova' like 
## pair-wise contrasts 
## hit list 
## plotMA, scatter plots, 

library(DESeq2)
library("pheatmap")
library("RColorBrewer")
library(Cairo)


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
  rownames(sampleDistMatrix) <-  rld$type
  colnames(sampleDistMatrix) <- NULL
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

plotScatter = function(dds,contrast,th=0.01){
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
  abline(-log10(2),1,col="black",lty=2,lwd=2)
  abline(log10(2),1,col="black",lty=2,lwd=2)
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

resTable = function(res,dds,th=0.01){
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



###################
###




#input = "R_test.tab"
name = basename(input)
name = sub("R_","",name)
lname = strsplit(name,"\\.")[[1]]
exclude = c("bam","tab","txt")
new = paste(lname[!lname%in%exclude],collapse = "_")
setup=read.table(input, comment.char="")

print(setup)
files  = setup[,1]

nn = gsub("_deseq","",new)
annoF = paste(nn,"anno.csv",sep="_")
#"uniq_filtered_MACS2_withInput_anno.csv"
print(annoF)
anno = read.table(annoF,sep=";",header=T,comment.char="",quote = "")
rownames(anno)=paste0("peak",1:nrow(anno))

ofInt = c("seqnames", "start", "end", "width", "annotation", "geneId" , "transcriptId" , "distanceToTSS")

anno = anno[,ofInt]

tmp=list()
for (i in files){
  print(i)
  tmp[[i]]=read.table(i,row.names = 1,header=F,sep="\t")
  colnames(tmp[[i]])=i
}
print("ok")
countTab=do.call("cbind",tmp)
countTab = countTab[-grep("__",rownames(countTab)),]

print("2")
setup$samp = ifelse(setup$V3=="sample",1,0)
colData = data.frame(type=as.factor(setup$V2))
rownames(colData)=colnames(countTab)

print("dds")
dds <- DESeqDataSetFromMatrix(countData = countTab,
                              colData = colData,
                              design = ~ type)



## plots based on all samples
runVisPlot(dds,paste(new,"all",sep="_"))
pdf(paste(paste(new,"all",sep="_"),"deseq_fig4.pdf",sep="_"))
plotRep(dds)
dev.off()

CairoPNG(paste(paste(new,"all",sep="_"),"deseq_fig4.png",sep="_"))
plotRep(dds)
dev.off()



## only samples


colData=data.frame(type=colData$type[which(setup$samp==1),drop=TRUE])
countTab=countTab[,setup$samp==1]
rownames(colData)=colnames(countTab)
dds <- DESeqDataSetFromMatrix(countData = countTab,
                              colData = colData,
                              design = ~ type)


runVisPlot(dds,paste(new,"samples",sep="_"))
pdf(paste(paste(new,"samples",sep="_"),"deseq_fig4.pdf",sep="_"))
plotRep(dds)
dev.off()

CairoPNG(paste(paste(new,"samples",sep="_"),"deseq_fig4.png",sep="_"))
plotRep(dds)
dev.off()


## "anova" like plus all pairwise contrasts

if (nlevels(colData$type)>2){
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- nbinomLRT(dds, reduced = ~ 1)
  res_anova <- as.data.frame(results(dds,name="Intercept"))
  res_ord_anova = res_anova[order(res_anova$padj),]
  for (n in levels(colData$type)){
    cc = as.matrix(rowMeans(counts(dds,normalized=TRUE)[,which(colData$type==n),drop=F]))
    print(head(cc))
    colnames(cc)=paste(n,"norm.counts",sep="_")
    res_ord_anova=cbind(res_ord_anova,cc[rownames(res_ord_anova),])
    colnames(res_ord_anova)[ncol(res_ord_anova)]=paste(n,"norm.counts",sep="_")
  }
  res_ord_anova=res_ord_anova[,!colnames(res_ord_anova)%in%c("log2FoldChange", "lfcSE","stat")]
  
  res_ord_anova = cbind(res_ord_anova,anno[rownames(res_ord_anova),])
  write.csv(res_ord_anova,file = paste(new,"anovaLike.csv",sep="_"))
  }
  
  

dds = DESeq(dds)
save(dds, file=paste(new,"dds.Rdata",sep="_"))


## contrasts
## hhh
pos = levels(colData$type)
contrL = list()
k=0
for (i in 1:length(pos)){
  for (j in (i+1):length(pos)){
    if (j<=(length(pos)) & j!=i){
      k=k+1
    contrL[[k]]=c("type", pos[i], pos[j] )
    }
  }
}

resL = list()
for (i in 1:length(contrL)){
  resL[[i]]= results(dds,contrast = contrL[[i]])
}

## plots
## ma
for (i in 1:length(resL)){
  pdf(paste0(paste(new,paste(contrL[[i]],collapse = "_"),sep="_"),"ma.pdf"))
  plotMA(resL[[i]],main=paste(contrL[[i]],collapse = "_"))
  dev.off()
  
  CairoPNG(paste0(paste(new,paste(contrL[[i]],collapse = "_"),sep="_"),"ma.png"))
  plotMA(resL[[i]],main=paste(contrL[[i]],collapse = "_"))
  dev.off()
  
  pdf(paste0(paste(new,paste(contrL[[i]],collapse = "_"),sep="_"),"scatter.pdf"))
  plotScatter(dds,contrast = contrL[[i]])
  dev.off()
  
  CairoPNG(paste0(paste(new,paste(contrL[[i]],collapse = "_"),sep="_"),"scatter.png"))
  plotScatter(dds,contrast = contrL[[i]])
  dev.off()
  
}

## table

for (i in 1:length(resL)){
resOrder = resTable(resL[[i]],dds)
resOrder = cbind(resOrder,anno[rownames(resOrder),])
write.csv(resOrder,file=paste(paste(new,paste(contrL[[i]],collapse = "_"),sep="_"),"csv",sep="."))
}



# annotate


## table on fold change

