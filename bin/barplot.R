args = commandArgs(trailingOnly=TRUE)


## first arg is type

type=args[1]

## rest files

files = args[2:length(args)]

### load libraries
library(wesanderson)
library(Cairo)


### define colors
mapcolors = wes_palette(n=4,name="GrandBudapest")[c(1:2,4,3)]
dupcolors = wes_palette(n=2,name="Royal1")
colors = c(c("lightblue","darkblue"),c("lightgreen","darkgreen"),c("pink","darkred"),c("yellow","darkorange"))[1:length(files)]

##
indata = matrix(NA,nrow=14,ncol=length(files))
colnames(indata)=files

for (f in files){
  indata[,f]=read.table(f,sep="\t",row.names = 1)[,1]
  rownames(indata) = rownames(read.table(f,sep="\t",row.names = 1))
} 


### split matrix in three
mapStat = indata[1:4,]
colnames(mapStat)=sapply(strsplit(colnames(mapStat),"\\."),"[[",1)
rownames(mapStat)=gsub(":","",rownames(mapStat))

lStst = indata[6:10,]
rownames(lStst)=gsub(":","",rownames(lStst))

dups = indata[c(2,5),]
rownames(dups)=gsub(":","",rownames(dups))
colnames(dups)=sapply(strsplit(colnames(dups),"\\."),"[[",1)

### create barplots
pdf(paste(type,"fig1.pdf",sep="_"))
barplot(mapStat,beside=T,col=mapcolors,main=type,las=2)
legend("topright",rownames(mapStat),fill=mapcolors)
dev.off()

CairoPNG(paste(type,"fig1.png",sep="_"))
barplot(mapStat,beside=T,col=mapcolors,main=type,las=2)
legend("topright",rownames(mapStat),fill=mapcolors)
dev.off()


pdf(paste(type,"fig2.pdf",sep="_"))
barplot(dups,beside=T,col=dupcolors,main=type,las=2)
legend("topright",rownames(dups),fill=dupcolors)
dev.off()

CairoPNG(paste(type,"fig2.png",sep="_"))
barplot(dups,beside=T,col=dupcolors,main=type,las=2)
legend("topright",rownames(dups),fill=dupcolors)
dev.off()


pdf(paste(type,"fig3.pdf",sep="_"))
par(mar=c(8,3,1,1))
barplot(t(lStst),beside=T,col=colors,las=2,main=type)
legend("topleft",colnames(lStst),fill=colors)
dev.off()

CairoPNG(paste(type,"fig3.png",sep="_"))
par(mar=c(8,3,1,1))
barplot(t(lStst),beside=T,col=colors,las=2,main=type)
legend("topleft",colnames(lStst),fill=colors)
dev.off()


