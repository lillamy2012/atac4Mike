
args = commandArgs(trailingOnly=TRUE)


source("../R/functions.R")
if (length(args)==0){
  stop("Provide tab filewith files to read")
} else if (length(args)>1) {
  print ("more than one input file, only first will be processed")
} else
  input=args[1]


### load libraries
library(wesanderson)
library(Cairo)

### define colors
#colors = c(wes_palette(n=4, name="Zissou")[c(2:1,3:4)],"darkblue","orange") ## need to change
mapcolors = wes_palette(n=4,name="GrandBudapest")[c(1:2,4,3)]
dupcolors = wes_palette(n=2,name="Royal1")

### in argument with files to use
#input = "../tables/R_uniq_filtered.bam.numbers.tab"
files=read.table(input)
print("ok")
fileF=read.table("../files.tab",comment.char="")
print(fileF)
colors=defineColors(fileF)
files  = files[,1]
files = checkInput(input=files,file=colors)
print(colors)
print(files)
colors = colors$col

#setwd("../stats")

###
name = basename(input)
name = sub("R_","",name)
lname = strsplit(name,"\\.")[[1]]
exclude = c("bam","tab","txt")
new = paste(lname[!lname%in%exclude],collapse = "_")

### read in data to matrix
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
pdf(paste(new,"fig1.pdf",sep="_"))
barplot(mapStat,beside=T,col=mapcolors,main=new,las=2)
legend("topright",rownames(mapStat),fill=mapcolors)
dev.off()

CairoPNG(paste(new,"fig1.png",sep="_"))
barplot(mapStat,beside=T,col=mapcolors,main=new,las=2)
legend("topright",rownames(mapStat),fill=mapcolors)
dev.off()


pdf(paste(new,"fig2.pdf",sep="_"))
barplot(dups,beside=T,col=dupcolors,main=new,las=2)
legend("topright",rownames(dups),fill=dupcolors)
dev.off()

CairoPNG(paste(new,"fig2.png",sep="_"))
barplot(dups,beside=T,col=dupcolors,main=new,las=2)
legend("topright",rownames(dups),fill=dupcolors)
dev.off()


pdf(paste(new,"fig3.pdf",sep="_"))
par(mar=c(8,3,1,1))
barplot(t(lStst),beside=T,col=colors,las=2,main=new)
legend("topleft",colnames(lStst),fill=colors)
dev.off()

CairoPNG(paste(new,"fig3.png",sep="_"))
par(mar=c(8,3,1,1))
barplot(t(lStst),beside=T,col=colors,las=2,main=new)
legend("topleft",colnames(lStst),fill=colors)
dev.off()


