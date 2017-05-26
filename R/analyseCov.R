library(wesanderson)
library(MASS)
library(Cairo)

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

### define colors
#colors = c(wes_palette(n=4, name="Zissou")[c(2:1,3:4)],"darkblue","orange") ## need to change
mapcolors = wes_palette(n=4,name="GrandBudapest")[c(1:2,4,3)]
dupcolors = wes_palette(n=2,name="Royal1")

### in argument with files to use
#input = "../tables/R_uniq_filtered.bam.numbers.tab"
files=read.table(input)
fileF=read.table("../files.tab",comment.char="")
colors=defineColors(fileF)
files  = files[,1]
files = checkInput(input=files,file=colors)
print(colors)
print(files)
colors = colors$col


###
name = basename(input)
name = sub("R_","",name)
lname = strsplit(name,"\\.")[[1]]
exclude = c("bam","tab","txt")
new = paste(lname[!lname%in%exclude],collapse = "_")


cov=list()
for (i in files){
  cov[[i]] = read.table(i)
  cov[[i]]$V1 = cov[[i]]$V1 ## 
  rownames(cov[[i]])=cov[[i]]$V1
  cov[[i]]=cov[[i]][,2,drop=FALSE]
}

## to matrix
nrows = max(sapply(cov,function(x)  max(as.numeric(rownames(x)))))
ncols = length(cov)
coverage = matrix(0,ncol=ncols,nrow=nrows+1)
colnames(coverage)=names(cov)
rownames(coverage)=0:nrows

for (i in names(cov)){
  coverage[rownames(cov[[i]]),i]=unlist(cov[[i]])
}

## expand and est disp etc
expand = list()

maxU = ifelse(nrow(coverage)>1000,1000,nrow(coverage))

for (i in files){
  expand[[i]] = rep(as.numeric(rownames(coverage)[1:maxU]),coverage[1:maxU,i])
}

summ = lapply(expand, summary)
vars = lapply(expand,var)
x = lapply(expand, function(x) fitdistr(sample(x,100000),"negative binomial"))
disps = sapply(x,function(y) y[[1]]["size"])

## create plots
#matplot(coverage,type="l",lty=1,lwd=3,col=colors)

maxS = ifelse(nrow(coverage)>200,200,nrow(coverage))
small = coverage[1:maxS,]
if(nrow(coverage)>200){
    extra = colSums(coverage[201:nrow(coverage),])
    small = rbind(small,extra)
}
pdf(paste(new,"fig1.pdf",sep="_"))
matplot((small),type="l",lty=1,lwd=3,col=colors,xlab="depth",ylab="",main=new)
legend("topright",colnames(small),fill=colors)
dev.off()
CairoPNG(paste(new,"fig1.png",sep="_"))
matplot((small),type="l",lty=1,lwd=3,col=colors,xlab="depth",ylab="",main=new)
legend("topright",colnames(small),fill=colors)
dev.off()


s_cum = apply(small,2,cumsum)
ydisp= seq(round(min(disps)),round(max(disps)),length.out = 11)
xcov = seq(round(min(s_cum)),round(max(s_cum)),length.out = 11)
lm_m = lm(ydisp~xcov)
new_scale = (disps-lm_m$coefficients[1])/lm_m$coefficients[2]

pdf(paste(new,"fig2.pdf",sep="_"))
matplot(s_cum,type="l",lty=1,lwd=3,col=colors,xlab="depth",ylab="",main=new)
legend("bottomright",colnames(small),fill=colors)
axis(4,at = xcov,labels = ydisp)
points(y=new_scale,x=sapply(x,function(y) y[[1]]["mu"]),col=colors,pch=19,cex=2)
dev.off()

CairoPNG(paste(new,"fig2.png",sep="_"))
matplot(s_cum,type="l",lty=1,lwd=3,col=colors,xlab="depth",ylab="",main=new)
legend("bottomright",colnames(small),fill=colors)
axis(4,at = xcov,labels = ydisp)
points(y=new_scale,x=sapply(x,function(y) y[[1]]["mu"]),col=colors,pch=19,cex=2)
dev.off()



