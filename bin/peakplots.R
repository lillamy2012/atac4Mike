#!/usr/bin/env Rscript


library(GenomicRanges)
library(wesanderson)
library(Cairo)
library(ChIPseeker)
library(UpSetR)
library(gtools)

w_col=wes_palette("Moonrise3",5)[c(1,5,3)]
tot_merged=read.csv("master_anno.csv",sep = ";", stringsAsFactors=F,comment="#")
psets = grep("narrowPeak",colnames(tot_merged))
pdf("peakOV.fig1.pdf")
upset(tot_merged,main.bar.color=w_col,sets.bar.color=w_col[1:2])
dev.off()

CairoPNG("peakOV.fig1.png")
upset(tot_merged,main.bar.color=w_col,sets.bar.color=w_col[1:2])
dev.off()



## length distr, of individual and of master
tot_merged$peakL = tot_merged$end-tot_merged$start
tot_merged$levels = factor(tot_merged$seqnames, levels=unique(mixedsort(tot_merged$seqnames)))
pdf("peakOV.fig2.pdf")
boxplot(tot_merged$peakL~tot_merged$levels,las=2,varwidth=T,col=wes_palette("Moonrise3",2)[2])
dev.off()

CairoPNG("peakOV.fig2.png")
boxplot(tot_merged$peakL~tot_merged$levels,las=2,varwidth=T,col=wes_palette("Moonrise3",2)[2])
dev.off()



tot_merged$comb=ifelse(rowSums(tot_merged[,psets])==2,2,tot_merged[,psets[2]])


### add individ
psets = grep("narrowPeak",colnames(tot_merged),value=T)
peaks=list()
for(i in psets){
  peaks[[i]] = readPeakFile(i,header=F)
  }
pdf("peakOV.fig3.pdf")
par(mfrow=c(1,2))
boxplot(tot_merged$peakL~tot_merged$comb,col=w_col,names=c(sapply(strsplit(psets,"_"),"[[",1),"both"),main="Master")
boxplot(lapply(peaks,function(x) y=width(x)),las=1,col=w_col[1:2],names=sapply(strsplit(psets,"_"),"[[",1),main="narrowPeak")
dev.off()

CairoPNG("peakOV.fig3.png")
par(mfrow=c(1,2))
boxplot(tot_merged$peakL~tot_merged$comb,col=w_col,names=c(sapply(strsplit(psets,"_"),"[[",1),"both"),main="Master")
boxplot(lapply(peaks,function(x) y=width(x)),las=1,col=w_col[1:2],names=sapply(strsplit(psets,"_"),"[[",1),main="narrowPeak")
dev.off()

