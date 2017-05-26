library(rtracklayer)
library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb = TxDb.Athaliana.BioMart.plantsmart28

peaks = readGFF("../full/uniq_filtered_MACS2_withInput.gff")
peakGR=as(peaks,"GRanges")
anno = annotatePeak(peakGR, tssRegion=c(-900, 900), 
             TxDb=txdb)