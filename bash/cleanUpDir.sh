#!/bin/bash


cd /lustre/scratch/users/elin.axelsson/atac4Mike/test
## pdf/png to plots
mkdir -p plots
mv *.pdf plots
mv *.png plots

## csv to outfiles
mkdir -p outfiles
mv *.csv outfiles

##gff to giff
mkdir -p gff
mv *.gff gff

##Rdata to Rdata
mkdir -p Rdata
mv *.Rdata Rdata

##tpm 
mkdir -p tpm
mv *tpm.tab tpm

##counts
mkdir -p counts
mv *counts.tab counts

## fpkm
mv *fpkm.tab tpm
 
## tables
mkdir -p tables
mv R_* tables

##files used for files
mkdir -p used_stats
mv *.tab used_stats

## remove python out
rm y_data.txt

## mv to better name
mv results stats_results

