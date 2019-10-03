######################################
# 10.02.2019
# AnnotationHub example
# BISC 481
######################################

## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# AnnotationHub
BiocManager::install("AnnotationHub")
# rtracklayer (might be required)
BiocManager::install("rtracklayer")
# Reference genome
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

## Initialization
library(AnnotationHub)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

## Information retreival
ah <- AnnotationHub()
ah
unique(ah$dataprovider)
unique(ah$speices)

ah <- subset(ah, species == "Homo sapiens")
ah
ah <- query(ah, c("H3K4me3", "Gm12878", "Roadmap"))
ah

## Genome sequence retreival
getFasta(ah[["AH29706"]], Hsapiens, width = 150, filename = "tmp.fa")
