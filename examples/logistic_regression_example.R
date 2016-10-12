######################################
# 12.10.2016
# Logistic regression on ChIP-seq data
# BISC 481
######################################

## Install packages
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("AnnotationHub")
biocLite("GenomicRanges")
biocLite("BSgenome.Mmusculus.UCSC.mm10")


## Initialization
library(DNAshapeR)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)


## Data retreive
seqLength <- 30 #500
sampleSize <- 2000 #42045
workingPath <- "/Users/lester/GitHub/BISC481/CTCF/"

# Bound (ChIP-seq)
ah <- AnnotationHub()
#unique(ah$dataprovider)
#query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH28451"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "bound_30.fa"))
          
# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
  tmpChrName <- sample(chrName, 1)
  tmpChrLength <- chrLength[tmpChrName]
  tmpGRanges <- GRanges( seqnames = tmpChrName, strand = "+",
                         IRanges(start = sample(1:tmpChrLength,1), width = seqLength) )
  if( length(findOverlaps(ctcfPeaks, tmpGRanges)) == 0 ){
    randomGr <- c( randomGr, tmpGRanges)
    print(length(randomGr))
  }else{
    print(paste(length(randomGr), "overlap"))
  }
}
randomGr

# Overlap checking
findOverlaps(ctcfPeaks, randomGr)

# Fasta file generation
getFasta(randomGr, Mmusculus, width = seqLength, 
         filename = paste0(workingPath, "unbound_30.fa"))


## Generate data from classifcation (assign 1 to bound and 0 to non-bound)
# bound
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape(paste0(workingPath,"ctcf.fa"))


## Encode feature vectors
featureType <- c("1-mer")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance)

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc
