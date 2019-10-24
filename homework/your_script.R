#Question 4
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DNAshapeR")

install.packages("caret")

install.packages("DNAshapeR")


######################################
# 10.02.2019
# Multiple Linear Regression (MLR) example
# BISC 481
######################################

## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret")

## Initialization
library(DNAshapeR)
library(caret)
workingPath <- "/Users/lizziliz/Downloads/BISC481-master-4/gcPBM/"

## Predict DNA shapes
fn_fasta <- paste0(workingPath, "Max.txt.fa")
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
head(featureVector)

## Build MLR model by using Caret
# Data preparation
fn_exp <- paste0(workingPath, "Max.txt")
exp_data <- read.table(fn_exp)
df <- data.frame(affinity=exp_data$V2, featureVector)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction without L2-regularized
model <- train (affinity~ ., data = df, trControl=trainControl, 
                method = "lm", preProcess=NULL)
summary(model)

# Prediction with L2-regularize
model2 <- train(affinity~., data = df, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2
result <- model2$results$Rsquared[1]
result

##R squared DATA
## MER+SEQUENCE
#max= 0.8646
#mad=0.8628
#myc=0.8549
##MER
#max=0.7855
#mad=0.7754
#myc=0.7785
#plotting
--
  plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm"),
axis.text.x = element_text(colour="black", size=12),
axis.text.y = element_text(colour="black", size=12),
axis.title.x = element_text(colour="black", size=12),
axis.title.y = element_text(colour="black", size=12),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
axis.line = element_line(colour = "black"),
axis.text = element_text(colour ="black"),
axis.ticks = element_line(colour = "black")
)

## Data preparation
SequenceShape <- c(0.8646,0.8628,0.8549)
Sequence <- c(0.7855,0.7754,0.7785)


## Ploting
ggplot() +
  geom_point(aes(x = Sequence, y = SequenceShape), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) 

#Question 5
## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# DNAshapeR
BiocManager::install("DNAshapeR")
# Caret
install.packages("caret")

## Install packages
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("AnnotationHub")
BiocManager::install("GenomicRanges")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

## Initialization
library(DNAshapeR)
library(AnnotationHub)
library(BSgenome.Mmusculus.UCSC.mm10)

## Data retreive
seqLength <- 50 #500
sampleSize <- 2000 #42045
workingPath <- "/Users/lizziliz/Downloads/BISC481-master-6/CTCF/"


# Bound (ChIP-seq)
ah <- AnnotationHub()
unique(ah$dataprovider)
query(ah, "H3K4me3")
ctcfPeaks <- ah[["AH46100"]]
seqlevelsStyle(ctcfPeaks) <- "UCSC"
getFasta( GR = sample(ctcfPeaks, sampleSize), BSgenome = Mmusculus, width = seqLength, 
          filename = paste0(workingPath, "bound_500.fa"))

# Unbound (random regions w/o overlapping)
chrName <- names(Mmusculus)[1:22]
chrLength <- seqlengths(Mmusculus)[1:22]
randomGr <- GRanges()

while ( length(randomGr) < sampleSize )  {
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
         filename = paste0(workingPath, "unbound_500.fa"))


# Extract sample sequences 
fn <- file.path("/Users", "lizziliz", "Downloads", "BISC481-master-4", "CTCF", "unbound_500.fa")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)

# Extract sample sequences --> this one worked better for me
fn <- file.path("/Users", "lizziliz", "Downloads", "BISC481-master-4", "CTCF", "bound_500.fa")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)

##Question 6
######################################
# 10.02.2019
# Logistic regression on ChIP-seq data
# BISC 481
######################################
##no change in workingpath; can simiply enter in data 


## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)
library(Biostrings)
workingPath <- "/Users/lizziliz/Downloads/BISC481-master-4/CTCF/"

## Generate data for the classifcation (assign Y to bound and N to non-bound)
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
pred <- getShape(paste0(workingPath, "ctcf.fa"))


## Encode feature vectors
featureType <- c("1-mer", "1-shape")
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
