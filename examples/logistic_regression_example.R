######################################
# 12.10.2016
# Logistic regression on ChIP-seq data
# BISC 481
######################################

## Install packages
install.packages("caret")
install.packages("e1071")
install.packages("ROCR")


## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)
workingPath <- "/Users/test/Downloads/BISC481/CTCF/"

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
