######################################
# 01.10.2016
# Logistic Regression (Classification) example
# BISC 481
######################################

## Install packages
# Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# DNAshapeR
biocLite("DNAshapeR")
# Caret
install.packages("caret")
install.packages("e1071")
install.packages("ROCR")

## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)


## Prepare data
# Predict DNA shapes
fn_fasta <- "/Users/lester/BISC481/examples/SELEX_ScrWT.fa"
pred <- getShape(fn_fasta)
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType)

# Transform binding affinity data to bound-unbound data
fn_exp <- "/Users/lester/BISC577/examples/SELEX_ScrWT.txt"
exp_data <- read.table(fn_exp)
exp_data$V2 <- ifelse(exp_data$V2 > 0.464, "Y", "N")

## Perform logistic regression
# Prepare data
df <- data.frame(isBound = exp_data$V2, featureVector)

# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)

# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

# Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs_new )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance)

# Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc
