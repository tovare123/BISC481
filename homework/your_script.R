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

Question #5
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
workingPath <- "/Users/lizziliz/Downloads/BISC481-master-4/CTCF/"

# Extract sample sequences
fn <- system.file("extdata", "CGRsample.fa", package = "DNAshapeR")

# Predict DNA shapes
pred <- getShape(fn)

# Generate ensemble plots
## Generate one at a time
plotShape(pred$MGW)
plotShape(pred$ProT)
plotShape(pred$Roll)
plotShape(pred$HelT)
