######################################
# 01.10.2016
# Housing example
# BISC 481
######################################

library(MASS)

# Load housing data from text file
df <- read.table("/Users/lester/BISC577/examples/house_price.txt")

# Define features and response value
x <- data.matrix( df[1:5] )
y <- data.matrix( df[6] )

# Calculate theta by Normal Equation
theta <- ginv(t(x) %*% x) %*% t(x) %*% y
theta

# Predict the response values for testing data
df2 <- read.table("/Users/lester/BISC577/examples/house_price_testing.txt")
x_testing <- data.matrix( df2[1:5] )
y_testing <- data.matrix( df2[6] )

y_prediction <- x_testing %*% theta
y_prediction

# Assess the model with R-squared
rs <- 1 - ( sum( y_testing - y_prediction ) ^ 2 ) / 
  sum( (y_testing - mean( y_testing ) ) ^ 2 ) 
rs
