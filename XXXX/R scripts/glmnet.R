library(glmnet)
library(tidyverse)
library(readr)
library(pROC)

set.seed(16)

# import csv files
omics2  <- read_csv("data/omics2.csv")
omics3  <- read_csv("data/omics3.csv")
omics1  <- read_csv("data/omics1.csv")
label   <- read_csv("data/label.csv")
test    <- read_csv("data/test_sample.csv")

# aggregate features and outcome into one dataframe. remove the sample ID from subsequent dataframes.
dat <- cbind(omics1, omics2[, -1], omics3[, -1], label[, -1])

# split data into train set and test set
dat_train <- dat[!(dat$Sample %in% test$Sample), ]
dat_test  <- dat[dat$Sample %in% test$Sample, ]

# obtain number of columns of aggregated dataset
n_col <- ncol(dat)

# convert dfs to right format for ridge regression
x_train <- as.matrix(dat_train[, 2:(n_col - 1)])
x_test  <- as.matrix(dat_test[, 2:(n_col - 1)])
y_train <- as.matrix(dat_train[, n_col])
y_test  <- as.matrix(dat_test[, n_col])

# run ridge regression on training set with cross-validation
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "binomial")

# run ridge regression on training set with lambda.1se
ridge <- glmnet(x_train, y_train, alpha = 0, family = "binomial", lambda = cv_ridge$lambda.min)

# make predictions on test set
test_ridge <- ridge %>% predict(newx = x_test)

# transform predictions into probabilities, then into classes
probabilities <- 1 / (1 + exp(-test_ridge))
pred_class <- ifelse(probabilities > 0.5, 1, 0)

# aggregate observed and predicted class in a new dataframe
result_ridge <- cbind(test, y_test, pred_class)
colnames(result_ridge) <- c("Sample", "Observed", "Predicted")

# write predictions file
write_csv(result_ridge, "results/result_ridge.csv")
