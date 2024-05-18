library(glmnet)
library(tidyverse)
library(readr)
library(ggplot2)
library(reshape2)

set.seed(900)

# import data .csv files
Meth <- read_csv("data/Meth_bal.csv")
RPPA <- read_csv("data/RPPA_bal.csv")
mRNA <- read_csv("data/mRNA_bal.csv")
tumor <- read_csv("data/outcome_bal.csv")
test_ID <- read_csv("data/test_sample.csv")

# concatenate the 3 omics + class into one dataframe (remove the sample columns that repeat)
dat_balanced <- cbind(Meth, RPPA[, -1], mRNA[, -1], tumor[, -1])

# separate the data into train/test sets
dat_train <- dat_balanced[!(dat_balanced$Sample %in% test_ID$Sample), ]
dat_test <- dat_balanced[(dat_balanced$Sample %in% test_ID$Sample), ]

# obtain number of columns of aggregated dataset
n_col <- as.numeric(ncol(dat_balanced))

# I wasn't sure whether the CV splits are done randomly or in an ordered way, so I shuffled the training data
shuffled_ids <- sample(nrow(dat_train))

# prepare data to be used as Ridge input. Separate outcome from predictors and turn data frames into matrices.
x_train <- as.matrix(dat_train[shuffled_ids, 2:(n_col - 1)])
x_test  <- as.matrix(dat_test[, 2:(n_col - 1)])
# for some reason R thinks n_col is a character despite working above, so I use the column number instead
y_train <- as.matrix(dat_train[shuffled_ids, 1392])
y_test  <- as.matrix(dat_test[, 1392])

# run ridge regression on training set with cross-validation
cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "multinomial")

# run ridge regression on training set using lambda.1se
ridge <- glmnet(x_train, y_train, alpha = 0, family = "multinomial", lambda = cv_ridge$lambda.1se)

# use model to predict test set
test_ridge <- ridge %>% predict(newx = x_test)

# obtain predicted classes
# because there are 3 tumor classes, I select the one with the highest probability for each observation
pred_class <- apply(test_ridge, 1, which.max)
pred_class <- pred_class - 1

# confusion matrix and accuracy score
cf_mat <- table(observed = y_test, predicted = pred_class)
cf_mat

accuracy <- sum(diag(cf_mat)) / sum(cf_mat)
accuracy
# accuracy is 0.5
