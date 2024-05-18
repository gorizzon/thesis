library(UBL)
library(tidyverse)

set.seed(90)

# import .csv files
Meth <- read_csv("raw_data/Meth.csv")
RPPA <- read_csv("raw_data/RPPA.csv")
mRNA <- read_csv("raw_data/mRNA.csv")
tumor <- read_csv("raw_data/outcome_class.csv")
test_ID <- read_csv("data/test_sample.csv")

##### Oversampling #####

# concatenate the 3 omics data sets (remove the sample columns that repeat)
dat <- cbind(Meth, RPPA[, -1], mRNA[, -1])

# two columns have the same name, so make them uniquely indentifiable
colnames(dat) <- make.unique(colnames(dat))

# marge the class data
dat <- right_join(dat, tumor, by = "Sample")
# recode outcome as factor
dat$Tumor <- as.factor(dat$Tumor)

# only training set can be oversampled, so subset the complete data set
dat_train <- dat[!(dat$Sample %in% test_ID$Sample), ]

# Apply SMOTE to create a balanced training set
train_balanced <- SmoteClassif(Tumor ~., dat_train[, -1], C.perc = "balance")

# Create unique IDs for the oversampled training set
Sample <- paste0("balanced_", c(1:nrow(train_balanced)), sep = "")
# Add the ID column to the oversampled data
train_balanced <- cbind(Sample, train_balanced)
# SMOTE created a balanced training set of 275 observations

# Divide the data back into the 3 separate omics data sets
Meth_train_bal <- train_balanced %>% select(colnames(Meth))
mRNA_train_bal <- train_balanced %>% select(colnames(mRNA))
RPPA_train_bal <- train_balanced %>% select(colnames(RPPA))

# Write the balanced training data sets as csv files
write_csv(Meth_train_bal, file = "data/Meth_train_bal.csv")
write_csv(mRNA_train_bal, file = "data/mRNA_train_bal.csv")
write_csv(RPPA_train_bal, file = "data/RPPA_train_bal.csv")

# Write the balanced class file (including only the sample and Tumor class)
write_csv(train_balanced[, c(1, 1392)], "data/outcome_train_bal.csv")


# Add the test data back in to recreate the complete data sets
dat_test <- dat[(dat$Sample %in% test_ID$Sample), ]
dat_balanced <- rbind(train_balanced, dat_test)

# Again divide the training + test balanced data back into separate omics
Meth_bal <- dat_balanced %>% select(colnames(Meth))
mRNA_bal <- dat_balanced %>% select(colnames(mRNA))
RPPA_bal <- dat_balanced %>% select(colnames(RPPA))

# Write csv files
write_csv(Meth_bal, file = "data/Meth_bal.csv")
write_csv(mRNA_bal, file = "data/mRNA_bal.csv")
write_csv(RPPA_bal, file = "data/RPPA_bal.csv")

# Write class csv file
write_csv(dat_balanced[, c(1, 1392)], "data/outcome_bal.csv")
