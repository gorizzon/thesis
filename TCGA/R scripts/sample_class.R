library(readr)
library(dplyr)

set.seed(100)

# read in csv file
outcome <- read_delim("data/raw_outcome.csv", 
                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

# process ID names to match the omics data sets
outcome$Sample <- gsub("\\-", "", outcome$Sample)

# recode Tumor with digits and remove unkown ("TX")
outcome$Tumor <- as.factor(outcome$Tumor)
outcome$Tumor <- recode(outcome$Tumor, T1 = 0, T2 = 1, T3 = 2, T4 = 3, TX = 99)

# keep only the cases for which omics data is available
outcome_subset <- outcome[outcome$Sample %in% Meth$Sample, ]
outcome_subset <- outcome_subset[outcome_subset$Sample != "TCGAB6A0RM", ]

# write the subsetted data frame as csv file
write_csv(outcome_subset, "data/outcome_class.csv")

# sample test set at random
test <- sample(outcome_subset$Sample, 70)
test <- as.data.frame(test)
# column name has to match
colnames(test) <- c("Sample")

# write the data frame as a csv file
write_csv(test, "data/test_sample.csv")

