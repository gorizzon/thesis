library(readr)
library(tidyverse)

set.seed(16)

###### Import csv files ######

omics1 <- read_csv("data/omics1.csv")
omics2 <- read_csv("data/omics2.csv")
omics3 <- read_csv("data/omics3.csv")

###### Data manipulation ######

# the omics2 file has the patient class in the second column, so we create a separate class data frame from the ID and class columns
# here label is: 0 = naive, 1 = resistant
label <- omics2[, 1:2]
colnames(label) <- c("Sample", "Label")

# reformat IDs so that they match in all data sets
omics3$Sample <- gsub("_Se1", "", omics3$Sample)
omics1$Sample <- gsub("_PL1", "", omics1$Sample)

# omics2 and omics3 have a smaller number of patients, so we intersect to find the patients that are in all data sets
# and then subset all data sets with these IDs. In the end there's 78 participants that have data for all 3 omics and a label
patients <- intersect(omics3$Sample, omics2$Sample_Name)

omics2 <- omics2[omics2$Sample_Name %in% patients, -2] %>%
  rename_with(~"Sample", Sample_Name)

omics3 <- omics3[omics3$Sample %in% patients, ]

omics1 <- omics1[omics1$Sample %in% patients, ]

label <- label[label$Sample %in% patients, ]

###### Norm to 0-1 ######

omics2_norm <- apply(omics2[, 2:ncol(omics2)], 2, function(x) (x - (min(x)))/max((x) - min(x)))
omics2 <- cbind(omics2[, 1], omics2_norm)

omics3_norm <- apply(omics3[, 2:ncol(omics3)], 2, function(x) (x - (min(x)))/max((x) - min(x)))
omics3 <- cbind(omics3[, 1], omics3_norm)

omics1_norm <- apply(omics1[, 2:ncol(omics1)], 2, function(x) (x - (min(x)))/max((x) - min(x)))
omics1 <- cbind(omics1[, 1], omics1_norm)

###### Export csv files ######

write_csv(omics2, file = "processed_data/omics2.csv")
write_csv(omics3, file = "processed_data/omics3.csv")
write_csv(omics1, file = "processed_data/omics1.csv")

write_csv(label, "processed_data/label.csv")

###### Randomly sample a test set ######

test <- sample(label$Sample, ceiling(nrow(label) * 0.2))
test <- as.data.frame(test)
colnames(test) <- c("Sample")

# export the list of test IDs
write_csv(test, "processed_data/test_sample.csv")

##### Export csv files with only training cases ######

omics1_train <- omics1[!(omics1$Sample %in% test$Sample), ]
omics2_train <- omics2[!(omics2$Sample %in% test$Sample), ]
omics3s_train  <- omics3[!(omics3$Sample %in% test$Sample), ]
label_train  <- label[!(label$Sample %in% test$Sample), ]

write_csv(omics1_train, file = "processed_data/omics1_train.csv")
write_csv(omics2_train, file = "processed_data/omics2_train.csv")
write_csv(omics3_train, file = "processed_data/omics3s_train.csv")

write_csv(label_train, "processed_data/label_train.csv")
