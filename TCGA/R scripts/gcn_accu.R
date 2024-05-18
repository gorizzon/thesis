library(readr)
library(tidyverse)

# import data .csv files
Meth <- read_csv("data/Meth_bal.csv")
RPPA <- read_csv("data/RPPA_bal.csv")
mRNA <- read_csv("data/mRNA_bal.csv")
tumor <- read_csv("data/outcome_bal.csv")
test_ID <- read_csv("data/test_sample.csv")

# import MoGCN predictions
gcn_pred_class <- read_csv("results/GCN_predicted_data.csv")

# concatenate the 3 omics + class into one dataframe (remove the sample columns that repeat)
dat_balanced <- cbind(Meth, RPPA[, -1], mRNA[, -1], tumor[, -1])

# obtain the test set
dat_test <- dat_balanced[(dat_balanced$Sample %in% test_ID$Sample), ]

# merge the MoGCN predictions to the test data
classes_df <- right_join(dat_test[, c(1, 1392)], gcn_pred_class, by = join_by(Sample))

# confusion matrix and accuracy score
cf_mat_gcn <- table(observed = classes_df$Tumor, predicted = classes_df$predict_label)
cf_mat_gcn

accuracy_gcn <- sum(diag(cf_mat_gcn)) / sum(cf_mat_gcn)
accuracy_gcn
# accuracy is 0.4714286