library(tidyverse)

## Function: Given a df, create a confusion matrix and calculate the AUC score ##
calculate_metrics <- function(df, observed_colname, predicted_colname) {
  
  # Create confusion matrix
  confusion_matrix <- table(df[[observed_colname]], df[[predicted_colname]])
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  # Calculate AUC score
  auc_score <- roc(df[[observed_colname]], df[[predicted_colname]])$auc
  
  # Return confusion matrix and AUC score
  list(confusion_matrix, accuracy, auc_score)
}

## Import the Ridge predicted classes ##

result_ridge <- read.csv("results/result_ridge.csv")
# the file has three columns: Sample, Observed, Predicted

## Import the GCN predicted classes ##

pred_file <- "results/GCN_predicted_data.csv"

# Import the results .csv file
pred_gcn <- read.csv(pred_file)

# Rename columns of the results df
colnames(pred_gcn) <- c("Sample", "Predicted")

# Merge the predicted values to the observed values (we use the Sample and Observed columns in the ridge df)
merged_gcn <- right_join(result_ridge[, 1:2], pred_gcn, by = "Sample")

## Results ##

performance_gcn <- calculate_metrics(merged_gcn, "Observed", "Predicted")
# Area under the curve: 0.5
# Accuracy 0.5
# $confusion_matrix
# 
#   0 1
# 0 2 6
# 1 2 6

performance_ridge <- calculate_metrics(result_ridge, "Observed", "Predicted")
# Accuracy 0.4375
# Area under the curve: 0.5625
# $confusion_matrix
#   0 1
# 0 3 5
# 1 4 4


#### Exploratory data analysis of the raw omics data ####

library(ggplot2)
library(prcomp)
library(RColorBrewer)

# run PCA
pca_omics1 <- prcomp(omics1[, -1], scale. = TRUE)
pca_omics2 <- prcomp(omics2[, -1], scale. = TRUE)
pca_omics3 <- prcomp(omics3[, -1], scale. = TRUE)

## Omics 1 ##

# extract the first two PCs
summary(pca_omics1) # first 2 PCs explain 0.35 of variance
pc1_omics1 <- pca_omics1$x[, 1]
pc2_omics1 <- pca_omics1$x[, 2]
pc_df_omics1 <- data.frame(PC1 = pc1_omics1, PC2 = pc2_omics1, Observed = as.factor(result_ridge$Observed))

# List of ggplot elements that repeat in all graphs
gg_theme <- list(
  stat_ellipse(linetype = 2, lwd = 0.5),
  labs(x = "First Principal Component", y = "Second Principal Component"),
  theme_linedraw(),
  scale_color_brewer(labels = c("Naïve", "Resistant"), palette="Accent"), 
  scale_fill_brewer(palette="Accent"),
  theme(legend.position = "bottom", 
        panel.grid.major = element_line(linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 2.8))
)

# plot the two PC by Observed class
ggplot(pc_df_omics1, aes(x = PC1, y = PC2, color = Observed)) +
  geom_point() +
  gg_theme

## Omics 2 ##

# extract the first two PCs
summary(pca_omics2) # first 2 PCs explain 0.23 of variance
pc1_omics2 <- pca_omics2$x[, 1]
pc2_omics2 <- pca_omics2$x[, 2]
pc_df_omics2 <- data.frame(PC1 = pc1_omics2, PC2 = pc2_omics2, Observed = as.factor(result_ridge$Observed))

# plot the two PC by Observed class
ggplot(pc_df_omics2, aes(x = PC1, y = PC2, color = Observed)) +
  geom_point() +
  gg_theme

## Omics 3 ##

# extract the first two PCs
summary(pca_omics3) # first 2 PCs explain 0.2694 of variance
pc1_omics3 <- pca_omics3$x[, 1]
pc2_omics3 <- pca_omics3$x[, 2]
pc_df_omics3 <- data.frame(PC1 = pc1_omics3, PC2 = pc2_omics3, Observed = as.factor(result_ridge$Observed))

# plot the two PC by Observed class
ggplot(pc_df_omics3, aes(x = PC1, y = PC2, color = Observed)) +
  geom_point() +
  gg_theme
