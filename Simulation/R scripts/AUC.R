library(tidyverse)
library(readr)
library(pROC)

set.seed(22)

####### Function: create confusion matrix and calculate auc score from a data frame of observed and predicted classes #######

# Given a df, create a confusion matrix and calculate the AUC score
calculate_metrics <- function(df, observed_colname, predicted_colname) {
  
  # Create confusion matrix
  confusion_matrix <- table(df[[observed_colname]], df[[predicted_colname]])
  
  # Calculate AUC score
  auc_score <- roc(df[[observed_colname]], df[[predicted_colname]])$auc
  
  # Return confusion matrix and AUC score
  list(confusion_matrix = confusion_matrix, auc_score = auc_score)
}

####### Collect GCN predictions #######

# Create an empty list of lists
all_data <- list()

# List all synthetic data folders
folder_list <- list.dirs("data", full.names = TRUE) %>%
  # Path is saved as first element, so remove it
  .[-1] %>%
  # Add / to make path name complete
  paste0("/", sep = "")

# Loop through each folder to import observed classes
for (i in 1:length(folder_list)) {
  
  # Create a new empty list per folder
  list_name <- paste0("iter_", LETTERS[i])
  all_data[[list_name]] <- list()
  
  # List all files that match the pattern "class_permXX.csv"
  files <- list.files(path = folder_list[i], pattern = "class_perm\\d{2}\\.csv", full.names = TRUE)
  
  # Loop through each file in the folder
  for (file in files) {
    
    # Format file name
    file_name <- basename(file) %>% gsub("\\.csv$", "", .)
    
    # Read csv file
    df <- read.csv(file)

    # Append the current df to its list
    all_data[[list_name]][[file_name]] <- df
  }
}

## Now import the GCN predicted classes ##

# List all results folders
folder_list <- list.dirs("results", full.names = TRUE) %>%
  # Path is saved as first element, so remove it
  .[-1] %>%
  # Keep only the iteration folders, remove innermost children
  .[seq(1, length(.), by = 16)]

# Loop through each folder to import predicted classes and merge them with the observed classes
for (i in 1:length(folder_list)) {
  
  # Recreate names of lists inside all_data
  iter_name <- paste0("iter_", LETTERS[i])
  
  # List subfolders within each simulation iteration
  subfolders <- list.dirs(folder_list[i], full.names = TRUE) %>%
    # remove parent and add /
    .[-1] %>% paste0("/", sep = "")
  
  # Loop through each subfolder to import predictions
  for (j in 1:length(subfolders)) {
    
    # Define file path
    pred_file <- paste0(subfolders[j], "GCN_predicted_data.csv", sep = "")
    
    # Import the results .csv file
    pred_df <- read.csv(pred_file)
    
    # Rename columns of the results df
    colnames(pred_df) <- c("Sample", "predicted_label")
    
    # Read dataframe with observed data from all_data
    obs_df <- all_data[[iter_name]][[j]]
    
    # Merge the pred_df to obs_df, right join to only keep test cases
    merged_df <- right_join(obs_df, pred_df, by = "Sample")
    
    # Update the dataframe in all_data
    all_data[[iter_name]][[j]] <- merged_df
  }
}

# all_data is a list of 20 lists, the 20 lists represent the iterations of the simulation.
# Within each of these 20 lists, there are 15 data frames with the Sample ID, Observed class and GCN Predicted class for the 15 combinations of sample size * number of parameters.

####### Calculate GCN metrics #######

performance_gcn <- all_data

# Iterate through each list in performance_gcn
for (list_name in names(performance_gcn)) {
  
  # Iterate through each dataframe in the list
  for (df_name in names(performance_gcn[[list_name]])) {
    
    # Get the dataframe
    df <- performance_gcn[[list_name]][[df_name]]
    
    # Calculate confusion matrix and AUC score
    metrics <- calculate_metrics(df, "observed_label", "predicted_label")
    
    # Replace the dataframe with confusion matrix and AUC score
    performance_gcn[[list_name]][[df_name]] <- metrics
  }
}

####### Calculate Ridge metrics #######

readRDS("analyses/ridge.rds")

performance_ridge <- net_aggr_data

n_perm <- 15

# Iterate through each list in performance_ridge
for (list_name in names(performance_ridge)) {
  
  # Iterate through each dataframe in the list
  for (i in 1:n_perm) {
    
    df <- performance_ridge[[list_name]][[i]]
    
    # Calculate confusion matrix and AUC score
    metrics <- calculate_metrics(df, "Observed", "Predicted")
    
    # Replace the dataframe with confusion matrix and AUC score
    performance_ridge[[list_name]][[i]] <- metrics
  }
}

# In the end performance_gcn and performance_ridge are both lists of 20 lists of 15 lists.
# Within each of the 20 lists there are 15 lists containing two items: the confusion matrix and ROC AUC for the specific simulation condition in that simulation iteration.

saveRDS(performance_gcn, file = "analyses/performance_gcn.rds")
saveRDS(performance_ridge, file = "analyses/performance_ridge.rds")
