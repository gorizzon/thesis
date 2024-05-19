library(glmnet)
library(tidyverse)
library(readr)

set.seed(22)

n_perm <- 15
n_iter <- LETTERS[1:20]

# Create 20 empty lists for each iteration
for (letter in n_iter) {
  list_name <- paste0("net_multi_", letter)
  assign(list_name, list())
}

# Function to perform ridge regression on the 18 test sets
net_multi <- function(data_path, n_perm, iter_list){
  
  for (i in 1:n_perm){
    
    # add a leading zero to single digit iterations
    iter <- formatC(i, width = 2, format = "d", flag = "0")
    
    # define path names to import data files
    oma <- paste(data_path, "oma_perm", iter, ".csv", sep = "")
    omb <- paste(data_path, "omb_perm", iter, ".csv", sep = "")
    omc <- paste(data_path, "omc_perm", iter, ".csv", sep = "")
    class <- paste(data_path, "class_perm", iter, ".csv", sep = "")
    test_ID <- paste(data_path, "test_perm", iter, ".csv", sep = "")
    
    # import data .csv files
    oma <- read_csv(oma)
    omb <- read_csv(omb)
    omc <- read_csv(omc)
    class <- read_csv(class)
    test_ID <- read_csv(test_ID)
    
    # rename column of test dataframe
    colnames(test_ID) <- "Sample_ID"
    
    # aggregate features and outcome into one dataframe
    dat <- cbind(oma, omb[, -1], omc[, -1], class[, -1])
    
    # split data into train set and test set
    train <- dat[!(dat$Sample_ID %in% test_ID$Sample_ID), ]
    test  <- dat[dat$Sample_ID %in% test_ID$Sample_ID, ]
    
    # obtain number of columns of aggregated dataset
    n_col <- ncol(dat)
    
    x_train <- as.matrix(train[, 2:(n_col - 1)])
    x_test  <- as.matrix(test[, 2:(n_col - 1)])
    y_train <- as.matrix(train[, n_col])
    y_test  <- as.matrix(test[, n_col])
    
    # run ridge regression on training set with cross-validation
    cv_ridge <- cv.glmnet(x_train, y_train, alpha = 0, family = "binomial")
    
    # run ridge regression on training set with lambda.1se
    ridge <- glmnet(x_train, y_train, alpha = 0, family = "binomial", lambda = cv_ridge$lambda.1se)
    
    # run regression on test set
    test_ridge <- ridge %>% predict(newx = x_test)
    
    # transform predictions into probabilities, then into classes
    probabilities <- 1 / (1 + exp(-test_ridge))
    pred_class <- ifelse(probabilities > 0.5, 1, 0)
    
    # aggregate observed and predicted class in a new dataframe
    result <- cbind(test_ID, y_test, pred_class)
    colnames(result) <- c("Sample", "Observed", "Predicted")
    
    # append dataframe to list
    iter_list[[i]] <- result
    
  }
  
  return(iter_list)
}

# Empty list to aggregate all the lists with ridge predictions
net_aggr_data <- list()

# The loop only works in this weird roundabout way where I need to create 20 empty lists but they remain empty after the loop.
# In the end net_aggr_data is a list of 20 lists, the 20 lists represent one iteration of the simulation.
# Within each of these 20 lists, there are 15 data frames with the Sample ID, Observed class and Ridge Predicted class for the 15 combinations of sample size * number of parameters.
for (letter in n_iter) {
  
  results_list <- paste0("net_", letter)
  path <- paste0("data/data_", letter, "/")
  empty_list <- paste0("net_multi_", letter)
  
  function_args <- list(path, n_perm, get(empty_list))
  ridge <- do.call(net_multi, function_args)
  
  net_aggr_data[[results_list]] <- ridge
}

# Export net_aggr_data
saveRDS(net_aggr_data, file = "analyses/ridge.rds")