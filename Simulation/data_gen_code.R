library(OmicsPLS)
library(PO2PLS)
library(readr)
library(tidyverse)

#set.seed(100)

# sample(1000, 19)
# 252 516 300 718 636 873  42 812 312 164 380 696 770 285 174  86 531 784 982

#set.seed(252)
#set.seed(516)
#set.seed(300)
#set.seed(718)
#set.seed(636)
#set.seed(873)
#set.seed(42)
#set.seed(812)
#set.seed(312)
# set.seed(164)
# set.seed(380)
# set.seed(696)
# set.seed(770)
# set.seed(285)
# set.seed(174)
# set.seed(86)
# set.seed(531)
# set.seed(784)
# set.seed(982)


gen_par3 <- function(p1, p2, p3, r, rx1, rx2, rx3, alpha = noise){
  prm <- generate_params(p1, p2, r, rx1, rx2, alpha)
  prm3 <- generate_params(p3, p2, r, rx3, rx2, alpha)
  list(W1 = prm$W, 
       W2 = prm$C,
       W3 = prm3$W,
       Wo1 = prm$Wo,
       Wo2 = prm$Co,
       Wo3 = prm3$Wo,
       SigT = prm$SigT,
       SigTo1 = prm$SigTo,
       SigTo2 = prm$SigUo,
       SigTo3 = prm3$SigTo,
       sig2E1 = prm$sig2E,
       sig2E2 = prm$sig2F,
       sig2E3 = prm3$sig2E,
       beta_T = sort(runif(r, -1, 1))
  )
}

gen_dat3 <- function(N, params, alpha_out = 0.1, sparse = FALSE){
  W1 <- params$W1
  W2 <- params$W2
  W3 <- params$W3
  Wo1 <- params$Wo1
  Wo2 <- params$Wo2
  Wo3 <- params$Wo3
  
  r <- ncol(W1)
  rx1 <- ncol(Wo1)
  rx2 <- ncol(Wo2)
  rx3 <- ncol(Wo3)
  p1 <- nrow(W1)
  p2 <- nrow(W2)
  p3 <- nrow(W3)
  
  if(sparse){
    W1[-(1:(p1/4)),] <- 0
    W2[-(1:(p2/4)),] <- 0
    W3[-(1:(p3/4)),] <- 0
  }
  # print(W1)
  
  SigT = params$SigT
  SigTo1 = params$SigTo1 + 1e-06 * SigT[1] * (params$SigTo1[1] == 0)
  SigTo2 = params$SigTo2 + 1e-06 * SigT[1] * (params$SigTo2[1] == 0)
  SigTo3 = params$SigTo3 + 1e-06 * SigT[1] * (params$SigTo3[1] == 0)
  
  #N <- 3/2 * N
  
  Tt <- matrix(rnorm(N * r), N, r) %*% chol(SigT)
  To1 <- matrix(rnorm(N * rx1), N, rx1) %*% chol(SigTo1)
  To2 <- matrix(rnorm(N * rx2), N, rx2) %*% chol(SigTo2)
  To3 <- matrix(rnorm(N * rx3), N, rx3) %*% chol(SigTo3)
  
  E1 <- matrix(rnorm(N * p1), N, p1) * sqrt(params$sig2E1)
  E2 <- matrix(rnorm(N * p2), N, p2) * sqrt(params$sig2E2)
  E3 <- matrix(rnorm(N * p3), N, p3) * sqrt(params$sig2E3)
  
  X1 <- Tt %*% t(W1) + To1 %*% t(Wo1) + E1
  X2 <- Tt %*% t(W2) + To2 %*% t(Wo2) + E2
  X3 <- Tt %*% t(W3) + To3 %*% t(Wo3) + E3
  
  outc <- Tt %*% params$beta_T
  outc <- outc + rnorm(N, 0, sd = sqrt( alpha_out/(1-alpha_out)*var(outc) ))
  
  outc_bin <- 1*(outc > 0)
  
  return(list(X1 = X1, X2 = X2, X3 = X3, outc = outc_bin))
  
}

# alpha_out is default at 0.1
create_data <- function(iteration, test_sample, N, p1, p2, p3, r, rx1, rx2, rx3, noise){
  
  prm <- gen_par3(p1, p2, p3, r, rx1, rx2, rx3, alpha = noise)
  
  dat <- gen_dat3(N, prm, sparse = TRUE)
  
  Sample_ID <- paste(iteration, "_ID", seq(from = 1, to = N), sep = "")
  
  write_csv(cbind(Sample_ID,as.data.frame(dat$X1)), paste("syn_data/oma_", iteration, ".csv", sep = ""))
  write_csv(cbind(Sample_ID,as.data.frame(dat$X2)), paste("syn_data/omb_", iteration, ".csv", sep = ""))
  write_csv(cbind(Sample_ID,as.data.frame(dat$X3)), paste("syn_data/omc_", iteration, ".csv", sep = ""))
  
  # column names match GCN
  ID_class <- as.data.frame(dat$outc) # %>% rename(observed_label = dat$outc)
  colnames(ID_class) <- c("observed_label")
  # print(colnames(ID_class))
  write_csv(cbind(Sample = Sample_ID, ID_class), paste("syn_data/class_", iteration, ".csv", sep = ""))

  test_ID <- paste(iteration, "_ID", test_sample, sep = "")
  
  write_csv(as.data.frame(test_ID), paste("syn_data/test_", iteration, ".csv", sep = ""))
  
}

# sample randomly half of test set from 0s and half from 1s
test_sample_100  <- sample(c(1:100), 20)
test_sample_250  <- sample(c(1:250), 50)
test_sample_500  <- sample(c(1:500), 100)
test_sample_1000 <- sample(c(1:1000), 200)
test_sample_5000 <- sample(c(1:5000), 1000)


create_data <- function(folder, iteration, test_sample, N, p1, p2, p3, r, rx1, rx2, rx3, noise){
  
  prm <- gen_par3(p1, p2, p3, r, rx1, rx2, rx3, alpha = noise)
  
  dat <- gen_dat3(N, prm, sparse = TRUE)
  
  Sample_ID <- paste(iteration, "_ID", seq(from = 1, to = N), sep = "")
  
  write_csv(cbind(Sample_ID,as.data.frame(dat$X1)), paste(folder, "oma_", iteration, ".csv", sep = ""))
  write_csv(cbind(Sample_ID,as.data.frame(dat$X2)), paste(folder, "omb_", iteration, ".csv", sep = ""))
  write_csv(cbind(Sample_ID,as.data.frame(dat$X3)), paste(folder, "omc_", iteration, ".csv", sep = ""))
  
  # column names match GCN
  ID_class <- as.data.frame(dat$outc) # %>% rename(observed_label = dat$outc)
  colnames(ID_class) <- c("observed_label")
  # print(colnames(ID_class))
  write_csv(cbind(Sample = Sample_ID, ID_class), paste(folder, "class_", iteration, ".csv", sep = ""))
  
  test_ID <- paste(iteration, "_ID", test_sample, sep = "")
  
  write_csv(as.data.frame(test_ID), paste(folder, "test_", iteration, ".csv", sep = ""))
  
}

create_iters <- function(seed_n, folder){
  
  set.seed(seed_n)
  
  create_data(folder, iteration = "perm01", test_sample = test_sample_100, N = 100, 
              p1 = 100, p2 = 100, p3 = 100, r = 5, rx1 = 2, rx2 = 2, rx3 = 1,
              noise = 0.7)
  
  create_data(folder, iteration = "perm02", test_sample = test_sample_250, N = 250, 
              p1 = 100, p2 = 100, p3 = 100, r = 5, rx1 = 2, rx2 = 2, rx3 = 1,
              noise = 0.7)
  
  create_data(folder, iteration = "perm03", test_sample = test_sample_500, N = 500, 
              p1 = 100, p2 = 100, p3 = 100, r = 5, rx1 = 2, rx2 = 2, rx3 = 1,
              noise = 0.7)
  
  create_data(folder, iteration = "perm04", test_sample = test_sample_1000, N = 1000, 
              p1 = 500, p2 = 500, p3 = 500, r = 5, rx1 = 2, rx2 = 2, rx3 = 1,
              noise = 0.7)
  
  create_data(folder, iteration = "perm05", test_sample = test_sample_5000, N = 5000, 
              p1 = 500, p2 = 500, p3 = 500, r = 5, rx1 = 2, rx2 = 2, rx3 = 1,
              noise = 0.7)
  
  
  
  create_data(folder, iteration = "perm06", test_sample = test_sample_100, N = 100, 
              p1 = 500, p2 = 500, p3 = 500, r = 25, rx1 = 9, rx2 = 8, rx3 = 8,
              noise = 0.7)
  
  create_data(folder, iteration = "perm07", test_sample = test_sample_250, N = 250, 
              p1 = 500, p2 = 500, p3 = 500, r = 25, rx1 = 9, rx2 = 8, rx3 = 8,
              noise = 0.7)
  
  create_data(folder, iteration = "perm08", test_sample = test_sample_500, N = 500, 
              p1 = 500, p2 = 500, p3 = 500, r = 25, rx1 = 9, rx2 = 8, rx3 = 8,
              noise = 0.7)
  
  create_data(folder, iteration = "perm09", test_sample = test_sample_1000, N = 1000, 
              p1 = 500, p2 = 500, p3 = 500, r = 25, rx1 = 9, rx2 = 8, rx3 = 8,
              noise = 0.7)
  
  create_data(folder, iteration = "perm10", test_sample = test_sample_5000, N = 5000, 
              p1 = 500, p2 = 500, p3 = 500, r = 25, rx1 = 9, rx2 = 8, rx3 = 8,
              noise = 0.7)
  
  
  
  create_data(folder, iteration = "perm11", test_sample = test_sample_100, N = 100, 
              p1 = 1000, p2 = 1000, p3 = 1000, r = 50, rx1 = 17, rx2 = 17, rx3 = 16,
              noise = 0.7)
  
  create_data(folder, iteration = "perm12", test_sample = test_sample_250, N = 250, 
              p1 = 1000, p2 = 1000, p3 = 1000, r = 50, rx1 = 17, rx2 = 17, rx3 = 16,
              noise = 0.7)
  
  create_data(folder, iteration = "perm13", test_sample = test_sample_500, N = 500, 
              p1 = 1000, p2 = 1000, p3 = 1000, r = 50, rx1 = 17, rx2 = 17, rx3 = 16,
              noise = 0.7)
  
  create_data(folder, iteration = "perm14", test_sample = test_sample_1000, N = 1000, 
              p1 = 1000, p2 = 1000, p3 = 1000, r = 50, rx1 = 17, rx2 = 17, rx3 = 16,
              noise = 0.7)
  
  create_data(folder, iteration = "perm15", test_sample = test_sample_5000, N = 5000, 
              p1 = 1000, p2 = 1000, p3 = 1000, r = 50, rx1 = 17, rx2 = 17, rx3 = 16,
              noise = 0.7)
  
}


create_iters(100, "data_A/")
create_iters(252, "data_B/")
create_iters(516, "data_C/")
create_iters(300, "data_D/")
create_iters(718, "data_E/")
create_iters(636, "data_F/")
create_iters(873, "data_G/")
create_iters(42, "data_H/")
create_iters(812, "data_I/")
create_iters(312, "data_J/")
create_iters(164, "data_K/")
create_iters(380, "data_L/")
create_iters(696, "data_M/")
create_iters(770, "data_N/")
create_iters(285, "data_O/")
create_iters(174, "data_P/")
create_iters(86, "data_Q/")
create_iters(531, "data_R/")
create_iters(784, "data_S/")
create_iters(982, "data_T/")

