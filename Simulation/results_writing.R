library(data.table)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(RColorBrewer)

#### Import the performance files ####

readRDS("analyses/performance_gcn.rds")
readRDS("analyses/performance_ridge.rds")

#### Replace confusion matrices with accuracy ####

gcn_perf_df <- performance_gcn

for (it in 1:20) {
  
  for (perm in 1:15) {
    
    permut <- gcn_perf_df[[it]][[perm]]
    conf_mat <- permut[[1]]
    
    if (ncol(conf_mat) == 1) {
      conf_mat <- cbind(conf_mat, 0)
    }
    
    TP <- conf_mat[1]
    TN <- conf_mat[4]
    
    accuracy <- (TP+TN)/sum(conf_mat)
    
    gcn_perf_df[[it]][[perm]][[1]] <- accuracy
    
  }
  
}

ridge_perf_df <- performance_ridge

for (it in 1:20) {
  
  for (perm in 1:15) {
    
    permut <- ridge_perf_df[[it]][[perm]]
    conf_mat <- permut[[1]]
    
    if (ncol(conf_mat) == 1) {
      conf_mat <- cbind(conf_mat, 0)
    }
    
    TP <- conf_mat[1]
    TN <- conf_mat[4]
    
    accuracy <- (TP+TN)/sum(conf_mat)
    
    ridge_perf_df[[it]][[perm]][[1]] <- accuracy
    
  }
  
}

#### Transform lists into long-form data frames ####

gcn_auc <- data.table::rbindlist(gcn_perf_df, idcol = TRUE) %>%
  rename(iteration = `.id`) %>%
  filter(row_number() %% 2 == 0)  %>%
  pivot_longer(cols = colnames(.)[2:16],
               names_to = "permutation",
               values_to = "performance") %>%
  mutate(performance = round(as.numeric(.$performance), 2))

gcn_accuracy <- data.table::rbindlist(gcn_perf_acc, idcol = TRUE) %>%
  rename(iteration = `.id`) %>%
  filter(row_number() %% 2 != 0)  %>%
  pivot_longer(cols = colnames(.)[2:16],
               names_to = "permutation",
               values_to = "performance") %>%
  mutate(performance = round(as.numeric(.$performance), 2))


ridge_auc <- data.table::rbindlist(ridge_perf_df, idcol = TRUE) %>%
  rename(iteration = `.id`) %>%
  filter(row_number() %% 2 == 0)  %>%
  pivot_longer(cols = colnames(.)[2:16],
               names_to = "permutation",
               values_to = "performance") %>%
  mutate(performance = round(as.numeric(.$performance), 2))

ridge_accuracy <- data.table::rbindlist(ridge_perf_df, idcol = TRUE) %>%
  rename(iteration = `.id`) %>%
  filter(row_number() %% 2 != 0)  %>%
  pivot_longer(cols = colnames(.)[2:16],
               names_to = "permutation",
               values_to = "performance") %>%
  mutate(performance = round(as.numeric(.$performance), 2))

#### Combine Ridge and GCN metrics into one long-format data frame ####

# Define the conditions of the 15 simulation scenarios
scenarios <- data.frame(N = rep(c(100, 250, 500, 1000, 5000), 3), Ft = rep(c(100, 500, 1000), each = 5))

all_auc <- cbind(gcn_auc, ridge_auc$performance) %>%
  rename(MoGCN = performance,
         Ridge = "ridge_auc$performance") %>%
  mutate(N = rep(rep(c(100, 250, 500, 1000, 5000), 3), 20),
         Ft = rep(rep(c(100, 500, 1000), each = 5), 20)) %>%
  pivot_longer(., cols = c(MoGCN, Ridge), names_to = "Model", values_to = "AUC")
  
#### Plot ROC AUC values ####

Ft_names <- as_labeller(c(`100` = "Number of Features = 100", `500` = "Number of Features = 500",`1000` = "Number of Features = 1000"))

gg_theme <- list(
  theme_classic(),
  scale_color_brewer(labels = c("MoGCN", "Ridge"), palette = "Accent"),
  scale_fill_brewer(palette = "Accent")
)

# Violin
ggplot(all_auc, aes(x = as.factor(N), y = AUC, color = Model, fill = Model)) +
  geom_violin(position = position_dodge(width = 0.5), alpha = 0.2, size = 0.2, trim = T) +
  geom_jitter(position = position_dodge(width = 0.5)) +
  facet_wrap(~Ft, scales = "free", ncol = 1, dir = "h", labeller = Ft_names) +
  labs(x = "Number of Observations", y = "ROC AUC") +
  ylim(0.38, 1.1) +
  theme_classic() +
  theme(legend.position = "bottom")

# Sina
ggplot(all_auc, aes(x = as.factor(N), y = AUC, color = Model, fill = Model)) +
  facet_wrap(~Ft, scales = "free", ncol = 1, dir = "h", labeller = Ft_names) +
  geom_violin(position = position_dodge(width = 0.7), alpha = 0.2, size = 0.1) +
  geom_sina(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.1) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", size = 0.1) +
  labs(x = "Number of Observations", y = "ROC AUC") +
  ylim(0.38, 1.1) +
  theme_classic() +
  theme(legend.position = "bottom")

# Box
ggplot(all_auc, aes(x = as.factor(N), y = AUC, color = Model)) +
  facet_wrap(~Ft, scales = "free", ncol = 1, dir = "h", labeller = Ft_names) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", size = 0.1) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "black", size = 0.1) +
  labs(x = "Number of Observations", y = "ROC AUC") +
  ylim(0.38, 1.1) +
  gg_theme +
  theme(axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2.8),
        legend.margin = margin(10,20, 0, 40))

#### Calculate average AUC per condition across 20 iterations ####

View(all_auc %>%
  group_by(Model, N, Ft) %>%
  summarize(auc = mean(AUC)))
