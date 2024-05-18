library(ggplot2)
library(prcomp)
library(RColorBrewer)

# import data
Meth_bal <- read_csv("data/Meth_bal.csv")
RPPA_bal <- read_csv("data/RPPA_bal.csv")
mRNA_bal <- read_csv("data/mRNA_bal.csv")

# run PCA
pca_mRNA <- prcomp(mRNA_bal[, -1], scale. = TRUE)
pca_RPPA <- prcomp(RPPA_bal[, -1], scale. = TRUE)
pca_Meth <- prcomp(Meth_bal[, -1], scale. = TRUE)

###### mRNA ######

# extract the first two PCs
summary(pca_mRNA) # explain 0.257 of variance
pc1_mRNA <- pca_mRNA$x[, 1]
pc2_mRNA <- pca_mRNA$x[, 2]
pc_df_mRNA <- data.frame(PC1 = pc1_mRNA, PC2 = pc2_mRNA, Tumor = dat_balanced$Tumor)

# create list to include all the ggplot extras that repeat
gg_theme <- list(
  stat_ellipse(linetype = 2, lwd = 0.5),
  labs(x = "First Principal Component", y = "Second Principal Component"),
  theme_linedraw(),
  scale_color_brewer(labels = c("T1", "T2", "T3"), palette="Accent"), 
  scale_fill_brewer(palette="Accent"),
  theme(legend.position = "bottom", 
        panel.grid.major = element_line(linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(vjust = -3),
        axis.title.y = element_text(vjust = 2.8))
)

# plot the two PC by Tumor class
ggplot(pc_df_mRNA, aes(x = PC1, y = PC2, color = Tumor)) +
  geom_point() +
  gg_theme

###### RPPA ######

# extract the first two PCs
summary(pca_RPPA) # explains 0.2615
pc1_RPPA <- pca_RPPA$x[, 1]
pc2_RPPA <- pca_RPPA$x[, 2]
pc_df_RPPA <- data.frame(PC1 = pc1_RPPA, PC2 = pc2_RPPA, Tumor = dat_balanced$Tumor)

# plot the two PC by Tumor class
ggplot(pc_df_RPPA, aes(x = PC1, y = PC2, color = Tumor)) +
  geom_point()  +
  gg_theme

###### Methylation ######

# extract the first two PCs
summary(pca_Meth) # explains 0.2696
pc1_Meth <- pca_Meth$x[, 1]
pc2_Meth <- pca_Meth$x[, 2]
pc_df_Meth <- data.frame(PC1 = pc1_Meth, PC2 = pc2_Meth, Tumor = dat_balanced$Tumor)

# plot the two PC by Tumor class
ggplot(pc_df_Meth, aes(x = PC1, y = PC2, color = Tumor)) +
  geom_point()  +
  gg_theme
