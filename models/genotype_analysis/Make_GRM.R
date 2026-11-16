####################################
#      Script: Make GRM RESPOT
# TRAA JUNE 2024

# subjects: genotyped RESPOT 2020-2022
####################################


####################################
#       Ready environment
####################################

###### Set working directory #######
setwd("~/RESPOT/GWAS_drone")
rm(list=ls())

set.seed(123)


library(dplyr)
library(AGHmatrix)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(ggrepel)


######################################################
#                Load data elements                  #
######################################################

load("data/genotypes/Geno_2020_to_2024_QC_imputed.RData")



######################################################
#            Make GRM                                #
######################################################

gt <- as.matrix(gt)

gt_ratio <- gt/4

G_VanRaden <- Gmatrix(gt_ratio, 
                      method="VanRaden", 
                      ploidy=4,
                      thresh.missing = 0.5,
                      ploidy.correction = TRUE,
                      ratio=TRUE)

# Calculate inbreeding coefficient as the mean of the diagonal
inbreeding <- (mean(diag(G_VanRaden))-1)/(4-1)

# Boxplot of the diagonal elements
png("plots/Boxplot_diagonal_G_matrix.png")
boxplot(diag(G_VanRaden), main = "Boxplot of G diagonal elements")
dev.off()

# Heatmap of the GRM
png("plots/G_matrix_heatmap.png", width = 12, height = 12, units = 'in', res = 300)
heatmap(G_VanRaden, labRow = NA, labCol = NA)
dev.off()

# Save the GRM
save(G_VanRaden, file = "data/genotypes/G_matrix_geno_2020_to_2024_QC_imputed.RData")


######################################################
#            PCA of the GRM                          #
######################################################

# Load the market segment annotation to color PCA plot by market segment
load("data/phenotypes/Market_segment_for_gt_clones.RData")

idx <- rownames(gt)

clone_market_ordered <- clone_market[order(match(clone_market$IdX, idx)), ]
save(clone_market_ordered, file = "data/phenotypes/Market_segment_for_gt_clones_ordered.RData")

G <- as.data.frame(G_VanRaden)

###### Principal component analysis of G matrix ######
pca_G <- prcomp(G)

# Using this PCA function (instead of Sergey's MDAtools package)
## sdev^2 = eigenvalues
## x = scores

eig <- (pca_G$sdev)^2 # eigenvalues
variance <- eig*100/sum(eig) # variance explained by each component
head(variance)
scores <- as.data.frame(pca_G$x) # scores - for the GRM scores and loadings are symmetrical in plotting space, but vary in scale because the GRM is symmetrical in rows and columns
head(scores[ ,1:4])
loadings <- as.data.frame(pca_G$rotation) #loadings
head(loadings[ ,1:4])

var_G <- data.frame("variance" = variance, 
                    "PC" = c(1:nrow(G)))

variance_plot <- ggplot(data = subset(var_G, PC <= 50), aes(x = PC, y = variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Principal component", y = "Variance explained [%]") +
  ggtitle("Variance explained by Principal Components") +
  theme_few()
plot_grid(variance_plot)

ggsave("plots/PCA_RESPOT_G_VanRaden_var_G.png", width = 5, height = 5)

##### PCA plot from scores ######

pc12 = ggplot(scores, aes(x = PC1, y = PC2, colour = clone_market_ordered$GRP_Overall), alpha = 0.5) +
  geom_point() + 
  theme_few() +
  labs(x = paste0("PC1 (", round(variance[1], 2)," %)"),
       y = paste0("PC2 (", round(variance[2], 2)," %)"),
       #title = "PC1 vs. PC2",
       col = "Market group")
pc13 = ggplot(scores, aes(x = PC1, y = PC3, colour = clone_market_ordered$GRP_Overall)) +
  geom_point() + 
  theme_few() +
  labs(x = paste0("PC1 (", round(variance[1], 2)," %)"),
       y = paste0("PC3 (", round(variance[3], 2)," %)"),
       #title = "PC1 vs. PC3"
       )

plot_grid(
  plot_grid(pc12 + theme(legend.position = "none"), 
            pc13 + theme(legend.position = "none"), nrow = 1),
  cowplot::get_legend(pc12), 
  rel_widths = c(3, 0.6)
)

ggsave("plots/PCA_RESPOT_G_VanRaden_PCA_G.png", width = 11, height = 5)

