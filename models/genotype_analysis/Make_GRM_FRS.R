##########################################

# Script: GWAS drone data - Make GRM for FRS

# Trine Aalborg
# Version: 2024-09-04

# RESPOT

##########################################


###################################
###### Set working directory ######
###################################

setwd("~/RESPOT/GWAS_drone")
rm(list=ls())

# Set pakcages

library(plyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(cowplot)
library(AGHmatrix)


###################################


###################################
###### Load data files       ######
###################################

load("data/genotypes/FRS_final_gt.RData")

# Load the market segment annotation to color PCA plot by market segment
load("data/phenotypes/Market_segment_for_gt_clones.RData")

###################################


# Remove clones with missing marker segment

Missing_GRP <- unique(clone_market[is.na(clone_market$GRP_Overall), "IdX"])

gt_FRS <- gt_FRS[!(rownames(gt_FRS) %in% Missing_GRP), ]

######################################################
#            Make GRM                                #
######################################################

gt <- as.matrix(gt_FRS)


gt_ratio <- gt/4

G_VanRaden <- Gmatrix(gt_ratio, 
                      method="VanRaden", 
                      ploidy=4,
                      thresh.missing = 0.5,
                      ploidy.correction = TRUE,
                      ratio=TRUE)

G_FRS <- G_VanRaden

# Calculate inbreeding coefficient as the mean of the diagonal
inbreeding <- (mean(diag(G_FRS))-1)/(4-1)

# Boxplot of the diagonal elements
png("plots/Boxplot_diagonal_G_matrix_FRS.png")
boxplot(diag(G_FRS), main = "Boxplot of G diagonal elements")
dev.off()

# Heatmap of the GRM
png("plots/G_matrix_heatmap_FRS.png", width = 12, height = 12, units = 'in', res = 300)
heatmap(G_FRS, labRow = NA, labCol = NA)
dev.off()

# Save the GRM
save(G_FRS, file = "data/genotypes/G_matrix_geno_2020_to_2024_QC_imputed_FRS.RData")


######################################################
#            PCA of the GRM                          #
######################################################


idx <- rownames(gt)

clone_market_ordered <- clone_market[order(match(clone_market$IdX, idx)), ]
clone_market_ordered <- clone_market_ordered[clone_market_ordered$IdX %in% idx, ]
save(clone_market_ordered, file = "data/phenotypes/Market_segment_for_gt_clones_ordered_FRS.RData")

G <- as.data.frame(G_FRS)

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
                    "PC" = c(1:nrow(gt)))

variance_plot <- ggplot(data = subset(var_G, PC <= 50), aes(x = PC, y = variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Principal component", y = "Variance explained [%]") +
  ggtitle("Variance explained by Principal Components - FRS") +
  theme_few()
plot_grid(variance_plot)

ggsave("plots/PCA_RESPOT_G_VanRaden_var_G_FRS.png", width = 5, height = 5)


# Calculate cumulative variance
var_G <- var_G %>%
  arrange(PC) %>%
  mutate(cumulative_variance = cumsum(variance))

# Find the number of PCs required to explain at least 95% of the variance
pc_95 <- min(var_G$PC[var_G$cumulative_variance >= 95])

# Create the cumulative variance plot
cumulative_variance_plot <- ggplot(data = subset(var_G, PC <= 100), aes(x = PC, y = cumulative_variance)) +
  geom_line(color = "skyblue", linewidth = 1) +
  geom_point(color = "skyblue", size = 2) +
  labs(x = "Principal Component", y = "Cumulative Variance Explained [%]") +
  ggtitle("Cumulative Variance Explained by Principal Components - FRS") +
  theme_minimal() +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
  geom_vline(xintercept = pc_95, linetype = "dashed", color = "red") +
  annotate("text", x = pc_95, y = 50, label = paste0("PCs needed: ", pc_95), color = "red", vjust = -1)

# Plot the cumulative variance
plot_grid(cumulative_variance_plot)

ggsave("plots/PCA_RESPOT_G_VanRaden_cum_var_FRS.png", width = 5, height = 5)



##### PCA plot from scores ######

pc12 = ggplot(scores, aes(x = PC1, y = PC2, colour = clone_market_ordered$GRP_Overall), alpha = 0.5) +
  geom_point() + 
  theme_few() +
  labs(x = paste0("PC1 (", round(variance[1], 2)," %)"),
       y = paste0("PC2 (", round(variance[2], 2)," %)"),
       title = "Emergence",
       col = "Market group")
pc13 = ggplot(scores, aes(x = PC1, y = PC3, colour = clone_market_ordered$GRP_Overall)) +
  geom_point() + 
  theme_few() +
  labs(x = paste0("PC1 (", round(variance[1], 2)," %)"),
       y = paste0("PC3 (", round(variance[3], 2)," %)"),
       title = " "
  )

plot_grid(
  plot_grid(pc12 + theme(legend.position = "none"), pc13 + theme(legend.position = "none"), nrow = 1),
  get_legend(pc12), 
  rel_widths = c(3, 0.6)
)

ggsave("plots/PCA_RESPOT_G_VanRaden_PCA_G_FRS.png", width = 11, height = 5)

##### Year geno plot of PC #####

ID_kept_yeargeno <- read.csv("~/RESPOT/GWAS_drone/data/original_danespo_files/genotypes/Helene_QC_imputed/ID_kept_yeargeno", sep="")
colnames(ID_kept_yeargeno) <- c("IdX", "Rank", "Year_geno")

clone_market_ordered <- left_join(clone_market_ordered, ID_kept_yeargeno, by = "IdX")

pcyear = ggplot(scores, aes(x = PC1, y = PC2, colour = as.factor(clone_market_ordered$Year))) +
  geom_point() + 
  theme_few() +
  labs(x = paste0("PC1 (", round(variance[1], 2)," %)"),
       y = paste0("PC2 (", round(variance[2], 2)," %)"),
       title = "Emergence",
       col = "Year genotyped")

pcyear

ggsave("plots/PCA_RESPOT_G_VanRaden_PCA_G_FRS_yeargeno.png", width = 6, height = 5)

library(ComplexHeatmap)
library(circlize)
metadata <- clone_market_ordered

grp_colors <- setNames(
  c(RColorBrewer::brewer.pal(length(unique(na.omit(metadata$GRP_Overall))), "Set3"), "grey"),
  c(unique(na.omit(metadata$GRP_Overall)), "NA")
)

# Replace NA values in GRP_Overall with "NA" as a string for annotation
metadata$GRP_Overall <- ifelse(is.na(metadata$GRP_Overall), "NA", metadata$GRP_Overall)

# Create the row annotation with the updated colors
row_annotation <- rowAnnotation(
  Market = metadata$GRP_Overall,
  col = list(Market = grp_colors), 
  show_legend = TRUE
)

# Define a white-to-purple color scale for the heatmap
heatmap_colors <- colorRamp2(
  breaks = c(min(G_FRS, na.rm = TRUE), max(G_FRS, na.rm = TRUE)),
  colors = c("white", "purple")
)

# Plot the heatmap with the side bar and custom colors
png("plots/G_matrix_heatmap_FRS_with_annotation.png", width = 12, height = 12, units = 'in', res = 300)
Heatmap(
  G_FRS,
  col = heatmap_colors,          # Apply the custom white-to-purple color scale
  cluster_rows = TRUE,           # Enable row clustering
  cluster_columns = TRUE,        # Enable column clustering
  show_row_names = FALSE,
  show_column_names = FALSE,
  right_annotation = row_annotation
)
dev.off()













