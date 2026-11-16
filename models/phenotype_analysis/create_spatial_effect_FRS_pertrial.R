##########################################

# Script: GWAS drone data - emergence 2023 - create spatial effect

# Trine Aalborg
# Version: 2024-09-03

# RESPOT

##########################################

###################################
###### Data explanation      ######
###################################

# We have two main traits recorded by drone

# 1) Emergence
#   - A combination trait of 3 individual traits: Plant count, Total Plant Biomass, Plant biomass uniformity 
#   - At every flight (time point), 3 grades are given:
#      1. time point: grades 9-8-7 (top score)
#      2. time point: grades 6-5-4
#      3. time point: grades 3-2-1
#   - The grades are given once and are a combination of the 3 indiv. traits
#   - This separates 3 groups, fast-mid-slow
#   - The first time point is given (flown) when 1/3 of plants are emerged
#   - At the second time point, grades are given to only newly emerged plants


# 2) Late blight resistance (SKT)
#   -  1 is susceptible, 9 is non-susceptible
#   -  Recorded in 6 time points

# Hence GWAS analysis should rely only on Emergence and SKT (late blight) as the other recorded phenotypes
# are already incorporated in the estimation of the total Emergence phenotype. 

# The Emergence1-4 are for 4 replicates, replicate 4 have very few entries and is the latest flight. 
# Based on data analysis, replicates 1-3 are identical for Emergence. 


###################################
###### Set working directory ######
###################################

rm(list = ls())
setwd("~/RESPOT/GWAS_drone")

# Set pakcages

library(plyr)
library(dplyr)
library(tidyr)
library(Matrix)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(cowplot)

###################################

###################################
###### Load data files       ######
###################################

load("data/phenotypes/FRS_final_pt.RData")

pt <- pt_FRS

###################################
# Recode Trial into a numeric code
pt$Trial <- as.numeric(as.factor(pt$Trial))

# Create unique codes for spatial effects                  
pt$spid <- (pt$Year * 100000000000) + (pt$Trial * 10000000) + (pt$Field * 100000) + (pt$Lane * 10000) + (pt$Row * 100) + (pt$Block + 1000)
length(unique(pt$spid)) # Check the unique spatial IDs

# Create Year-Field-Trial ID
pt$Field3 <- paste(pt$Year, pt$Field, pt$Trial, pt$Lane, sep = "_")

unique_lanes_per_field3 <- pt %>%
  group_by(Field3) %>%
  summarise(n_unique_lanes = n_distinct(Lane)) %>%
  ungroup()
all(unique_lanes_per_field3$n_unique_lanes == 1) # as there are some trials with multiple lanes, the spatial grid should be computed for each lane in each trial. 

# Initialize empty list to collect all results
spatial_list <- list()

# Loop through each unique Field3
for (f in unique(pt$Field3)) {
  
  sub_df <- pt[pt$Field3 == f, ]
  
  df_temp <- data.frame(x = sub_df$Block + 1, 
                        y = sub_df$Row + 1, 
                        spid = sub_df$spid)
  
  # Create matrix dimensions (with buffer)
  nrow_mat <- max(df_temp$y) + 2
  ncol_mat <- max(df_temp$x) + 2
  mat <- matrix(rep(NA, nrow_mat * ncol_mat), ncol = ncol_mat)
  
  # Fill matrix with spid values
  for(i in 1:nrow(df_temp)) {
    mat[df_temp$y[i], df_temp$x[i]] <- as.character(df_temp$spid[i])
  }
  
  # Fill empty cells with random unique dummy IDs
  mat[is.na(mat)] <- sample(1e5:1.5e5, size = sum(is.na(mat)), replace = FALSE)
  
  # Extract neighbors
  df_temp$east      <- mat[df_temp$y     + nrow(mat) * (df_temp$x - 2)]
  df_temp$west      <- mat[df_temp$y     + nrow(mat) * (df_temp$x)]
  df_temp$south     <- mat[df_temp$y + 1 + nrow(mat) * (df_temp$x - 1)]
  df_temp$north     <- mat[df_temp$y - 1 + nrow(mat) * (df_temp$x - 1)]
  df_temp$northeast <- mat[df_temp$y - 1 + nrow(mat) * (df_temp$x - 2)]
  df_temp$northwest <- mat[df_temp$y - 1 + nrow(mat) * (df_temp$x)]
  df_temp$southeast <- mat[df_temp$y + 1 + nrow(mat) * (df_temp$x - 2)]
  df_temp$southwest <- mat[df_temp$y + 1 + nrow(mat) * (df_temp$x)]
  
  # Keep only relevant columns
  df_temp <- df_temp[, c("spid", "east", "west", "south", "north", 
                         "northeast", "northwest", "southeast", "southwest")]
  
  spatial_list[[f]] <- df_temp
}

# Combine all spatial neighbor data
df_spatial_neighbors <- do.call(rbind, spatial_list)

# Merge with pt
pt_sp <- left_join(pt, df_spatial_neighbors, by = "spid")

# Track original row order using a stable ID
pt_sp$rowid <- seq_len(nrow(pt_sp))
pt_sp$spid_char <- as.character(pt_sp$spid)
 
###################################

#pt_sp <- pt_sp[,-c(32:33)] # Remove x and y coordinates
#pt_sp <- pt_sp[, c(2:15, 30:32, 1, 33:40, 16:29, 41:42)] # Reorder columns

Y <- pt_sp
str(Y)

###################################
# Create incidence matrices for spatial effects
Y$spid <- as.factor(as.character(Y$spid))
Zspid <- model.matrix(~Y$spid - 1)

Y$west <- as.factor(as.character(Y$west))
Zwest <- model.matrix(~Y$west - 1)

Y$east <- as.factor(as.character(Y$east))
Zeast <- model.matrix(~Y$east - 1)

Y$south <- as.factor(as.character(Y$south))
Zsouth <- model.matrix(~Y$south - 1)

Y$north <- as.factor(as.character(Y$north))
Znorth <- model.matrix(~Y$north - 1)

Y$northwest <- as.factor(as.character(Y$northwest))
Znorthwest <- model.matrix(~Y$northwest - 1)

Y$northeast <- as.factor(as.character(Y$northeast))
Znortheast <- model.matrix(~Y$northeast - 1)

Y$southwest <- as.factor(as.character(Y$southwest))
Zsouthwest <- model.matrix(~Y$southwest - 1)

Y$southeast <- as.factor(as.character(Y$southeast))
Zsoutheast <- model.matrix(~Y$southeast - 1)

# Adjust column names for each matrix
colnames(Zspid) <- substring(colnames(Zspid), 7)
colnames(Zwest) <- substring(colnames(Zwest), 7)
colnames(Zeast) <- substring(colnames(Zeast), 7)
colnames(Zsouth) <- substring(colnames(Zsouth), 8)
colnames(Znorth) <- substring(colnames(Znorth), 8)
colnames(Znorthwest) <- substring(colnames(Znorthwest), 12)
colnames(Znortheast) <- substring(colnames(Znortheast), 12)
colnames(Zsouthwest) <- substring(colnames(Zsouthwest), 12)
colnames(Zsoutheast) <- substring(colnames(Zsoutheast), 12)

###################################
# Combine the incidence matrices into a full matrix
temp <- cbind(Zspid, Zwest, Zsouth, Znorth, Znorthwest, Znortheast, Zsouthwest, Zeast, Zsoutheast)
dim(temp) # Check the dimensions of the matrix

# Create summarized version of the incidence matrix by summing duplicates
temp2 <- sapply(unique(colnames(temp)), function(x) rowSums(temp[, colnames(temp) == x, drop = FALSE]))


# Ensure row order matches original by using rowid
Y_ordered <- Y[order(Y$rowid), ]
rownames(temp2) <- Y_ordered$spid_char  # Explicit assignment of spid_char


# Calculate the spatial similarity matrix
RR <- temp2 %*% t(temp2)
RR_scaled <- RR / mean(diag(RR)) # Scale the matrix

# Assign correct spid names to rows and columns
rownames(RR_scaled) <- rownames(temp2)
colnames(RR_scaled) <- rownames(temp2)

# Save the final spatial similarity matrix and the data frame with spatial effects
save(RR_scaled, file = "data/phenotypes/spatial_matrix_scaled_FRSpertrial.RData")
save(pt_sp, file = "data/phenotypes/FRS_final_pt_spatialpertrial.RData")


###################################
###### PCA of SP-matrix      ######
###################################


pca_RR <- prcomp(RR_scaled)

# Using this PCA function (instead of Sergey's MDAtools package)
## sdev^2 = eigenvalues
## x = scores

eig <- (pca_RR$sdev)^2 # eigenvalues
variance <- eig*100/sum(eig) # variance explained by each component
head(variance)
scores <- as.data.frame(pca_RR$x) # scores - for the GRM scores and loadings are symmetrical in plotting space, but vary in scale because the GRM is symmetrical in rows and columns
head(scores[ ,1:4])
loadings <- as.data.frame(pca_RR$rotation) #loadings
head(loadings[ ,1:4])

var_RR <- data.frame("variance" = variance, 
                     "PC" = c(1:nrow(RR_scaled)))

variance_plot <- ggplot(data = subset(var_RR, PC <= 200), aes(x = PC, y = variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Principal component", y = "Variance explained [%]") +
  ggtitle("Variance explained by Principal Components - FRS") +
  theme_few()
plot_grid(variance_plot)

ggsave("plots/PCA_RESPOT_var_spatialmatrix_FRS_pertrial.png", width = 5, height = 5)

png("plots/heatmap_RESPOT_spatialmatrix_FRSpertrial.png", width = 500, height = 500)
heatmap(RR_scaled, 
        Rowv = NA,  # No row clustering
        Colv = NA  # No column clustering
        )
dev.off()

png("plots/heatmap_RESPOT_spatialmatrix_FRSpertrial_1_500.png", width = 500, height = 500)
heatmap(RR_scaled[1:500,1:500], 
        Rowv = NA,  # No row clustering
        Colv = NA  # No column clustering
)
dev.off()
