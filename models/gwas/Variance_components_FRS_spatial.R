##########################################

# Script: GWAS drone data - test GWAS sommer - FRS

# Trine Aalborg
# Version: 2024-09-18

# RESPOT

##########################################


###################################
###### Set working directory ######
###################################

rm(list = ls())
setwd("~/RESPOT/GWAS_drone")

# Set pakcages

library(plyr)
library(dplyr)
library(sommer)
library(GWASpoly)

###################################


###################################
###### Load data files       ######
###################################

genofile <- "results/GWASpoly_genofile_FRS_spatial.csv"

phenofile <- "results/GWASpoly_phenofile_FRS_spatial_BLUE.csv" 

pt <- read.csv("results/GWASpoly_phenofile_FRS_spatial_BLUE.csv", check.names = F)
pt$IdX <- as.character(pt$IdX)

##############################################################
###### Correct phenotypes for spatial and fixed effects ######
##############################################################


###### Make GRM ######

data <- read.GWASpoly(ploidy = 4, 
                      pheno.file = phenofile,
                      geno.file = genofile,
                      format = "numeric",
                      n.traits = 1,
                      delim = ",")

# Create GRM without Leave-One-Chromosome-Out
data.loco <- set.K(data,
                   LOCO=FALSE, # make a single GRM, without LOCO
                   n.core=2)

GRM <- data.frame(data.loco@K$all)

colnames(GRM) <- gsub("^X", "", colnames(GRM))

GRM <- as.matrix(GRM)

pca_G <- prcomp(GRM)

scores <- as.data.frame(pca_G$x)

PC1_3 <- scores[,1:3]

PC1_3$IdX <- rownames(GRM)

pt <- left_join(pt, PC1_3, by = "IdX")


# fit the LMM model to correct for spatial and fixed effects at the same time
mod <- mmer(BLUE ~ PC1 + PC2,
            random = ~ vsr(IdX, Gu = GRM),# + IdX, # replicates are accommodated in GWASpoly automatically without specifying a line effect (and no additional random variance components can be added manually)
            data = pt, verbose = FALSE)

mod_summary <- summary(mod)


# Extract the coefficients of the fixed effects and the variance component of the spatial effect for downstream analysis
vc <- mod_summary$varcomp[, "VarComp"]
se <- mod_summary$varcomp[, "VarCompSE"]
varcomp_df <- data.frame(Variance = vc, SE = se, Component = c("Additive", "Res"))

save(varcomp_df, file = "results/GWASpoly_sommer_varcomp_FRS_spatial_BLUE.RData")

h2 <- varcomp_df[varcomp_df$Component == "Additive", 1] / sum(varcomp_df[,1])

save(h2, file = "results/GWASpoly_sommer_h2_FRS_spatial_BLUE.RData")
