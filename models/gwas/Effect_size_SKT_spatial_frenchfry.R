##########################################

# Script: GWAS drone data - SKT

# Trine Aalborg
# Version: 2024-09-24

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
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(cowplot)
library(GWASpoly)
library(AGHmatrix)
library(qqman)

###################################


###################################
###### Load data files       ######
###################################

load("results/GWAS_SKT_result_spatial_BLUE_frenchfry.RData")

geno_original <- read.csv("results/GWASpoly_genofile_SKT_spatial.csv") 

load("results/GWAS_RESPOT_SKT_allchr_bonf_sigSNP_spatial_BLUE_frenchfry.RData")

pt <- read.csv("results/GWASpoly_phenofile_SKT_spatial_BLUE_frenchfry.csv")

load("data/genotypes/SKT_final_gt.RData")
# Remove IdX not present in pt (those with missing market segment)
gt <- gt_SKT[rownames(gt_SKT) %in% as.character(unique(pt$IdX)), ]


###################################


###################################
###### Calculate MAF         ######
###################################


geno <- geno_original

### Calculate MAF ###

# Calculate the total allele count per row (after replacing NAs)
allele_sum <- rowSums(geno[, 4:ncol(geno)], na.rm = TRUE)

# Number of individuals is now the same for each row since NAs are replaced
num_individuals <- ncol(geno[, 4:ncol(geno)]) 

# Total number of alleles is 4 times the number of individuals (tetraploid system)
total_alleles <- 4 * num_individuals

# Calculate the alternative allele frequency
alt_allele_freq <- allele_sum / total_alleles

# Calculate the reference allele frequency (1 - alternative allele frequency)
ref_allele_freq <- 1 - alt_allele_freq

# Calculate the Minor Allele Frequency (MAF) as the smaller of the two
maf <- pmin(alt_allele_freq, ref_allele_freq)

# Add MAF as a new column to the data frame
geno$MAF <- maf

maf_frame <- geno[, c("SNP", "chr", "pos", "MAF")]
colnames(maf_frame) <- c("Marker", "chr", "pos", "MAF")

###################################



###################################
###### Extract effect sizes  ######
###################################

beta <- data.frame(data.loco.scan@effects)
beta$Marker <- rownames(beta)
colnames(beta) <- c("beta", "Marker")

# Add correlation coefficient to data frame
sig <- left_join(snp_select_bonferroni_SKT, beta, by = "Marker")

# Add MAF to data frame
sig <- left_join(sig, maf_frame[, c("Marker", "MAF")], by = "Marker")


###################################



