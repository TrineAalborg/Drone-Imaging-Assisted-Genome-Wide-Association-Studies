####################################
#      Script: Prep market segment, trial and year groups for phenotypes
# TRAA NOVEMBER 2024

# subjects: genotyped RESPOT 2020-2024
####################################


####################################
#       Ready environment
####################################

###### Set working directory #######
setwd("/home/bio.aau.dk/vm82rp/RESPOT/GWAS_drone")
rm(list=ls())

set.seed(123)

###### Load packages ######

library(dplyr)
library(readxl)
library(ggplot2)
library(ggthemes)
library(ggrepel)

###############################################
#          Load data elements                 #
###############################################



# Final genotype matrix post filtering
load("data/genotypes/Geno_2020_to_2024_QC_imputed.RData")

# List of clone parents, and market segment
clones <- read_excel("data/original_danespo_files/genotypes/SNP Genotypes incl. 2024.xlsx", sheet = 8,
                         col_types = c("numeric", "text", "text", "text", "numeric", "text", "text", "numeric", "numeric"))
colnames(clones) <- c("Lbnr", "SOR", "IdX", "Remarks", "StamNrRef", "Mother", "Father", "GRP", "Year")
clones$SOR <- toupper(clones$SOR)

clones <- clones[!(clones$Year == 2019), ]

# Remove duplicate entries
clones <- clones %>%
  distinct(IdX, .keep_all = TRUE)

###############################################
#          Filter to genotyped IDX            #
###############################################

gt_clones <- rownames(gt)

clones_sub <- clones[clones$IdX %in% gt_clones, -c(1,4:5,9)]
clones_sub <- distinct(clones_sub)

###############################################
#          Market segment groups              #
###############################################

market <- read_excel("data/original_danespo_files/Marketsegment_Managementtype.xlsx", 
                     col_types = c("numeric", "text", "text"))


market$GRP <- as.factor(market$GRP)
clones_sub$GRP <- as.factor(clones_sub$GRP)

clones_sub <- left_join(clones_sub, market, by = "GRP")

###############################################
#          Reduce to genotyped clones         #
###############################################

gt_market <- data.frame("IdX" = gt_clones)
gt_market <- left_join(gt_market, clones_sub, by = "IdX")


##### Plot of market segment distributions in the genotyped batches ######

ggplot(gt_market, aes(x = GRP_Overall, fill = Target)) +
  geom_bar(position = "stack", color = "black") +
  labs(x = "Market segment", y = "No. clones") +
  theme(legend.position = "right") +
  theme_few()

ggsave(filename = "plots/Market_segment_histogram_genotypedIDX.png", height = 5, width = 7)

##### Save the market segment annotations for phenotype data ######
clone_market <- gt_market
save(clone_market, file = "data/phenotypes/Market_segment_for_gt_clones.RData")



