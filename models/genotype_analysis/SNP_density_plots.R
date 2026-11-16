####################################
#      Script: Plot SNP density
# TRAA JUNE 2024

# subjects: genotyped RESPOT 2020-2023
####################################


####################################
#       Ready environment
####################################

###### Set working directory #######
setwd("~/RESPOT/GWAS_drone")
rm(list=ls())

set.seed(123)

###### Load packages ######

library(dplyr)
library(CMplot)

###############################################
#          Load data elements                 #
###############################################

# Read the QC'ed genotype file
load("data/genotypes/Geno_2020_to_2024_QC_imputed.RData")

# Read the SNP positional file
map <- read.table("data/original_danespo_files/genotypes/map.txt", header = T)

# Replace "-" with "." in the probe name (as is done when loading genotypes, side effect)
map$Probe <- gsub("-", ".", map$Probe)


###############################################
#          Visualize SNP distribution         #
###############################################

# SNP density plot on chromosome according to position - SNP density post filtering

# Order according to 1) chromosome and 2) position.
map <- map[order(map$Chr, map$Position),]

# Subset the map to only include the SNPs left after filtering
map_sub <- map[map$Probe %in% colnames(gt), ]

colnames(map_sub) <- c("SNP", "chr", "pos")

save(map_sub, file = "data/genotypes/SNP_map_Geno_2020_to_2024_QC_imputed.RData")

############### OBS ## OBS ## OBS ###############
# Count of the number of SNPs in on each chromosome after filtering
count <- table(map_sub$chr)
count

# Exclude scaffolds (chromsome 0 and 13?)

map_sub_red <- map_sub[map_sub$chr != 0 & map_sub$chr != 13, ]
############### OBS ## OBS ## OBS ###############

### SNP density plot ###

setwd("/home/bio.aau.dk/vm82rp/RESPOT/GWAS_drone/plots/") # setwd for where the plot is stored automatically

CMplot(map_sub_red,plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",file.name="",dpi=300,
       main="SNP density, 1 Mb window",file.output=TRUE,verbose=TRUE,width=9,height=6)


