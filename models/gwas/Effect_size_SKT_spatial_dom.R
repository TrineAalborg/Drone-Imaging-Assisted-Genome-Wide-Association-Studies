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

load("results/GWAS_SKT_result_spatial_BLUE_dom.RData")

geno_original <- read.csv("results/GWASpoly_genofile_SKT_spatial.csv")

load("results/GWAS_RESPOT_SKT_allchr_bonf_sigSNP_spatial_BLUE_dom.RData")

pt <- read.csv("results/GWASpoly_phenofile_SKT_spatial_BLUE_dom.csv")

load("data/genotypes/SKT_final_gt.RData")
# Remove IdX not present in pt (those with missing market segment)
gt <- gt_SKT[rownames(gt_SKT) %in% as.character(unique(pt$IdX)), ]

load("results/GWASpoly_sommer_varcomp_SKT_spatial_BLUE.RData")


###################################


###################################
###### Calculate MAF         ######
###################################

### Impute with mode (like GWASpoly) ###

geno <- geno_original

## Define a function to replace NAs with the mode in each row
#replace_na_with_mode <- function(row_data) {
#  # Find the mode (most common value) in the row (excluding NAs)
#  mode_value <- as.numeric(names(sort(table(row_data), decreasing = TRUE)[1]))
#  
#  # Replace NAs with the mode
#  row_data[is.na(row_data)] <- mode_value
#  
#  return(row_data)
#}
#
## Apply this replacement for each row in the columns 4:970
#geno[, 4:970] <- t(apply(geno[, 4:970], 1, replace_na_with_mode))

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

map <- data.frame(data.loco.scan@map)

beta_alt <- data.frame(data.loco.scan@effects$BLUE$`1-dom-alt`, data.loco.scan@scores$BLUE$`1-dom-alt`)
beta_alt$Marker <- map$Marker
colnames(beta_alt) <- c("beta", "score", "Marker")
beta_alt <- beta_alt[!(is.na(beta_alt$beta)),]

beta_ref <- data.frame(data.loco.scan@effects$BLUE$`1-dom-ref`, data.loco.scan@scores$BLUE$`1-dom-ref`)
beta_ref$Marker <- map$Marker
colnames(beta_ref) <- c("beta", "score", "Marker")
beta_ref <- beta_ref[!(is.na(beta_ref$beta)),]

beta <- rbind(beta_alt, beta_ref)
beta <- left_join(beta, map, by = "Marker")
beta$Chrom <- as.numeric(as.character(beta$Chrom))

beta <- beta %>% arrange(Chrom, Position)

# Add correlation coefficient to data frame
sig <- left_join(snp_select_bonferroni_SKT, beta, by = c("Marker","Chrom", "Position", "Ref", "Alt"))

# Add MAF to data frame
maf_frame <- maf_frame[,c("Marker", "MAF")]
sig <- left_join(sig, maf_frame[, c("Marker", "MAF")], by = "Marker")


###################################




###################################
###### SNP effect            ######
###################################

# The genetic variance explained by the SNP (in the unit of the phenotype, as the coefficient beta is)

sig$SNP_effect <- 2 * (sig$beta)^2 * sig$MAF * (1 - sig$MAF) # definition for (Park et al., 2010)



# To get it in % of the phenotypic variance, we need to calculate that. 
# I get the variance components of the sommer model, as the Vp is the sum of all variance components (all not specified in the model is in the Ve)

Vp <- sum(varcomp_df$Variance)

sig$Perc_Vp <- (sig$SNP_effect/Vp)*100


# To get it in % of the additive genetic variance, we need to additive genetic variance component. Though calculated as part of the GWASpoly, it is not extracted.
# I recalculate it using a mixed model on the correct genotypes using the full GRM 

Va <- varcomp_df[varcomp_df$Component == "Additive", "Variance"]

sig$Perc_Va <- (sig$SNP_effect/Va)*100

save(sig, file = "results/GWAS_RESPOT_SKT_allchr_bonf_sigSNP_SNPeffect_spatial_BLUE_dom.RData")

library(writexl)
write_xlsx(sig, "plots/GWAS_RESPOT_SKT_allchr_bonf_sigSNP_SNPeffect_spatial_BLUE_dom.xlsx")


###################################


###################################


###################################
###### SNP dosage plot       ######
###################################

# Plot the frequency of each dosage for the QTL SNPs against the phenotype

pt <- read.csv("results/GWASpoly_phenofile_SKT_spatial.csv")

pt <- pt[,c("IdX", "AUDPC_adj")]
pt$IdX <- as.character(pt$IdX)

# Calculate average phenotype per individual
avg_AUDPC <- aggregate(AUDPC_adj ~ IdX, data = pt, FUN = mean)

QTL.snps <- snp_select_bonferroni_SKT %>% group_by(as.factor(Chrom)) %>% filter(Pcorr == min(Pcorr)) %>% summarise(Marker = Marker[1]) 

SNP.1 <- geno[geno$SNP == QTL.snps$Marker[1],-ncol(geno)]

SNP.1 <- SNP.1[,-c(1:3)]

SNP.1 <- as.data.frame(t(as.matrix(SNP.1)))
colnames(SNP.1) <- "gt"
SNP.1$IdX <- rownames(SNP.1)
SNP.1$IdX <- gsub("^X", "", rownames(SNP.1))

pt_SNP.1 <- left_join(avg_AUDPC,SNP.1, by = "IdX")
pt_SNP.1 <- na.omit(pt_SNP.1)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.1$gt_round <- round(pt_SNP.1$gt)

# Assuming your data frame is named df
ggplot(pt_SNP.1, aes(x = as.factor(gt_round), y = AUDPC_adj)) +
  geom_boxplot() +
  labs(#title = "Box Plot of SKT by dosage",
       x = "Genotype",
       y = "AUDPC") +
  theme_few()
ggsave("plots/GWAS_RESPOT_allchr_SKT_spatial_BLUE_dosage_SNP1_dom.png", height = 3, width = 3)



SNP.2 <- geno[geno$SNP == QTL.snps$Marker[2],-ncol(geno)]

SNP.2 <- SNP.2[,-c(1:3)]

SNP.2 <- as.data.frame(t(as.matrix(SNP.2)))
colnames(SNP.2) <- "gt"
SNP.2$IdX <- rownames(SNP.2)
SNP.2$IdX <- gsub("^X", "", rownames(SNP.2))

pt_SNP.2 <- left_join(avg_AUDPC,SNP.2, by = "IdX")
pt_SNP.2 <- na.omit(pt_SNP.2)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.2$gt_round <- round(pt_SNP.2$gt)

# Assuming your data frame is named df
ggplot(pt_SNP.2, aes(x = as.factor(gt_round), y = AUDPC_adj)) +
  geom_boxplot() +
  labs(#title = "Box Plot of SKT by dosage",
       x = "Genotype",
       y = "AUDPC") +
  theme_few()
ggsave("plots/GWAS_RESPOT_allchr_SKT_spatial_BLUE_dosage_SNP2_dom.png", height = 3, width = 3)


SNP.3 <- geno[geno$SNP == QTL.snps$Marker[3],-ncol(geno)]

SNP.3 <- SNP.3[,-c(1:3)]

SNP.3 <- as.data.frame(t(as.matrix(SNP.3)))
colnames(SNP.3) <- "gt"
SNP.3$IdX <- rownames(SNP.3)
SNP.3$IdX <- gsub("^X", "", rownames(SNP.3))

pt_SNP.3 <- left_join(avg_AUDPC,SNP.3, by = "IdX")
pt_SNP.3 <- na.omit(pt_SNP.3)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.3$gt_round <- round(pt_SNP.3$gt)

# Mean phenotype for each allele dosage
mean_phenotype_SNP.3 <- pt_SNP.3 %>%
  group_by(gt_round) %>%
  summarise(
    AUDPC_adj = mean(AUDPC_adj, na.rm = TRUE),
    count = n()
  )

# Assuming your data frame is named df
p3 <- ggplot(pt_SNP.3, aes(x = as.factor(gt_round), y = AUDPC_adj)) +
  geom_boxplot() +
  labs(title = "QTL3",
       x = "Genotype",
       y = "AUDPC") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + 
  theme_few() + geom_text(data = mean_phenotype_SNP.3,
                          aes(x = as.factor(gt_round), y = Inf, label = round(AUDPC_adj, 2)),
                          vjust = 1.5, size = 3, color = "black")+
  # Add count values below the boxplots
  geom_text(data = mean_phenotype_SNP.3,
            aes(x = as.factor(gt_round), y = -Inf, label = count),
            vjust = -1, size = 3, color = "black")

p3
ggsave("plots/GWAS_RESPOT_allchr_SKT_spatial_BLUE_dosage_SNP3_dom.png", height = 3, width = 3)
save(p3, file = "plots/GWAS_RESPOT_allchr_SKT_spatial_BLUE_dosage_SNP3_dom.RData")


##########
