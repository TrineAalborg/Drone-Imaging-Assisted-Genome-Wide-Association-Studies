##########################################

# Script: GWAS drone data - FRS

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

load("results/GWAS_FRS_result_spatial_BLUE.RData")

geno_original <- read.csv("results/GWASpoly_genofile_FRS_spatial.csv")

load("results/GWAS_RESPOT_FRS_allchr_bonf_sigSNP_spatial_BLUE_fixchr2_chr5_2marker.RData")

qtl8 <- snp_select_bonferroni_FRS

load("results/GWAS_RESPOT_FRS_allchr_bonf_sigSNP_spatial_BLUE.RData")

snp_select_bonferroni_FRS <- rbind(snp_select_bonferroni_FRS, qtl8)

pt <- read.csv("results/GWASpoly_phenofile_FRS_spatial_BLUE.csv")

load("data/genotypes/FRS_final_gt.RData")
# Remove IdX not present in pt (those with missing market segment)
gt <- gt_FRS[rownames(gt_FRS) %in% as.character(unique(pt$IdX)), ]

load("results/GWASpoly_sommer_varcomp_FRS_spatial_BLUE.RData")


###################################


###################################
###### Calculate MAF         ######
###################################

geno <- geno_original

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
sig <- left_join(snp_select_bonferroni_FRS, beta, by = "Marker")

# Add MAF to data frame
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

save(sig, file = "results/GWAS_RESPOT_FRS_allchr_bonf_sigSNP_SNPeffect_spatial_BLUE.RData")

library(writexl)
write_xlsx(sig, "plots/GWAS_RESPOT_FRS_allchr_bonf_sigSNP_SNPeffect_spatial_BLUE.xlsx")

###################################


###################################
###### SNP dosage plot       ######
###################################

# Plot the frequency of each dosage for the QTL SNPs against the phenotype

pt <- read.csv("results/GWASpoly_phenofile_FRS_spatial.csv")

pt <- pt[,c("IdX", "Emergence1")]
pt$IdX <- as.character(pt$IdX)

#QTL.snps <- snp_select_bonferroni_FRS %>% group_by(as.factor(Chrom)) %>% filter(p == min(p)) %>% summarise(Marker = Marker[1]) 
QTL.snps <- snp_select_bonferroni_FRS[order(snp_select_bonferroni_FRS$Pcorr), ] # use as there is two QTLs on chrom 5

SNP.1 <- geno[geno$SNP == "AX.339298605",-ncol(geno)]

SNP.1 <- SNP.1[,-c(1:3)]

SNP.1 <- as.data.frame(t(as.matrix(SNP.1)))
colnames(SNP.1) <- "gt"
SNP.1$IdX <- rownames(SNP.1)
SNP.1$IdX <- gsub("^X", "", rownames(SNP.1))

pt_SNP.1 <- left_join(pt,SNP.1, by = "IdX")
pt_SNP.1 <- na.omit(pt_SNP.1)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.1$gt_round <- round(pt_SNP.1$gt)

# Mean phenotype for each allele dosage
mean_phenotype_SNP.1 <- pt_SNP.1 %>%
  group_by(gt_round) %>%
  summarise(
    Emergence1 = mean(Emergence1, na.rm = TRUE),
    count = n()
  )

# Assuming your data frame is named df
p1 <- ggplot(pt_SNP.1, aes(x = as.factor(gt_round), y = Emergence1)) +
  geom_boxplot() +
  labs(#title = "Box Plot of Emergence by dosage",
       x = "Genotype",
       y = "Emergence") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + 
  theme_few() + geom_text(data = mean_phenotype_SNP.1,
                          aes(x = as.factor(gt_round), y = Inf, label = round(Emergence1, 2)),
                          vjust = 1.5, size = 3, color = "black")+
  # Add count values below the boxplots
  geom_text(data = mean_phenotype_SNP.1,
            aes(x = as.factor(gt_round), y = -Inf, label = count),
            vjust = -1, size = 3, color = "black")
p1
ggsave("plots/GWAS_RESPOT_allchr_FRS_spatial_BLUE_dosage_SNP1.png", height = 3, width = 3)





########################

SNP.2 <- geno[geno$SNP == "AX.339335105",-ncol(geno)] # QTL2 is the third lowest p-value

SNP.2 <- SNP.2[,-c(1:3)]

SNP.2 <- as.data.frame(t(as.matrix(SNP.2)))
colnames(SNP.2) <- "gt"
SNP.2$IdX <- rownames(SNP.2)
SNP.2$IdX <- gsub("^X", "", rownames(SNP.2))

pt_SNP.2 <- left_join(pt,SNP.2, by = "IdX")
pt_SNP.2 <- na.omit(pt_SNP.2)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.2$gt_round <- round(pt_SNP.2$gt)

# Mean phenotype for each allele dosage
mean_phenotype_SNP.2 <- pt_SNP.2 %>%
  group_by(gt_round) %>%
  summarise(
    Emergence1 = mean(Emergence1, na.rm = TRUE),
    count = n()
  )

# Assuming your data frame is named df
p2 <- ggplot(pt_SNP.2, aes(x = as.factor(gt_round), y = Emergence1)) +
  geom_boxplot() +
  labs(#title = "Box Plot of Emergence by dosage",
       x = "Genotype",
       y = "Emergence") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + 
  theme_few() + geom_text(data = mean_phenotype_SNP.2,
                          aes(x = as.factor(gt_round), y = Inf, label = round(Emergence1, 2)),
                          vjust = 1.5, size = 3, color = "black")+
  # Add count values below the boxplots
  geom_text(data = mean_phenotype_SNP.2,
            aes(x = as.factor(gt_round), y = -Inf, label = count),
            vjust = -1, size = 3, color = "black")
p2
ggsave("plots/GWAS_RESPOT_allchr_FRS_spatial_BLUE_dosage_SNP2.png", height = 3, width = 3)


########################

SNP.3 <- geno[geno$SNP == "AX.339337515",-ncol(geno)] # QTL3 is the second lowest p-value

SNP.3 <- SNP.3[,-c(1:3)]

SNP.3 <- as.data.frame(t(as.matrix(SNP.3)))
colnames(SNP.3) <- "gt"
SNP.3$IdX <- rownames(SNP.3)
SNP.3$IdX <- gsub("^X", "", rownames(SNP.3))

pt_SNP.3 <- left_join(pt,SNP.3, by = "IdX")
pt_SNP.3 <- na.omit(pt_SNP.3)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.3$gt_round <- round(pt_SNP.3$gt)

# Mean phenotype for each allele dosage
mean_phenotype_SNP.3 <- pt_SNP.3 %>%
  group_by(gt_round) %>%
  summarise(
    Emergence1 = mean(Emergence1, na.rm = TRUE),
    count = n()
  )

# Assuming your data frame is named df
p3 <- ggplot(pt_SNP.3, aes(x = as.factor(gt_round), y = Emergence1)) +
  geom_boxplot() +
  labs(#title = "Box Plot of Emergence by dosage",
    x = "Genotype",
    y = "Emergence") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + 
  theme_few() + geom_text(data = mean_phenotype_SNP.3,
                          aes(x = as.factor(gt_round), y = Inf, label = round(Emergence1, 2)),
                          vjust = 1.5, size = 3, color = "black")+
  # Add count values below the boxplots
  geom_text(data = mean_phenotype_SNP.3,
            aes(x = as.factor(gt_round), y = -Inf, label = count),
            vjust = -1, size = 3, color = "black")
p3
ggsave("plots/GWAS_RESPOT_allchr_FRS_spatial_BLUE_dosage_SNP3.png", height = 3, width = 3)



########################

SNP.4 <- geno[geno$SNP == "AX.339312238",-ncol(geno)] # QTL3 is the second lowest p-value

SNP.4 <- SNP.4[,-c(1:3)]

SNP.4 <- as.data.frame(t(as.matrix(SNP.4)))
colnames(SNP.4) <- "gt"
SNP.4$IdX <- rownames(SNP.4)
SNP.4$IdX <- gsub("^X", "", rownames(SNP.4))

pt_SNP.4 <- left_join(pt,SNP.4, by = "IdX")
pt_SNP.4 <- na.omit(pt_SNP.4)

library(ggplot2)
library(ggthemes)

# Imputation has generated non-integer genotypes. For simplicity, these are rounded to closest whole number
pt_SNP.4$gt_round <- round(pt_SNP.4$gt)

# Mean phenotype for each allele dosage
mean_phenotype_SNP.4 <- pt_SNP.4 %>%
  group_by(gt_round) %>%
  summarise(
    Emergence1 = mean(Emergence1, na.rm = TRUE),
    count = n()
  )

# Assuming your data frame is named df
p4 <- ggplot(pt_SNP.4, aes(x = as.factor(gt_round), y = Emergence1)) +
  geom_boxplot() +
  labs(#title = "Box Plot of Emergence by dosage",
    x = "Genotype",
    y = "Emergence") +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + 
  theme_few() + geom_text(data = mean_phenotype_SNP.4,
                          aes(x = as.factor(gt_round), y = Inf, label = round(Emergence1, 2)),
                          vjust = 1.5, size = 3, color = "black")+
  # Add count values below the boxplots
  geom_text(data = mean_phenotype_SNP.4,
            aes(x = as.factor(gt_round), y = -Inf, label = count),
            vjust = -1, size = 3, color = "black")
p4
ggsave("plots/GWAS_RESPOT_allchr_FRS_spatial_BLUE_dosage_SNP4.png", height = 3, width = 3)




##########

ggsave("plots/GWAS_RESPOT_allchr_FRS_spatial_BLUE_dosage_SNP1SNP2SNP3SNP4.png", 
       plot_grid(p1 + ggtitle("QTL5"), p2  + ggtitle("QTL6"), p3 + ggtitle("QTL7"), p4 + ggtitle("QTL8"), ncol = 4), 
       height = 2.5, width = 10, dpi = 300)
