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
library(ggplot2)
#library(ggrepel)
#library(ggthemes)
library(cowplot)
library(GWASpoly)
library(AGHmatrix)
library(qqman)

###################################


###################################
###### Load data files       ######
###################################


load("data/phenotypes/spatial_matrix_scaled_FRS.RData")

load("data/phenotypes/FRS_final_pt_spatial.RData")
pt_FRS <- pt_sp

ID_kept_yeargeno <- read.csv("~/RESPOT/GWAS_drone/data/original_danespo_files/genotypes/Helene_QC_imputed/ID_kept_yeargeno", sep="")
colnames(ID_kept_yeargeno) <- c("IdX", "Rank", "Year_geno")

pt_FRS <- left_join(pt_FRS, ID_kept_yeargeno, by = "IdX")

load("data/genotypes/FRS_final_gt.RData")

load("data/genotypes/SNP_map_Geno_2020_to_2024_QC_imputed.RData")



###################################

# All SNPs are kept, but this is a check
SNP_map <- map_sub[map_sub$SNP %in% colnames(gt_FRS), ]

# Order the markers according to the ordered map and not the 19BetSol probe order
SNP <- SNP_map$SNP

gt_FRS_ordered <- gt_FRS[ ,SNP]


##################################
#   Standardize
##################################

# Based on the phenotype distribution plots, standardization of the FRS phenotype across 2023 and 2024 is not needed, they follow similar distributions

##################################



###################################
###### GWAS prep             ######
###################################

# Prepare geno.file for GWASpoly
gt_t <- data.frame(t(gt_FRS_ordered))

gt_t$SNP <- rownames(gt_t)

# Remove the "X" prefix from the column names (is automatically introduced in conversion of numerical rownames to column names)
colnames(gt_t) <- sub("^X", "", colnames(gt_t))

geno.file <- left_join(SNP_map, gt_t, by = "SNP")

# Remove IdX with missing marker segment

Missing_GRP <- unique(pt_FRS[is.na(pt_FRS$GRP_Overall), "IdX"])

geno.file2 <- geno.file[ , setdiff(names(geno.file), Missing_GRP)]

#write.csv(geno.file2, "results/GWASpoly_genofile_FRS_spatial.csv", row.names = F)

# Prepare pheno.file for GWASpoly

pheno.file <- pt_FRS[ ,c("IdX", "Emergence1", "TFL", "GRP_Overall", "spid", "Year_geno")]

# Remove IdX with missing marker segment

pheno.file2 <- pheno.file[!(pheno.file$IdX %in% Missing_GRP), ]

#write.csv(pheno.file2, "results/GWASpoly_phenofile_FRS_spatial.csv", row.names = F)


###################################

##############################################################
###### Calculate BLUEs for the genotype effect          ######
##############################################################

library(sommer)

# Remove spid without marker segment in the spatial matrix

pheno.file2$IdX <- as.factor(as.character(pheno.file2$IdX))
pheno.file2$IdX <- relevel(pheno.file2$IdX, ref = "1")
pheno.file2$TFL <- as.factor(as.character(pheno.file2$TFL))
pheno.file2$GRP_Overall <- as.factor(as.character(pheno.file2$GRP_Overall))
pheno.file2$spid <- as.factor(as.character(pheno.file2$spid))
pheno.file2$Year_geno <- as.factor(as.character(pheno.file2$Year_geno))

# spid with missing marker segment
spid <- pheno.file[pheno.file$IdX %in% Missing_GRP, "spid"]

RR_red <- RR_scaled[!rownames(RR_scaled) %in% spid, !colnames(RR_scaled) %in% spid] 

# fit the LMM model to correct for spatial and fixed effects at the same time
mod <- mmer(Emergence1 ~ 1 + IdX,
            random = ~ vsr(spid, Gu = RR_red) + TFL, # no random line effect, as replicates are accommodated in GWASpoly automatically
            data = pheno.file2)

mod_summary <- summary(mod)


# Extract the coefficients of the fixed effects and the variance component of the spatial effect for downstream analysis
vc <- mod_summary$varcomp[, "VarComp"]
se <- mod_summary$varcomp[, "VarCompSE"]
varcomp_df <- data.frame(Variance = vc, SE = se)


# Extract the BLUEs for genotype to include in GWAS

BLUE <- mod$Beta
BLUE$Effect <- as.character(BLUE$Effect)
# BLUE of IdX1 is the intercept (reference level) - add the intercept to all other coefficients, to include that in analysis
BLUE$Effect[1] <- "IdX1"
BLUE <- BLUE[,2:3]
BLUE$Estimate[2:nrow(BLUE)] <- BLUE$Estimate[2:nrow(BLUE)] + BLUE$Estimate[1]

BLUE$Effect <- gsub("IdX", "", BLUE$Effect)
colnames(BLUE) <- c("IdX", "BLUE")
BLUE <- BLUE[2:nrow(BLUE),]

write.csv(BLUE, "results/GWASpoly_phenofile_FRS_spatial_BLUE_fixchr5_1marker.csv", row.names = F)



# Fix chr91 SNP (include in phenofile)

snp <- as.data.frame(geno.file2[geno.file2$SNP == "AX.339335105", 4:ncol(geno.file2)])

snp_transposed <- as.data.frame(t(snp))

snp_transposed <- data.frame(Variable = rownames(snp_transposed), Value = snp_transposed[,1], row.names = NULL)

colnames(snp_transposed) <- c("IdX", "AX.339335105")

# replace missing values with the most frequent genotype (same imputation as GWASpoly) - a.k.a. the mode genotype

mode <- as.numeric(names(sort(table(snp_transposed$AX.339335105), decreasing = TRUE)[1]))
snp_transposed$AX.339335105[is.na(snp_transposed$AX.339335105)] <- mode

# save the new phenotype file

pt <- left_join(BLUE, snp_transposed, by = "IdX")



# remove the SNP from the geno.file2 and save it

geno.file2 <- geno.file2[!geno.file2$SNP == "AX.339335105", ]



write.csv(geno.file2, "results/GWASpoly_genofile_FRS_spatial_BLUE_fixchr5_1marker.csv", row.names = F)

write.csv(pt, "results/GWASpoly_phenofile_FRS_spatial_BLUE_fixchr5_1marker.csv", row.names = F)

##################################



###################################
###### GWAS                  ######
###################################

# Read data as GWASpoly.data element

genofile <- "results/GWASpoly_genofile_FRS_spatial_BLUE_fixchr5_1marker.csv"

phenofile <- "results/GWASpoly_phenofile_FRS_spatial_BLUE_fixchr5_1marker.csv" # Update to the corrected phenotypes (residuals of LMM model)

data <- read.GWASpoly(ploidy = 4, 
                      pheno.file = phenofile,
                      geno.file = genofile,
                      format = "numeric",
                      n.traits = 1,
                      delim = ",")

# Create Leave-One-Chromosome_Out GRMs 

data.loco <- set.K(data,
                   LOCO=TRUE,
                   n.core=4)


# Test markers for significance

N <- nrow(gt_FRS) #Population size
params <- set.params(geno.freq = 1 - 5/N,
                     fixed = c("AX.339335105"),
                     fixed.type = c("factor"),
                     n.PC = 2) # updated as part of revision from 3 to 2

data.loco.scan <- GWASpoly(data=data.loco,
                           models=c("additive"),
                           traits=c("BLUE"),
                           params=params,
                           n.core=4)

save(data.loco.scan, file = "results/GWAS_FRS_result_spatial_BLUE_fixchr5_1marker.RData")

qq.plot(data.loco.scan,trait="BLUE") + ggtitle(label="Emergence1_BLUE") + theme(legend.position = "none")

ggsave("plots/GWAS_FRS_qq_loco_spatial_BLUE_fixchr5_1marker.png", height = 8, width = 8)


data2 <- set.threshold(data.loco.scan,method="Bonferroni",level=0.05)

#data3 <- set.threshold(data.loco.scan,method="M.eff",level=0.05) # Is not much different from the Bonferroni in this case, but a bit less conservative

#p <- manhattan.plot(data2,traits="BLUE")
#p <- p + 
#  geom_point(shape = 19, aes(color = as.factor(data2@map$Chrom)), show.legend = FALSE) +  # Full points and no legend
#  scale_color_manual(values = rep(c("darkgray", "chartreuse3"), length.out = length(unique(data2@map$Chrom)))) +  # Alternating colors
#  theme_classic() +  # Classic theme
#  theme(axis.text.x = element_text(vjust = 0.5),  # Adjust axis text
#        legend.position = "none")  
#
#p
#
#ggsave("plots/GWAS_FRS_manhattan_bonferroni_spatial_BLUE_fixchr5_1marker.png", height = 8, width = 12)

###################################



###################################
###### LD plot               ######
###################################

p <- LD.plot(data2, max.loci=1000)

p + xlim(0,30) 

ggsave("plots/GWAS_FRS_LDplot_spatial_BLUE_fixchr5_1marker.png", height = 8, width = 12)

###################################


######################################
###### Genomic inflation factor ######
######################################

results <- data.frame(data.loco.scan@scores)
results$Marker <- rownames(results)
results$p <- 10^(-results$additive)

map <- data.frame(data.loco.scan@map)

###### QQ plot of -log10(p) to evaluate for bias ######
# also using "qqman" package



png(filename = "plots/GWAS_RESPOT_allchr_FRS_QQprelambda_spatial_BLUE_fixchr5_1marker.png", width = 400, height = 400)
qq(results$p, 
   main = "Q-Q plot of GWAS p-values - FRS, pre-lambda_spatial_BLUE_fixchr5_1marker")
dev.off()

###### Calculate genomic inflation factor from p-values ######

chisq <- qchisq(results$p, 1, lower.tail = FALSE) # convert p-values to chi-square values, DF = 1
lambda.gc <- median(chisq)/qchisq(0.5,1) # calculate the genomic inflation factor as the obs median divided by the expected median (chi-square of 50th percentile, DF = 1)

save(lambda.gc, file = "results/lambda_gc_FRS_spatial_BLUE_fixchr5_1marker.RData")

# If it is > 1, there is an indication of systematic bias in the results, that needs to be corrected. 
# lambda in [1, 1.05] = OK, lambda in [1.1, 1.2] = suggest some level of inflation (bias sources, pop. strat.), lambda in [> 1.2] = there is indication of substatial biases, that  could lead to false positive findings. ALL IS DEPENDANT ON THE SPECIFIC CASE
# (Kristen et al., 2018) uses the genomic inflation factor to correct p-value results. 

### Correct p-value using genomic inflation factor

corr.chisq <- chisq/lambda.gc # correct the chi-square values 
results$Pcorr <- pchisq(corr.chisq, df = 1, lower.tail = FALSE) # convert the corrected chi-square values back to p-values and store in the GWAS_UU2 result matrix

###### QQ plot of -log(pcorr) to control bias correction ######

png(filename = "plots/GWAS_RESPOT_allchr_FRS_QQpostlambda_spatial_BLUE_fixchr5_1marker.png", width = 400, height = 400)
qq(results$Pcorr, 
   main = "Q-Q plot of GWAS p-values - FRS, post-lambda_spatial_BLUE_fixchr5_1marker")
dev.off()

####################################


######################################################
###### Manhattan plot of corrected GWAS results ######
######################################################

results_corr <- left_join(results, map, by = "Marker")
results_corr$Chrom <- as.numeric(as.character(results_corr$Chrom))

# Joined QQ plot

png(filename = "plots/GWAS_RESPOT_allchr_FRS_QQtotal_spatial_BLUE_fixchr5_1marker.png", width = 1000, height = 500)
par(mfrow = c(1,2))
qq(results_corr$p, main = "Q-Q plot of GWAS p-values - FRS, pre-lambda_spatial_BLUE_fixchr5_1marker", ylim = c(0, -log10(min(results_corr$p))))  
qq(results_corr$Pcorr, main = "Q-Q plot of GWAS p-values - FRS, post-lambda_spatial_BLUE_fixchr5_1marker", ylim = c(0, -log10(min(results_corr$p))))
dev.off()


bonferroni <- 0.05 / nrow(results)

#Meff <- data3@threshold

results_corr <- results_corr[!(results_corr$Chrom == 0 | results_corr$Chrom == 13), ]

png(file = "plots/GWAS_RESPOT_allchr_FRS_manhattan_bonferroni_Pcorr_spatial_BLUE_fixchr5_1marker.png", width = 1000, height = 500)
qqman::manhattan(results_corr, chr = "Chrom", bp = "Position", snp = "Marker", p = "Pcorr", suggestiveline = -log10(bonferroni), genomewideline = FALSE, col = c("darkgrey", "chartreuse3"), main = "Manhattan plot of FRS_Emergence1 GWAS, bonferroni corr_spatial_BLUE_fixchr5_1marker")
dev.off()

#png(file = "plots/GWAS_RESPOT_allchr_FRS_manhattan_Meff_Pcorr_spatial_BLUE.png", width = 1000, height = 500)
#qqman::manhattan(results_corr, chr = "Chrom", bp = "Position", snp = "Marker", p = "Pcorr", suggestiveline = Meff, genomewideline = FALSE, col = c("darkgrey", "chartreuse3"), main = "Manhattan plot of FRS_Emergence1 GWAS, Meff corr_spatial_BLUE")
#dev.off()

#############################################


#############################################
###### Select significant SNPs for SKT ######
#############################################

# Bonferroni correction
snp_select_bonferroni_FRS <- results_corr[results_corr$Pcorr < bonferroni, ] 

save(snp_select_bonferroni_FRS, file = "results/GWAS_RESPOT_FRS_allchr_bonf_sigSNP_spatial_BLUE_fixchr5_1marker.RData")

#############################################
























