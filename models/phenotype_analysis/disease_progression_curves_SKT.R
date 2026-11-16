##########################################

# Script: GWAS drone data - late blight 2022-2023 - Disease progression curves

# Trine Aalborg
# Version: 2024-09-02

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
library(readxl)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(reshape2)
library(agricolae)

###################################


###################################
###### Load data files       ######
###################################

load("data/phenotypes/Late_Blight_2022.RData")
load("data/phenotypes/Late_Blight_2023.RData")
load("data/phenotypes/Late_Blight_2024.RData")

###################################

###################################
###### Recode 1-9 to 9-1     ######
###################################

# It is recoded to create a AUDPC instead of a AUHPC

# 2022
Late_blight_2022[, 8:11] <- lapply(Late_blight_2022[, 8:11], function(column) {
  10 - column
})

# 2023
Late_blight_2023[, 9:14] <- lapply(Late_blight_2023[, 9:14], function(column) {
  10 - column
})

# 2024
Late_blight_2024[, 10:13] <- lapply(Late_blight_2024[, 10:13], function(column) {
  10 - column
})

###################################

###################################
###### Calculate AUDPC       ######
###################################

# Function to calculate AUDPC using agricolae
calculate_audpc <- function(scores, times) {
  audpc(scores, times, type = "relative") # Use relative AUDPC
}

# 2022
df_22 <- Late_blight_2022[, c(1, 8:11)]
colnames(df_22) <- c("Lbnr", "Time1", "Time2", "Time3", "Time4")

# Dates of the time points
date_points_22 <- as.Date(c("2022-07-18", "2022-08-05", "2022-08-09", "2022-08-11"))

# Convert dates to numeric values: days since the first date
numeric_time_points_22 <- as.numeric(date_points_22 - date_points_22[1])

# Apply the relative AUDPC calculation
df_22$AUDPC <- apply(df_22[, -1], 1, function(scores) calculate_audpc(scores, numeric_time_points_22))

# Join to full late blight data
Late_blight_2022_AUDPC <- left_join(Late_blight_2022, df_22, by = "Lbnr")


# 2023
df_23 <- Late_blight_2023[, c(1, 9:14)]
colnames(df_23) <- c("Lbnr", "Time1", "Time2", "Time3", "Time4", "Time5", "Time6")

# Dates of the time points
date_points_23 <- as.Date(c("2023-07-12", "2023-07-31", "2023-08-02", "2023-08-03", "2023-08-04", "2023-08-11"))

# Convert dates to numeric values: days since the first date
numeric_time_points_23 <- as.numeric(date_points_23 - date_points_23[1])

# Apply the relative AUDPC calculation
df_23$AUDPC <- apply(df_23[, -1], 1, function(scores) calculate_audpc(scores, numeric_time_points_23))

# Join to full late blight data
Late_blight_2023_AUDPC <- left_join(Late_blight_2023, df_23, by = "Lbnr")


# 2024
df_24 <- Late_blight_2024[, c(1, 10:13)]
colnames(df_24) <- c("Lbnr", "Time1", "Time2", "Time3", "Time4")

# Dates of the time points
date_points_24 <- as.Date(c("2024-07-18", "2024-07-25", "2024-07-29", "2024-07-31"))

# Convert dates to numeric values: days since the first date
numeric_time_points_24 <- as.numeric(date_points_24 - date_points_24[1])

# Apply the relative AUDPC calculation
df_24$AUDPC <- apply(df_24[, -1], 1, function(scores) calculate_audpc(scores, numeric_time_points_24))

# Join to full late blight data
Late_blight_2024_AUDPC <- left_join(Late_blight_2024, df_24, by = "Lbnr")




###################################
###### Create output files   ######
###################################

save(Late_blight_2022_AUDPC, file = "data/phenotypes/Late_Blight_2022_AUDPC.RData")
save(Late_blight_2023_AUDPC, file = "data/phenotypes/Late_Blight_2023_AUDPC.RData")
save(Late_blight_2024_AUDPC, file = "data/phenotypes/Late_Blight_2024_AUDPC.RData")

###################################



