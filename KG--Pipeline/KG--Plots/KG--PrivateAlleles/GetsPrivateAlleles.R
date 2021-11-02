### The BEGINNING ~~~~~
##
# ~ Gets numbers of private alleles based on a .geno file| By George Pacheco


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(tidyverse, extrafont)


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Loads data ~
Geno <- read.table(gzfile("OnlyKherANDKsp.geno.gz"), header = FALSE)


# Expands Geno by adding AlleleState ~
Geno$KherAlleleState <- rep(NA, nrow(Geno))
Geno$KspESPAlleleState <- rep(NA, nrow(Geno))
Geno$QuickNumbers <- rep(NA, nrow(Geno))
Geno$QuickNumbers2 <- rep(NA, nrow(Geno))


# Fills KherAlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "KherAlleleState"] <-
    ifelse(all(Geno[ROW, c(5:8, 11)] == 2) & all(Geno[ROW, c(9:10, 12:22)] == 0) | 
           all(Geno[ROW, c(5:8, 11)] == 0) & all(Geno[ROW, c(9:10, 12:22)] == 2), "HomoFixedKspESP-Kher",
           "Other")}

# Gets numbers for KspESPAlleleState ~
table(Geno$KherAlleleState)


# Fills KspESPAlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "KspESPAlleleState"] <-
   ifelse(all(Geno[ROW, c(5:8, 11)] == 0) | all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 2), "KspESP_REF",
   ifelse(all(Geno[ROW, c(5:8, 11)] == 2) | all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 0), "KspESP_ALT",
   "Other"))}


# Fills KherAlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "KherAlleleState"] <-
    ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) & all(Geno[ROW, c(5:8, 11)] == 1) |
           all(Geno[ROW, c(9:10, 12:22)] == 2) & all(Geno[ROW, c(5:8, 11)] == 1), "HomoFixedKher-HetKspESP",
    "Other")}


# Gets numbers for KherAlleleState ~
table(Geno$KherAlleleState)


# Fills KherAlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "QuickNumbers"] <-
    ifelse(all(Geno[ROW, c(9:10, 12:22)] == 2) | all(Geno[ROW, c(9:10, 12:22)] == 0), "HomoFixedKher",       
    "Other")}


# Fills KherAlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "QuickNumbers2"] <-
    ifelse(all(Geno[ROW, c(5:8, 11)] == 0) | all(Geno[ROW, c(5:8, 11)] == 2) &
           all(Geno[ROW, c(9:10, 12:22)] == 0) | all(Geno[ROW, c(9:10, 12:22)] == 2),
    "HomoFixedKher-KspESP",       
    "Other")}


# Gets numbers for KspESPAlleleState ~
table(Geno$KspESPAlleleState)


# Gets numbers for KherAlleleState ~
table(Geno$QuickNumbers)


# Gets numbers for KherAlleleState ~
table(Geno$QuickNumbers2)



#
##
### The END ~~~~~