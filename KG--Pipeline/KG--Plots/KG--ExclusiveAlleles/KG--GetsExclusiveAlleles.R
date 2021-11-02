### The BEGINNING ~~~~~
##
# ~ Gets numbers of exclusive alleles based on a .geno file| By George Pacheco.


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


# Expands Geno by adding the AlleleState columns ~
Geno$HomoFixedKher_AlleleState <- rep(NA, nrow(Geno))
Geno$HomoFixedKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$HomoFixed_KherKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$HomoFixed_Kher_Het_KspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$ExclusiveToKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$ExclusiveToKher_AlleleState <- rep(NA, nrow(Geno))


## 1) How many of the 5,688 SNPs are homozygous and fixed (either REF or ALT) across all 13 Kher individuals?


# Fills HomoFixedKher_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixedKher_AlleleState"] <-
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) | 
          all(Geno[ROW, c(9:10, 12:22)] == 2), "HomoFixed_Kher",
          "Other")}

# Gets numbers for HomoFixedKher_AlleleState ~
table(Geno$HomoFixedKher_AlleleState)


## 2) How many of the 5,688 SNPs are homozygous and fixed (either REF or ALT) across all 5 KspESP individuals?


# Fills HomoFixedKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixedKspESP_AlleleState"] <-
  ifelse(all(Geno[ROW, c(5:8, 11)] == 0) | 
         all(Geno[ROW, c(5:8, 11)] == 2), "HomoFixed_KspESP",
         "Other")}

# Gets numbers for HomoFixedKspESP_AlleleState ~
table(Geno$HomoFixedKspESP_AlleleState)


## 3) How many of the 4,976 fixed in Kher is also fixed (either REF or ALT) across all 5 KspESP individuals?

# Fills HomoFixed_KherKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixed_KherKspESP_AlleleState"] <-
    ifelse(all(Geno[ROW, c(5:8, 11)] == 0) & all(Geno[ROW, c(9:10, 12:22)] == 2) | 
           all(Geno[ROW, c(5:8, 11)] == 2) & all(Geno[ROW, c(9:10, 12:22)] == 0), "HomoFixed_KherKspESP",
           "Other")}

# Gets numbers for HomoFixed_KherKspESP_AlleleState ~
table(Geno$HomoFixed_KherKspESP_AlleleState)


## 4) How many of the 4,976 fixed in Kher (either REF or ALT) is heterozygous in KspESP?

# Fills HomoFixed_Kher_Het_KspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixed_Kher_Het_KspESP_AlleleState"] <-
    ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) & all(Geno[ROW, c(5:8, 11)] == 1) |
           all(Geno[ROW, c(9:10, 12:22)] == 2) & all(Geno[ROW, c(5:8, 11)] == 1), "HomoFixed_Kher_Het_KspESP",
           "Other")}


# Gets numbers for HomoFixed_Kher_Het_KspESP_AlleleState ~
table(Geno$HomoFixed_Kher_Het_KspESP_AlleleState)


## 5) How many SNPs contain alleles that are exclusive to KspESP?


# Fills ExclusiveToKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "ExclusiveToKspESP_AlleleState"] <-
   ifelse(all(Geno[ROW, c(5:8, 11)] == 0) | all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 2), "ExclusiveToKspESP",
   ifelse(all(Geno[ROW, c(5:8, 11)] == 2) | all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 0), "ExclusiveToKspESP",
   "Other"))}


# Gets numbers for ExclusiveToKspESP_AlleleState ~
table(Geno$ExclusiveToKspESP_AlleleState)


## 6) How many SNPs contain alleles that are exclusive to KspESP?


# Fills ExclusiveToKher_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "ExclusiveToKher_AlleleState"] <-
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) | all(Geno[ROW, c(9:10, 12:22)] == 1) & all(Geno[ROW, c(5:8, 11)] == 2), "ExclusiveToKher",
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 2) | all(Geno[ROW, c(9:10, 12:22)] == 1) & all(Geno[ROW, c(5:8, 11)] == 0), "ExclusiveToKher",
   "Other"))}


# Gets numbers for ExclusiveToKher_AlleleState ~
table(Geno$ExclusiveToKher_AlleleState)


#
##
### The END ~~~~~