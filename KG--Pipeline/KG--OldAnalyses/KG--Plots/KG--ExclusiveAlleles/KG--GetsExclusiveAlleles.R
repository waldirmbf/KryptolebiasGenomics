### The BEGINNING ~~~~~
##
# ~ Gets numbers of exclusive alleles based on a .geno file | By George Pacheco.


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


## 1) How many of the 5,688 SNPs are homozygous and fixed (either REF or ALT) across all 13 Kher individuals?


# Fills HomoFixedKher_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixedKher_AlleleState"] <-
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) | 
          all(Geno[ROW, c(9:10, 12:22)] == 2), "HomoFixed_Kher",
          "Other")}


# Gets numbers for HomoFixedKher_AlleleState ~
table(Geno$HomoFixedKher_AlleleState)


# Creates HomoFixedKher dataframe ~
HomoFixedKher_df <- Geno[ grep("HomoFixed_Kher", Geno$HomoFixedKher_AlleleState), ]


# Expands HomoFixedKher_df by adding the AlleleState columns ~
HomoFixedKher_df$HomoFixed_KherKspESP_AlleleState <- rep(NA, nrow(HomoFixedKher_df))
HomoFixedKher_df$HomoFixed_Kher_Het_KspESP_AlleleState <- rep(NA, nrow(HomoFixedKher_df))


## 2) How many of the 4,976 SNPS that are fixed in Kher are also fixed in KspESP?

# Fills HomoFixed_KherKspESP_AlleleState ~
for(ROW in 1:nrow(HomoFixedKher_df)){
  HomoFixedKher_df[ROW, "HomoFixed_KherKspESP_AlleleState"] <-
    ifelse(all(HomoFixedKher_df[ROW, c(5:8, 11)] == 0) & all(HomoFixedKher_df[ROW, c(9:10, 12:22)] == 2) | 
           all(HomoFixedKher_df[ROW, c(5:8, 11)] == 2) & all(HomoFixedKher_df[ROW, c(9:10, 12:22)] == 0), "HomoFixed_KherKspESP",
           "Other")}


# Gets numbers for HomoFixed_KherKspESP_AlleleState ~
table(HomoFixedKher_df$HomoFixed_KherKspESP_AlleleState)


## 3) How many of the 4,976 SNPs fixed in Kher are heterozygous in KspESP?

# Fills HomoFixed_Kher_Het_KspESP_AlleleState ~
for(ROW in 1:nrow(HomoFixedKher_df)){
  HomoFixedKher_df[ROW, "HomoFixed_Kher_Het_KspESP_AlleleState"] <-
    ifelse(all(HomoFixedKher_df[ROW, c(9:10, 12:22)] == 0) & all(HomoFixedKher_df[ROW, c(5:8, 11)] == 1) |
           all(HomoFixedKher_df[ROW, c(9:10, 12:22)] == 2) & all(HomoFixedKher_df[ROW, c(5:8, 11)] == 1), "HomoFixed_Kher_Het_KspESP",
           "Other")}


# Gets numbers for HomoFixed_Kher_Het_KspESP_AlleleState ~
table(HomoFixedKher_df$HomoFixed_Kher_Het_KspESP_AlleleState)


# Expands Geno by adding further AlleleState columns ~
Geno$ExclusiveToKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$ExclusiveToKher_AlleleState <- rep(NA, nrow(Geno))
Geno$HomoFixedKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$HetKspESP_AlleleState <- rep(NA, nrow(Geno))
Geno$HetKher_AlleleState <- rep(NA, nrow(Geno))
Geno$HetKspESP_Kher_AlleleState <- rep(NA, nrow(Geno))


## 4) How many of the 5,688 SNPs are homozygous and fixed (either REF or ALT) across all 5 KspESP individuals?


# Fills HomoFixedKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "HomoFixedKspESP_AlleleState"] <-
   ifelse(all(Geno[ROW, c(5:8, 11)] == 0) | 
         all(Geno[ROW, c(5:8, 11)] == 2), "HomoFixed_KspESP",
         "Other")}


# Gets numbers for HomoFixedKspESP_AlleleState ~
table(Geno$HomoFixedKspESP_AlleleState)


## 5) How many SNPs contain alleles that are exclusive to KspESP?


# Fills ExclusiveToKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
  Geno[ROW, "ExclusiveToKspESP_AlleleState"] <-
   ifelse(all(Geno[ROW, c(5:8, 11)] == 0) & all(Geno[ROW, c(9:10, 12:22)] == 2) |
          all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 2), "ExclusiveToKspESP_REF",
   ifelse(all(Geno[ROW, c(5:8, 11)] == 2) & all(Geno[ROW, c(9:10, 12:22)] == 0) |
          all(Geno[ROW, c(5:8, 11)] == 1) & all(Geno[ROW, c(9:10, 12:22)] == 0), "ExclusiveToKspESP_ALT",
   "Other"))}


# Gets numbers for ExclusiveToKspESP_AlleleState ~
table(Geno$ExclusiveToKspESP_AlleleState)


## 6) How many SNPs contain alleles that are exclusive to Kher?


# Fills ExclusiveToKher_AlleleState ~
for(ROW in 1:nrow(Geno)){
 Geno[ROW, "ExclusiveToKher_AlleleState"] <-
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 0) & all(Geno[ROW, c(5:8, 11)] == 2) |
          all(Geno[ROW, c(9:10, 12:22)] == 1) & all(Geno[ROW, c(5:8, 11)] == 2), "ExclusiveToKher_REF",
   ifelse(all(Geno[ROW, c(9:10, 12:22)] == 2) & all(Geno[ROW, c(5:8, 11)] == 0) |
          all(Geno[ROW, c(9:10, 12:22)] == 1) & all(Geno[ROW, c(5:8, 11)] == 0), "ExclusiveToKher_ALT",
   "Other"))}


# Gets numbers for ExclusiveToKspESP_AlleleState ~
table(Geno$ExclusiveToKher_AlleleState)


## 7) How many SNPs are heterozygous in KspESP?


# Fills HetKspESP_AlleleState ~
for(ROW in 1:nrow(Geno)){
 Geno[ROW, "HetKspESP_AlleleState"] <-
   ifelse(all(Geno[ROW, c(5:8, 11)] == 1), "HetKspESP",
          "Other")}


# Gets numbers for HetKspESP_AlleleState ~
table(Geno$HetKspESP_AlleleState)


## 8) How many SNPs are heterozygous in Kher?


# Fills HetKher_AlleleState ~
for(ROW in 1:nrow(Geno)){
   Geno[ROW, "HetKher_AlleleState"] <-
    ifelse(all(Geno[ROW, c(9:10, 12:22)] == 1), "HetKher",
           "Other")}


# Gets numbers for HetKher_AlleleState ~
table(Geno$HetKher_AlleleState)


## 9) How many SNPs are heterozygous in KspESP & Kher?


# Fills HetKspESP_Kher_AlleleState ~
for(ROW in 1:nrow(Geno)){
 Geno[ROW, "HetKspESP_Kher_AlleleState"] <-
    ifelse(all(Geno[ROW, c(5:8, 11)] == 1) &
           all(Geno[ROW, c(9:10, 12:22)] == 1), "HetKspESP_Kher",
           "Other")}


# Gets numbers for HetKspESP_Kher_AlleleState ~
table(Geno$HetKspESP_Kher_AlleleState)


#
##
### The END ~~~~~