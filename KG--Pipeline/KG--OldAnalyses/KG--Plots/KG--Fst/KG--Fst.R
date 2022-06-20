
##### Package installation ####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("SeqArray")
install.packages("pacman")

##### Loading required packages ####
library(gdsfmt)
library(SNPRelate)
library(SeqArray)

##### Data Input ####
setwd("C:/Users/Waldir Miron/Desktop/IntrogressionAnalysis_tests/FST")

KryptosVCF<- system.file("extdata", "ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_VcfwithAD.vcf.gz", package="SNPRelate")

KryptosVCF <- "ES-Article--AllSamples_NoKbraNoKgra_WithWGSs_SNPs_VcfwithAD.vcf.gz"

# Parsing VCF as a ".gds" object
snpgdsVCF2GDS(KryptosVCF, "KryptosVCF.gds", method="biallelic.only")
snpgdsSummary("KryptosVCF.gds")
genofile <- snpgdsOpen("KryptosVCF.gds")
genofile



##### Fst analysis considering KherSOUTH and KherCENTRAL different species ####
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
## Input population codes per sample ##
pop_code <- scan("KryptosIDs.txt", what=character())
pop_code


##### KspES vs KherSouth #####

flag <- pop_code %in% c("KspES", "KherSOUTH")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KspESvsKherSOUTH_FST <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
               method="W&C84")

# Weir and Cockerham weighted Fst estimate
KspESvsKherSOUTH_FST$Fst

# Weir and Cockerham mean Fst estimate
KspESvsKherSOUTH_FST$MeanFst


#### KspES vs KherCentral ####

flag <- pop_code %in% c("KspES", "KherCENTRAL")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KspESvsKherCENTRAL_FST <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                                  method="W&C84")
# Weir and Cockerham weighted Fst estimate
KspESvsKherCENTRAL_FST$Fst

# Weir and Cockerham mean Fst estimate
KspESvsKherCENTRAL_FST$MeanFst

#### KspES vs Kmar ####

flag <- pop_code %in% c("KspES", "Kmar")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KspESvsKmar <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                                    method="W&C84")

# Weir and Cockerham weighted Fst estimate
KspESvsKmar$Fst



# Weir and Cockerham mean Fst estimate
KspESvsKmar$MeanFst

#### KspES vs Koce ####

flag <- pop_code %in% c("KspES", "Koce")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KspESvsKoce <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                         method="W&C84")

# Weir and Cockerham weighted Fst estimate
KspESvsKoce$Fst

# Weir and Cockerham mean Fst estimate
KspESvsKoce$MeanFst

#### KherSouth vs KherCentral ####

flag <- pop_code %in% c("KherSOUTH", "KherCENTRAL")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KherSOUTHvsKherCENTRAL <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                         method="W&C84")

# Weir and Cockerham weighted Fst estimate
KherSOUTHvsKherCENTRAL$Fst
# Weir and Cockerham mean Fst estimate
KherSOUTHvsKherCENTRAL$MeanFst

#### KherSouth vs Kmar ####

flag <- pop_code %in% c("KherSOUTH", "Kmar")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KherSOUTHvsKmar <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                                    method="W&C84")
# Weir and Cockerham weighted Fst estimate
KherSOUTHvsKmar$Fst
# Weir and Cockerham mean Fst estimate
KherSOUTHvsKmar$MeanFst

#### KherSouth vs Koce ####

flag <- pop_code %in% c("KherSOUTH", "Koce")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KherSOUTHvsKoce <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                             method="W&C84")
# Weir and Cockerham weighted Fst estimate
KherSOUTHvsKoce$Fst
# Weir and Cockerham mean Fst estimate
KherSOUTHvsKoce$MeanFst

#### KherCentral vs Kmar ####

flag <- pop_code %in% c("KherCENTRAL", "Kmar")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KherCENTRALvsKmar <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                             method="W&C84")
# Weir and Cockerham weighted Fst estimate
KherCENTRALvsKmar$Fst
# Weir and Cockerham mean Fst estimate
KherCENTRALvsKmar$MeanFst

#### KherCentral vs Koce ####

flag <- pop_code %in% c("KherCENTRAL", "Koce")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KherCENTRALvsKoce <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                               method="W&C84")
# Weir and Cockerham weighted Fst estimate
KherCENTRALvsKoce$Fst
# Weir and Cockerham mean Fst estimate
KherCENTRALvsKoce$MeanFst

#### Kmar vs Koce ####

flag <- pop_code %in% c("Kmar", "Koce")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KmarvsKoce <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                               method="W&C84")
# Weir and Cockerham weighted Fst estimate
KmarvsKoce$Fst
# Weir and Cockerham mean Fst estimate
KmarvsKoce$MeanFst

##### Fst analysis considering KherSOUTH and KherCENTRAL the same species ####
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
## Input population codes per sample ##
pop_code <- scan("KryptosIDs_4species.txt", what=character())
pop_code
##### KspES vs Kher #####

flag <- pop_code %in% c("KspES", "Kher")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KspESvsKher_FST <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                                  method="W&C84")
# Weir and Cockerham weighted Fst estimate
KspESvsKher_FST$Fst
# Weir and Cockerham mean Fst estimate
KspESvsKher_FST$MeanFst



##### Kher vs Kmar #####

flag <- pop_code %in% c("Kher", "Kmar")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KhervsKmar_FST <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                            method="W&C84")
# Weir and Cockerham weighted Fst estimate
KhervsKmar_FST$Fst
# Weir and Cockerham mean Fst estimate
KhervsKmar_FST$MeanFst
##### Kher vs Koce #####

flag <- pop_code %in% c("Kher", "Koce")
flag
sample.selelection <- sample.id[flag]
sample.selelection 
population.selection <- pop_code[flag]
population.selection 

KhervsKoce_FST <- snpgdsFst(genofile, sample.id=sample.selelection, population=as.factor(population.selection), autosome.only = FALSE,
                             method="W&C84")
# Weir and Cockerham weighted Fst estimate
KhervsKoce_FST$Fst
# Weir and Cockerham mean Fst estimate
KhervsKoce_FST$MeanFst

######## Heatmap plot ##########
##### Considering KherSOUTH and KherCENTRAL different species #####

pacman::p_load(pheatmap, tidyverse, reshape2)
setwd("C:/Users/Waldir Miron/Desktop/IntrogressionAnalysis_tests/FST")
data <- read.table("FstR_5species.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
data

# Adds column names 
colnames(data) <- c("Pop1", "Pop2","Weighted")
# Melts data ~
pops = union(data$Pop1, data$Pop2)
n = length(pops)

# Creates Fst-Sites matrix ~
Fst <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
for (i in 1:nrow(data)) {
  Fst[data[i, "Pop1"], data[i, "Pop2"]] = data[i, "Weighted"]
  Fst[data[i, "Pop2"], data[i, "Pop1"]] = data[i, "Weighted"]}

Fst

# Writes Fst-Sites matrix ~
write.table(Fst, "Fst-Matrix.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Creates & Saves the heatmap ~
pheatmap(Fst, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", border_color = "black", cellwidth = 40, cellheight = 40,
         treeheight_row = 0, treeheight_col = 0, angle_col = "45", fontsize = 8, filename = "pairwiseWeightedFst_5species_R.pdf")
?pheatmap

# Creates & Saves the heatmap with only one diagonal ~
Fst[upper.tri(Fst)] <- NA
# Creates & Saves the heatmap ~
pheatmap(Fst, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", border_color = "black", cellwidth = 40, cellheight = 40,
         treeheight_row = 0, treeheight_col = 0, angle_col = "45", fontsize = 8,cluster_rows=F, cluster_cols=F, na_col = "white", filename = "pairwiseWeightedFst_5speciesOneDiagonal_R.pdf")

##### Considering KherSOUTH and KherCENTRAL the same species #####

pacman::p_load(pheatmap, tidyverse, reshape2)
setwd("C:/Users/Waldir Miron/Desktop/IntrogressionAnalysis_tests/FST")
data <- read.table("FstR_4species.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
data

# Adds column names 
colnames(data) <- c("Pop1", "Pop2","Weighted")
# Melts data ~
pops = union(data$Pop1, data$Pop2)
n = length(pops)

# Creates Fst-Sites matrix ~
Fst <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))
for (i in 1:nrow(data)) {
  Fst[data[i, "Pop1"], data[i, "Pop2"]] = data[i, "Weighted"]
  Fst[data[i, "Pop2"], data[i, "Pop1"]] = data[i, "Weighted"]}

Fst

# Writes Fst-Sites matrix ~
write.table(Fst, "Fst-Matrix_4species.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Creates & Saves the heatmap ~
pheatmap(Fst, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", border_color = "black", cellwidth = 40, cellheight = 40,
         treeheight_row = 0, treeheight_col = 0, angle_col = "45", fontsize = 8, filename = "pairwiseWeightedFst_4species_R.pdf")


# Creates & Saves the heatmap with only one diagonal ~
Fst[upper.tri(Fst)] <- NA
# Creates & Saves the heatmap ~
pheatmap(Fst, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", border_color = "black", cellwidth = 40, cellheight = 40,
         treeheight_row = 0, treeheight_col = 0, angle_col = "45", fontsize = 8,cluster_rows=F, cluster_cols=F, na_col = "white", filename = "pairwiseWeightedFst_4speciesOneDiagonal_R.pdf")

