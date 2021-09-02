### The BEGINNING ~~~~~
##
# ~ Estimates Ploidy from GBS Data | Gompert & Mock, 2017 - Molecular Ecology


# Cleans the environment ~
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(gbs2ploidy, devtools, tidyverse)


# Loads SNPs data ~
KryptosIDsInfo <- read.csv("KryptoIDsInfo.csv")


# Loads more data ~
SNPs_data <- as.matrix(read.table("OnlyGBS_SNPs_VcfwithAD.txt", header = FALSE))


# Subsets data ~
a <- seq(1, 54, 2)
b <- seq(2, 54, 2)
cov1 <- SNPs_data[ ,a]
cov2 <- SNPs_data[ ,b]


# Gets individual heterozigosity ~
H <- apply(is.na(cov1) == FALSE, 2, mean)


# Gets coverage depth for heterozygous sites ~
D <- apply(cov1 + cov2, 2, mean, na.rm = TRUE)


# Generates data frame ~
GeneralData <- data.frame(D, H)


# Creates Plot Panel ~
HetvsDepthPlot <-
 ggplot(GeneralData, aes(x = D, y = H, color = KryptosIDsInfo$Species)) +
  geom_smooth(method = "lm", se = FALSE,  color = "darkred", fill = "blue") +
  geom_point(size = 4) +
  scale_color_discrete() +
  labs(x = "Mean Depth", y = ("Heterozygosity")) + 
  theme_bw() +
  theme(panel.grid = element_blank())

## 3. Estimate Allelic Proportions ~

# 3.1. Diploids, Triploids & Tetraploids ~


# Gets ploidy probabilities ~
propOutAllPlody <- estprops(cov1 = cov1, cov2 = cov2, props = c(0.25, 0.33, 0.5, 0.66, 0.75), mcmc.nchain = 3, mcmc.steps = 100, 
                            mcmc.burnin = 10, mcmc.thin = 5)
head(propOutAllPlody)


# Expands propOutAllPlody by adding the Species columns ~
for( i in seq_along(propOutAllPlody)){
  propOutAllPlody[[i]] <- as.data.frame(propOutAllPlody[[i]])
  propOutAllPlody[[i]]$Percentage <- colnames(propOutAllPlody[[i]])
  propOutAllPlody[[i]]$Probs <- rownames(propOutAllPlody[[i]])
  propOutAllPlody[[i]]$SampleID <- rep(KryptosIDsInfo$ID[i], nrow(propOutAllPlody[[i]]))
  propOutAllPlody[[i]]$Species <- rep(KryptosIDsInfo$Species[i], nrow(propOutAllPlody[[i]]))}

propOutAllPlodyUp <- do.call("rbind", propOutAllPlody)
rownames(propOutAllPlodyUp) <- NULL

head(propOutAllPlodyUp)

# Converts DF from wide into long ~
propOutAllPlodyUp <- gather(propOutAllPlodyUp, Percentage, Value, "2.5%", "25%", "50%", "75%", "97.5%")


# Mean Of Each Species
Results = data.frame(lapply(propOutAllPlodyUp, function(x) {aggregate(x, list(propOutAllPlodyUp$Value), FUN = mean)$x}))

# Creates Plots ~

 ggplot(propOutAllPlodyUp, aes(x = Value, fill = Species), colour = "000000") +
  geom_density(alpha = .15, size = .3) +
  scale_x_continuous("Global Depth (X)",
                     #breaks = c(5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 55000, 60000, 65000, 70000),
                     #labels = c("5K", "10K", "15K","20K", "25K", "30K", "35K", "40K", "45K", "50K", "55K", "60K", "65K","70K"),
                     #limits = c(0, 71000),
                     expand = c(0,0)) +
  scale_y_continuous("Density",
                     #breaks = c(0.000025, 0.00005, 0.000075, 0.0001, 0.000125), 
                     #labels = c("2.5e-05", "5e-05", "7.5e-05", "1e-04", "1.25e-04"), 
                     #limits = c(0, 0.000145),
                     expand = c(0,0)) +
  theme(panel.background = element_rect(fill = '#ffffff'),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text = element_text(size = 9, color = "#000000"),
        axis.ticks = element_line(size = .3, color = "#000000"),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title.x = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, face = "bold", color = "#000000", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.position = "none")








###############################################


      ###### 3.2. Diploids and triploids ######
propOutDiploidTriplod<-estprops(cov1=cov1,cov2=cov2,props=c(0.33, 0.5, 0.66),mcmc.nchain=3,mcmc.steps=100,mcmc.burnin=10,mcmc.thin=5)
propOutDiploidTriplod
      ###### 3.3. Diploids and tretraploids ######
propOutDiploidTetraploid<-estprops(cov1=cov1,cov2=cov2,props=c(0.25, 0.5, 0.75),mcmc.nchain=3,mcmc.steps=100,mcmc.burnin=10,mcmc.thin=5)
propOutDiploidTetraploid
      ######  3.4.  Plot allelic proportions and 95% confidence ########
#### Here I only run for all potential allele frequencies. To change for specific ploidy combinations, change variable 'propOut'
## For all individuals
for(i in 1:27){
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$ID[[i]]))
}
dev.off()

## For KsES
par(mfrow=c(2,3))
for(i in 1:7){
  if(i==5||i==6){
    next
  } 
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$ID[[i]]))
}
dev.off()

## For Koce
par(mfrow=c(3,3))
for(i in 8:16){
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$ID[[i]]))
}
dev.off()

## For KherSouth
par(mfrow=c(4,4))

for(i in 5:27){
  if(i>=7 && i<=16){
    next
  } 
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$ID[[i]]))
}
dev.off()


## For KherCENTRAL
par(mfrow=c(2,1))
for(i in 28:29){
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$id[[i]]))
}
dev.off()

## For Kmar
par(mfrow=c(2,2))
for(i in 30:33){
  plot(propOutAllPlody[[i]][,3],ylim=c(0,1),axes=FALSE,xlab="Allelic ratio",ylab="Posterior probability")
  axis(1,at=1:5,c("1:3","1:2","1:1","2:1","3:1"))
  axis(2)
  box()
  segments(1:5,propOutAllPlody[[i]][,1],1:5,propOutAllPlody[[i]][,5])
  title(main=paste(KryptosIDsInfo$id[[i]]))
}

######  3.Estimating ploidy without ploidy information  #### 
      ######  3.1. For two potential cytotypes  ######
TwoCyotypes_noPloidy<-estploidy(alphas=propOutDiploidTriplod,het=H,depth=D,train=FALSE,nclasses=2,ids=KryptosIDsInfo[,1],pcs=1:2)

TwoCyotypes_noPloidy$pcscrs

TwoCyotypes_noPloidy$pp

##PCA
ggbiplot(prcomp(TwoCyotypes_noPloidy$pcscrs), ellipse = TRUE, groups = KryptosIDsInfo$Species)

## Write probabilities file for 2 cytotypes without known ploidy
#write.csv2(TwoCyotypes_noPloidy$pp,"PPs_2cytotypes.txt")

## Write PCA scores for 2 cytotypes without known ploidy
#write.csv2(TwoCyotypes_noPloidy$pcscrs,"PCscores_2cytotypes.txt")

      ######  3.2. For three potential cytotypes  ######
ThreeCyotypes_noPloidy<-estploidy(alphas=propOutAllPlody,het=H,depth=D,train=FALSE,nclasses=3,ids=KryptosIDsInfo[,1],pcs=1:2)

ThreeCyotypes_noPloidy$pcwghts

ThreeCyotypes_noPloidy$pp

##PCA

ggbiplot(prcomp(ThreeCyotypes_noPloidy$pcscrs),  groups = KryptosIDsInfo$Species)+
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    text = element_text(family = 'Times'),
    legend.position = 'right',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


?prcomp()

## Write probabilities file for 3 cytotypes without known ploidy
#write.csv2(ThreeCyotypes_noPloidy$pp,"PPs_3cytotypes.txt")

## Write PCA scores for 3 cytotypes without known ploidy
#write.csv2(ThreeCyotypes_noPloidy$pcscrs,"PCscores_3cytotypes.txt")

      ######  3.3. For four potential cytotypes ######
FourCyotypes_noPloidy<-estploidy(alphas=propOutAllPlody,het=H,depth=D,train=FALSE,nclasses=4,ids=KryptosIDsInfo[,1],pcs=1:2)

FourCyotypes_noPloidy$pcwghts

FourCyotypes_noPloidy$pp

##PCA
ggbiplot(prcomp(FourCyotypes_noPloidy$pcscrs), ellipse = TRUE, groups = KryptosIDsInfo$Species)

## Write probabilities file for 4 cytotypes without known ploidy
#write.csv2(FourCyotypes_noPloidy$pp,"PPs_4cytotypes.txt")

## Write PCA scores for 4 cytotypes without known ploidy
#write.csv2(FourCyotypes_noPloidy$pcscrs,"PCscores_4cytotypes.txt")

######  4.Estimating ploidy with ploidy information ###### 
      ##### 4.1. For two potential cytotypes  ######

truep<-KryptosIDsInfo[[3]]
truep
truep<-KryptosIDsInfo[[3]]
truep
trn<-sort(sample(1:33,16,replace=FALSE))
trn
truep[-trn]<-NA
str(truep)
TwoCyotypes_Ploidy<-estploidy(alphas=propOut,het=H,depth=D,train=TRUE,pl=KryptosIDsInfo$Ploidy, set=trn, nclasses=2,ids=KryptosIDsInfo[,1],pcs=1:2)
TwoCyotypes_Ploidy


## Write probabilities file for 2 cytotypes with known ploidy
write.csv2(TwoCyotypes_withPloidy$pp,"PPs_2cytotypes.txt")

## Write PCA scores for 2 cytotypes with known ploidy
write.csv2(TwoCyotypes_withPloidy$pcscrs,"PCscores_2cytotypes.txt")
