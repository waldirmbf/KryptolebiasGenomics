### The BEGINNING ~~~~~
##
# ~ Plots FPGP--PopGenEstimates | By Marie-Christine RUFENER & George PACHECO


# Cleans the environment ~ 
rm(list=ls())


# Sets working directory ~
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Loads required packages ~
pacman::p_load(scales, extrafont, dplyr, grid, lubridate, cowplot, egg, tidyverse, stringr, reshape)


# Load helper function ~
source("utilities.R")


# Imports extra fonts ~
loadfonts(device = "win", quiet = TRUE)


# Loads datasets ~
PopGen <- read.table("FPG--PopGenEstimates.txt", sep = "\t", header = FALSE); head(PopGen)
Hets <- read.table("FPG--Heterozygosity.txt", sep = "\t", header = FALSE); head(Hets)


# Adds column names to datasets ~
colnames(PopGen) <- c("Population", "NSites", "Nucleotide_Diversity", "Watson_Theta", "Tajima_D"); head(PopGen)
colnames(Hets) <- c("Sample_ID", "Population", "Het", "DataType"); head(Hets)


# Corrects Population names in Hets:
levels(Hets$Population <- sub("FeralUT", "SaltLakeCity", Hets$Population))
levels(Hets$Population <- sub("FeralVA", "Virginia", Hets$Population))


# Expands Hets by adding BioStatus ~
PopGen$BioStatus <- ifelse(PopGen$Population %in% c("Torshavn","Ejde","Sumba","LjosAir","Kunoy","Nolsoy","Crete",
                                                    "Sardinia","Vernelle","WadiHidan","PigeonIsland","Trincomalee"), "Remote_Localities_Within_Natural_Range",
                    ifelse(PopGen$Population %in% c("Guimaraes","Lisbon","Barcelona","Berlin","Cambridge",
                                                    "Colombo","Copenhagen","London","Prague","Jihlava","Abadeh",
                                                    "Isfahan","Lahijan","Nowshahr","Tehran","TelAviv"), "Urban_Localities_Within_Natural_Range",
                    ifelse(PopGen$Population %in% c("SaltLakeCity","Denver","Virginia","TlaxcalaDeXicohtencatl",
                                                    "MexicoCity","Monterrey","SanCristobalDeLasCasas","Santiago",
                                                    "Salvador","Tatui","Johannesburg","Nairobi","Perth"), "Localities_Outside_Natural_Range",
                    ifelse(PopGen$Population %in% c("TelAvivColony","Wattala", "Wellawatte"), "Captive_Populations", "Outgroups")))); head(PopGen)


# Expands Hets by adding BioStatus ~
Hets$BioStatus <- ifelse(Hets$Population %in% c("Torshavn","Ejde","Sumba","LjosAir","Kunoy","Nolsoy","Crete",
                                                "Sardinia","Vernelle","WadiHidan","PigeonIsland","Trincomalee"), "Remote_Localities_Within_Natural_Range",
                  ifelse(Hets$Population %in% c("Guimaraes","Lisbon","Barcelona","Berlin","Cambridge",
                                                "Colombo","Copenhagen","London","Prague","Jihlava","Abadeh",
                                                "Isfahan","Lahijan","Nowshahr","Tehran","TelAviv"), "Urban_Localities_Within_Natural_Range",
                  ifelse(Hets$Population %in% c("SaltLakeCity","Denver","Virginia","TlaxcalaDeXicohtencatl",
                                                "MexicoCity","Monterrey","SanCristobalDeLasCasas","Santiago",
                                                "Salvador","Tatui","Johannesburg","Nairobi","Perth"), "Localities_Outside_Natural_Range",
                  ifelse(Hets$Population %in% c("TelAvivColony","Wattala", "Wellawatte"), "Captive_Populations", "Outgroups")))); head(Hets)


# Tidies up data frames ~
levels(PopGen$Population)
PopGen$Population <- as.factor(gsub(" ", "", PopGen$Population))


# Sets BioStatus from character -> factor (better for data manipulation & necessary for plotting)
PopGen$BioStatus <- as.factor(PopGen$BioStatus)
Hets$BioStatus <- as.factor(Hets$BioStatus)


# Removes unwanted populations ~
UnwantedPops <- c("Virginia", "Ejde", "Sumba", "LjosAir", "Kunoy", "Nolsoy", "Cambridge", "Jihlava",
                  "Wattala", "Wellawatte", "Srisoria", "Cpalumbus", "Crupestris")
PopGen <- filter(PopGen, !Population %in% UnwantedPops)
Hets <- filter(Hets, !Population %in% UnwantedPops)


# Converts DF from wide into long ~
PopGenUp <- gather(PopGen, Estimate, Value, "Nucleotide_Diversity", "Watson_Theta", "Tajima_D")


# Adds data ID column to each DF ~
PopGenUp$ID <- factor(paste("PopGen"))
Hets$ID <- factor(paste("Hets"))


# Binds the 2 DFs based on common columns ~
fulldf <- mybind(PopGenUp, Hets)


# Includes label for empty factor level (related to PHS) ~
idx <- which(fulldf$ID == "Hets")
fulldf[idx,"Estimate"] <- rep("PHS", length(idx))
fulldf$Estimate <- factor(fulldf$Estimate)


# Reorders factor levels ~
fulldf$Estimate <-
  factor(fulldf$Estimate, ordered = T, levels = c("PHS",
                                                  "Nucleotide_Diversity",
                                                  "Watson_Theta",
                                                  "Tajima_D"))

# Corrects the facet labels ~
ylabel <- c("Nucleotide_Diversity" = "Nucelotide Diversity",
            "Tajima_D"= "Tajima's D",
            "Watson_Theta" = "Watson's Theta",
            "PHS"= "Heterozygosity")


# Corrects population names ~
levels(fulldf$Population <- sub("Torshavn", "Tórshavn", fulldf$Population))
levels(fulldf$Population <- sub("WadiHidan", "Wadi Hidan", fulldf$Population))
levels(fulldf$Population <- sub("Tatui", "Tatuí", fulldf$Population))
levels(fulldf$Population <- sub("PigeonIsland", "Pigeon Island", fulldf$Population))
levels(fulldf$Population <- sub("Guimaraes", "Guimarães", fulldf$Population))
levels(fulldf$Population <- sub("SaltLakeCity", "Salt Lake City", fulldf$Population))
levels(fulldf$Population <- sub("MexicoCity", "Mexico City", fulldf$Population))
levels(fulldf$Population <- sub("TlaxcalaDeXicohtencatl", "Tlaxcala de Xicohténcatl", fulldf$Population))
levels(fulldf$Population <- sub("SanCristobalDeLasCasas", "San Cristóbal de las Casas", fulldf$Population))
levels(fulldf$Population <- sub("TelAvivColony", "Tel Aviv Colony", fulldf$Population))
levels(fulldf$Population <- sub("TelAviv", "Tel Aviv", fulldf$Population))


# Reorders populations ~
fulldf$Population <- factor(fulldf$Population, ordered = T,
                            levels = c("Tórshavn", "Crete", "Sardinia", "Vernelle", "Wadi Hidan", "Pigeon Island", "Trincomalee",
                                       "Lisbon", "Guimarães", "Barcelona", "London", "Berlin", "Copenhagen", "Prague", "Tel Aviv", 
                                       "Abadeh", "Tehran", "Lahijan", "Nowshahr", "Isfahan", "Colombo", "Denver", "Salt Lake City", 
                                       "Tlaxcala de Xicohténcatl", "Mexico City", "Monterrey","San Cristóbal de las Casas", "Santiago", 
                                       "Salvador", "Tatuí", "Johannesburg", "Nairobi", "Perth", "Tel Aviv Colony"))

# Reorders BioStatus ~
fulldf$BioStatus <- factor(fulldf$BioStatus, ordered = T,
                            levels = c("Remote_Localities_Within_Natural_Range",
                                       "Urban_Localities_Within_Natural_Range",
                                       "Localities_Outside_Natural_Range",
                                       "Captive_Populations",
                                       "Outgroups"))


# Creates the panel ~
PopGennEstimates <- 
 ggplot() +
  geom_boxplot(data = subset(fulldf, ID == "Hets"),
               aes(x = Population, y = Het, fill = BioStatus), show.legend = FALSE, outlier.colour = "black", outlier.shape = 21, outlier.size = 1.85, width = .3, lwd = .1) +
  geom_point(data = subset(fulldf, ID =="PopGen"),
             aes(x = Population, y = Value, fill = BioStatus), colour = "black", shape = 21, size = 3.5, alpha = .9) +
  facet_grid(Estimate ~. , scales = "free", labeller = labeller(Estimate = ylabel)) +
  scale_fill_manual(values = c("#44AA99", "#F0E442", "#E69F00", "#56B4E9"),
                    labels = gsub("_", " ", levels(fulldf$BioStatus))) +
  scale_colour_manual(values = c("#44AA99", "#F0E442", "#E69F00", "#56B4E9")) +
  theme(panel.background = element_rect(fill = "#ffffff"),
        panel.grid.major.x = element_line(color = "#ededed", linetype = "dashed", size = .00005),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(colour = "#000000", size = .3),
        axis.title = element_blank(),
        axis.text.x = element_text(colour="#000000", size = 16, face = "bold", family = "Helvetica", angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(color="#000000", size = 16, family = "Helvetica"),
        axis.ticks.x = element_line(color="#000000", size = .3),
        axis.ticks.y = element_line(color="#000000", size = .3),
        strip.background.y = element_rect(colour = "#000000", fill = "#d6d6d6", size = 0.3),
        strip.text = element_text(colour = "#000000", size = 12, face = "bold", family = "Georgia"),
        legend.position = "top",
        legend.margin = margin(t = 0, b = 0, r = 0, l = 0),
        legend.box.margin = margin(t = 10, b = 20, r = 0, l = 0),
        legend.key = element_rect(fill = NA),
        legend.background =element_blank()) +
  guides(fill = guide_legend(title = "Biological Status", title.theme = element_text(size = 16, face = "bold", family = "Helvetica"),
                             label.theme = element_text(size = 14, family = "Helvetica"),
                             override.aes = list(size = 5, alpha = .9)), colour = "none")


# Saves the panel ~
ggsave(PopGennEstimates, file = "FPG--PopGenEstimates.pdf", device = cairo_pdf, width = 12, height = 8, scale = 1.5, dpi = 600)


#
##
### The END ~~~~~