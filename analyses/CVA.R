# Sexual dimorphism in bat wing morphology – variation among foraging styles

# Authors: Dominique G. Maucieri[1], Austin J. Ashbaugh[2] and Jessica M. Theodor[2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Department of Biology, University of Calgary, Calgary, Alberta, T2N 1N4, Canada 
#
# Corresponding Author: Dominique G. Maucieri, dominiquemaucieri@uvic.ca

# Script to conduct CVA analyses

##############################

## load packages 
library(rgl)
library(geomorph)
library(ggplot2)
library(dplyr)
library(Morpho)
library(dplyr)

## Set your working directory
# Make sure that this contains the "Maucieri_etal_2021_CJZ" folder
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

load("data/PCA.Rdata")
load("data/classifier.Rdata")
load("data/filelist.Rdata")
load("data/PCSCORES.Rdata")

theme_DGM <- function () { 
  theme_classic(base_size=12, base_family="Times New Roman") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=18, face="bold"), 
                                                                     axis.text = element_text(size=16, face="plain"), legend.text=element_text(size=16, face="plain"), legend.title = element_text(size=18, face="bold"))}

##############################

#making vectors for different classifiers
GUILDS <- 1:length(classifier$Number)  
SEX <- 1:length(classifier$Number)  
DIET <- 1:length(classifier$Number)  

for (i in 1:length(classifier$Number)) {
  
  xx <- filelist[i]
  xx <-gsub(".txt", "", paste(xx))
  
  yy <- subset(classifier, Number == xx)
  
  GUILDS[i] <- yy$SexandDiet
  
  SEX[i] <- yy$Sex
  
  DIET[i] <- yy$Diet
  
}

#Just Diet CVA
DIETCVA <- CVA(wingpca$x[,c(1:13)], DIET, plot = TRUE)
DIETCVA

DIET_CVAscores <- as.data.frame(DIETCVA$CVscores)

CVAscores_DIET <- as.data.frame(matrix(data = NA, ncol = 4, nrow = length(GUILDS)))
colnames(CVAscores_DIET) <- c("ID", "CV1", "CV2", "Diet")
CVAscores_DIET$ID <- PCSCORES$ID
CVAscores_DIET$CV1 <- DIET_CVAscores$`CV 1`
CVAscores_DIET$CV2 <- DIET_CVAscores$`CV 2`
CVAscores_DIET$Diet <- DIET

CVA_hull <- CVAscores_DIET %>% group_by(Diet) %>%  slice(chull(CV1, CV2))
CVA_hull$Diet <- factor(CVA_hull$Diet, levels = c("g", "a", "f"))
CVAscores_DIET$Diet <- factor(CVAscores_DIET$Diet, levels = c("g", "a", "f"))

diet_colors_CVA <- c("#440154FF", "#21908CFF", "#FDE725FF")

Figure_S4 <- ggplot(CVAscores_DIET, aes(x = CV1, y = CV2)) + labs(x="CV1", y="CV2", title="", colour="Diet", size=2)+
  geom_point(aes(shape=Diet,colour=Diet),size=3)+scale_shape_manual(values=c(2,0,1))+ scale_color_manual(values = diet_colors_CVA,
  breaks = c("g", "a", "f"), labels = c("Gleaner", "Aerial", "Frugivore")) +
  aes(fill = factor(Diet)) + geom_polygon(data = CVA_hull, alpha = 0.2) +scale_fill_manual(values = diet_colors_CVA) +
  guides(shape = FALSE, fill = FALSE, colour = guide_legend(override.aes = list(shape = c(17,15,19), fill = "white"))) + theme_DGM()+ 
  labs(color="Foraging Guild")
Figure_S4

# tiff("analyses/figures/Figure_S4.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_S4
# dev.off()

diet_colors_gsc <- c("black","black","black")

Figure_2 <- ggplot(CVAscores_DIET, aes(x = CV1, y = CV2)) + labs(x="CV1", y="CV2", title="", size=2)+
  geom_point(aes(shape=Diet,colour=Diet),size=3)+scale_shape_manual(values=c(2,0,1))+ scale_color_manual(values = diet_colors_gsc,
  breaks = c("g", "a", "f"), labels = c("Gleaner", "Aerial", "Frugivore")) +
  aes(fill = factor(Diet)) + geom_polygon(data = CVA_hull, alpha = 0.2, color = "black") + scale_fill_grey(start = 0, end = .9) +
  guides(fill = FALSE, shape = FALSE, color = guide_legend(override.aes = list(shape = c(2,0,1), fill = "white"))) + theme_DGM()+ 
  labs(color="Foraging Guild")
Figure_2

# tiff("analyses/figures/Figure_2.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_2
# dev.off()


# Figure 2 and Figure S4 sample sizes

CVA_Diet_samplesizes <- CVAscores_DIET %>% group_by(Diet) %>% summarise(sample_size = length(Diet))



