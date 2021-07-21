# Sexual dimorphism in bat wing morphology – variation among foraging styles

# Authors: Dominique G. Maucieri[1], Austin J. Ashbaugh[2] and Jessica M. Theodor[2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Department of Biology, University of Calgary, Calgary, Alberta, T2N 1N4, Canada
#
# Corresponding Author: Dominique G. Maucieri, dominiquemaucieri@uvic.ca

# Script to conduct multivariate analyses

##############################

## load packages
library(rgl)
library(geomorph)
library(ggplot2)
library(ggfortify)
library(dispRity)
library(dplyr)

## Set your working directory
# Make sure that this contains the "Maucieri_etal_2021_CJZ" folder
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

load("data/PCA.Rdata")
load("data/classifier.Rdata")
load("data/filelist.Rdata")
load("data/landmarks.Rdata")

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

Procrustes <- gpagen(coords, Proj = TRUE)
Centroids <- Procrustes$Csize
Procrustes_Coords <- Procrustes$coords
gdf_allometry<-geomorph.data.frame(coords=Procrustes_Coords,Csize=Centroids)
reg.res <- procD.lm(coords~Csize, data = gdf_allometry,RRPP=TRUE,logsz=TRUE,iter=999)
corrected_coords_wings <- reg.res$GM$residuals + replicate(length(Centroids), mshape(Procrustes_Coords))
coords_wings<-corrected_coords_wings

Raw_PCSCORES <- as.data.frame(wingpca$x)
PCSCORES <- as.data.frame(matrix(data = NA, ncol = 6, nrow = length(Raw_PCSCORES$PC1)))
colnames(PCSCORES) <- c("ID", "PC1", "PC2", "Guild", "Sex", "Diet")
PCSCORES$ID <- rownames(Raw_PCSCORES)
PCSCORES$PC1 <- Raw_PCSCORES$PC1
PCSCORES$PC2 <- Raw_PCSCORES$PC2
PCSCORES$Guild <- GUILDS
PCSCORES$Sex <- SEX
PCSCORES$Diet <- DIET

#Save PC scores to .RDATA file
# save(PCSCORES, file = "data/PCSCORES.RData")

#PC 1 --> 14 explains 95% of the morphological variation
Multivariate_Analysis <- geomorph.data.frame(fine = wingpca$x[,c(1:13)], sex = PCSCORES$Sex, diet = PCSCORES$Diet)

Procrustes_ANOVA <- procD.lm ((fine ~ sex*diet), RRPP=TRUE,iter=999, data = Multivariate_Analysis)
summary(Procrustes_ANOVA)

disparity_data<-cbind(wingpca$x[,c(1:13)],diet = PCSCORES$Diet)
disparity_data_df <- as.data.frame(disparity_data)
cols.num <- c(1:13)
disparity_data_df[cols.num] <- sapply(disparity_data_df[cols.num],as.numeric)

Diet_list<-list(
  a=which(disparity_data_df$diet == "a"),
  f=which(disparity_data_df$diet == "f"),
  g=which(disparity_data_df$diet == "g"))

#Run  disparity analyses
#NOTE - When running the disparity analysis you will get differing values as the dispRity.per.group() function uses bootstrapping
# which is a randomized function, so we have set our seed so you can get the same results we did in our study. 
set.seed(794)
disparity_bats <- dispRity.per.group(disparity_data_df[,c(1:13)],Diet_list,metric = pairwise.dist)
summary(disparity_bats)
# plot(disparity_bats)

#Run disparity t-tests w bonferroni correction
bats_post_HOC<-test.dispRity(disparity_bats, test = t.test, comparisons = "pairwise",
                             concatenate = FALSE, correction = "bonferroni")
bats_post_HOC[[3]]

#Plot for Diet
twoD_coords_wings<- two.d.array(coords_wings)
twoD_coords_wings_wDiet <- cbind(twoD_coords_wings, DIET)
colnames(twoD_coords_wings_wDiet) <- c(rep(c("X", "Y"), times = 19),"Diet")
twoD_coords_wings_wDiet <- as.data.frame(twoD_coords_wings_wDiet)
twoD_coords_wings_wDiet$Diet <- factor(twoD_coords_wings_wDiet$Diet, levels = c("g", "a", "f"))

diet_colors <- c("#440154FF", "#21908CFF", "#FDE725FF")

Figure_S3 <- autoplot(wingpca,data=twoD_coords_wings_wDiet,colour="Diet",shape="Diet",frame = TRUE)+geom_point(aes(shape=DIET,colour=DIET),size=3)+
  scale_shape_manual(values=c(2,0,1))+ scale_fill_manual(values = diet_colors) + scale_color_manual(values=diet_colors, breaks = c("g", "a", "f"), labels = c("Gleaner", "Aerial", "Frugivore")) +
  theme_DGM()+ guides(shape = FALSE, fill = FALSE, colour = guide_legend(override.aes = list(shape = c(2,0,1), fill = "white"))) + 
  labs(color="Foraging Guild")
Figure_S3

# tiff("analyses/figures/Figure_S3.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_S3
# dev.off()


diet_colors_gsc <- c("black","black","black")

Figure_1 <- autoplot(wingpca,data=twoD_coords_wings_wDiet,colour="Diet",shape="Diet",frame = TRUE)+geom_point(aes(shape=DIET,colour=DIET),size=3)+
  scale_shape_manual(values=c(2,0,1))+scale_color_manual(values=diet_colors_gsc, breaks = c("g", "a", "f"), labels = c("Gleaner", "Aerial", "Frugivore")) + 
  scale_fill_grey(start = 0.5, end = .54) + theme_DGM()+ guides(fill = FALSE, shape = FALSE, color = guide_legend(override.aes = list(shape = c(2,0,1), fill = "white")))+ 
  labs(color="Foraging Guild")
Figure_1

# tiff("analyses/figures/Figure_1.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_1
# dev.off()

#Figure 1 and Figure S3 sample sizes
Diet.df <- as.data.frame(twoD_coords_wings_wDiet)
Diet_vector <- Diet.df[,37:39]
Diet_samplesizes <- Diet_vector %>% group_by(Diet) %>% summarise(sample_size = length(Diet))



