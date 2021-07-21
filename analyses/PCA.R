# Sexual dimorphism in bat wing morphology – variation among foraging styles

# Authors: Dominique G. Maucieri[1], Austin J. Ashbaugh[2] and Jessica M. Theodor[2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Department of Biology, University of Calgary, Calgary, Alberta, T2N 1N4, Canada 
#
# Corresponding Author: Dominique G. Maucieri, dominiquemaucieri@uvic.ca

# Script to conduct PCA analyses

##############################

## load packages 
library(rgl)
library(geomorph)
library(ggplot2)
library(ggfortify)
library(viridis)
library(dplyr)

## Set your working directory
# Make sure that this contains the "Maucieri_etal_2021_CJZ" folder
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

load("data/landmarks.Rdata")
load("data/classifier.Rdata")
load("data/filelist.Rdata")

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


#Crusties
Procrustes <- gpagen(coords, Proj = TRUE)
plot(Procrustes)

Centroids <- Procrustes$Csize

Procrustes_Coords <- Procrustes$coords
x <- two.d.array(Procrustes_Coords)

#Examining effect of allometry 

#Make a gdf for use in the static allometry model
gdf_allometry<-geomorph.data.frame(coords=Procrustes_Coords,Csize=Centroids)

#Allometry
reg.res <- procD.lm(coords~Csize, data = gdf_allometry,RRPP=TRUE,logsz=TRUE,iter=999)
plot(reg.res,type="regression",predictor=log(Centroids))

#Check results of static allometric regression
summary(reg.res)
par(mfrow=c(2,2))
plot(reg.res)
par(mfrow=c(1,1))

#Extract residual shape coordinates bc significant result
corrected_coords_wings <- reg.res$GM$residuals + replicate(length(Centroids), mshape(Procrustes_Coords))

#PLot and rename coords
coords_wings<-corrected_coords_wings
plotAllSpecimens(corrected_coords_wings)

#to save procrustes figure
# tiff("analyses/figures/Procrustes.tiff", width=12, height = 10, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# plotAllSpecimens(corrected_coords_wings)
# dev.off()

#Check specimens
mshape_wings<-mshape(corrected_coords_wings)

#MPlot all specimens
par(mfrow=c(3,3))
for (i in 1:205){
  plotRefToTarget(mshape_wings,corrected_coords_wings[,,i],mag=1.5,method="TPS")
  mtext(as.character(i),side=3)
}
dev.off()

#Reformat two-s-array into 3D array like uncorrected coordinates
twoD_coords_wings<- two.d.array(coords_wings)
wingpca<-prcomp(twoD_coords_wings)
summary(wingpca)

#scree plot of varation
barplot(wingpca$sdev)

#Corrected
coords_wings_wGuides <- cbind(twoD_coords_wings, GUILDS)
colnames(coords_wings_wGuides) <- c(rep(c("X", "Y"), times = 19),"Guild")
Guild.df <- as.data.frame(coords_wings_wGuides)
Guild.df$Guild <- factor(Guild.df$Guild, levels = c("Female Gleaner", "Male Gleaner", "Female Insectivore", "Male Insectivore", "Female Frugivore", "Male Frugivore"))



#PLot for removed centroid correction
Figure_S2 <- autoplot(wingpca,data=Guild.df,colour="Guild",shape="Guild",frame = TRUE)+geom_point(aes(shape=GUILDS,colour=GUILDS),size=3)+
  scale_shape_manual(values=c(0,1,2,3,4,5), breaks = c("Female Gleaner", "Male Gleaner", "Female Insectivore", "Male Insectivore", "Female Frugivore", "Male Frugivore"), labels = c("Female Gleaner", "Male Gleaner", "Female Aerial", "Male Aerial", "Female Frugivore", "Male Frugivore")) +
  theme_DGM() + scale_fill_viridis(discrete=TRUE, breaks = c("Female Gleaner", "Male Gleaner", "Female Insectivore", "Male Insectivore", "Female Frugivore", "Male Frugivore"), labels = c("Female Gleaner", "Male Gleaner", "Female Aerial", "Male Aerial", "Female Frugivore", "Male Frugivore"))+ 
  scale_color_viridis(discrete=TRUE, breaks = c("Female Gleaner", "Male Gleaner", "Female Insectivore", "Male Insectivore", "Female Frugivore", "Male Frugivore"), labels = c("Female Gleaner", "Male Gleaner", "Female Aerial", "Male Aerial", "Female Frugivore", "Male Frugivore"))+ 
  labs(color="Sex and Foraging Guild", shape="Sex and Foraging Guild", fill="Sex and Foraging Guild") + theme(legend.title = element_text(size=16, face="bold"))
Figure_S2

#Save figure which was Figure S2 in the manuscript
# tiff("analyses/figures/Figure_S2.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_S2
# dev.off()

#Sample Sizes for Figure S2:

Guild_vector <- Guild.df[,37:39]
guild_samplesizes <- Guild_vector %>% group_by(Guild) %>% summarise(sample_size = length(Guild))


#Warp grid generation
gm_pca_wings<-gm.prcomp(coords_wings)
par(mfrow = c(2,2))
plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp1$min, mag=-0.5,method="TPS")
plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp1$max, mag=-0.5,method="TPS")
plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp2$min, mag=-0.5,method="TPS")
plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp2$max, mag=-0.5,method="TPS")
dev.off()

#Save the warp grids which were added to Figure 1, S2, and S3

# tiff("analyses/figures/PC1_min.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp1$min, mag=-0.5,method="TPS")
# dev.off()
# 
# tiff("analyses/figures/PC1_max.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp1$max, mag=-0.5,method="TPS")
# dev.off()
# 
# tiff("analyses/figures/PC2_min.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp2$min, mag=-0.5,method="TPS")
# dev.off()
# 
# tiff("analyses/figures/PC2_max.tiff", width=6, height = 4, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# plotRefToTarget(mshape_wings, gm_pca_wings$shapes$shapes.comp2$max, mag=-0.5,method="TPS")
# dev.off()

#Save PCA results to .RDATA file
# save(wingpca, file = "data/PCA.RData")


