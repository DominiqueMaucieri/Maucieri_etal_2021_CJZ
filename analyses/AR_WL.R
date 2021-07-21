# Sexual dimorphism in bat wing morphology – variation among foraging styles

# Authors: Dominique G. Maucieri[1], Austin J. Ashbaugh[2] and Jessica M. Theodor[2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Department of Biology, University of Calgary, Calgary, Alberta, T2N 1N4, Canada
#
# Corresponding Author: Dominique G. Maucieri, dominiquemaucieri@uvic.ca

# Script to conduct aspect ratio and wing loading analyses

##############################

## load packages
library(ggplot2)
library(dplyr)
library(geomorph)
library(dplyr)

## Set your working directory
# Make sure that this contains the "Maucieri_etal_2021_CJZ" folder
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

load("data/classifier.Rdata")
load("data/landmarks.Rdata")

theme_DGM <- function () { 
  theme_classic(base_size=12, base_family="Times New Roman") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                                                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=18, face="bold"), 
                                                                     axis.text = element_text(size=16, face="plain"), legend.text=element_text(size=16, face="plain"), legend.title = element_text(size=18, face="bold"))}

##############################

indvARWL <- classifier

for (i in 1:length(indvARWL$Species)) {
  if(indvARWL$Diet[i] == "f"){indvARWL$Diet[i] <- "Frugivore"} else {
    if(indvARWL$Diet[i] == "a"){indvARWL$Diet[i] <- "Aerial"} else {
      if(indvARWL$Diet[i] == "g"){indvARWL$Diet[i] <- "Gleaner"} else {
          indvARWL$Diet[i] <- "BLANK"}}}}

indvARWL$Wing_Loading <- ((indvARWL$Bat_Mass / 1000) * 9.81) / (((indvARWL$Wing_Area)*2)/10000)

indvARWL$Aspect_Ratio <- (((indvARWL$Wing_Length)*2)^2)/((indvARWL$Wing_Area)*2)

indvARWL$ID_Number <- paste(indvARWL$Species, indvARWL$Sex, sep = " ")

indvARWLnopreg <- subset(indvARWL, Reproductive_Status != "P")

MaleBats <- subset(indvARWL, Sex=="M")
FemaleBats <- subset(indvARWL, Sex=="F")
Maleaerial <- subset(indvARWL, Sex=="M" & Diet=="a")
Femaleaerial <- subset(indvARWL, Sex=="F" & Diet=="a")
Malegleaner <- subset(indvARWL, Sex=="M" & Diet=="g")
Femalegleaner <- subset(indvARWL, Sex=="F" & Diet=="g")
Malefrug <- subset(indvARWL, Sex=="M" & Diet=="f")
Femalefrug <- subset(indvARWL, Sex=="F" & Diet=="f")

MaleBats_nopreg <- subset(indvARWLnopreg, Sex=="M")
FemaleBats_nopreg <- subset(indvARWLnopreg, Sex=="F")
Maleaerial_nopreg <- subset(indvARWLnopreg, Sex=="M" & Diet=="a")
Femaleaerial_nopreg <- subset(indvARWLnopreg, Sex=="F" & Diet=="a")
Malegleaner_nopreg <- subset(indvARWLnopreg, Sex=="M" & Diet=="g")
Femalegleaner_nopreg <- subset(indvARWLnopreg, Sex=="F" & Diet=="g")
Malefrug_nopreg <- subset(indvARWLnopreg, Sex=="M" & Diet=="f")
Femalefrug_nopreg <- subset(indvARWLnopreg, Sex=="F" & Diet=="f")


means_wpreg <- indvARWL %>% group_by(Sex, Diet) %>% dplyr::summarise(Mean_WL= mean(Wing_Loading), SE_WL= (sd(Wing_Loading))/(sqrt(length(Wing_Loading))), 
               Mean_AR= mean(Aspect_Ratio), SE_AR= (sd(Aspect_Ratio))/(sqrt(length(Aspect_Ratio))), sample_size = length(Aspect_Ratio))

means_nopreg <- indvARWLnopreg %>% group_by(Sex, Diet) %>% dplyr::summarise(Mean_WL= mean(Wing_Loading), SE_WL= (sd(Wing_Loading))/(sqrt(length(Wing_Loading))), 
                Mean_AR= mean(Aspect_Ratio), SE_AR= (sd(Aspect_Ratio))/(sqrt(length(Aspect_Ratio))), sample_size = length(Wing_Loading))


indvARWL$Status <- "Blank"
indvARWL$Centroid_Size <- "Blank"

Procrustes <- gpagen(coords, Proj = TRUE)
Centroids <- Procrustes$Csize

Centroids_df <- as.data.frame(Centroids)
Centroids_df$Number <- rownames(Centroids_df)

for(i in 1:length(indvARWL$Species)){
  if(indvARWL$Reproductive_Status[i] == "P"){indvARWL$Status[i] <- "Preg"}else{
    indvARWL$Status[i] <- "Not_Preg"}
  
  zz <- subset(Centroids_df, Number == indvARWL$Number[i])
  indvARWL$Centroid_Size[i] <- zz$Centroids
}

WingLoad <- indvARWL

WingLoadFemale <- subset(WingLoad, Sex == "F")

lm_Prego <- t.test(log(WingLoadFemale$Wing_Loading) ~ WingLoadFemale$Status)
lm_Prego
hist(log(WingLoadFemale$Wing_Loading))

WingLoadFemale %>% group_by(Status) %>% summarise(mean = mean(Wing_Loading), standard_error = ((sd(Wing_Loading))/(sqrt(length(Wing_Loading)))))

#Aspect Ratio

hist(indvARWL$Aspect_Ratio)

AR_aov <- aov(Aspect_Ratio ~ Sex * Diet + (Species %in% Diet), data = indvARWL)
summary(AR_aov)
# par(mfrow=c(2,2))
# plot(AR_aov)

WL_aov <- aov(log(Wing_Loading) ~ Sex * Diet + (Species %in% Diet), data = indvARWLnopreg)
summary(WL_aov)
# par(mfrow=c(2,2))
# plot(WL_aov)


indvARWL$Diet <- factor(indvARWL$Diet, levels = c("Gleaner", "Aerial", "Frugivore"))

indvARWLnopreg$Diet <- factor(indvARWLnopreg$Diet, levels = c("Gleaner", "Aerial", "Frugivore"))


Figure_3 <- ggplot(indvARWL, aes(x = Diet, y = Aspect_Ratio, fill = Sex)) +labs(x="Foraging guild", y="Aspect ratio", title="", fill="Sex", size=3) +
  geom_boxplot()+ theme_DGM()+scale_fill_manual(values = c("#FFFFFF", "#A9A9A9"))
Figure_3

# tiff("analyses/figures/Figure_3.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_3
# dev.off()

Figure_4 <- ggplot(indvARWLnopreg, aes(x = Diet, y = Wing_Loading, fill = Sex)) +labs(x="Foraging guild", y=expression(bold(paste("Wing loading (N/", "m"^2, ")"))), title="", fill="Sex", size=3) + 
  geom_boxplot()+ theme_DGM()+scale_fill_manual(values = c("#FFFFFF", "#A9A9A9"))
Figure_4

# tiff("analyses/figures/Figure_4.tiff", width=9, height = 6, units = "in", pointsize = "14", compression = "none", type = "cairo", res = 300)
# Figure_4
# dev.off()






