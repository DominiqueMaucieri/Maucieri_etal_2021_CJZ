# Sexual dimorphism in bat wing morphology – variation among foraging styles

# Authors: Dominique G. Maucieri[1], Austin J. Ashbaugh[2] and Jessica M. Theodor[2]
# Institution: [1] Department of Biology, University of Victoria, Victoria, British Columbia, V8P 5C2, Canada
# [2] Department of Biology, University of Calgary, Calgary, Alberta, T2N 1N4, Canada
#
# Corresponding Author: Dominique G. Maucieri, dominiquemaucieri@uvic.ca

# Script to conduct permutation tests

##############################

## load packages
library(dplyr)

## Set your working directory
# Make sure that this contains the "Maucieri_etal_2021_CJZ" folder
setwd("C:/Users/...") # If on a PC
setwd("/Users/...") # If on a Mac

load("data/classifier.Rdata")

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


perm_matrix <- as.data.frame(matrix(data = NA, nrow = 1000, ncol = 12))
colnames(perm_matrix) <- c("teststat1", "teststat2", "teststat3", "teststat4", "teststat5", "teststat6", "teststat7", "teststat8", "teststat9", "teststat10", "teststat11", "teststat12")

permARWL <- indvARWL
permARWLnopreg <- indvARWLnopreg

for(i in 1:1000){
  
  set.seed(57 + i)
  
  permARWL$Aspect_Ratio_Trial <- NA
  
  permARWL$Aspect_Ratio_Trial <- sample(indvARWL$Aspect_Ratio, replace = FALSE)
  
  perm_AR_sexanddiet <- permARWL %>% group_by(SexandDiet) %>% summarise(meanAR = mean(Aspect_Ratio_Trial))
  
  F.Frug.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Female Frugivore")$meanAR
  M.Frug.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Male Frugivore")$meanAR
  F.Glean.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Female Gleaner")$meanAR
  M.Glean.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Male Gleaner")$meanAR
  F.Insect.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Female Insectivore")$meanAR
  M.Insect.perm.ar.sexanddiet <- subset(perm_AR_sexanddiet, SexandDiet == "Male Insectivore")$meanAR
  
  perm_matrix$teststat1[i] <- (((F.Glean.perm.ar.sexanddiet/M.Glean.perm.ar.sexanddiet)+(F.Insect.perm.ar.sexanddiet/M.Insect.perm.ar.sexanddiet))/2) - 
    (F.Frug.perm.ar.sexanddiet/M.Frug.perm.ar.sexanddiet)
  perm_matrix$teststat2[i] <- ((F.Glean.perm.ar.sexanddiet/M.Glean.perm.ar.sexanddiet)-(F.Insect.perm.ar.sexanddiet/M.Insect.perm.ar.sexanddiet))
  perm_matrix$teststat3[i] <- (F.Glean.perm.ar.sexanddiet - M.Glean.perm.ar.sexanddiet)
  perm_matrix$teststat4[i] <- (F.Insect.perm.ar.sexanddiet - M.Insect.perm.ar.sexanddiet)
  
  
  perm_AR_diet <- permARWL %>% group_by(Diet) %>% summarise(meanAR = mean(Aspect_Ratio_Trial))
  
  Frug.perm.ar.diet <- subset(perm_AR_diet, Diet == "Frugivore")$meanAR
  Glean.perm.ar.diet <- subset(perm_AR_diet, Diet == "Gleaner")$meanAR
  Insect.perm.ar.diet <- subset(perm_AR_diet, Diet == "Aerial")$meanAR
  
  perm_matrix$teststat5[i] <- (((Insect.perm.ar.diet + Glean.perm.ar.diet)/2)-(Frug.perm.ar.diet))
  perm_matrix$teststat6[i] <- (Glean.perm.ar.diet) - (Insect.perm.ar.diet)
  
  
  permARWLnopreg$Wing_Loading_Trial <- NA
  
  permARWLnopreg$Wing_Loading_Trial <- sample(indvARWLnopreg$Wing_Loading, replace = FALSE)
  
  perm_WL_sexanddiet <- permARWLnopreg %>% group_by(SexandDiet) %>% summarise(meanWL = mean(Wing_Loading_Trial))
  
  F.Frug.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Female Frugivore")$meanWL
  M.Frug.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Male Frugivore")$meanWL
  F.Glean.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Female Gleaner")$meanWL
  M.Glean.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Male Gleaner")$meanWL
  F.Insect.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Female Insectivore")$meanWL
  M.Insect.perm.wl.sexanddiet <- subset(perm_WL_sexanddiet, SexandDiet == "Male Insectivore")$meanWL
  
  perm_matrix$teststat7[i] <- (((F.Glean.perm.wl.sexanddiet/M.Glean.perm.wl.sexanddiet)+(F.Insect.perm.wl.sexanddiet/M.Insect.perm.wl.sexanddiet))/2) - 
    (F.Frug.perm.wl.sexanddiet/M.Frug.perm.wl.sexanddiet)
  perm_matrix$teststat8[i] <- ((F.Glean.perm.wl.sexanddiet/M.Glean.perm.wl.sexanddiet)-(F.Insect.perm.wl.sexanddiet/M.Insect.perm.wl.sexanddiet))
  perm_matrix$teststat9[i] <- (F.Glean.perm.wl.sexanddiet - M.Glean.perm.wl.sexanddiet)
  perm_matrix$teststat10[i] <- (F.Insect.perm.wl.sexanddiet - M.Insect.perm.wl.sexanddiet)
  
  
  perm_WL_diet <- permARWLnopreg %>% group_by(Diet) %>% summarise(meanWL = mean(Wing_Loading_Trial))
  
  Frug.perm.wl.diet <- subset(perm_WL_diet, Diet == "Frugivore")$meanWL
  Glean.perm.wl.diet <- subset(perm_WL_diet, Diet == "Gleaner")$meanWL
  Insect.perm.wl.diet <- subset(perm_WL_diet, Diet == "Aerial")$meanWL
  
  perm_matrix$teststat11[i] <- (((Insect.perm.wl.diet + Glean.perm.wl.diet)/2)-(Frug.perm.wl.diet))
  perm_matrix$teststat12[i] <- (Glean.perm.wl.diet) - (Insect.perm.wl.diet)
  
}


##AR
means.real.ar.sexanddiet <- indvARWL %>% group_by(SexandDiet) %>% summarise(meanAR = mean(Aspect_Ratio))

F.Frug.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Female Frugivore")$meanAR
M.Frug.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Male Frugivore")$meanAR
F.Glean.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Female Gleaner")$meanAR
M.Glean.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Male Gleaner")$meanAR
F.Insect.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Female Insectivore")$meanAR
M.Insect.real.ar.sexanddiet <- subset(means.real.ar.sexanddiet, SexandDiet == "Male Insectivore")$meanAR

teststat1.real.ar.sexanddiet <- (((F.Glean.real.ar.sexanddiet/M.Glean.real.ar.sexanddiet)+(F.Insect.real.ar.sexanddiet/M.Insect.real.ar.sexanddiet))/2) - (F.Frug.real.ar.sexanddiet/M.Frug.real.ar.sexanddiet)
teststat1.real.ar.sexanddiet
(length(perm_matrix$teststat1[perm_matrix$teststat1 < teststat1.real.ar.sexanddiet])*2)/1000
hist(perm_matrix$teststat1)
abline(v= teststat1.real.ar.sexanddiet, col="red")


teststat2.real.ar.sexanddiet <- ((F.Glean.real.ar.sexanddiet/M.Glean.real.ar.sexanddiet)-(F.Insect.real.ar.sexanddiet/M.Insect.real.ar.sexanddiet))
teststat2.real.ar.sexanddiet
(length(perm_matrix$teststat2[perm_matrix$teststat2 > teststat2.real.ar.sexanddiet])*2)/1000
hist(perm_matrix$teststat2)
abline(v= teststat2.real.ar.sexanddiet, col="red")


teststat3.real.ar.sexanddiet <- (F.Glean.real.ar.sexanddiet - M.Glean.real.ar.sexanddiet)
teststat3.real.ar.sexanddiet
(length(perm_matrix$teststat3[perm_matrix$teststat3 > teststat3.real.ar.sexanddiet])*2)/1000
hist(perm_matrix$teststat3)
abline(v= teststat3.real.ar.sexanddiet, col="red")


teststat4.real.ar.sexanddiet <- (F.Insect.real.ar.sexanddiet - M.Insect.real.ar.sexanddiet)
teststat4.real.ar.sexanddiet
(length(perm_matrix$teststat4[perm_matrix$teststat4 > teststat4.real.ar.sexanddiet])*2)/1000
hist(perm_matrix$teststat4)
abline(v= teststat4.real.ar.sexanddiet, col="red")


means.real.ar.diet <- indvARWL %>% group_by(Diet) %>% summarise(meanAR = mean(Aspect_Ratio))

Frug.real.ar.diet <- subset(means.real.ar.diet, Diet == "Frugivore")$meanAR
Glean.real.ar.diet <- subset(means.real.ar.diet, Diet == "Gleaner")$meanAR
Insect.real.ar.diet <- subset(means.real.ar.diet, Diet == "Aerial")$meanAR

teststat5.real <- (((Insect.real.ar.diet + Glean.real.ar.diet)/2)-(Frug.real.ar.diet))
teststat5.real
(length(perm_matrix$teststat5[perm_matrix$teststat5 < teststat5.real])*2)/1000
hist(perm_matrix$teststat5)
abline(v= teststat5.real, col="red")

teststat6.real <- (Insect.real.ar.diet - Glean.real.ar.diet)
teststat6.real
(length(perm_matrix$teststat6[perm_matrix$teststat6 > teststat6.real])*2)/1000
hist(perm_matrix$teststat6)
abline(v= teststat6.real, col="red")


#WL
means.real.wl.sexanddiet <- indvARWLnopreg %>% group_by(SexandDiet) %>% summarise(meanWL = mean(Wing_Loading))

F.Frug.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Female Frugivore")$meanWL
M.Frug.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Male Frugivore")$meanWL
F.Glean.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Female Gleaner")$meanWL
M.Glean.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Male Gleaner")$meanWL
F.Insect.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Female Insectivore")$meanWL
M.Insect.real.wl.sexanddiet <- subset(means.real.wl.sexanddiet, SexandDiet == "Male Insectivore")$meanWL


teststat7.real.wl.sexanddiet <- (((F.Glean.real.wl.sexanddiet/M.Glean.real.wl.sexanddiet)+(F.Insect.real.wl.sexanddiet/M.Insect.real.wl.sexanddiet))/2) - (F.Frug.real.wl.sexanddiet/M.Frug.real.wl.sexanddiet)
teststat7.real.wl.sexanddiet
(length(perm_matrix$teststat7[perm_matrix$teststat7 < teststat7.real.wl.sexanddiet])*2)/1000
hist(perm_matrix$teststat7)
abline(v= teststat7.real.wl.sexanddiet, col="red")


teststat8.real.wl.sexanddiet <- ((F.Glean.real.wl.sexanddiet/M.Glean.real.wl.sexanddiet)-(F.Insect.real.wl.sexanddiet/M.Insect.real.wl.sexanddiet))
teststat8.real.wl.sexanddiet
(length(perm_matrix$teststat8[perm_matrix$teststat8 > teststat8.real.wl.sexanddiet])*2)/1000
hist(perm_matrix$teststat8)
abline(v= teststat8.real.wl.sexanddiet, col="red")


teststat9.real.wl.sexanddiet <- (F.Glean.real.wl.sexanddiet - M.Glean.real.wl.sexanddiet)
teststat9.real.wl.sexanddiet
(length(perm_matrix$teststat9[perm_matrix$teststat9 > teststat9.real.wl.sexanddiet])*2)/1000
hist(perm_matrix$teststat9)
abline(v= teststat9.real.wl.sexanddiet, col="red")


teststat10.real.wl.sexanddiet <- (F.Insect.real.wl.sexanddiet - M.Insect.real.wl.sexanddiet)
teststat10.real.wl.sexanddiet
(length(perm_matrix$teststat10[perm_matrix$teststat10 > teststat10.real.wl.sexanddiet])*2)/1000
hist(perm_matrix$teststat10)
abline(v= teststat10.real.wl.sexanddiet, col="red")


means.real.wl.diet <- indvARWLnopreg %>% group_by(Diet) %>% summarise(meanWL = mean(Wing_Loading))

Frug.real.wl.diet <- subset(means.real.wl.diet, Diet == "Frugivore")$meanWL
Glean.real.wl.diet <- subset(means.real.wl.diet, Diet == "Gleaner")$meanWL
Insect.real.wl.diet <- subset(means.real.wl.diet, Diet == "Aerial")$meanWL

teststat11.real <- (((Insect.real.wl.diet + Glean.real.wl.diet)/2)-(Frug.real.wl.diet))
teststat11.real
(length(perm_matrix$teststat11[perm_matrix$teststat11 < teststat11.real])*2)/1000
hist(perm_matrix$teststat11)
abline(v= teststat11.real, col="red")


teststat12.real <- (Insect.real.wl.diet - Glean.real.wl.diet)
teststat12.real
(length(perm_matrix$teststat12[perm_matrix$teststat12 > teststat12.real])*2)/1000
hist(perm_matrix$teststat12)
abline(v= teststat12.real, col="red")


