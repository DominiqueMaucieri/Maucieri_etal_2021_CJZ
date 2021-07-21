
****

Data and R code accompanying:

<b>Sexual dimorphism in bat wing morphology – variation among foraging styles</b>

DOI: 10.1139/cjz-2021-0035

Authors: [Dominique G. Maucieri](https://dominiquemaucieri.com), [Austin J. Ashbaugh](https://austinashbaugh.wixsite.com/austinjashbaugh) and [Jessica M. Theodor](https://people.ucalgary.ca/~jtheodor/)

****

[analyses folder](analyses/)

In the order they were conducted:

* [PCA.R](analyses/PCA.R) script to conduct PCA analysis of wing shape landmarks
* [Multivariate.R](analyses/Multivariate.R) script for multivariate analysis of wing shape landmark data
* [CVA.R](analyses/CVA.R) script for CVA analysis of wing shape landmarks
* [AR_WL.R](analyses/AR_WL.R) script for analysis of aspect ratio and wing loading data
* [Permutation.R](analyses/Permutation.R) script for permutation test examining aspect ratio and wing loading data
* <b>[figures folder](analyses/figures)</b> figures used in the manuscript and supplemental materials which were produced using the above scripts

[data folder](data/)

* [landmarks.RData](data/landmarks.RData)
	* ```coords```: coordinates for all 19 landmarks for each individual bat
	
* [classifier.RData](data/classifier.RData)
	* ```Number```: Bat IDs
	* ```Country```: Country that each bat was captured in
	* ```Bat_Mass```: Mass of each bat (g)
	* ```Species```: Species of each bat
	* ```Diet```: Dominant foraging style by each bat
		* ```f``` = frugivore
		* ```a``` = aerial insectivore
		* ```g``` = gleaning insectivore
	* ```Sex```: Sex of each bat
		* ```F``` = Female
		* ```M``` = Male
	* ```Reproductive_Status```: Reproductive status of each bat
		* ```NR``` = Not in a reproductive stage (Males and Females)
		* ```TD``` = Testes have descended (Males)
		* ```NOP``` = Not obviously pregnant (Females)
		* ```P``` = Pregnant (Females)
		* ```L``` = Lactating (Females)
		* ```PL``` = Post lactation (Females)
	* ```Sex_w_pregnant```: Sex of each bat with whether females are pregnant or not pregnant
		* ```FP``` = Female (pregnant)
		* ```FNP``` = Female (not pregnant)
		* ```M``` = Male
	* ```Wing_Length```: Length of bat wing (cm)
	* ```Wing_Area```: Area of bat wing (cm^2)
	* ```SexandDiet```: Sex and diet grouping variable

* [filelist.RData](data/filelist.RData)
	* ```filelist```: list of filenames that link to the order of landmark coordinates in ```coords``` data
	
* [PCA.RData](data/PCA.RData)
	* ```wingpca```: PCA output from [PCA.R](analyses/PCA.R)
	
* [PCSCORES.RData](data/PCSCORES.RData)
	* ```ID```: Bat IDs
	* ```PC1```: Principle component 1 scores for each bat
	* ```PC2```: Principle component 2 scores for each bat
	* ```Guild```: Sex and diet grouping variable
	* ```Sex```: Sex of each bat
		* ```F``` = Female
		* ```M``` = Male
	* ```Diet```: Dominant foraging style by each bat
		* ```f``` = frugivore
		* ```a``` = aerial insectivore
		* ```g``` = gleaning insectivore
