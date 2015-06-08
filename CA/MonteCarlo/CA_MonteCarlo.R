#--------------------------------------------------------------------------------
# Name:         CA_MonteCarlo.r
# Purpose:      MonteCarlo Cellular Automata Simulation test file. Stand-alone script. 
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      09/20/2013
# Copyright:    (c) 2013 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

#install packages
#install.packages(c("rgdal","raster","lubridate","rgrass7","optparse", "plotrix"))

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
setwd("~/My Dropbox/SHARED_SERVER")

#Path to folders in which you want to save all your vector & raster files
fOutput <- 'Output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myfunctionsMC_CA.R')

##Load all required libraries
print('Loading required libraries...')
load.packages()

###Let's set all simulation parameters:

##Number of Monte Carlo runs
cat('\nHow many MonteCarlo Simulation Runs?\n')
NRuns <- scan(n=1)
#NRuns <- 100

##First Year of Simulation 
cat('\nWhat Is The First Year Of Simulation?\n')  
start_time <- scan(n=1)  #(input from terminal console)
#start_time <- 2003

##Last Year of Simulation 
cat('\nWhat Is The Last Year Of Simulation?\n')
end_time <- scan(n=1)  #(input from terminal console)
#end_time <- 2020

##Set the age at which colonies start producing first alates
cat('\nAt What Age Do Colonies Generate First Alates?\n')
ColAge_swarmers <- scan(n=1) #(input from terminal console)
#ColAge_swarmers <- 4

##Check for max density of source colonies, based on the user input 
##This module returns a list (to access list elements use $ sign)
cat('\nType The Maximum Colony Density --colonies per hectare (1ha = 100 sq. meters)\n')
maxdensity <- scan(n=1)
#maxdensity <- 7

##Set the maximum life expectancy for a colony
cat('\nWhat Is The Maximum Life Expectancy Of a Colony (In Years)?\n')
MaxAge <- scan(n=1)  #(input from terminal console)
#MaxAge <- 20

##Set survival probability of termite alates 
cat('\nWhat Is The Survival Probability of Termite Alates (Type Fraction)?\n')
Survival <- scan(n=1)  #(input from terminal console)
#Survival <- 0.01

##Set male-female ratio of colony reproductives 
cat('\nWhat Is Male-Female Ratio of Reprodutives in The Colony (Type Fraction)?\n')
SexRatio  <- scan(n=1)  #(input from terminal console)
#SexRatio <- 0.5

#Read look-up table for number of new colonies
nCol_table <- read.table('./NewCol_table.csv', header = T, stringsAsFactors = F, sep=',')
nCol_table <- nCol_table[nCol_table$ratio == SexRatio, ]

##choose number of pixels to be randomly infested across the study area
#n_inf_pixel <- 20

grain <- 100 #cell resolution in METERS...each cell corresponds to a single termite colony
KernelSize <- 31 #cells =  (e.g. KernelSize = 10 --->  5 * 100 = 500 x 500 m window)

#Let's generate a Euclidean Distance Kernel (in an automatic manner)
#and transform it into a Probability Kernel (by choosing a reference prob. distr.)
#then choose an average distance (e.g. 200meters)
probKernel <- generate.Kernel(grain, KernelSize, FUN = "exp", alpha = 200)   #use "exp" for exponential kernel and "gauss" for gaussian kernel 								 
#probKernel <- generate.Kernel(grain, KernelSize, FUN = "tnorm", avg = 200, sd = 200)   #or truncated normal distribution

unsuitable_habitat <- 'habitat_block.shp'
##Read background non-suitable habitat layer
habitat_block <- read.NSHabitat(unsuitable_habitat, grain, PolyToRaster = TRUE)

##OR
#unsuitable_habitat <- 'habitat_block.img'
#habitat_block <- read.NSHabitat(unsuitable_habitat, grain)

#study area extent
ext <- list(Left = 0, Right = grain * KernelSize * 10, Bottom = 0, Top = grain * KernelSize * 10)																 
## Let's define the extent of study area
## !!define projection based on unsuitable habitat
init <- generate.Area(ext, filename = './coords.txt')  #use n_inf_pixel alternatively to filename = ... 

##Create a database to store occupied areas over time & across all MC simulation runs
area.dataset <- as.data.frame(matrix(NA, ncol=1, nrow=NRuns))
colnames(area.dataset) <- start_time
rownames(area.dataset) <- paste('Run',seq(NRuns),sep='')

Tot_area <- cellStats(init$rr > 0, stat='sum') * grain  #in units of the grain (meters)
Tot_area_km2 <- Tot_area/1000000  
area.dataset[,1] <- Tot_area_km2 #the starting year the area is always the same in all MonteCarlo model runs

##LOOP for each MC simulation run
for (run in seq(NRuns)){
		
	simArea <- init$rr 	
	Age.lst <- init$age 
	Colonies.lst <- init$cl
	Alates <- simArea
	Alates[] <- 0
	
	print(paste('Run',run))
		
	for (year in seq(start_time,end_time)){

		print(paste('Year',year))

		if (year == start_time) {
				
			writeRaster(simArea > 0, filename=paste('./',fOutput,'/Run',run,'_',year,'.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)
			
			#plot(simArea > 0, main = paste("Run ", run, "\n", year, sep=''))
			
		}else{	
			
			##Increase age by 1 in those cells with at least one colony inside
			Age.lst <- lapply(Age.lst, FUN=function(x){ if (!is.null(x)) x <- x + 1 })
			
			##Remove colonies over MaxAge...
			Age.lst <- lapply(Age.lst, FUN=max.age_fun)
			
			##Whenever a single element within age list has no more colonies inside (died b/c of age),
			##update the corresponding list element in colony list and set it == 0 
			Colonies.lst <- mapply(FUN=function(x,y){ if(is.null(x)) y <- 0 else y; return(list(y)) } , Age.lst, Colonies.lst) 
			
			if ( !any(unlist(Colonies.lst) > 0) ) {
			
				print(paste('No colonies survived at year', year))
				simArea[] <- unlist(Colonies.lst)
				writeRaster(simArea > 0, filename=paste('./',fOutput,'/Run',run,'_',year,'.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)
				##Estimate the approx. area covered by current colonies in Km^2
				Tot_area <- cellStats(simArea > 0, stat='sum') * grain  #in units of the grain (meters)
				Tot_area_km2 <- Tot_area/1000000  
				area.dataset[run,year-start_time+1] <- Tot_area_km2 
				colnames(area.dataset)[year-start_time+1] <- year
				maxyear <- year
				end_time <- maxyear
				##EXIT the inner loop and go to the next MC simulation run
				break
			}
			
			##Generate total number of offspring per cell by using the age of each existing colony
			Alates.lst <- lapply(Age.lst, FUN=function(x){if (!is.null(x)) alates.gen(x, scenario='optimistic')})
			Alates.count <- lapply(Alates.lst, FUN=function(x){if(is.null(x))x <- 0 else sum(x)})
			
			Alates[] <- unlist(Alates.count)
			Alates <- round(focal(Alates, w=as.matrix(probKernel), fun=sum))
			
			##remove alates from unsuitable habitat
			Alates[habitat_block[] == 1] <- 0
			
			NewCol <- Alates
			##use look-up table now to decide how many new colonies based on alates number in each cell
			NewCol[] <- sapply(Alates[], FUN=function(x) newCol.gen(x, tab=nCol_table))

			NewCol.lst <- sapply(rasterToPoints(NewCol)[,3], FUN=list)	
			
			out <- newCol.addAge(l1 = Age.lst, l2 = Colonies.lst, l3 = NewCol.lst)		
			Age.lst <- out$age
			Colonies.lst <- out$cl
			
			rm(NewCol.lst)
			
			simArea[] <- unlist(Colonies.lst)
			writeRaster(simArea > 0, filename=paste('./',fOutput,'/Run',run,'_',year,'.img',sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)
			
			#plot(simArea > 0, main = paste("Run ", run, "\n", year, sep=''))
			
			##Estimate the approx. area covered by current colonies in Km^2
			##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
			Tot_area <- cellStats(simArea > 0, stat='sum') * grain  #in units of the grain (meters)
			Tot_area_km2 <- Tot_area/1000000  
			area.dataset[run,year-start_time+1] <- Tot_area_km2 #the starting year the area is always the same in all MonteCarlo model runs
			colnames(area.dataset)[year-start_time+1] <- year	
			
			##Remove objects not needed from the memory
			rm(Tot_area)
		}
	}
	save.image('CA_MC.RData')
}

##Print on screen a summary of statistics across all MC runs 
summary.stats()

##Ask whether the user wants to compute and save an 'occupancy envelope' or not 
envelope.raster()






