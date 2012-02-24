#--------------------------------------------------------------------------------
# Name:         IBM_MonteCarlo.r
# Purpose:      MonteCarlo Simulation. Stand-alone script. 
# Author:       Francesco Tonini
# Email: 	    f_tonini@hotmail.com
# Created:      11/10/2011
# Copyright:    (c) 2011 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.14.0 64-bit version(http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

##Define a working directory
##If you are in a Mac OSX/Linux environment make sure to specify the appropriate path
##Always use either / or \\ to specify the path 
setwd ('~/Desktop/Temp')

#Path to a folder in which you want to save all your raster files
workdir_Raster <- '~/Desktop/Temp/Raster/'

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myfunctionsMC.R')

###Let's set all simulation parameters:

##Number of Monte Carlo runs
cat('\nHow many MonteCarlo Simulation Runs?\n')
NRuns <- scan(n=1)

##First Year of Simulation (input from terminal console)
cat('\nWhat Is The First Year Of Simulation?\n')
start_time <- scan(n=1)

##Last Year of Simulation (input from terminal console)
cat('\nWhat Is The Last Year Of Simulation?\n')
end_time <- scan(n=1)
while(end_time < start_time) {cat('\n!!Last Year Of Simulation Cannot Be Inferior To The First Year!!\n'); end_time <- scan(n=1)}

##Set the age at which colonies start producing first swarmers
cat('\nAt What Age Do Colonies Generate First Swarmers?\n')
ColAge_swarmers <- scan(n=1)

##Set the maximum life expectancy for a colony (input from terminal console)
cat('\nWhat Is The Maximum Life Expectancy Of a Colony (In Years)?\n')
MaxAge <- scan(n=1)

##Apply a radius representing the max attraction distance btw Males & Females
cat('\nType Max Attraction Distance Between Males & Females (In meters)?\n')
radius <- scan(n=1)
while(radius == 0) {cat('Attraction Distance Must Be > 0\n'); radius <- scan(n=1)}
 
##Load all required libraries
print('Loading Required Libraries...')
load.packages()

##Read background non-suitable habitat layer
habitat_block <- read.NSHabitat()

##Read the Input File w/ Starting Points (colonies)
##By default it assigns a random age to all input colonies
starting_colonies <- read.file(random.age = 'random')

##Check habitat suitability and remove colonies falling within non-suitable habitat
starting_colonies <- habitat.survival(starting_colonies)

##Define a Simulation Extent
extent <- simulation.extent()

##Check for max density of source colonies, based on the user input 
##This module returns a list (to access list elements use $ sign)
cat('\nType The Maximum Colony Density --colonies per hectare (1ha = 100 sq. meters)\n')
maxdensity <- scan(n=1)
lst <- max.density.source()

##lst is a list with outputs returned by the max.density.source() module
starting_colonies <- lst$starting_colonies
maxdensity <- lst$maxdensity 

##Create a database to store occupied areas over time & across all MC simulation runs
area.dataset <- as.data.frame(matrix(0,ncol=1,nrow=NRuns))
colnames(area.dataset) <- start_time
rownames(area.dataset) <- paste('Run',seq(NRuns),sep='')

##Create a Tk window element as a container
root <- tktoplevel()

##Define the label of the first progress bar
l1 <- tk2label(root)
##Define the length of the first progress bar
pb1 <- tk2progress(root, length = 300)
##Configure the first progress bar's max value (min is always = 0)
##Max will be the number of MC simulation runs - 1 (since index starts at 0)
tkconfigure(pb1, value=0, maximum = NRuns-1)

##Define the label of the second progress bar
l2 <- tk2label(root)
##Define the length of the second progress bar
pb2 <- tk2progress(root, length = 300)
##Configure the second progress bar's max value (min is always = 0)
##Max will be the number of subdivisions - 1 (since index starts at 0)
tkconfigure(pb2, value=0, maximum = 2)

##Pack all labels and progress bars inside the Tk container
tkpack(l1)
tkpack(pb1)
tkpack(l2)
tkpack(pb2)

##Refresh the Tk window
tcl('update')

##LOOP for each MC simulation run
for (run in seq(NRuns)){
	
	##Configure the MC Runs progress bar by updating its label and current value
	tkconfigure(l1, text = paste('Run', run))
    tkconfigure(pb1, value = run - 1)
	
	##Let's now create an object called colonies (we do not delete starting_colonies because it can still be used)
	colonies <- starting_colonies
	
	##LOOP for each year of simulation
	for (year in seq(start_time,end_time)){
		
		tkconfigure(l2, text = paste('Year',year))
        tkconfigure(pb2, value = 0)
		tcl('update')
		
		if (year == start_time) {
			
			##Estimate the approx. area covered by current colonies in Km^2
			##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
			Tot_area <- raster.save(colonies,year,writeRas=TRUE)
			area.dataset[run,year-start_time+1] <- Tot_area
			colnames(area.dataset)[year-start_time+1] <- year
			
			##Remove objects not needed from the memory
			rm(Tot_area)
			
			##Update progress bar (move progress segment to the end of the bar)
			tkconfigure(pb2, value = 2)
			tcl('update')
			
		}else{
		
			##Increase Age and Timecount by 1
			colonies <- age.increase()
			
			##Only keep current colonies (remove the past year)
			colonies <- colonies[which(colonies$Timecount == year),]
			
			##Colonies older than MaxAge die (i.e. removed from dataset)
			if (length(which(colonies$Age > MaxAge)) > 0) colonies = colonies[-which(colonies$Age > MaxAge),]
			
			##Generate offspring for each existing colony
			##Change the scenario to 'pessimistic' if colonies generate more individuals at a younger age
			##(other scenarios can be defined within the script 'myfunctionsMC.R' and can be manually changed)
			out <- offspring.generate(survival_prob = 0.01, male_prob = 0.5, scenario = 'optimistic', dist.mean = 200)
			pop <- out$pop
			flag <- out$flag
			
			##Configure the progress bar by updating its value 
			##(the green bar will move 1 segment forward)
			tkconfigure(pb2, value = 1)
			tcl('update')
			
			##flag = 1 if ONE or MORE individuals crossed the simulation boundaries
			##If flag is not NULL print a warning on screen, EXIT (break) the loop and move to the next MC simulation run 
			if (!is.null(flag)) {
			
				cat('\n')
				cat('WARNING: ONE OR MORE INDIVIDUALS LAY OUTSIDE THE SIMULATION EXTENT!!\n')
				cat(paste('The Current & All Following Simulation Runs Will Be Terminated Before Year: ',year,'\n',sep=''))
				cat('\n')
				maxyear <- year - 1
				
				##Overwrite end_time with a value corresponding to the current year - 1
				##so that in the next MC run the simulation will stop at the same year
				##This is done to be able to compare across MC runs for same years
				end_time <- maxyear
				
				##EXIT the inner loop and go to the next MC simulation run
				break
			
			}else{ maxyear <- NULL}
			
			##Only if there are AT LEAST two swarmers check for Male-Female within a 'radius' buffer
			if (nrow(pop) > 0){
				
				##Check for Male-Female within a 'radius' buffer
				new_colonies <- new.colonies.create()
							
				##Store new colonies in the main colony dataset
				colonies <- new.colonies.stack()
				
				##Configure the progress bar by updating its value 
				##(the green bar will move 1 segment forward)
				tkconfigure(pb2, value = 2)
				tcl('update')
					
				##Uncomment one of the following modules when having only ONE source of invasion
				##and you are interested in keeping ONLY colonies laying on the fringe:	
				##colonies <- convex.hull(colonies) 
					
				##Check habitat suitability of current colonies
				colonies <- habitat.survival(colonies)
				
				##If there are no colonies left exit the LOOP, set the new end_time to the current year
				##and EXIT the LOOP
				if (nrow(colonies) == 0) {
					print(paste('No Colonies Survived At Year', year))
					Tot_area <- raster.save(colonies,year,writeRas=TRUE)
					area.dataset[run,year-start_time+1] <- Tot_area
					colnames(area.dataset)[year-start_time+1] <- year
					maxyear <- year
					end_time <- maxyear
					##EXIT the inner loop and go to the next MC simulation run
					break
				}
			
				##Remove objects not needed from the memory 
				rm(list=c('pop','new_colonies'))
			}	
			
			##Estimate the approx. area covered by current colonies in Km^2
			##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
			Tot_area <- raster.save(colonies,year,writeRas=TRUE)
			area.dataset[run,year-start_time+1] <- Tot_area
			colnames(area.dataset)[year-start_time+1] <- year
			
			##Remove objects not needed from the memory
			rm(Tot_area)
		
		}
		
	}
	
	save.image()
}	

cat('\nSimulation Is Over!\n')

##Print on screen a summary of statistics across all MC runs 
summary.stats()

##Ask whether the user wants to compute and save an 'occupancy envelope' or not 
envelope.raster()







