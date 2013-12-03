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

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
mainDir <- 'C:/Temp'

##Let's set the working directory
setwd(file.path(mainDir))

#Path to folders in which you want to save all your vector & raster files
workdir_Raster <- 'Raster'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(file.path(mainDir, workdir_Raster), showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myfunctionsMC.R')

##Load all required libraries
print('Loading required libraries...')
load.packages()

###Let's set all simulation parameters:

##Number of Monte Carlo runs
cat('\nHow many MonteCarlo simulation runs?\n')
NRuns <- scan(n=1)

##First Year of Simulation (input from terminal console)
cat('\nWhat is the first year of simulation?\n')
start_time <- scan(n=1)

##Last Year of Simulation (input from terminal console)
cat('\nWhat is the last year of simulation?\n')
end_time <- scan(n=1)
while(end_time < start_time) {
	tkmessageBox(title = "Warning", message = 'The last year of simulation cannot be inferior to the first year!...Please type again'
				, icon = "warning", type = "ok")
	end_time <- scan(n=1)
}

##Set the age at which colonies start producing first swarmers
cat('\nAt what age do colonies generate the first swarmers?\n')
ColAge_swarmers <- scan(n=1)

##Set the maximum life expectancy for a colony (input from terminal console)
cat('\nWhat is the maximum life expectancy of a colony (in years)?\n')
MaxAge <- scan(n=1)

##Apply a radius representing the max attraction distance btw Males & Females
cat('\nType max attraction distance between males & females (in meters)?\n')
radius <- scan(n=1)
while(radius == 0) {
	tkmessageBox(title = "Warning", message = 'The attraction distance must be > 0!...Please type again'
				, icon = "warning", type = "ok")
	radius <- scan(n=1)
}
 
##Read background non-suitable habitat layer
habitat_block <- read.NSHabitat()

##Read the Input File w/ Starting Points (colonies)
##By default it assigns a random age to all input colonies
starting_colonies <- read.file(random.age = 'random')

##Check habitat suitability and remove colonies falling within non-suitable habitat
##starting_colonies <- habitat.survival(starting_colonies)

##Define a Simulation Extent
extent <- simulation.extent()

##Check for max density of source colonies, based on the user input 
cat('\nType the maximum colony density in colonies per hectare (1ha = 100 sq. meters)\n')
maxdensity <- scan(n=1)
lst <- max.density.source() ##This module returns a list (to access list elements use $ sign)
starting_colonies <- lst$starting_colonies
maxdensity <- lst$maxdensity 

##Create a database to store occupied areas over time & across all MC simulation runs
area.dataset <- as.data.frame(matrix(0,ncol=1,nrow=NRuns))
colnames(area.dataset) <- start_time
rownames(area.dataset) <- paste('Run',seq(NRuns),sep='')

##Save shapefile of non duplicated starting colonies
##year <- 2003
##shp.save(year)

##If we only have ONE source of invasion, define an empty vector in which we will store the 
##mean Eucl. dist. of all existing colonies from it, at each time step
if (nrow(starting_colonies) == 1) {
	avgdistSource <- as.data.frame(matrix(0,ncol=1,nrow=NRuns))
	colnames(avgdistSource) <- start_time
	rownames(avgdistSource) <- paste('Run',seq(NRuns),sep='')
}

##Create a Tk window element as a container
##root <- tktoplevel()
##tktitle(root) <- 'MonteCarlo Simulation'

##Define the label of the first progress bar
##l1 <- tk2label(root)
##Define the length of the first progress bar
##pb1 <- tk2progress(root, length = 300)
##Configure the first progress bar's max value (min is always = 0)
##Max will be the number of MC simulation runs - 1 (since index starts at 0)
##tkconfigure(pb1, value=0, maximum = NRuns-1)

##Define the label of the second progress bar
##l2 <- tk2label(root)
##Define the length of the second progress bar
##pb2 <- tk2progress(root, length = 300)
##Configure the second progress bar's max value (min is always = 0)
##Max will be the number of subdivisions - 1 (since index starts at 0)
##tkconfigure(pb2, value=0, maximum = 2)

##Pack all labels and progress bars inside the Tk container
##tkpack(l1)
##tkpack(pb1)
##tkpack(l2)
##tkpack(pb2)

##Refresh the Tk window
##tcl('update')

##LOOP for each MC simulation run
for (run in seq(NRuns)){
	
	##Configure the MC Runs progress bar by updating its label and current value
	##tkconfigure(l1, text = paste('Run', run))
    ##tkconfigure(pb1, value = run - 1)
	
	##Let's now create an object called colonies (we do not delete starting_colonies because it can still be used)
	colonies <- starting_colonies
	
	##print RUN on screen
	print(paste('Run',run))
	
	##LOOP for each year of simulation
	for (year in seq(start_time,end_time)){
		
		##tkconfigure(l2, text = paste('Year',year))
        ##tkconfigure(pb2, value = 0)
		##tcl('update')
		
		##print YEAR on screen
		print(paste('Year',year))
		
		if (year == start_time) {
			
			##Estimate the approx. area covered by current colonies in Km^2
			##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
			Tot_area <- raster.save(colonies, year, writeRas=TRUE)
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
			if (any(colonies$Age >= ColAge_swarmers)) {
				
				out <- offspring.generate(survival_prob = 0.01, male_prob = 0.5, scenario = 'optimistic', dist.mean = 200) 
				pop <- out$pop
				flag <- out$flag
				
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
				
				}
			
			}else{ out <- NULL }
			
			##Configure the progress bar by updating its value 
			##(the green bar will move 1 segment forward)
			tkconfigure(pb2, value = 1)
			tcl('update')
			
			##Only if there are AT LEAST two swarmers check for Male-Female within a 'radius' buffer
			if (!is.null(out)){
				
				##Check for Male-Female within a 'radius' buffer
				new_colonies <- new.colonies.create()
						
				##Store new colonies in the main colony dataset
				if (nrow(new_colonies) > 0) {
			
					colonies <- new.colonies.stack()
					
					##Check habitat suitability of current colonies
					colonies <- habitat.survival(colonies)
					
					##Uncomment one of the following modules when having only ONE source of invasion
					##and you are interested in keeping ONLY colonies laying on the fringe:		
					##colonies <- convex.hull(colonies) 
					
					##If there are no colonies left exit the LOOP, set the new end_time to the current year
					##and EXIT the LOOP
					if (nrow(colonies) == 0) {
						print(paste('No colonies survived at year', year))
						Tot_area <- raster.save(colonies,year,writeRas=TRUE)
						area.dataset[run,year-start_time+1] <- Tot_area
						colnames(area.dataset)[year-start_time+1] <- year
						maxyear <- year
						end_time <- maxyear
						##EXIT the inner loop and go to the next MC simulation run
						break
					}
			
				}
			
				##Remove objects not needed from the memory 
				rm(list=c('pop', 'new_colonies'))
			}	
			
			##Configure the progress bar by updating its value 
			##(the green bar will move 1 segment forward)
			tkconfigure(pb2, value = 2)
			tcl('update')
		
			##Estimate the approx. area covered by current colonies in Km^2
			##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
			Tot_area <- raster.save(colonies,year,writeRas=TRUE)
			area.dataset[run,year-start_time+1] <- Tot_area
			colnames(area.dataset)[year-start_time+1] <- year
			
			##If we are starting with ONE source of invasion, calculate the mean Eucl. dist. of all colonies 
			##from the source
			if (nrow(starting_colonies) == 1) {
				avgdistSource[run,year-start_time+1] <- dist.source(colonies)
				colnames(avgdistSource)[year-start_time+1] <- year
			}
			
			##Remove objects not needed from the memory
			rm(Tot_area)
		
		}
		
	}	
	
	save.image()
}	

#Final message of simulation is over
tkmessageBox(title = "Message", message = 'The simulation is over!', icon = "info", type = "ok")	

##Suppress the progress bar
tkdestroy(root)

##Print on screen a summary of statistics across all MC runs 
summary.stats()

##Ask whether the user wants to compute and save an 'occupancy envelope' or not 
envelope.raster()
	







