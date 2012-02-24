#--------------------------------------------------------------------------------
# Name:         IBM_Simple.r
# Purpose:      Spatially-explicit IBM stochastic simulation. Stand-alone script. 
# Author:       Francesco Tonini
# Email: 	    f_tonini@hotmail.com
# Created:      11/10/2011
# Copyright:    (c) 2011 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.14.0 64-bit version(http://www.r-project.org/)
#--------------------------------------------------------------------------------------

##Define a working directory
##If you are in a Mac OSX/Linux environment make sure to specify the appropriate path
##Always use either / or \\ to specify the path 
setwd ('~/Desktop/Temp')

#Path to folders in which you want to save all your vector & raster files
workdir_Raster <- '~/Desktop/Temp/Raster/'
workdir_Shapefiles <- '~/Desktop/Temp/Shapefiles'

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myfunctions.R')

###Let's set all simulation parameters:

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
	
##Define an empty vector in which we will store the total area covered at each time step
area <- rep(0, (end_time - start_time) + 1)

##Create a Tk window element as a container
root <- tktoplevel()

##Define the label of the progress bar
l <- tk2label(root)
##Define the length of the progress bar
pb <- tk2progress(root, length = 300)
##Configure the second progress bar's max value (min is always = 0)
##Max will be the number of subdivisions - 1 (since index starts at 0)
tkconfigure(pb, value=0, maximum = 2)

##Pack all labels and progress bars inside the Tk container
tkpack(l)
tkpack(pb)

##Refresh the Tk window
tcl('update')

##Let's store all starting colonies in a new dataset
##This way we do overwrite it and can still use it later
colonies <- starting_colonies

##LOOP for each year of simulation
for (year in seq(start_time,end_time)){
	
	##Configure the progress bar by updating its label and current value
	tkconfigure(l, text = paste('Year', year))
    tkconfigure(pb, value = 0)
	tcl('update')

	if (year == start_time) {
		
		##Save current colonies in a shapefile (.shp) file
		##Using our shp.save() module
		shp.save(workdir_Shapefiles, start_time)
		
		##Estimate the approx. area covered by current colonies in Km^2
		##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
		Tot_area <- raster.save(colonies, start_time, writeRas=TRUE)
		area[1] <- Tot_area
			
		##Remove objects not needed from the memory
		rm(Tot_area)
			
		##Update progress bar (move progress segment to the end of the bar)
		tkconfigure(pb, value = 2)
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
		##(other scenarios can be defined within the script 'myfunctions.R' and can be manually changed)
		pop <- offspring.generate(survival_prob = 0.01, male_prob = 0.5, scenario = 'optimistic', dist.mean = 200)
		
		##Configure the progress bar by updating its value 
		##(the green bar will move 1 segment forward)
		tkconfigure(pb, value = 1)
		tcl('update')
		
		##Only if there are AT LEAST two swarmers check for Male-Female within a 'radius' buffer
		if (nrow(pop) > 0){
			
			##Check for Male-Female within a 'radius' buffer
			new_colonies <- new.colonies.create()
						
			##Store new colonies in the main colony dataset
			colonies <- new.colonies.stack()
			
			##Configure the progress bar by updating its value 
			##(the green bar will move 1 segment forward)
			tkconfigure(pb, value = 2)
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
				Tot_area <- raster.save(colonies,year,writeRas=FALSE)
				area[(year-start_time)+1] <- Tot_area
				##EXIT the inner loop
				break
			}
		
			##Remove objects not needed from the memory 
			rm(list=c('pop', 'new_colonies'))
		}	
		
		##Save current colonies in a .shp file
		shp.save(workdir_Shapefiles, year)
		
		##Estimate the approx. area covered by current colonies in Km^2 (Point to Raster operation is done to estimate the area)
		##Also, if writeRas=TRUE a raster file (.grd by default) is produced
		Tot_area <- raster.save(colonies,year,writeRas=TRUE)
		area[(year-start_time)+1] <- Tot_area
		
		##Remove objects not needed from the memory
		rm(Tot_area)

	}

}
	
cat('\nSimulation Is Over!\n')

##Print on screen a summary of statistics
summary.stats()











