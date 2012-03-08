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

##Define the main working directory
##Make sure to specify the appropriate path using either / or \\ to specify the path 
mainDir <- 'C:/Temp'

##Let's set the working directory
setwd(file.path(mainDir))

#Path to folders in which you want to save all your vector & raster files
workdir_Raster <- 'Raster'
workdir_Shapefiles <- 'Shapefiles'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
dir.create(file.path(mainDir, workdir_Raster), showWarnings = FALSE)
dir.create(file.path(mainDir, workdir_Shapefiles), showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('myfunctions.R')

##Load all required libraries
print('Loading required libraries...')
load.packages()

###Let's set all simulation parameters:

##First Year of Simulation (input from terminal console)
cat('\nWhat is the first year of simulation?\n')
start_time <- scan(n=1)

##Last Year of Simulation (input from terminal console)
cat('\nWhat Is The Last Year Of Simulation?\n')
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
starting_colonies <- habitat.survival(starting_colonies)

##Define a Simulation Extent
extent <- simulation.extent()

##Check for max density of source colonies, based on the user input 
##This module returns a list (to access list elements use $ sign)
cat('\nType the maximum colony density in colonies per hectare (1ha = 100 sq. meters)\n')
maxdensity <- scan(n=1)
lst <- max.density.source()

##lst is a list with outputs returned by the max.density.source() module
starting_colonies <- lst$starting_colonies
maxdensity <- lst$maxdensity 
	
##Define an empty vector in which we will store the total area covered at each time step
area <- rep(0, (end_time - start_time) + 1)

##If we only have ONE source of invasion, define an empty vector in which we will store the 
##mean Eucl. dist. of all existing colonies from it, at each time step
if (nrow(starting_colonies) == 1) avgdistSource <- rep(0, (end_time - start_time) + 1) 

##Create a Tk window element as a container
root <- tktoplevel()
tktitle(root) <- 'Single-Run Simulation'

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
		shp.save(start_time)
		
		##Estimate the approx. area covered by current colonies in Km^2
		##If writeRas=TRUE a raster file (.img or as defined in source file) is produced
		Tot_area <- raster.save(colonies, start_time, writeRas=TRUE)
		area[(year-start_time)+1] <- Tot_area
		
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
		if (any(colonies$Age >= ColAge_swarmers)) pop <- offspring.generate(survival_prob = 0.01, male_prob = 0.5, scenario = 'optimistic', dist.mean = 200) else pop <- NULL
		
		##Configure the progress bar by updating its value 
		##(the green bar will move 1 segment forward)
		tkconfigure(pb, value = 1)
		tcl('update')
		
		##Only if there are AT LEAST two swarmers check for Male-Female within a 'radius' buffer
		if (!is.null(pop)){
			
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
					Tot_area <- raster.save(colonies,year,writeRas=FALSE)
					area[(year-start_time)+1] <- Tot_area
					##EXIT the inner loop
					break
				}
			
			}
								
			##Remove objects not needed from the memory 
			rm(list=c('pop', 'new_colonies'))	
		}

		##Configure the progress bar by updating its value 
		##(the green bar will move 1 segment forward)
		tkconfigure(pb, value = 2)
		tcl('update')
		
		##Save current colonies in a .shp file
		shp.save(year)
		
		##Estimate the approx. area covered by current colonies in Km^2 (Point to Raster operation is done to estimate the area)
		##Also, if writeRas=TRUE a raster file (.grd by default) is produced
		Tot_area <- raster.save(colonies,year,writeRas=TRUE)
		area[(year-start_time)+1] <- Tot_area
		
		##If we are starting with ONE source of invasion, calculate the mean Eucl. dist. of all colonies 
		##from the source
		if (nrow(starting_colonies) == 1) avgdistSource[(year-start_time)+1] <- dist.source(colonies)
		
		##Remove objects not needed from the memory
		rm(Tot_area)

	}

}

tkmessageBox(title = "Message", message = 'The simulation is over!', icon = "info", type = "ok")	

##Suppress the progress bar
tkdestroy(root)

##Print on screen a summary of statistics
summary.stats()



