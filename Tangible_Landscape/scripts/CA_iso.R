#--------------------------------------------------------------------------------
# Name:         CA_iso.r
# Purpose:      Lattice-based simulation of the spread of invasive termite N. corniger over a heterogeneous landscape. 
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      09/20/2013
# Copyright:    (c) 2013 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

#install packages
#install.packages(c("rgdal","raster","lubridate","rgrass7","optparse","plotrix"))

#load packages:
suppressPackageStartupMessages(library(raster))    #Raster operation and I/O. Depends R (≥ 2.15.0)
suppressPackageStartupMessages(library(rgdal))     #Geospatial data abstraction library. Depends R (≥ 2.14.0)
suppressPackageStartupMessages(library(lubridate)) #Make dealing with dates a little easier. Depends R (≥ 3.0.0)
suppressPackageStartupMessages(library(rgrass7))   #Interface Between GRASS 7 GIS and R. Depends R (≥ 2.12)
suppressPackageStartupMessages(library(optparse))  #Parse args from command line
suppressPackageStartupMessages(library(plotrix))   #Add text annotations to plot

##Define the main working directory based on the current script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
base_name <- dirname(sub(file_arg, "", initial_options[grep(file_arg, initial_options)]))
setwd(paste(sep="/", base_name, ".."))

#Path to folders in which you want to save all your vector & raster files
#fOutput <- 'output'

##Create a physical copy of the subdirectory folder(s) where you save your output
##If the directory already exists it gives a warning, BUT we can suppress it using showWarnings = FALSE
#dir.create(fOutput, showWarnings = FALSE)

##Use an external source file w/ all modules (functions) used within this script. 
##Use FULL PATH if source file is not in the same folder w/ this script
source('./scripts/myfunctions_CA.r')

###Input simulation parameters: #####
option_list = list(
  make_option(c("-h","--habitat"), action="store", default=NA, type='character', help="input suitable habitat raster map (mask)"),
  make_option(c("-src","--sources"), action="store", default=NA, type='character', help="initial sources of infection shapefile"),
  make_option(c("-img","--image"), action="store", default=NA, type='character', help="background satellite raster image for plotting"),
  make_option(c("-s","--start"), action="store", default=NA, type='integer', help="start year"),
  make_option(c("-e","--end"), action="store", default=NA, type='integer', help="end year"),
  make_option(c("-t","--tab"), action="store", default=NA, type='character', help="look-up table for mating probability"),
  make_option(c("-a","--age"), action="store", default=4, type='integer', help="age of fist alates production"),
  make_option(c("-md","--maxd"), action="store", default=1, type='integer', help="max colony density per grid cell"),
  make_option(c("-ma","--maxage"), action="store", default=20, type='integer', help="max colony age"),
  make_option(c("-sv","--surv"), action="store", default=0.01, type='numeric', help="alates survival rate"),
  make_option(c("-sx","--sex"), action="store", default=0.5, type='numeric', help="sex ratio (male prevalence)"),
  make_option(c("-kd","--kdist"), action="store", default=NA, type='numeric', help="max distance (meters) covered by the dispersal kernel window"),
  make_option(c("-kt","--ktype"), action="store", default='exp', type='character', help="dispersal kernel distribution"),
  make_option(c("-o","--output"), action="store", default=NA, type='character', help="basename for output GRASS raster maps"),
)

opt = parse_args(OptionParser(option_list=option_list))

##Input raster --> SUITABLE HABITAT MASK:
##Read background non-suitable habitat layer
habitat_block <- readRAST(opt$habitat)
habitat_block <- raster(habitat_block)  #transform 'sp' obj to 'raster' obj
#habitat_block <- raster('./layers/habitat_block')

##Initial sources of infestation:
src_pnt <- readVECT(opt$sources)
#src_pnt <- readOGR('./layers', layer = 'init_colonies')

##background satellite image for plotting
bkr_img <- raster(paste('./layers/', opt$image, sep='')) 


##Start-End date: 
start <- opt$start
end <- opt$end
if (start > end) stop('start date must precede end date!!')

#build time series for simulation steps:
dd_start <- as.POSIXlt(as.Date(paste(start,'-01-01',sep='')))
dd_end <- as.POSIXlt(as.Date(paste(end,'-01-01',sep='')))
tstep <- as.character(seq(dd_start, dd_end, 'years'))

#create formatting expression for padding zeros depending on total number of steps
formatting_str <- paste("%0", floor( log10( length(tstep) ) ) + 1, "d", sep='')

##Age at which colonies start producing first alates
ColAge_swarmers <- opt$age
#maximum colony density
maxdensity <- opt$maxd
##Set the maximum life expectancy for a colony
MaxAge <- opt$maxage
##Set survival probability of termite alates 
Survival <- opt$surv
##Set male-female ratio of colony reproductives 
SexRatio  <- opt$sex  

##Read look-up table for number of new colonies
nCol_table <- read.table(paste(opt$tab,'.csv', sep=''), header = T, stringsAsFactors = F, sep=',')
nCol_table <- nCol_table[nCol_table$ratio == SexRatio, ]

##Dispersal kernel:
#max cut-off distance for the 'moving window'
if(is.na(opt$kdist)) stop('max kernel distance parameter needs to be specified!')
if(opt$kdist < res(habitat_block)[1]) stop('max kernel distance cannot be less than spatial resolution of study area!')
if(opt$kdist > (extent(habitat_block)[2]-extent(habitat_block)[1])) warning('dispersal kernel size larger than study area: please consider reducing it!')

if(!opt$ktype %in% c('exp', 'gauss')) step('kernel distribution type must be equal to "exp" or "gauss"')
probKernel <- generate.Kernel(habitat_block, opt$kdist, FUN=opt$ktype, alpha=200) 

#create raster from initial infestation points
colonies_rast <- rasterize(src_pnt@coords[,1:2], habitat_block, background=0, fun='count')
#list of colony counts per grid cell
colonies_lst <- sapply(colonies_rast[], FUN=list)
#cell ID for pixels with colonies
idx_col <- cellFromXY(colonies_rast, src_pnt@coords[,1:2])
#list of all cell IDs for simulation raster
idx <- sapply(1:ncell(colonies_rast), FUN=list)
#list of colony ages for each cell
age <- lapply(idx, FUN = function(x) if (x %in% idx_col) src_pnt@data$Age[x == idx_col] else NULL)

#open window screen
windows(width = 10, height = 10, xpos=350, ypos=50, buffered = FALSE)
#quartz()  #use this on Mac OSX
#x11()     #use this on Linux (not tested!)

#plot background image
plot(bkr_img)
#plot coordinates for plotting text:
xpos <- (bbox(I_rast)[1,2] + bbox(I_rast)[1,1]) / 2
ypos <- bbox(I_rast)[2,2] - 150

simArea <- init$rr 	
Age.lst <- init$age 
Colonies.lst <- init$cl
Alates <- simArea
Alates[] <- 0
	
##LOOP for each year
for (tt in tstep){
  
  #split date string for raster time stamp
  split_date <- unlist(strsplit(tt, '-'))

	if (tt == tstep[1]) {
		
	  if(!any(unlist(colonies_lst) > 0)) stop('Simulation ended. There are no more termite colonies on the landscape!')
    
	  ##CALCULATE OUTPUT TO PLOT:
	  # 1) values as % colonized (infected) hectare
	  colonies_rast[] <- ifelse(colonies_rast[] == 0, NA, round(colonies_rast[]/7, 1))
    
    # 2) number of colonies per hectare
    colonies_rast[] <- ifelse(colonies_rast[] == 0, NA, colonies_rast[])
    
	  # 3) values as 0 (non infected) and 1 (infected) cell
	  #colonies_rast[] <- ifelse(colonies_rast[] > 0, 1, 0) 
	  #colonies_rast[] <- ifelse(colonies_rast[] > 0, 1, NA) 
    
	  #PLOT: overlay current plot on background image
	  #bks <- seq(0, 100, by=10)
	  #bks <- seq(0, mx, length = 10)
	  bks <- c(0, 0.25, 0.5, 0.75, 1)
	  #colors <- c("yellow","gold","orange","red")
	  image(colonies_rast, breaks=bks, col=rev(heat.colors(length(bks)-1, alpha=1)), axes=F, ann=F, useRaster=T,add=T, xaxs = "i", yaxs = "i")
	  boxed.labels(xpos, ypos, tt, bg="white", border=NA, font=2)
    
	  #WRITE TO FILE:
	  colonies_rast_sp <- as(colonies_rast, 'SpatialGridDataFrame')
	  writeRAST(colonies_rast_sp, vname=paste(opt$output, '_', sprintf(formatting_str, 0), sep=''), overwrite=TRUE) #write to GRASS raster file
	  execGRASS('r.timestamp', map=paste(opt$output, '_', sprintf(formatting_str, cnt), sep=''), date=paste(split_date[3], months_names[as.numeric(split_date[2])], split_date[1]))
	  
	  #writeRaster(colonies_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='FLT4S', overwrite=TRUE) # % infected as output
	  #writeRaster(colonies_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='INT1U', overwrite=TRUE) # nbr. colonies as output
	  #writeRaster(colonies_rast, filename=paste('./', fOutput, '/', opt$output, '_', sprintf(formatting_str, 0), sep=''), format='HFA', datatype='LOG1S', overwrite=TRUE)  # 0=non infected 1=infected output
			
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

##Print on screen a summary of statistics across all MC runs 
summary.stats()

##Ask whether the user wants to compute and save an 'occupancy envelope' or not 
envelope.raster()






