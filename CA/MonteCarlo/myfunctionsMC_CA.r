#--------------------------------------------------------------------------------
# Name:         myfunctionsMC_CA.r
# Purpose:      Modules (functions) called by the main script
# Author:       Francesco Tonini
# Email:        f_tonini@hotmail.com
# Created:      09/20/2013
# Copyright:    (c) 2013 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.15 64-bit version(http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

##CALL PACKAGES MODULE
##This module is used to call all required packages
load.packages <- function()
{
	
	#setRepositories(ind=1:2)
	
	pkg <- c("rgdal","raster","msm")
	w <- which(pkg %in% row.names(installed.packages()) == FALSE)
	if (length(w) > 0) install.packages(pkg)[w] 
	#install.packages(c("rgdal","raster"))
	
	update.packages(pkg, ask = FALSE, dependencies = c('Suggests'))
	
	library(raster)			#Raster operation and I/O
	library(rgdal)		
	#library(msm)			#for Truncated Normal Distribution
	
	cat('\nAll libraries have been installed/loaded!\n')
}

##HABITAT/BACKGROUND LAYER MODULE
##This module is used to read a background landscape/habitat layer
read.NSHabitat <- function(layer, grain, PolyToRaster = FALSE)
{

	cat('Reading unsuitable habitat layer...\n')
	
	##Define all accepted file formats (OGR & GDAL)
	OGR_ext <- c('shp')
	GDAL_ext <- c('grd','asc','sdat','rst','nc','tif','envi','bil','img')
	
	##Strip the extension from the file name
	file.ext <- unlist(strsplit(layer, '\\.'))[2]
	
	##Strip the name from the file name
	file.name <- unlist(strsplit(basename(layer), '\\.'))[1]
		
	##If the file extension is .shp read the vector file using
	##the readOGR() function from the rgdal package
	##Otherwise the routine assumes it is a raster file and uses
	##the readGDAL() function
	if (any(file.ext %in% OGR_ext)) { 
	
		#For 'rgdal' library all supported OGR data formats are listed under ogrDrivers()
		habitat_block <- readOGR(dsn=dirname(layer), layer = file.name)
		
		##Check if the background layer is georeferenced
		if (!is.na(proj4string(habitat_block))){
		
			##If the coord. system of the layer is Geographic ("LAT-LON")
			##Ask the user to define parameters of a Projected coord. system
			if(substring(proj4string(habitat_block),7,13) == 'longlat'){
				
				cat('UTM zone (e.g. 17):\n')
				UTM_zone <- scan(n=1,what='')
				cat('Ellipsoid (e.g. GRS80):\n')
				Ellips <- scan(n=1,what='')
				cat('Datum (e.g. NAD83):\n')
				Datum <- scan(n=1,what='')
				
				CRS <- paste(' +proj=utm +zone=',UTM_zone,' +ellps=',Ellips,' +datum=',Datum,' +units=m +no_defs +towgs84=0,0,0',sep='')
				
				#To Change Projection and/or datum use the spTransform() function within package 'rgdal'
				habitat_block <- spTransform(habitat_block,CRS=CRS)
			
			}
			
		}else{
			
			cat('UTM zone (e.g. 17):\n')
			UTM_zone <- scan(n=1,what='')
			cat('Ellipsoid (e.g. GRS80):\n')
			Ellips <- scan(n=1,what='')
			cat('Datum (e.g. NAD83):\n')
			Datum <- scan(n=1,what='')
				
			CRS <- paste(' +proj=utm +zone=',UTM_zone,' +ellps=',Ellips,' +datum=',Datum,' +units=m +no_defs +towgs84=0,0,0',sep='')
		
			proj4string(habitat_block) <- CRS
		}
		
		if (PolyToRaster == TRUE){
		
			#get bounding box
			left <- habitat_block@bbox[1,1]
			right <- habitat_block@bbox[1,2]
			bottom <- habitat_block@bbox[2,1]
			top <- habitat_block@bbox[2,2]
			
			#create an empty raster container
			r <- raster(xmn=left,xmx=right,ymn=bottom,ymx=top)
			
			h <- (r@extent@xmax - r@extent@xmin) / grain
			if (!is.integer(h)) new.xmx <- r@extent@xmin + round(h) * grain; r@extent@xmax <- new.xmx
			v <- (r@extent@ymax - r@extent@ymin) / grain
			if (!is.integer(v)) new.ymx <- r@extent@ymin + round(v) * grain; r@extent@ymax <- new.ymx
			
			res(r) <- grain
			habitat_block <- rasterize(habitat_block, r, field=1)
			
		}
		
	}else{
	
		#If the extension is not .shp, read as raster file (if possible) gdalDrivers()
		habitat_block <- raster(layer)
		
		##Check if the background layer is georeferenced
		if (!is.na(proj4string(habitat_block))){
		
			##If the coord. system of the layer is Geographic ("LAT-LON")
			##Ask the user to define parameters of a Projected coord. system
			if(substring(proj4string(habitat_block),8,14) == 'longlat'){
				
				cat('UTM zone (e.g. 17):\n')
				UTM_zone <- scan(n=1,what='')
				cat('Ellipsoid (e.g. GRS80):\n')
				Ellips <- scan(n=1,what='')
				cat('Datum (e.g. NAD83):\n')
				Datum <- scan(n=1,what='')
				
				CRS <- paste(' +proj=utm +zone=',UTM_zone,' +ellps=',Ellips,' +datum=',Datum,' +units=m +no_defs +towgs84=0,0,0',sep='')
				
				#To Change Projection and/or datum use the spTransform() function within package 'rgdal'
				habitat_block <- spTransform(habitat_block,CRS=CRS)
			
			}
			
		}else{
			
			cat('UTM zone (e.g. 17):\n')
			UTM_zone <- scan(n=1,what='')
			cat('Ellipsoid (e.g. GRS80):\n')
			Ellips <- scan(n=1,what='')
			cat('Datum (e.g. NAD83):\n')
			Datum <- scan(n=1,what='')
				
			CRS <- paste(' +proj=utm +zone=',UTM_zone,' +ellps=',Ellips,' +datum=',Datum,' +units=m +no_defs +towgs84=0,0,0',sep='')
		
			proj4string(habitat_block) <- CRS
		}
		
		#if resolution of raster habitat file is different than our simulation grain
		#it has to be adjusted
		if (any(res(habitat_block) != grain)){
			
			#get bounding box
			left <- habitat_block@extent@xmin
			right <- habitat_block@extent@xmax
			bottom <- habitat_block@extent@ymin
			top <- habitat_block@extent@ymax
			
			#create an empty raster container
			r <- raster(xmn=left,xmx=right,ymn=bottom,ymx=top)
			
			h <- (r@extent@xmax - r@extent@xmin) / grain
			if (!is.integer(h)) new.xmx <- r@extent@xmin + round(h) * grain; r@extent@xmax <- new.xmx
			v <- (r@extent@ymax - r@extent@ymin) / grain
			if (!is.integer(v)) new.ymx <- r@extent@ymin + round(v) * grain; r@extent@ymax <- new.ymx
			
			res(r) <- grain
			habitat_block <- resample(habitat_block, r, method='ngb')
			writeRaster(habitat_block, filename='./habitat_block_resample',format='HFA', datatype='LOG1S', overwrite=TRUE)
		
		}
		
	}	
	
	if (PolyToRaster == TRUE) writeRaster(habitat_block, filename="./habitat_block", format='HFA', datatype='LOG1S', overwrite=TRUE)
		
	cat('Done!\n')
	
	return(habitat_block)
}

generate.Kernel <- function(grain, KernelSize, xy = c(0,0), FUN, alpha){

	
	if (KernelSize %% 2 == 0) {
		warning('You specified a moving window with an even number of cells! Size would be increased to nearest odd number of cells!')
		KernelSize <- KernelSize + 1
	}
	
	r <- raster(xmn= -(grain * KernelSize/2), xmx= grain * KernelSize/2, ymn= -(grain * KernelSize/2), ymx= grain * KernelSize/2)
	res(r) <- grain 
	distKernel <- distanceFromPoints(r, xy) 
	
	if (FUN == 'exp'){
		c <- 1
		distKernel[] <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(distKernel[]/alpha)^c  ) 
		distKernel[] <- distKernel[] / sum(distKernel[])
	}else if (FUN == 'gauss'){
		c <- 2
		distKernel[] <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(distKernel[]/alpha)^c  ) 
		distKernel[] <- distKernel[] / sum(distKernel[])
	}

	return(distKernel)

}

generate.Area <- function(ext, n_inf_pixel = NULL, filename = NULL){
	
	out <- list()
	if (ext$Right - ext$Left < grain * KernelSize | ext$Top - ext$Bottom < grain * KernelSize) stop("The extent of the chosen study area is less than the chosen distKernel size!\nPlease choose larger dimensions!")
	
	if (!is.null(filename)) {
	
		cat('\nReading Input File...\n')
	
		##Define all accepted file formats
		Extensions <- c('txt','csv')
		
		##Strip the extension from the file name
		file.ext <- unlist(strsplit(filename, '\\.'))[2]
				
		if (file.ext == 'txt') {
			starting_colonies <- as.matrix(read.delim(file = filename, header = TRUE, stringsAsFactors = FALSE)) #For TAB-delimited files
		}else if (file.ext == 'csv'){
			starting_colonies <- as.matrix(read.table(file = filename, header = TRUE, stringsAsFactors = FALSE, sep=',')) #For Comma-delimited files
		}
		
		##Grab the header row and change it to upper-case, as a default
		header <- toupper(colnames(starting_colonies))
		correct_labels <- c('LAT','LON','LNG','X','Y')
		
		##Check for label inconsistencies and/or errors:
				
		##Store the two coordinate values into 2 different variables (regardless of whether they are lat-lon or x-y)
		FirstCoord <- starting_colonies[,substring(header,1,3) == 'LAT' | substring(header,1,3) == 'Y']
		SecondCoord <- starting_colonies[,substring(header,1,3) == 'LON' | substring(header,1,3) == 'LNG' | substring(header,1,3) == 'X']
		
		##Build a matrix with the 2 coordinates
		MatrCoord <- cbind(as.numeric(FirstCoord),as.numeric(SecondCoord)) 
		
		Coord1Order <- which(substring(header,1,3) == 'LAT' | substring(header,1,3) == 'Y')
		Coord2Order <- which(substring(header,1,3) == 'LON' | substring(header,1,3) == 'LNG' | substring(header,1,3) == 'X')
		
		##If the first header label is either LAT or Y, switch columns in the input file to have first LON/LNG/X and then LAT/Y
		if (Coord1Order < Coord2Order) MatrCoord <- cbind(as.numeric(SecondCoord),as.numeric(FirstCoord)) 
		
		starting_colonies <- MatrCoord
		
		##Remove all duplicate points from the uploaded file
		if(any(duplicated(MatrCoord))) {
		
			exclude <- duplicated(MatrCoord)
			starting_colonies <- starting_colonies[-which(exclude),]
			cat(paste('\nRemoved', sum(exclude),'Duplicate Records!'))
			cat('\n')
		}
		
		##Call columns 'X' and 'Y'
		colnames(starting_colonies) = c('X','Y')
		
		##If coord. system is geographic project coords according to the background habitat-layer
		if ( any(substring(header,1,3) == "LAT") & (any(substring(header,1,3) == "LON") | any(substring(header,1,3) == "LNG") ) ){
			
			cat('Detected geographic coordinates!\n')
			cat('Projecting coordinates...\n')
			
			##Use the project() function from 'rgdal' package to project geog.coordinates...use coord system of the habitat backgroud layer
			Matr_proj <- project(starting_colonies, proj4string(habitat_block), inv = FALSE)
				
			starting_colonies[,1] <- Matr_proj[,1]
			starting_colonies[,2] <- Matr_proj[,2]
			cat('Coordinates projected according to the background layer!\n')
		
		}
		
		cat('Creating colonies dataset...\n')
		
		##Convert the input matrix to a dataframe
		starting_colonies <- as.data.frame(starting_colonies)
		
		##Rasterize points
		r <- raster()
		extent(r) <- c(ext$Left, ext$Right, ext$Bottom, ext$Top)
		res(r) <- grain		
		z <- rasterize(starting_colonies[,1:2], r, background=0, fun='count')
		z[z > maxdensity] <- maxdensity
		
		cat('Dataset created!\n')
	
	}else{
	
		z <- raster()
		extent(z) <- c(ext$Left, ext$Right, ext$Bottom, ext$Top)
		res(z) <- grain
		z[] <- 0
		
		a <- rasterToPoints(z)
		cond <- apply(a >= ((ext$Right - ext$Left) / 2) - 1000 & a <= ((ext$Right - ext$Left) / 2) + 1000, 1, sum)
		if (n_inf_pixel > length(cond)) stop('number of initial infested pixels cannot exceed available number of cells!') 
		ss <- sample(which(cond == 2), n_inf_pixel)
		z[ss] <- sample(1:maxdensity, n_inf_pixel, replace=T)
		
	}
	
	out$cl <- sapply(rasterToPoints(z)[,3], FUN=list)
	out$age <- lapply(out$cl, FUN=function(x){if(x > 0) sample(0:MaxAge, x, replace=T)})
	out$rr <- z
	
	return(out)

}

habitat.survival <- function(a, z){
	
	out <- list()
	w <- habitat_block[] == 1
	z[w] <- 0

	out$rr <- z	
	out$cl <- sapply(rasterToPoints(z)[,3], FUN=list)
	out$age <- mapply(FUN=function(x,y){if(x == 0) y <- NULL; return(list(y))}, out$cl, a)

	return(out)

}

max.age_fun <- function(x){
		
	if (!is.null(x) & any(x > 20)) {
		x <- x[-which(x > MaxAge)] 
		if(length(x) == 0) x <- NULL
	}else{x}

	return(x)
}

alates.gen <- function(x, scenario='optimistic'){
	
		if (scenario == 'optimistic'){
		
			z <- ifelse(x < ColAge_swarmers, 0, 100000 * Survival)
			z <- ifelse(x >= ColAge_swarmers & x < 10, 1000 * Survival, z)
			z <- ifelse(x >= 10 & x < 15, 10000 * Survival, z)
			
		}else if (scenario == 'pessimistic'){
		
			z <- ifelse(x < ColAge_swarmers, 0, 100000 * Survival)
			z <- ifelse(x >= ColAge_swarmers & x < ColAge_swarmers * 2, 10000 * Survival, z)
			z <- ifelse(x >= ColAge_swarmers * 2 & x < ColAge_swarmers * 3, 50000 * Survival, z)
		
		}
	
	return(z)
}

newCol.gen <- function(x, tab){
	if (!is.na(x) & x > 1){
		idx <- sample(1:1000, 1)
		int <- findInterval(x, tab[,1])
		sub <- tab[(int - 1000 + 1):int,3]
		x <- sub[idx]
	}else{x <- 0}
	return(x)
}

newCol.addAge <- function(l1,l2,l3){

	out <- list()

	for (i in 1:length(l1)){
		
		if(l3[[i]] == 0) next
		if(l2[[i]] < maxdensity) {
			if (l2[[i]] + l3[[i]] > maxdensity) {
				l1[[i]] <- append(l1[[i]], rep(0, maxdensity - l2[[i]])) 
				l2[[i]] <- length(l1[[i]]) 
			}else{
				l1[[i]] <- append(l1[[i]], rep(0, l3[[i]]))
				l2[[i]] <- length(l1[[i]]) 
			}
		}

	}

	out$age <- l1
	out$cl <- l2

	return(out)

}

##SUMMARY STATISTICS MODULE
##This module is used to calculate required stastistics on simulation results
summary.stats <- function()
{
		
	#Now let's have a look at the final stats
	simYears <- seq(start_time,end_time)
	
	#Store the information about the area covered after each time step 
	Area_dataset <- data.frame(Year = simYears)
		
	SD <- sapply(area.dataset, sd)
	MEAN <- sapply(area.dataset, mean)
	Area_dataset$MEAN <- round(MEAN,2)
	Area_dataset$SD <- round(SD,2)
	Area_dataset$CV <- round(SD / MEAN, 2)
	
	if (nrow(Area_dataset) > 1) {
		
		Area_dataset$areaSpeed <- 0
		Area_dataset$percGrowth <- 0
		
		for (i in 2:nrow(Area_dataset)){
					
			Area_dataset$areaSpeed[i] <- round(MEAN[i] - MEAN[i-1], 2) #Km^2 / year				
			Growth <- (MEAN[i] - MEAN[i-1]) / MEAN[i-1]
			Area_dataset$percGrowth[i] <- round(Growth * 100)
			
			if (n_inf_pixel == 1) {
				
				Area_dataset$EuclDistSpeed[i] <- round(MEAN2[i] - MEAN2[i-1], 2)  #Km / year	
				GrowthEucl <- (MEAN2[i] - MEAN2[i-1]) / MEAN2[i-1]
				if (!is.finite(GrowthEucl)) GrowthEucl <- 0
				Area_dataset$EuclDistGrowth[i] <- round(GrowthEucl * 100)
			}
			
		}
		
		colnames(Area_dataset) <- colnames(Area_dataset) <- c('Years','Avg. Area (Km^2)', 'St.Dev. Area (Km^2)',
										'CV Area', 'Avg. Km^2/year', 'Avg. Areal Growth (%)')
		var1 <- c('Years', Area_dataset[,1])
		var2 <- c('Avg_Area_Km^2', Area_dataset[,2])
		var3 <- c('St.Dev_Area_Km^2', Area_dataset[,3])
		var4 <- c('CV_Area', Area_dataset[,4])
		var5 <- c('Avg.Km^2/year', Area_dataset[,5])
		var6 <- c('Avg.Areal_Growth(%)', Area_dataset[,6])
		
		myArray <- c(var1,var2,var3,var4,var5,var6)	
	
		
	}else{
	
		colnames(Area_dataset) <- c('Years','Avg. Area (Km^2)', 'St.Dev. Area (Km^2)','CV Area')
		var1 <- c('Years', Area_dataset[,1])
		var2 <- c('Avg_Area_Km^2', Area_dataset[,2])
		var3 <- c('St.Dev_Area_Km^2', Area_dataset[,3])
		var4 <- c('CV_Area', Area_dataset[,4])
		
		myArray <- c(var1,var2,var3,var4)

	}
	
	#This part is run only in the case of an invasion starting from ONE source
	#So that we can use the average euclidean distance from the center (invasion point)
	#when comparing the expasion rate of simulation Vs. theoretical uniform distribution
	###if (n_inf_pixel == 1)###
	#Radial Increase general formula: [sqrt(A1)/pi - sqrt(A0)/pi ] / t1- t0	
	
	#Final message of simulation is over
	print('Final statistics saved to main folder!')

	write.table(Area_dataset, './Final_Stats.csv', row.names=F, sep=',')

		
}

##OCCUPANCY ENVELOPE MODULE
##This module is used to compute occupancy envelopes across MonteCarlo simulation runs
##and save raster files accordingly
envelope.raster <- function()
{
	
	#Create a variable with all the simulation years, repeated based on the nbr. of MonteCarlo runs
	#etiqYear <- as.character(rep(seq(start_time,end_time),each=NRuns))
	
	#Create a variable with all the names of the raster files created during the MonteCarlo runs
	#etiqRun <- paste('Run',seq(NRuns),'_',etiqYear,sep='')
		
	#LOOP through each simulation year
	for (yrs in start_time:end_time){
		
		RasDir <- paste('./', fOutput, '/', sep='')
		#Store full path and name of all raster files corresponding to the current year of the LOOP 
		etiq <- paste(RasDir,'Run',seq(NRuns),'_',yrs,'.img',sep='')
				
		#Use the stack() function of the 'raster' package to stack all layers stored in 'etiq'
		#Note: our rasters only have 1 band...in case you have multispectral images you can add that
		Ras_stack <- stack(etiq)
				
		##OVERLAY FUNCTION ('raster package')
		#Overlay all rasters of the stack and sum values for each pixel/cell
		Overlap <- overlay(Ras_stack, fun=function(x){sum(x)})
				
		#Based on what the user chose as occupancy envelope, calculate the final raster and save it
		#We only want to keep pixels that are occupied AT LEAST ONCE across all MonteCarlo runs
		MaskValue <- NRuns * 0
		Overlap1 <- Overlap
		Overlap1[Overlap1 == MaskValue] <- NA
		Overlap1[Overlap1 > MaskValue] <- 1
		writeRaster(Overlap1,filename=paste(RasDir,'Occupancy',yrs,'_0','.img',sep=''),
					format='HFA', datatype='LOG1S', overwrite=TRUE)
				
		#We only want to keep pixels that are occupied AT LEAST 25% of all MonteCarlo runs
		#MaskValue <- NRuns * 0.25
		#Overlap2 <- Overlap
		#Overlap2[Overlap2 < MaskValue] <- NA
		#Overlap2[Overlap2 >= MaskValue] <- 1
		#writeRaster(Overlap2,filename=paste(RasDir,'Occupancy',yrs,'_25','.img',sep=''),
		#			format='HFA', datatype='LOG1S', overwrite=TRUE)
				
		#We only want to keep pixels that are occupied AT LEAST 50% of all MonteCarlo runs
		MaskValue <- NRuns * 0.5
		Overlap3 <- Overlap
		Overlap3[Overlap3 < MaskValue] <- NA
		Overlap3[Overlap3 >= MaskValue] <- 1
		writeRaster(Overlap3,filename=paste(RasDir,'Occupancy',yrs,'_50','.img',sep=''),
					format='HFA', datatype='LOG1S', overwrite=TRUE)
			
		#We only want to keep pixels that are occupied AT LEAST 75% of all MonteCarlo runs
		#MaskValue <- NRuns * 0.75
		#Overlap4 <- Overlap
		#Overlap4[Overlap4 < MaskValue] <- NA
		#Overlap4[Overlap4 >= MaskValue] <- 1
		#writeRaster(Overlap4,filename=paste(RasDir,'Occupancy',yrs,'_75','.img',sep=''),
		#			format='HFA', datatype='LOG1S', overwrite=TRUE)
				
		#We only want to keep pixels that are occupied across ALL MonteCarlo runs
		MaskValue <- NRuns 
		Overlap5 <- Overlap
		Overlap5[Overlap5 < MaskValue] <- NA
		Overlap5[Overlap5 == MaskValue] <- 1
		writeRaster(Overlap5,filename=paste(RasDir,'Occupancy',yrs,'_100','.img',sep=''),
					format='HFA', datatype='LOG1S', overwrite=TRUE)
			
		
	}
	
	##Delete all single raster files (comment this line in case you want to keep all of them)
	do.call(file.remove,list(list.files(RasDir, pattern='Run',full.names=TRUE)))
	
	##Now you have all raster files ready to be visualized in any GIS software, showing areas/pixels occupied
	##In more than X% of the MonteCarlo simulation runs
	print('occupancy raster has been saved to the Raster folder')
	print('Done')
		
	#END

}

