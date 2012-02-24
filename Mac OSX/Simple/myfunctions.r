#----------------------------------------------------------------------------
# Name:         myfunctions.r
# Purpose:      Modules (functions) called by the main script StandAloneIBM.r
# Author:       Francesco Tonini
# Email: 	    f_tonini@hotmail.com
# Created:      11/10/2011
# Copyright:    (c) 2011 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R-2.14.0 64-bit version(http://www.r-project.org/)
#-----------------------------------------------------------------------------------------

##CALL PACKAGES MODULE
##This module is used to call all required packages
load.packages <- function()
{
	
	setRepositories(ind=1:2)	
	#"msm","RANN","FNN","SDMTools","spdep","maptools"
	pkg <- c("tcltk2","spatstat","sp","rgdal","CircStats","raster")
	w <- which(pkg %in% row.names(installed.packages()) == FALSE)
	if (length(w) > 0) install.packages(pkg)[w]
	#install.packages(c('tcltk2','spatstat','sp','CircStats','raster','rgdal'))
	
	#Load (call) specif packages from the existing library collection into the current R session
	#library(msm) 	  		#for the truncated Normal distribution
	library(tcltk2)      	#for progress bar
	library(spatstat)  		#for quadratcount & disc buffers
	library(sp) 			#A package that provides classes and methods for spatial data
	library(rgdal) 			#rgdal: Bindings for the Geospatial Data Abstraction Library. Projections etc..
	library(CircStats) 		#for radians angle implementations (Keep this one otherwise not script will not work
							#in R CMD..not sure why?!) 
	library(raster)			#Raster operation and I/O
	#library(maptools)		#Shapefiles operation and I/O	
	
	cat('\nAll Libraries Have Been Loaded Successfully!\n')
}

##HABITAT/BACKGROUND LAYER MODULE
##This module is used to read a background landscape/habitat layer
read.NSHabitat <- function()
{
	##Read typed layer from terminal console
	cat('\nType the Name and Extension (e.g. .shp,.img,.grd, etc.) of the Non-Suitable Habitat Layer\n')
	layer_name <- scan(n=1,what='')

	cat('Reading Layer...\n')
	
	##Define all accepted file formats (OGR & GDAL)
	OGR_ext <- c('shp')
	GDAL_ext <- c('grd','asc','sdat','rst','nc','tif','envi','bil','img')
	
	##Strip the extension from the file name
	file.ext <- unlist(strsplit(layer_name, '\\.'))[2]
	
	##If there is no extension ask the user to include it
	while(is.na(file.ext)){
		cat('You Forgot To Add The File Extension!...Please Type Again:\n')
		layer_name <- scan(n=1,what='')
		file.ext <- unlist(strsplit(layer_name, '\\.'))[2]
	}
	
	##Strip the name from the file name
	file.name <- unlist(strsplit(layer_name, '\\.'))[1]
	
	##If the file extension is not one of the formats defined above
	##give error message
	while(!file.ext %in% c(OGR_ext,GDAL_ext)){
		cat('The extension of the background layer MUST be one of the following:\n')
		cat(c(OGR_ext,GDAL_ext))
		stop('Change format and restart the simulation')
	}
		
	##If the file extension is .shp read the vector file using
	##the readOGR() function from the rgdal package
	##Otherwise the routine assumes it is a raster file and uses
	##the readGDAL() function
	if (any(file.ext %in% OGR_ext)) { 
	
		#For 'rgdal' library all supported OGR data formats are listed under ogrDrivers()
		habitat_block <- readOGR('.', layer = file.name)
	
	}else{
	
		#If the extension is not .shp, read as raster file (if possible) gdalDrivers()
		habitat_block <- readGDAL(layer_name)
	}	
	
	##Check if the background layer is georeferenced
	if (!is.na(proj4string(habitat_block))){
	
		##If the coord. system of the layer is Geographic ("LAT-LON")
		##Ask the user to define parameters of a Projected coord. system
		if(substring(proj4string(habitat_block),8,14) == 'longlat'){
			
			cat('\nIn order to run the simulation the background layer must be projected\n')
			cat('Please Define The Following Projection Parameters:\n')
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
		
		cat('\nBackground layer not georeferenced...\n')
		cat('\nIn order to run the simulation the background layer must be projected\n')
		cat('Please Define The Following Projection Parameters:\n')
		cat('UTM zone (e.g. 17):\n')
		UTM_zone <- scan(n=1,what='')
		cat('Ellipsoid (e.g. GRS80):\n')
		Ellips <- scan(n=1,what='')
		cat('Datum (e.g. NAD83):\n')
		Datum <- scan(n=1,what='')
			
		CRS <- paste(' +proj=utm +zone=',UTM_zone,' +ellps=',Ellips,' +datum=',Datum,' +units=m +no_defs +towgs84=0,0,0',sep='')
	
		proj4string(habitat_block) <- CRS
	}
		
	cat('\nDone!\n')
	
	return(habitat_block)
}

##READ & FORMAT CHECK MODULE
##This module is used to read an input file containing 2D coordinates
##check that labels are well formatted, check the coordinate system,
##removes duplicate records, and builds a final dataset of starting colonies
read.file <- function(random.age)
{
	
	##Read typed layer from terminal console
	cat('\nType the Name and Extension (e.g. coords.txt) of the file \n')
	input <- scan(n=1,what='')
	
	cat('\nReading Input File...')
	
	##Define all accepted file formats
	Extensions <- c('txt','csv')
	
	##Strip the extension from the file name
	file.ext <- unlist(strsplit(input, '\\.'))[2]
	
	##If there is no extension ask the user to include it
	while(is.na(file.ext)){
		cat('You Forgot To Add The File Extension!...Please Type Again:\n')
		input <- scan(n=1,what='')
		file.ext <- unlist(strsplit(input, '\\.'))[2]
	}
		
	##If the file extension is not one of the formats defined above
	##give error message
	while(!file.ext %in% Extensions){
		cat('The extension of the background layer MUST be one of the following:\n')
		cat(Extensions)
		stop('Change format and restart the simulation')
	}
	
	if (file.ext == 'txt') {
		starting_colonies <- as.matrix(read.delim(file = input, header = TRUE, stringsAsFactors = FALSE)) #For TAB-delimited files
	}else if (file.ext == 'csv'){
		starting_colonies <- as.matrix(read.csv(file = input, header = TRUE, stringsAsFactors = FALSE)) #For Comma-delimited files
	}
	
	##Grab the header row and change it to upper-case, as a default
	header <- toupper(colnames(starting_colonies))
	correct_labels <- c('LAT','LON','LNG','X','Y')
	
	##Check for label inconsistencies and/or errors:
	
	##If both header labels do not correspond to any of the ones defined in 'correct_labels' give error message
	if( sum(substring(header,1,3) %in% correct_labels) != 2 ) stop('<<ERROR: Wrong Coordinate Labels in Uploaded File!>>')
	
	##If header labels are duplicates give error message
	if (sum(substring(header,1,3) == 'LAT') == 2 | sum(substring(header,1,3) == 'LON') == 2 | sum(substring(header,1,3) == 'LNG') == 2 
		| sum(header == 'X') == 2 | sum(header == 'Y') == 2 ) stop('<<ERROR: Incorrect Labels! Labels Cannot Be Identical!>>')
	
	##If the header labels do not correspond to LAT-LON/LON-LAT/LAT-LNG/LNG-LAT/X-Y/Y-X give error message
	if ((any(substring(header,1,3) == 'LAT') & any(header == 'Y')) | (any(substring(header,1,3) == 'LAT') & any(header == 'X')) 
		| (any(substring(header,1,3) == 'LON') & any(header == 'X')) | (any(substring(header,1,3) == 'LON') & any(header == 'Y'))
		| (any(substring(header,1,3) == 'LNG') & any(header == 'X')) | (any(substring(header,1,3) == 'LNG') & any(header == 'Y'))){
		
		stop('<<ERROR: Wrong Coordinate Labels in Uploaded File!>>')	
	}
	
	##Store the two coordinate values into 2 different variables (regardless of whether they are lat-lon or x-y)
	FirstCoord <- starting_colonies[,substring(header,1,3) == 'LAT' | substring(header,1,3) == 'Y']
    SecondCoord <- starting_colonies[,substring(header,1,3) == 'LON' | substring(header,1,3) == 'LNG' | substring(header,1,3) == 'X']
	
	##Build a matrix with the 2 coordinates
	MatrCoord <- cbind(as.numeric(FirstCoord),as.numeric(SecondCoord)) 
	
	##Remove all duplicate points from the uploaded file
	if(any(duplicated(MatrCoord))) {
	
		exclude <- duplicated(MatrCoord)
		starting_colonies <- starting_colonies[-which(exclude),]
		cat(paste('\nRemoved', sum(exclude),'Duplicate Records!'))
		cat('\n')
	}
	
	Coord1Order <- which(substring(header,1,3) == 'LAT' | substring(header,1,3) == 'Y')
	Coord2Order <- which(substring(header,1,3) == 'LON' | substring(header,1,3) == 'LNG' | substring(header,1,3) == 'X')
	
	##If the first header label is either LAT or Y, switch columns in the input file to have first LON/LNG/X and then LAT/Y
	if (Coord1Order < Coord2Order) MatrCoord <- cbind(as.numeric(SecondCoord),as.numeric(FirstCoord)) 
	
	starting_colonies <- MatrCoord
	
	##Call columns 'X' and 'Y'
	colnames(starting_colonies) = c('X','Y')
	
	##If coord. system is geographic project coords according to the background habitat-layer
	if ( any(substring(header,1,3) == "LAT") & (any(substring(header,1,3) == "LON") | any(substring(header,1,3) == "LNG") ) ){
		
		cat('Detected Geographic Coordinates!\n')
		cat('Projecting Coordinates...\n')
		
		##Use the project() function from 'rgdal' package to project geog.coordinates...use coord system of the habitat backgroud layer
		Matr_proj <- project(starting_colonies, proj4string(habitat_block), inv = FALSE)
			
		starting_colonies[,1] <- Matr_proj[,1]
		starting_colonies[,2] <- Matr_proj[,2]
		cat('Coordinates Projected According to the Background Layer!\n')
	
	}
	
	cat('Creating Colonies Dataset...\n')
	
	##Convert the input matrix to a dataframe
	starting_colonies <- as.data.frame(starting_colonies)
	
	##If there is ONLY ONE source of invasion then create a 'distance from the source' field
	##that will keep track of the expansion (in meters) from the original source
	if (nrow(starting_colonies) == 1) starting_colonies$Dist_source <- 0
	
	##Add a Timecount variable to the dataframe (initialized using the start_time value input by the user)
	starting_colonies$Timecount <- start_time
	##Add a colony ID variable to the dataframe
	starting_colonies$ID_col <- seq(from = 1, to = nrow(starting_colonies))
	
	##Check what age to be assigned to each colony and add an Age variable to the dataframe
	if (random.age == 'random') starting_colonies$Age <- sample(0:MaxAge,nrow(starting_colonies),replace=TRUE) else starting_colonies$Age <- as.integer(random.age)
	
	cat('Dataset Created!\n')
	
	return(starting_colonies)
}

##HABITAT SURVIVAL MODULE
##This module is used to check the survival of colonies across the simulation area. 
##Colonies falling within non-suitable areas are eliminated from the dataset
habitat.survival <- function(input)
{
	cat('\n')
	print('Checking Habitat Suitability...')
	
	SpatPntDataFrame <- input
	
	# Turn input data into a SpatialPointsDataFrame
	coordinates(SpatPntDataFrame) <- c('X','Y')
	
	# tell R that our coordinates are in the same reference system
	# as the polygon data -- BUT ONLY BECAUSE WE KNOW THIS IS THE CASE!
	proj4string(SpatPntDataFrame) <- proj4string(habitat_block)
	
	# combine is.na() with over() to do the containment test; note that we
	# need to 'denote' the layer to either a SpatialPolygons or SpatialGrid object first
	if(class(habitat_block)[1] == 'SpatialPolygonsDataFrame') {
		inside.layer <- which(!is.na(over(SpatPntDataFrame, habitat_block))) 
	}else{
		inside.layer <- which(!is.na(over(SpatPntDataFrame, habitat_block)))
	}
	
	if (length(inside.layer) > 0) {
	
		input <- input[-inside.layer,]
		input$ID_col <- seq(1,nrow(input))
		print('Removed Colonies From Non-suitable Habitat!')
	}
	
	print('Check Completed!')
	cat('\n')
	
	return(input)
}

##SIMULATION EXTENT MODULE
##This module is used to define the area extent for the simulation
simulation.extent <- function(cellsize=100) #default = 100 meters
{
	cat('\nDefine The Simulation Extent...\n')
	
	##Create a list (we will save the desired outputs inside) 
	extent <- list()
	
	cat('How Would You Like to Define The Simulation Extent?\n(Option 2 Recommended If Non-Suitable Habitat Layer Is Big!)\n')
	choice <- menu(c('Use The Whole Extent Of The Habitat-Suitability Layer','Specify a Distance From The Center of The Study Area'))
	if (choice == 1) {
		
		print('The simulation will stop as soon as a simulation boundary is reached!')
		
		if(class(habitat_block)[1] == 'SpatialPolygonsDataFrame') {
			Top <- habitat_block@bbox[2,2] #@bbox is an attribute of spatialdataframes 
			Left <- habitat_block@bbox[1,1]
			Bottom <- habitat_block@bbox[2,1]
			Right <- habitat_block@bbox[1,2]
		}else{
			Top <- extent(habitat_block)@ymax
			Left <- extent(habitat_block)@xmin
			Bottom <- extent(habitat_block)@ymin
			Right <- extent(habitat_block)@xmax
		}
		
		##Calculate the centroid (center) of the colony point cloud
		centroid.X <- mean(starting_colonies$X)
		centroid.Y <- mean(starting_colonies$Y)
		
		##Calculate the distance of bbox boundaries from the centroid
		DistanceX <- Right - centroid.X
		DistanceY <- Top - centroid.Y
		
		##Method to automatically find the smallest possible bounding box possible to be used as simulation extent:
		##use proportional increments to expand a grid around
		##the centroid until it reaches a distance not smaller than the farthest point (colony). Do it for X and Y coords
		
		##Increment in positive horizontal direction
		Delta1 <- round(DistanceX / cellsize)
		##Save the computed right simulation extent inside the list
		extent$Right <- (centroid.X + cellsize/2) + cellsize * Delta1
		##Increment in negative horizontal direction
		Delta2 <- round(-DistanceX / cellsize)
		##Save the computed left simulation extent inside the list
		extent$Left <- (centroid.X - cellsize/2) + cellsize * Delta2
		##Increment in positive vertical direction
		Delta3 <- round(DistanceY / cellsize)
		##Save the computed top simulation extent inside the list
		extent$Top <- (centroid.Y + cellsize/2) + cellsize * Delta3
		##Increment in negative vertical direction
		Delta4 <- round(-DistanceY / cellsize)
		##Save the computed bottom simulation extent inside the list
		extent$Bottom <- (centroid.Y - cellsize/2) + cellsize * Delta4
		
		##Calculate the nbr of columns of the simulation grid
		ncols <- length(seq(extent$Left,extent$Right,by=cellsize)) - 1
		##Calculate the nbr of rows of the simulation grid
		nrows <- length(seq(extent$Bottom,extent$Top,by=cellsize)) - 1
		
					
	}else{
		
		cat('Specify a BIG-ENOUGH extent...the simulation will stop as soon as a simulation boundary is reached!\n')
		centroid.X <- mean(starting_colonies$X)
		centroid.Y <- mean(starting_colonies$Y)
		
		cat('\nHorizontal Distance From The Center of The Point Cloud (in Km):\n')
		DistanceX <- scan(n=1)
		if (DistanceX <= 0) {cat('Error: The Distance Must Be Greater Than 0\n'); DistanceX <- scan(n=1)}
		cat('\nVertical Distance From The Center of The Point Cloud (in Km):\n')
		DistanceY <- scan(n=1)
		if (DistanceY <= 0) {cat('Error: The distance must be greater than 0\n'); DistanceY <- scan(n=1)}
		
		##Transform Km to meters
		DistanceX <- DistanceX * 1000
		DistanceY <- DistanceY * 1000
		
		##Save simulation extent
		Top <- centroid.Y + DistanceY
		Bottom <- centroid.Y - DistanceY
		Left <- centroid.X - DistanceX
		Right <- centroid.X + DistanceX
		
		#Method to automatically find the smallest envelope box possible
		#by using proportional increments and expanding the original cell around
		#the source point
		Delta1 <- round(DistanceX / cellsize)
		extent$Right <- (centroid.X + cellsize/2) + cellsize * Delta1
		Delta2 <- round(-DistanceX / cellsize)
		extent$Left <- (centroid.X - cellsize/2) + cellsize * Delta2
		Delta3 <- round(DistanceY / cellsize)
		extent$Top <- (centroid.Y + cellsize/2) + cellsize * Delta3
		Delta4 <- round(-DistanceY / cellsize)
		extent$Bottom <- (centroid.Y - cellsize/2) + cellsize * Delta4
		
		ncols <- length(seq(extent$Left,extent$Right,by=cellsize)) - 1
		nrows <- length(seq(extent$Bottom,extent$Top,by=cellsize)) - 1
		
		X <- starting_colonies$X
		Y <- starting_colonies$Y
		
		##Define the simulation extent as a mask
		msk <- data.frame(x = c(extent$Left,extent$Right,extent$Right,extent$Left), 
					  y = c(extent$Bottom,extent$Bottom,extent$Top,extent$Top))
		
		##Use the point.in.polygon() function ('sp' package) to check
		##which points fall inside the simulation extent mask
		idx <- point.in.polygon(point.x = starting_colonies$X, 
				  point.y = starting_colonies$Y, pol.x = msk[,1], pol.y = msk[,2])
		
		##Until one or more points fall outside the mask, ask the user to increase the distance from the centroid so that all points are included
		while (any(idx==0)) {
			
			cat('WARNING: one or more colonies lay outside the simulation extent you defined!\n')
			cat('Please increase the size:\n')
			cat('\nHorizontal Distance From The Center of The Point Cloud (in Km):\n')
			DistanceX <- scan(n=1)
			if (DistanceX <= 0) {cat('Error: The Distance Must Be Greater Than 0\n'); DistanceX <- scan(n=1)}
			cat('\nVertical Distance From The Center of The Point Cloud (in Km):\n')
			DistanceY <- scan(n=1)
			if (DistanceY <= 0) {cat('Error: The distance must be greater than 0\n'); DistanceY <- scan(n=1)}
			
			##Transform Km to meters
			DistanceX <- DistanceX * 1000
			DistanceY <- DistanceY * 1000
			
			##Save simulation extent
			Top <- centroid.Y + DistanceY
			Bottom <- centroid.Y - DistanceY
			Left <- centroid.X - DistanceX
			Right <- centroid.X + DistanceX
			
			##Method to automatically find the smallest possible bounding box possible to be used as simulation extent:
			##use proportional increments to expand a grid around
			##the centroid until it reaches a distance not smaller than the farthest point (colony). Do it for X and Y coords
			
			##Increment in positive horizontal direction
			Delta1 <- round(DistanceX / cellsize)
			##Save the computed right simulation extent inside the list
			extent$Right <- (centroid.X + cellsize/2) + cellsize * Delta1
			##Increment in negative horizontal direction
			Delta2 <- round(-DistanceX / cellsize)
			##Save the computed left simulation extent inside the list
			extent$Left <- (centroid.X - cellsize/2) + cellsize * Delta2
			##Increment in positive vertical direction
			Delta3 <- round(DistanceY / cellsize)
			##Save the computed top simulation extent inside the list
			extent$Top <- (centroid.Y + cellsize/2) + cellsize * Delta3
			##Increment in negative vertical direction
			Delta4 <- round(-DistanceY / cellsize)
			##Save the computed bottom simulation extent inside the list
			extent$Bottom <- (centroid.Y - cellsize/2) + cellsize * Delta4
			
			##Calculate the nbr of columns of the simulation grid
			ncols <- length(seq(extent$Left,extent$Right,by=cellsize)) - 1
			##Calculate the nbr of rows of the simulation grid
			nrows <- length(seq(extent$Bottom,extent$Top,by=cellsize)) - 1
			
			X <- starting_colonies$X
			Y <- starting_colonies$Y
			
			##Define the simulation extent as a mask
			msk <- data.frame(x = c(extent$Left,extent$Right,extent$Right,extent$Left), 
						  y = c(extent$Bottom,extent$Bottom,extent$Top,extent$Top))
			
			##Use the point.in.polygon() function ('sp' package) to check
			##which points fall inside the simulation extent mask
			idx <- point.in.polygon(point.x = starting_colonies$X, 
					  point.y = starting_colonies$Y, pol.x = msk[,1], pol.y = msk[,2])
			
		}
		
	}
	
	cat('\n')
	cat('Simulation Extent Defined:\n')
	cat(paste('Left:',round(extent$Left),'\n'))
	cat(paste('Right:',round(extent$Right),'\n'))
	cat(paste('Bottom:',round(extent$Bottom),'\n'))
	cat(paste('Top:',round(extent$Top),'\n'))
	cat(paste('N.Columns:',ncols,'\n'))
	cat(paste('N.Rows:',nrows,'\n'))
	cat('\n')
	
	return(extent)
	
}

##INITIAL MAX DENSITY MODULE
##This module is used to make sure that all starting points do not exceed a set maximum density 
max.density.source <- function(cellsize=100) #default = 100 meters
{	
	##Initialize an empty list in which to save all desired outputs
	lst <- list()
	
	cat('\n')
	print('Checking Density of Source Colonies...')
	cat('\n')
	
	##Using the simulation extent define grid breaks (based on cellsize resolution)	
	xbreaks <- seq(extent$Left,extent$Right,by=cellsize)
	ybreaks <- seq(extent$Bottom,extent$Top,by=cellsize)
	
	X <- starting_colonies$X
	Y <- starting_colonies$Y
	
	##Matrix with XY coord.
	Matr <- cbind(X,Y)	
	
	##2D spatial window (simulation extent)
	W <- c(extent$Left,extent$Right,extent$Bottom,extent$Top)
	
	##Coerce any reasonable kind of data to a point pattern
	##(an object of class "ppp") for use by the 'spatstat' package)
	pp <- as.ppp(Matr,W)	
	#plot(pp, pch='+')
	
	##Use quadratcount() function ('spatstat' package) to calculate
	##nbr of points fallins within each grid cell (quadrat)
	qX <- quadratcount(pp, xbreaks = xbreaks, ybreaks = ybreaks)
	#plot(qX, add=TRUE, col='red', cex=1.5, lty=2)
	
	##Invert the order of ybreaks to match the default order by
	##which the quadratcount() function goes
	ybreaks <- sort(ybreaks,decreasing=T)	
	
	##The [i,j] entry in the contingency table is the point count for the quadrat with coordinates 
	##(xbreaks[i],xbreaks[i+1]) by (ybreaks[i], ybreaks[i+1]). 
	
	##Check if any quadrat has a number of points exceeding the max density as by user input
	if(any(qX > maxdensity)){
		
		cat('=============================================================================\n')
		print('WARNING: The Starting Colonies Uploaded From File Have A Higher Density...')
		print('Maybe Not All Uploaded Points Are Colonies!! In That Case Pick Option 2')
		cat('=============================================================================\n')
		cat('\n')
		print('What Would You Like To Do?')
		choice <- menu(c('I Am Sure All Uploaded Points Are Different Colonies. Keep All Points And Find Automatically the Maximum Density',
						paste('Remove Points To Match The Density Of',maxdensity,'Colonies /ha'))) 
		
		if (choice == 1){
			
			NewMaxDen <- max(qX)
			cat('\n')
			cat(paste('The Max. Observed Density Is:',NewMaxDen,'Colonies /ha...\nWhat Would You Like To Do?')) 
			choice2 <- menu(c('Go Ahead With This Max. Density',paste('Remove Points To Match the Density Of',maxdensity,'Colonies /ha')))
			
			if (choice2 == 2) choice <- choice2 else maxdensity <- NewMaxDen
			
		}
		
		##If option 2 is picked
		if (choice == 2){
			
			##Initialize an empty vector in which we will save
			##the IDs of all points to be deleted
			Id_del <- c()
			
			columns <- length(xbreaks) - 1
			rows <- length(ybreaks) - 1
			
			##Check which grid cells have a density higher than the max density
			CellIdx <- which(qX > maxdensity)
			
			##Store the density of all grid cells w/ density higher
			##than max density
			if (length(CellIdx) == 1) freq <- qX[[CellIdx]] else freq <- qX[CellIdx] 
			
			##Loop through each one of the above cells			
			for (i in 1:length(CellIdx)){
		
				##Spot the positions of the cell within the simulation grid:
				
				##Column number (block nbr on X axis)
				blockx <- ceiling(CellIdx[[i]] / rows)
				
				##Row number (block nbr on Y axis)
				if (CellIdx[[i]] == 1) blocky <- 1 else blocky <- CellIdx[[i]] - (blockx-1)*rows
				
				##Identify points falling within the current cell				
				pnt_within <- starting_colonies$X > xbreaks[blockx] & starting_colonies$X <= xbreaks[blockx+1] & 
								starting_colonies$Y <= ybreaks[blocky] & starting_colonies$Y > ybreaks[blocky+1]
				
				##Extract those points from the dataset
				pnt_extract <- starting_colonies[pnt_within,]
				
				##Create a sequence of position IDs based on the nbr of extracted points
				pos_IDs <- seq(1:nrow(pnt_extract))
				
				##Random sample of position IDs to be deleted
				##The sample size is equal to the nbr of points exceeding the max density
				del_pos_IDs <- sample(pos_IDs, size = freq[i] - maxdensity)
				
				##Store the colony IDs of the ones randomly eliminated
				Id_del_temp <- pnt_extract$ID_col[del_pos_IDs]
				
				##Save and stack eliminated colony IDs
				Id_del <- c(Id_del,Id_del_temp)
			
			}
			
			##After the LOOP, remove colonies from the original dataset based on the eliminated IDs
			starting_colonies <- starting_colonies[-match(Id_del,starting_colonies$ID_col),] 
			
			##Overwrite the old colony IDs w/ new ones based on the nbr of colonies not eliminated
			starting_colonies$ID_col <- seq(1,nrow(starting_colonies))
			
		}
				
	}
	
	##Save the maximum density and the new starting colony dataset
	##To the empty list
	lst$maxdensity <- maxdensity
	lst$starting_colonies <- starting_colonies
	
	return(lst)	
}

##SAVE SHAPEFILE MODULE
##This module is used to save a point dataset to a shapefile (.shp)
shp.save <- function(folder, year)
{
	simulationShp <- colonies
	
	##Turn input point dataset into a SpatialPointsDataFrame
	coordinates(simulationShp) <- c('X','Y')
	
	##Define projected coord. according to the background/habitat layer
	proj4string(simulationShp) <- proj4string(habitat_block)
	
	##Write/Save the shapefile
	writeOGR(simulationShp, folder, layer = paste('Sim',year,sep=''), driver='ESRI Shapefile',overwrite_layer=TRUE)
}

##RASTER-AREA MODULE
##This module is used to convert current points to raster, calculate total occupied area, and save output raster
raster.save <- function(input, year, writeRas=TRUE, cellsize=100)  #by function default a raster is written in output
{
	
	##Using the simulation extent define grid breaks (based on cellsize resolution)	
	xbreaks <- seq(extent$Left,extent$Right, by=cellsize)
	ybreaks <- seq(extent$Bottom,extent$Top, by=cellsize)
	
	ncols <- length(xbreaks) - 1
	nrows <- length(ybreaks) - 1
	
	##Matrix w/ XY coordinates
	coords <- cbind(input$X,input$Y)
	
	##Store projection type from habitat layer
	CRS <- proj4string(habitat_block)
	
	if (nrow(input) > 0){
	
		##Create a raster using the raster() function of the 'raster' package
		RastLayer <- raster(xmn=extent$Left, xmx=extent$Right, ymx=extent$Top, ymn=extent$Bottom, ncol=ncols, nrow=nrows, crs=CRS)
		
		##Convert points-to-raster
		##A Point-To-Raster operation is done to estimate the overall area
		PntToRaster <- rasterize(coords, RastLayer) #by default (1-presence NA-absence) #
		
		if (writeRas == TRUE) writeRaster(PntToRaster,filename=paste(workdir_Raster,'Sim',year,'.img',sep=''),
										format='HFA', datatype='LOG1S', overwrite=TRUE)
		
		Tot_area_ha <- sum(!is.na(PntToRaster@data@values))
		Tot_area_km2 <- Tot_area_ha/100
	
	}else{
		
		##Create a raster using the raster() function of the 'raster' package
		RastLayer <- raster(matrix(NA),xmn=extent$Left, xmx=extent$Right, ymx=extent$Top, ymn=extent$Bottom, ncol=ncols, nrow=nrows, crs=CRS)
		
		##Convert points-to-raster
		##A Point-To-Raster operation is done to estimate the overall area
		if (writeRas == TRUE) writeRaster(RastLayer,filename=paste(workdir_Raster,'Sim',year,'.img',sep=''),
										format='HFA', datatype='LOG1S', overwrite=TRUE)
		
		Tot_area_ha <- 0
		Tot_area_km2 <- 0
	}
	
	return(Tot_area_km2)
}	

##AGE & TIME INCREASE MODULE
##This module is used to increase colony age and timestep by 1
age.increase <- function()
{
	
	##Isolate colonies from the previous time step
	past_colonies <- colonies[colonies$Timecount == year - 1,]
	
	##Overwrite the current year to the old one
	past_colonies$Timecount <- year 
	
	##Add +1 to the age of each colony
	past_colonies$Age <- past_colonies$Age + 1
	
	##After increasing age and timestep, store it as current colonies
	present_colonies <- past_colonies
	
	#Stack present colonies to the previous ones
	colonies <- rbind(colonies,present_colonies)
		
	return(colonies)

}

##SWARMER GENERATION MODULE
##This module is used to generate swarmers from each colony depending on its age
offspring.generate <- function(survival_prob, male_prob, scenario, dist.mean)
{
	
	if (any(colonies$Age >= ColAge_swarmers)){
	
		if (scenario == 'optimistic'){
				
			alates <- ifelse(colonies$Age < ColAge_swarmers, 0, 100000)
			alates <- ifelse(colonies$Age >= ColAge_swarmers & colonies$Age < 10, 1000, alates)
			alates <- ifelse(colonies$Age >= ColAge_swarmers & colonies$Age < 15, 10000, alates)
			
		}else if (scenario == 'pessimistic'){
			
			#TODO:: CHANGE IF NEEDED
			alates <- ifelse(colonies$Age < ColAge_swarmers, 0, 100000)
			alates <- ifelse(colonies$Age >= ColAge_swarmers & colonies$Age < ColAge_swarmers * 2, 10000, alates)
			alates <- ifelse(colonies$Age >= ColAge_swarmers & colonies$Age < ColAge_swarmers * 3, 50000, alates)
		}
		
		indiv <- alates * survival_prob

		X_temp <- rep(colonies$X, indiv)
		Y_temp <- rep(colonies$Y, indiv)

		#Initialize two coordinate variables based on the established (or existing) colonies... now as vectors of the entire data frame size
		distance <- rexp(sum(indiv), rate = 1/dist.mean)
		theta <- runif(sum(indiv), 0, 2 * pi)
		C <- cos(theta)
		S <- sin(theta)
		
		#XY coords (meters) using polar coordinate transformations
		X <- X_temp + S * distance
		Y <- Y_temp + C * distance
		
		pop <- data.frame(X,Y)  
		pop$Sex <- rbinom(sum(indiv),1,male_prob)
		pop$ID <- 1:nrow(pop)
		
		pop$X <- round(pop$X,2)
		pop$Y <- round(pop$Y,2)
	}
	
	#Define outline of desired mask
	msk <- data.frame(x = c(extent$Left,extent$Right,extent$Right,extent$Left), 
					  y = c(extent$Bottom,extent$Bottom,extent$Top,extent$Top))
			
	idx <- point.in.polygon(point.x = pop$X, point.y = pop$Y, pol.x = msk[,1], pol.y = msk[,2])		
	
	if(any(idx==0) & length(indiv) > 0) {
		
		cat('\n')
		cat('STOP: ONE OR MORE INDIVIDUALS LAY OUTSIDE THE SIMULATION EXTENT!!\n')
		cat(paste('The Simulation Has Been Terminated Before Completing Year: ',year,'\n',sep=''))
		stop()
		cat('\n')
		
	}
	
	return(pop)
}

##NEW COLONIES CREATION MODULE
##This module is used to generate new colonies within each grid cell, controlling for the maximum allowed density
new.colonies.create <- function(cellsize=100)
{
	#Initialize a data.frame in which I will store all my final new colonies
	#ID1 and ID2 correspond to individual IDs of a male-female pair, and X and Y are the coord. of the new colony 
	#formed by ID1 and ID2. X and Y are determined by summing X coord. of ID1 and ID2 and dividing by 2...same for Y coord.
	All_new_colonies <- data.frame(ID1=0,ID2=0,X=0,Y=0)
	
	#Using the simulation extend (defined in the main script) 
	#define grid intervals, depending on cellsize (grid resolution)	
	xbreaks <- seq(extent$Left,extent$Right,by=cellsize)
	ybreaks <- seq(extent$Bottom,extent$Top,by=cellsize)
		
	#Matrix w/ coord. of each individual
	Matr <- cbind(pop$X,pop$Y)
	
	#Matrix w/ coord. of each existing colony at the current time step
	MatrCol <- cbind(colonies$X,colonies$Y)	
	
	#2D spatial window (simulation extent)
	W <- c(extent$Left,extent$Right,extent$Bottom,extent$Top)
	
	#Use as.ppp to coerce any reasonable kind of data to a point pattern (an object of class 'ppp') 
	pp <- as.ppp(Matr,W) #in pp we have individuals
	ppCol <- as.ppp(MatrCol,W) #in ppCol we have all existing colonies
		
	#plot(pp, pch='+')
	#plot(ppCol, pch='+')
	
	#Use the function quadratcount (spatstat package) to compute number of points within each quadrat:
	
	#qX stores the # of individuals (frequencies) falling within each cell/quadrat
	qX <- quadratcount(pp, xbreaks = xbreaks, ybreaks = ybreaks)
	#plot(qX, add=TRUE, col='red', cex=1.5, lty=2)
	
	#qXCol stores the # of existing colonies (frequencies) falling within each square/quadrat
	qXCol <- quadratcount(ppCol, xbreaks = xbreaks, ybreaks = ybreaks)
	#plot(qXCol, add=TRUE, col='red', cex=1.5, lty=2)
	
	#Sort ybreaks decreasing to mimic the qX matrix attributes (i.e. first cell is the top-left in a 2D space)
	ybreaks <- sort(ybreaks,decreasing=T)
	
	#Store #of rows and columns of the spatial window
	columns <- length(xbreaks) - 1
	rows <- length(ybreaks) - 1
	
	#CELL.ID_MaxDen spots which cells have already the max possible density (i.e. no more new colonies will establish)
	CELL.ID_MaxDen <- which(qXCol == maxdensity)
	
	#CELL.ID_2Indiv spots which cells have at least 2 individuals (1 indiv. alone cannot create a new colony)
	#That way I ignore all other cells...
	CELL.ID_GT2Indiv <- which(qX >= 2)
	
	#Ignore (remove) cells already at maximum carrying capacity from the list of cells previously selected
	if ( any(CELL.ID_MaxDen) ) {
		commonCells <- match(CELL.ID_MaxDen,CELL.ID_GT2Indiv)
		#Make sure to remove any NA values from the list and keep only the common cells to be removed
		remove_cells <- commonCells[!is.na(commonCells)]
		CELL.ID_GT2Indiv <- CELL.ID_GT2Indiv[-remove_cells]	
	}	

	##MAIN LOOP: loop through each cell --containing at least 2 individuals---
	for (cellIdx in 1:length(CELL.ID_GT2Indiv)){  
		
		ColDensity <- qXCol[CELL.ID_GT2Indiv[cellIdx]]
		
		##Now we want to spot the position (row and column) of the current cell within the big matrix:
		
		#blockx indentifies the column 
		blockx <- ceiling(CELL.ID_GT2Indiv[cellIdx] / rows)
		
		#blocky indentifies the row
		blocky <- CELL.ID_GT2Indiv[cellIdx] - (blockx-1) * rows
		
		#Define the spatial window (extent) corresponding to the current cell		
		Left <- xbreaks[blockx] 
		Right <- xbreaks[blockx+1] 
		Bottom <- ybreaks[blocky+1]
		Top <- ybreaks[blocky]
		
		#Define outline of desired mask
		msk <- data.frame(x = c(Left,Right,Right,Left,Left), y = c(Bottom,Bottom,Top,Top,Bottom))
			
		idx <- point.in.polygon(point.x = pop$X, point.y = pop$Y, pol.x = msk[,1], pol.y = msk[,2])		
	
		#Store the selected points in a variable called select
		select <- which(idx > 0)
		indiv_selection <- pop[select,]
		
		#Add a variable called Flag and set if to FALSE (meaning that none of the selected individuals has paired
		#w/ a partner yet). The flag will be changed to TRUE as soon as they will find a partner...
		indiv_selection$Flag <- FALSE 
		
		#Initialize a data.frame where I will store all new colonies created inside the current cell!
		new_colonies <- data.frame(ID1=0,ID2=0,X=0,Y=0)
			
		##SECOND LOOP: loop through each individual selected within the current cell
		for (indiv in 1:nrow(indiv_selection)){
			
			#Go ahead only if the current individual has not found a partner yet (flag = FALSE), otherwise skip this individual 			
			if (indiv_selection$Flag[indiv] == FALSE){
				
				#Use disc function from the 'spatstat' package
				Buffer <- disc(radius, centre=c(indiv_selection$X[indiv], indiv_selection$Y[indiv])) 
				
				#Find points falling within buffer(other than the current individual--center of the buffer)
				BufferPoly = data.frame(x = Buffer$bdry[[1]]$x, y = Buffer$bdry[[1]]$y)
				idx <- point.in.polygon(point.x = indiv_selection[-indiv,]$X, point.y = indiv_selection[-indiv,]$Y, pol.x = BufferPoly[,1], pol.y = BufferPoly[,2])		
				
				#Go ahead only if there is AT LEAST one individual within the buffer
				if ( !any(idx > 0) ) next
				
				#Check wheater any individual within the buffer has already found a partner or not
				pos <- match(indiv_selection[-indiv,]$ID[which(idx > 0)],indiv_selection$ID)
				flag_check <- indiv_selection$Flag[pos] == FALSE
				
				#Go ahead only if, out of the potential partners, there is AT LEAST an individual that has not paired already 
				if ( sum(flag_check) == 0 ) next
				
				#Check if the individuals within the buffer can be a potential partner (must be different sex)
				Partner_heterosex <- indiv_selection$Sex[indiv] != indiv_selection$Sex[pos]
				
				#Go ahead only if, out of the potential partners, there is AT LEAST a heterosexual
				if ( sum(Partner_heterosex) == 0 ) next
					
				#Set the flag to TRUE to indicate that the current individual has found a partner
				indiv_selection$Flag[indiv] <- TRUE 
					
				pnt.inBuffID <- indiv_selection[-indiv,]$ID[which(idx > 0)]
						
				#Select IDs of heterosex individuals within the buffer
				ID_heterosex <- pnt.inBuffID[Partner_heterosex]							
						
				#Calculate the euclidean distance(s) btw the current individual and the heterosexual individuals
				eucdist <- sqrt( (indiv_selection$X[indiv] - indiv_selection$X[match(ID_heterosex,indiv_selection$ID)])^2 + 
								(indiv_selection$Y[indiv] - indiv_selection$Y[match(ID_heterosex,indiv_selection$ID)])^2 )
							
				#Select the CLOSEST heterosex individual
				ID_Selected <- ID_heterosex[which(eucdist == min(eucdist))]
						
				#Set the Flag of the selected partner to TRUE						
				indiv_selection$Flag[which(indiv_selection$ID == ID_Selected)] <- TRUE
						
				#Save a new record with IDs of partners and X-Y coord. of the new colony 
				#(as an assumption create the new point in the middle btw the two partners)
				new_colonies_temp <- data.frame(ID1=indiv_selection$ID[indiv], ID2=ID_Selected, 
									X=( indiv_selection$X[indiv] + indiv_selection$X[which(indiv_selection$ID == ID_Selected)] ) / 2, 
									Y=( indiv_selection$Y[indiv] + indiv_selection$Y[which(indiv_selection$ID == ID_Selected)] ) / 2)
							
						
				#Stack the new colony record to the main new colonies data.frame for the current cell						
				new_colonies <- rbind(new_colonies,new_colonies_temp) ##Can be improved without rbind
						
				#Only if the user defines a maximum density EXIT (break) the LOOP as soon as the max density is reached!!		
				#If the user does not define a max density the program will run much slower!
				if ( nrow(new_colonies) - 1 == maxdensity - ColDensity ) break
		
		    }
		}
		
		#Remove the first line of the data.frame that was set to all 0's.
		new_colonies <- new_colonies[-1,]
		
		#If new colonies data.frame has at least one record then store it into the main data.frame with all other
		#new colonies
		if (nrow(new_colonies) > 0) All_new_colonies <- rbind(All_new_colonies,new_colonies)			
		rm(new_colonies)
	}
	
	All_new_colonies <- All_new_colonies[-1,]
	
	return(All_new_colonies)
	
}

##FUNCTION TO STACK OLD & NEW COLONIES IN CURRENT TIMESTEP
##This module is used to stack the newly generated colonies to all existing ones
new.colonies.stack <- function()
{
	#If there is only ONE source of invasion in the model, compute the euclidean distance 
	#between each new created colony and the original source.
	if (nrow(starting_colonies) == 1) {
	new_colonies$Dist_source <- round(sqrt( (new_colonies$X - starting_colonies$X)^2 + (new_colonies$Y - starting_colonies$Y)^2 ),2)
	}
	#Remove the first two columns, since we are not interested in keeping IDs of pairs
	new_colonies <- new_colonies[,-c(1,2)]
	#Here I create a new time variable and add it as a reference to the new_colonies dataframe
	new_colonies$Timecount <- rep(year,nrow(new_colonies))
	#Here I create a colony ID variable and add it as a reference to the new_colonies dataframe
	new_colonies$ID_col <- seq(1,nrow(new_colonies))
	#Age field
	new_colonies$Age <- rep(0,nrow(new_colonies))
	colonies <- rbind(colonies,new_colonies)
	colonies$ID_col <- seq(1,nrow(colonies))
	
	return(colonies)

}

##CONVEX HULL MODULE
##This module is used to spot and keep only colonies laying on the fring of the invasion
convex.hull <- function(input)
{	

	Matr <- cbind(input$X,input$Y)	
	hpts <- chull(Matr)
	
	colonies <- input[hpts,] 
	colonies$ID_col <- seq(1,nrow(colonies))
	
	return(colonies)	

}

##SUMMARY STATISTICS MODULE
##This module is used to calculate required stastistics on simulation results
summary.stats <- function()
{
	#Now let's have a look at the final stats
	simYears <- seq(start_time,end_time)
	
	#Store the information about the area covered after each time step 
	Area_dataset <- data.frame(Year = simYears, Area = area)
	colnames(Area_dataset) <- c('Years','Area (Km^2)')
	
	#Print on console the dataset
	cat('\n===========================================================\n')
	cat('The Approx. Area (Km^2) Covered After Each Time Step Is:\n')
	print(Area_dataset)
	cat('\n=============================================================\n')
}








