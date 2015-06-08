#--------------------------------------------------------------------------------
# Name:         myfunctions_CA.r
# Purpose:      Modules (functions) called by the main script
# Author:       Francesco Tonini
# Email:        ftonini84@gmail.com
# Created:      09/20/2013
# Copyright:    (c) 2013 by Francesco Tonini
# License:    	GNU General Public License (GPL)
# Software:     Tested successfully using R version 3.0.2 (http://www.r-project.org/)
#-----------------------------------------------------------------------------------------


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

kernel2D <- function(x, d){
  
  rs <- res(x)
  if(d <= res(x)[1] | d <= res(x)[2]) stop('d must be >= grain of study area')
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  probKernel <- matrix(NA, ncol=nx, nrow=ny)
  probKernel[ceiling(ny/2), ceiling(nx/2)] <- 1
  
  return(probKernel)
  
}

generate.Kernel <- function(x, d, FUN, alpha){
  
  rs <- res(x)
  
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  
  m <- matrix(ncol=nx, nrow=ny)
  
  xr <- (nx * rs[1]) / 2
  yr <- (ny * rs[2]) / 2 
  r <- raster(m, xmn=-xr[1], xmx=xr[1], ymn=-yr[1], ymx=yr[1], crs='+proj=utm +zone=1')
  r[ceiling(ny/2), ceiling(nx/2)] <- 1
   
  dist <- as.matrix(distance(r)) 
  
  if (FUN == 'exp'){
    c <- 1
    m <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(dist/alpha)^c  ) 
  }else if (FUN == 'gauss'){
    c <- 2
    m <- ( c / 2 * alpha * gamma(1/c) ) * exp( -abs(dist/alpha)^c  ) 
  }
  
  #sum of weights should add up to 1
  m/sum(m)
  
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

