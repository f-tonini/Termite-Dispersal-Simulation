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

alates.gen <- function(x, scenario='baseline'){
  
  if (scenario == 'baseline'){
    
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
