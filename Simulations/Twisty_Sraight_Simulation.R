#------------------------------------------------------------------------------------------------#
#----Simulation study comparing twisty sampling design to alternative straight line sampling.----#
#----Data is simmulated for twisty sampling design of transects using actual transects.----------#
#----Transects were imported as shapefiles. Alternative sampling design is simulated.------------#
#----Script lasted edited by Matthew Farr (1/26/17)----------------------------------------------#
#------------------------------------------------------------------------------------------------#

#----------#
#-Set Seed-#
#----------#

set.seed(1985)

#-----------------------#
#-Set Working Directory-#
#-----------------------#

setwd("C:/Users/farrm/Documents/GitHub/CDSM/Simulations")

#---------------#
#-Load Libaries-#
#---------------#

library(rgdal)
library(sp)
library(dplyr)
library(tidyr)
library(jagsUI)

#--------------------------#
#-Create Sampling Boundary-#
#--------------------------#

#Easting
xlim <- c(715304, 752393)

#Northing
ylim <- c(9831970, 9857296)

#-------------------------------#
#-Simulate Set Parameter Values-#
#-------------------------------#

#Number of groups
N <- 1000

#Half normal detection parameter
sigma <- 300

#Mid point of each distance class
midpt <- seq(12.5, 650, 25)

#Index for distance class
nG <- length(midpt)

#Width of distance class
v <- 25

#Transect half width
B <- 650

#----------------------------#
#-Import Transect Shapefiles-#
#----------------------------#

#Directory for transects by site shapefile
d.dir <- "C:/Users/farrm/OneDrive/Hyena Project/datasets/NewDatasets/Rscripts/ExcelFiles/Site"

#Transects by site
Site1 <- readOGR(dsn = d.dir, layer = "Site1")
Site2 <- readOGR(dsn = d.dir, layer = "Site2")
Site3 <- readOGR(dsn = d.dir, layer = "Site3")
Site4 <- readOGR(dsn = d.dir, layer = "Site4")
Site5 <- readOGR(dsn = d.dir, layer = "Site5")
Site6 <- readOGR(dsn = d.dir, layer = "Site6")
Site7 <- readOGR(dsn = d.dir, layer = "Site7")
Site8 <- readOGR(dsn = d.dir, layer = "Site8")
Site9 <- readOGR(dsn = d.dir, layer = "Site9")
Site10 <- readOGR(dsn = d.dir, layer = "Site10")
Site11 <- readOGR(dsn = d.dir, layer = "Site11")
Site12 <- readOGR(dsn = d.dir, layer = "Site12")
Site13 <- readOGR(dsn = d.dir, layer = "Site13")
Site14 <- readOGR(dsn = d.dir, layer = "Site14")
Site15 <- readOGR(dsn = d.dir, layer = "Site15")
Site16 <- readOGR(dsn = d.dir, layer = "Site16")
Site17 <- readOGR(dsn = d.dir, layer = "Site17")

#-----------------------------------#
#-Sample Coordinates from Transects-#
#-----------------------------------#

s1p <- spsample(Site1, 30, type = "regular")
s2p <- spsample(Site2, 30, type = "regular")
s3p <- spsample(Site3, 30, type = "regular")
s4p <- spsample(Site4, 30, type = "regular")
s5p <- spsample(Site5, 30, type = "regular")
s6p <- spsample(Site6, 30, type = "regular")
s7p <- spsample(Site7, 30, type = "regular")
s8p <- spsample(Site8, 30, type = "regular")
s9p <- spsample(Site9, 30, type = "regular")
s10p <- spsample(Site10, 30, type = "regular")
s11p <- spsample(Site11, 30, type = "regular")
s12p <- spsample(Site12, 30, type = "regular")
s13p <- spsample(Site13, 30, type = "regular")
s14p <- spsample(Site14, 30, type = "regular")
s15p <- spsample(Site15, 30, type = "regular")
s16p <- spsample(Site16, 30, type = "regular")
s17p <- spsample(Site17, 30, type = "regular")

#----------------------#
#-Alternative Sampling-#
#----------------------#

#sampling area
mdE <- 733848.5
mdN <- 9844633

Et <- mdE - (13000/2)
Nt <- mdN + (12650/2)

Ep <- seq((Et + 650), (Et + (13000 - 650)), 1300)
Np <- seq(Nt, (Nt - 12650), -253)

#---------------#
#-Model Offsets-#
#---------------#

#Search area of each site
A.site <- as.vector(c(11.6542, 11.9619, 12.4702, 12.5182, 10.7843, 10.2384, 10.7495, 
                      12.0545, 9.0114, 11.2589, 10.4075, 9.7834, 11.8226, 10.5295,
                      11.5376, 14.8511, 14.0352))

#----------------------#
#-BUGS Model File Path-#
#----------------------#

modFile <- "C:/Users/farrm/Documents/GitHub/CDSM/Simulations/CDSM_model.txt"

#------------------#
#-Begin Iterations-#
#------------------#

#Number of iterations
niter <- 120	
subniter <- 25

#Starting iteration
iter <- 1

#Output vectors
bias <- array(NA, dim = c(5, 6, subniter, niter), 
              dimnames = list(c("Groups In", "Abundance In", "Groups", "Abundance", "Sigma"), 
                              c("Twst True", "Twst Est", "Twst Bias", "Str True", "Str Est", "Str Bias"),
                              NULL, NULL))
overlapcor <- array(NA, dim = c(niter, subniter))

system.time(while(iter <= niter){
  
  print (c(iter, "iter"))
  
  subiter <- 1       

  #-----------------------------------#
  #-Simulate Dynamic Parameter Values-#
  #-----------------------------------#
  
  #Simulate coordinates of groups
  u1 <- runif(N, xlim[1], xlim[2]) 
  u2 <- runif(N, ylim[1], ylim[2])
  
  #Group size
  lambda.group <- 2
  cs <- rpois(N, lambda.group) + 1
  
  #Abundance
  Ntotal <- sum(cs)
  
  #---------------------#
  #--BEGIN SUBSAMPLING--#
  #---------------------#
  
  while(subiter <= subniter){
    
    print (c(subiter, "subiter"))
    
    #----------------------------------------------#
    #-Combine Site Coordinates for Twisty Sampling-#
    #----------------------------------------------#
    
    #Easting
    X <- c(s1p@coords[,1], s2p@coords[,1], s3p@coords[,1], s4p@coords[,1], 
           s5p@coords[,1], s6p@coords[,1], s7p@coords[,1], s8p@coords[,1],
           s9p@coords[,1], s10p@coords[,1], s11p@coords[,1], s12p@coords[,1],
           s13p@coords[,1], s14p@coords[,1], s15p@coords[,1], s16p@coords[,1],
           s17p@coords[,1])
    
    #Northing
    Y <- c(s1p@coords[,2], s2p@coords[,2], s3p@coords[,2], s4p@coords[,2], 
           s5p@coords[,2], s6p@coords[,2], s7p@coords[,2], s8p@coords[,2],
           s9p@coords[,2], s10p@coords[,2], s11p@coords[,2], s12p@coords[,2],
           s13p@coords[,2], s14p@coords[,2], s15p@coords[,2], s16p@coords[,2],
           s17p@coords[,2])
    
    #-------------------#
    #-Initialize Values-#
    #-------------------#
    
    #Index for sites
    nsites <- 17
    
    #Index for transect points
    J <- length(X)
    
    #ID for sites
    si <- seq(0, J, (J/nsites))
    
    #ID for distance class
    di <- seq(0,650,25)
    
    #Distance class
    dclass <- rep(NA, N)
    
    #Minimum distance value
    dst <- rep(NA, N)
    
    #ID for nearest site
    q <- rep(NA, N)
    
    #Site
    site <- rep(NA, N)
    
    #Distance value to each transect point
    d <- array(NA, dim = c(N, J))
    
    #ID for groups less than 650 meters
    y <- rep(NA, N)
    
    #Index recorder
    index <- rep(NA, N)
    
    #---------------#
    #-Simulate Data-#
    #---------------#
    
    for(i in 1:N){
      for(j in 1:J){
        
        #Distance from each group to each point on the transect
        d[i,j] <- sqrt((u1[i] - X[j])^2 + (u2[i] - Y[j])^2)
      }
      
      #Distance to nearest point on the transect
      dst[i] <- min(d[i,])
      
      #Index of which point in 1:J is the nearest
      q[i] <- which.min(d[i,])
      
      for(j in 1:nsites){
        
        #Determine the site for each group
        if(si[j] < q[i] && q[i] <= si[j+1])
          site[i] <- j
      }
      
      #Index of which observation are within 650 meters of transect
      if(dst[i] < 650)
        y[i] <- 1
      index[i] <- i
    }
    
    #--------------#
    #-Harvest Data-#
    #--------------#
    
    #Dataframe that includes information on all groups
    Dtot <- cbind(y, index, u1, u2, site, cs)
    
    #Dataframe containing only groups within 650 meters to transect
    Din <- Dtot[complete.cases(Dtot),]
    
    #Number of groups within 650 meters
    Nin <- length(Din[,1])
    
    #Abundance within 650 meters
    Nintotal <- sum(Din[,6])
    
    #-----------------#
    #-Initialize Data-#
    #-----------------#
    
    #Remove groups not within 650 meters
    index <- index[y==1]
    index <- index[!is.na(index)]
    
    #Detection Probability
    p <- NULL
    
    #Number of captured ("detected") groups
    ncap <- rep(NA, Nin)
    
    #Distance Class
    dclass <- rep(NA, Nin)
    
    #---------------#
    #-Simulate Data-#
    #---------------#
    
    for(i in 1:Nin){
      
      #Detection probability using half-normal distance function
      p[i] <- exp(-dst[index[i]] * dst[index[i]] / (2 * sigma * sigma))
      
      #Simulate number of groups detected
      ncap[i] <- rbinom(1, 1, p[i])
      
      for(k in 1:nG){
        
        #Determine distance class for each group
        if(di[k] < dst[index[i]] && dst[index[i]] <= di[k+1])
          dclass[i] <- k
      }
    }
    
    #--------------#
    #-Harvest Data-#
    #--------------#
    
    #Add distance class, detection probability, and detection index to dataframe
    Din <- cbind(Din[,2:6],dclass, p, ncap)
    
    #Undetected groups as NAs
    for(i in 1:Nin){
      if(Din[i,8] == 0)
        Din[i,8] <- NA
    }
    
    #Dataframe of detected inidividuals
    Dcap <- Din[complete.cases(Din),]
    
    #Create observed number of groups per site
    y.new <- table(Dcap[,4])
    y.new <- as.data.frame(y.new)
    colnames(y.new) <- c("site", "freq")
    y.new$site <- as.integer(y.new$site)
    y.new <- tbl_df(y.new)
    
    #Add in sites with no detections
    miss <- y.new %>% expand(site = 1:nsites)
    miss$freq <- rep(0, length(miss))
    
    #Add missing sites into observed groups per site
    yobs <- full_join(y.new, miss, by = "site")
    yobs <- yobs %>% arrange(site)
    yobs <- as.numeric(yobs$freq.x)
    yobs[is.na(yobs)] <- 0
    
    #Site index for observed number of groups
    site <- Dcap[,4]
    
    #Distance class index for observed number of groups
    dclass <- Dcap[,6]
    
    #Number of observations
    nobs <- sum(yobs)
    
    #Group size
    gs <- Dcap[,5]
    
    #-----------------------#
    #-Generate Overlap Data-#
    #-----------------------#
    
    #Initialize overlap array
    overlap <- array(0, dim = c(Nin,J))
    
    for(i in 1:Nin){
      for(j in 1:J){
        
        #Harvest potential overlap data for groups within 650 meters
        if(ncap[i] == 1 && d[index[i],j] < 650)
          overlap[i,j] <- d[index[i],j]
      }
    }
    
    #Matirx of group ID (col 1) and transect point ID (col 2) within 650 meters
    OVA <- which(!overlap == 0, arr.ind = TRUE)
    OVA <- OVA[order(OVA[,1]),]
    
    #Initialize number of overlaps per site
    OVAsite <- NULL
    
    for(i in 1:length(OVA[,1])){
      for(j in 1:nsites){
        
        #Corresponding site ID for each group ID
        if(si[j] < OVA[i,2] && OVA[i,2] <= si[j+1])
          OVAsite[i] <- j
      }
    }
    
    #Combine site ID with dataframe
    OVA <- data.frame(OVA[,1], OVAsite)
    
    #Deletes duplicate values
    OVA <- unique(OVA)
    
    #Removes groups not seen in 2 sites
    OVA <- subset(OVA, duplicated(OVA[,1]) | duplicated(OVA[,1], fromLast = TRUE))
    
    #Determines number of overlap per site
    OVA <- group_by(OVA, OVAsite)%>%
      summarize(n_distinct(OVA...1.))
    colnames(OVA) <- c("site", "overlaps")
    
    #Adds sites with no overlap
    miss <- OVA %>% expand(site = 1:nsites)
    miss$"overlaps" <- rep(0, length(miss))
    OVA <- full_join(OVA, miss, by = "site")
    OVA <- OVA %>% arrange(site)
    OVA[is.na(OVA)] <- 0
    
    #-------------------#
    #-Compile BUGS data-#
    #-------------------#
    
    #Input data
    realD <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                      nobs = nobs, dclass = dclass, nsites = nsites, gs = gs, offset = A.site)
    
    #Initial values
    N.in <- yobs + 1
    
    inits <- function(){list(N = N.in, sigma = runif(17, 50, 350))} 
    
    #Parameters to monitor
    params<-c('sigma', 'Nin', 'Nintotal', 'Nreal', 'Nrealtotal')
    
    #MCMC settings
    
    nc <- 3
    ni <- 12000
    nb <- 2000
    nt <- 4
    
    #----------------#
    #-Run BUGS Model-#
    #----------------#
    
    realM <- jags(data = realD, inits = inits, parameters.to.save = params, model.file = modFile, 
                  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    #----------------------------------------#
    #-Save and Remove Data for Next Sampling-#
    #----------------------------------------#
    
    realVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal)
    
    rm(X, Y, nsites, J, si, di, dclass, dst, q, 
       site, d, y, index, Dtot, Din, Nin, Nintotal,
       p, ncap, Dcap, y.new, miss, yobs, nobs, gs,
       N.in, inits)
    
    #----------------------#
    #-Alternative Sampling-#
    #----------------------#
    
    X <- rep(Ep, rep(length(Np), length(Ep)))
    Y <- rep(Np, length(Ep))
    
    #-------------------#
    #-Initialize Values-#
    #-------------------#
    
    #Index for sites
    nsites <- 10
    
    #Index for transect points
    J <- length(X)
    
    #ID for sites
    si <- seq(0, J, (J/nsites))
    
    #ID for distance class
    di <- seq(0,650,25)
    
    #Minimum distance value
    dst <- rep(NA, N)
    
    #ID for nearest site
    q <- rep(NA, N)
    
    #Site
    site <- rep(NA, N)
    
    #Distance value to each transect point
    d <- array(NA, dim = c(N, J))
    
    #ID for groups less than 650 meters
    y <- rep(NA, N)
    
    #Index recorder
    index <- rep(NA, N)
    
    #---------------#
    #-Simulate Data-#
    #---------------#
    
    for(i in 1:N){
      for(j in 1:J){
        
        #Distance from each group to each point on the transect
        d[i,j] <- sqrt((u1[i] - X[j])^2 + (u2[i] - Y[j])^2)
      }
      
      #Distance to nearest point on the transect
      dst[i] <- min(d[i,])
      
      #Index of which point in 1:J is the nearest
      q[i] <- which.min(d[i,])
      
      for(j in 1:nsites){
        
        #Determine the site for each group
        if(si[j] < q[i] && q[i] <= si[j+1])
          site[i] <- j
      }
      
      #Index of which observation are within 650 meters of transect
      if(dst[i] < 650)
        y[i] <- 1
      index[i] <- i
    }
    
    #--------------#
    #-Harvest Data-#
    #--------------#
    
    #Dataframe that includes information on all groups
    Dtot <- cbind(y, index, u1, u2, site, cs)
    
    #Dataframe containing only groups within 650 meters to transect
    Din <- Dtot[complete.cases(Dtot),]
    
    #Number of groups within 650 meters
    Nin <- length(Din[,1])
    
    #Abundance within 650 meters
    Nintotal <- sum(Din[,6])
    
    #-----------------#
    #-Initialize Data-#
    #-----------------#
    
    #Remove groups not within 650 meters
    index <- index[y==1]
    index <- index[!is.na(index)]
    
    #Detection Probability
    p <- NULL
    
    #Number of captured ("detected") groups
    ncap <- rep(NA, Nin)
    
    #Distance Class
    dclass <- rep(NA, Nin)
    
    #---------------#
    #-Simulate Data-#
    #---------------#
    
    for(i in 1:Nin){
      
      #Detection probability using half-normal distance function
      p[i] <- exp(-dst[index[i]] * dst[index[i]] / (2 * sigma * sigma))
      
      #Simulate number of groups detected
      ncap[i] <- rbinom(1, 1, p[i])
      
      for(k in 1:nG){
        
        #Determine distance class for each group
        if(di[k] < dst[index[i]] && dst[index[i]] <= di[k+1])
          dclass[i] <- k
      }
    }
    
    #--------------#
    #-Harvest Data-#
    #--------------#
    
    #Add distance class, detection probability, and detection index to dataframe
    Din <- cbind(Din[,2:6],dclass, p, ncap)
    
    #Undetected groups as NAs
    for(i in 1:Nin){
      if(Din[i,8] == 0)
        Din[i,8] <- NA
    }
    
    #Dataframe of detected inidividuals
    Dcap <- Din[complete.cases(Din),]
    
    #Create observed number of groups per site
    y.new <- table(Dcap[,4])
    y.new <- as.data.frame(y.new)
    colnames(y.new) <- c("site", "freq")
    y.new$site <- as.integer(y.new$site)
    y.new <- tbl_df(y.new)
    
    #Add in sites with no detections
    miss <- y.new %>% expand(site = 1:nsites)
    miss$freq <- rep(0, length(miss))
    
    #Add missing sites into observed groups per site
    yobs <- full_join(y.new, miss, by = "site")
    yobs <- yobs %>% arrange(site)
    yobs <- as.numeric(yobs$freq.x)
    yobs[is.na(yobs)] <- 0
    
    #Site index for observed number of groups
    site <- Dcap[,4]
    
    #Distance class index for observed number of groups
    dclass <- Dcap[,6]
    
    #Number of observations
    nobs <- sum(yobs)
    
    #Group size
    gs <- Dcap[,5]
    
    #-------------------#
    #-Compile BUGS data-#
    #-------------------#
    
    #Input data
    altD <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                     nobs = nobs, dclass = dclass, nsites = nsites, gs = gs, offset = rep(1, nsites))
    
    #Initial values
    N.in <- yobs + 1
    
    inits <- function(){list(N = N.in, sigma = runif(10, 50, 350))}
    
    #----------------#
    #-Run BUGS Model-#
    #----------------#
    
    altM <- jags(data = altD, inits = inits, parameters.to.save = params, model.file = modFile, 
                 n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    #----------------------------------------#
    #-Save and Remove Data for Next Sampling-#
    #----------------------------------------#
    
    altVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal)
    
    rm(X, Y, nsites, J, si, di, dclass, dst, q, 
       site, d, y, index, Dtot, Din, Nin, Nintotal,
       p, ncap, Dcap, y.new, miss, yobs, nobs, gs,
       N.in, inits)
    
    #----------------#
    #-Bias Estimates-#
    #----------------#
      
    bias[1, 1, subiter, iter] <- realVals[[3]]
    bias[1, 2, subiter, iter] <- realM$mean$Nin
    bias[1, 3, subiter, iter] <- (abs(mean((realM$sims.list$Nin - realVals[[3]])/realVals[[3]])) * 100)
    bias[1, 4, subiter, iter] <- altVals[[3]]
    bias[1, 5, subiter, iter] <- altM$mean$Nin
    bias[1, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nin - altVals[[3]])/altVals[[3]])) * 100)
      
    bias[2, 1, subiter, iter] <- realVals[[4]]
    bias[2, 2, subiter, iter] <- realM$mean$Nintotal
    bias[2, 3, subiter, iter] <- (abs(mean((realM$sims.list$Nintotal - realVals[[4]])/realVals[[4]])) * 100)
    bias[2, 4, subiter, iter] <- altVals[[4]]
    bias[2, 5, subiter, iter] <- altM$mean$Nintotal
    bias[2, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nintotal - altVals[[4]])/altVals[[4]])) * 100)
      
    bias[3, 1, subiter, iter] <- N
    bias[3, 2, subiter, iter] <- realM$mean$Nreal
    bias[3, 3, subiter, iter] <- (abs(mean((realM$sims.list$Nreal - N)/N)) * 100)
    bias[3, 4, subiter, iter] <- N
    bias[3, 5, subiter, iter] <- altM$mean$Nreal
    bias[3, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nreal - N)/N)) * 100)
      
    bias[4, 1, subiter, iter] <- Ntotal
    bias[4, 2, subiter, iter] <- realM$mean$Nrealtotal
    bias[4, 3, subiter, iter] <- (abs(mean((realM$sims.list$Nrealtotal - Ntotal)/Ntotal)) * 100)
    bias[4, 4, subiter, iter] <- Ntotal
    bias[4, 5, subiter, iter] <- altM$mean$Nrealtotal
    bias[4, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nrealtotal - Ntotal)/Ntotal)) * 100)
      
    bias[5, 1, subiter, iter] <- sigma
    bias[5, 2, subiter, iter] <- mean(realM$mean$sigma)
    bias[5, 3, subiter, iter] <- (abs(mean((rowMeans(realM$sims.list$sigma) - sigma)/sigma)) * 100)
    bias[5, 4, subiter, iter] <- sigma
    bias[5, 5, subiter, iter] <- mean(altM$mean$sigma)
    bias[5, 6, subiter, iter] <- (abs(mean((rowMeans(altM$sims.list$sigma) - sigma)/sigma)) * 100)
    
    #---------------------#
    #-Overlap Correlation-#
    #---------------------#
    
    overlapcor[iter, subiter] <- cor(OVA[,2], (abs(colMeans((realM$sims.list$sigma - sigma)/sigma)) * 100))
    
    #--------------------#
    #-Next sub iteration-#
    #--------------------#
    
    subiter <- subiter + 1
    
  } #End sub loop
  
  #----------------#
  #-Next iteration-#
  #----------------#

  iter <- iter + 1
  
} #End main loop
)

summary.bias <- array(NA, dim = c(niter, 2))
for(i in 1:niter){
  for(j in 1:2){
  summary.bias[i,1] <- mean(bias[5,3,,i])
  summary.bias[i,2] <- mean(bias[5,6,,i])
  }
}

sim.save <- list(bias, overlapcor)
dput(sim.save, file ="C:/Users/farrm/Documents/GitHub/CDSM/Simulations/sim.save.R")
