#-------------------------------------------------------------------------------------------------#
#----Simulation study comparing winding sampling design to alternative straight line sampling.----#
#----Data is simmulated for twisty sampling design of transects using actual transects.-----------#
#----Transects were imported as shapefiles. Alternative sampling design is simulated.-------------#
#----Script created by Matthew Farr --------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

#----------#
#-Set Seed-#
#----------#

set.seed(1985)

#-----------------------#
#-Set Working Directory-#
#-----------------------#

setwd("./SupportingInformation")

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
d.dir <- "./Transects"

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

modFile <- "HMSDS_model.txt"

#------------------#
#-Begin Iterations-#
#------------------#

#Number of iterations
niter <- 100	
subniter <- 10

#Starting iteration
iter <- 1

#Output vectors
bias <- array(NA, dim = c(5, 6, subniter, niter), 
              dimnames = list(c("Groups In", "Abundance In", "Groups", "Abundance", "Sigma"), 
                              c("Twst True", "Twst Est", "Twst Bias", "Str True", "Str Est", "Str Bias"),
                              NULL, NULL))

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
    twistyD <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                      nobs = nobs, dclass = dclass, nsites = nsites, gs = gs, offset = A.site)
    
    #Initial values
    N.in <- yobs + 1
    
    inits <- function(){list(N = N.in, sigma = runif(17, 50, 350))} 
    
    #Parameters to monitor
    params<-c('sigma', 'Nin', 'Nintotal', 'Ntwisty', 'Ntwistytotal')
    
    #MCMC settings
    
    nc <- 3
    ni <- 12000
    nb <- 2000
    nt <- 4
    
    #----------------#
    #-Run BUGS Model-#
    #----------------#
    
    twistyM <- jags(data = twistyD, inits = inits, parameters.to.save = params, model.file = modFile, 
                  n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, parallel = TRUE)
    
    #----------------------------------------#
    #-Save and Remove Data for Next Sampling-#
    #----------------------------------------#
    
    twistyVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal, mean(p))
    
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
    
    altVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal, mean(p))
    
    rm(X, Y, nsites, J, si, di, dclass, dst, q, 
       site, d, y, index, Dtot, Din, Nin, Nintotal,
       p, ncap, Dcap, y.new, miss, yobs, nobs, gs,
       N.in, inits)
    
    #----------------#
    #-Bias Estimates-#
    #----------------#
      
    bias[1, 1, subiter, iter] <- twistyVals[[3]]
    bias[1, 2, subiter, iter] <- twistyM$mean$Nin
    bias[1, 3, subiter, iter] <- (abs(mean((twistyM$sims.list$Nin - twistyVals[[3]])/twistyVals[[3]])) * 100)
    bias[1, 4, subiter, iter] <- altVals[[3]]
    bias[1, 5, subiter, iter] <- altM$mean$Nin
    bias[1, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nin - altVals[[3]])/altVals[[3]])) * 100)
      
    bias[2, 1, subiter, iter] <- twistyVals[[4]]
    bias[2, 2, subiter, iter] <- twistyM$mean$Nintotal
    bias[2, 3, subiter, iter] <- (abs(mean((twistyM$sims.list$Nintotal - twistyVals[[4]])/twistyVals[[4]])) * 100)
    bias[2, 4, subiter, iter] <- altVals[[4]]
    bias[2, 5, subiter, iter] <- altM$mean$Nintotal
    bias[2, 6, subiter, iter] <- (abs(mean((altM$sims.list$Nintotal - altVals[[4]])/altVals[[4]])) * 100)
      
    bias[3, 1, subiter, iter] <- N
    bias[3, 2, subiter, iter] <- twistyM$mean$Ntwisty
    bias[3, 3, subiter, iter] <- (abs(mean((twistyM$sims.list$Ntwisty - N)/N)) * 100)
    bias[3, 4, subiter, iter] <- N
    bias[3, 5, subiter, iter] <- altM$mean$Ntwisty
    bias[3, 6, subiter, iter] <- (abs(mean((altM$sims.list$Ntwisty - N)/N)) * 100)
      
    bias[4, 1, subiter, iter] <- Ntotal
    bias[4, 2, subiter, iter] <- twistyM$mean$Ntwistytotal
    bias[4, 3, subiter, iter] <- (abs(mean((twistyM$sims.list$Ntwistytotal - Ntotal)/Ntotal)) * 100)
    bias[4, 4, subiter, iter] <- Ntotal
    bias[4, 5, subiter, iter] <- altM$mean$Ntwistytotal
    bias[4, 6, subiter, iter] <- (abs(mean((altM$sims.list$Ntwistytotal - Ntotal)/Ntotal)) * 100)
      
    bias[5, 1, subiter, iter] <- sigma
    bias[5, 2, subiter, iter] <- mean(twistyM$mean$sigma)
    bias[5, 3, subiter, iter] <- (abs(mean((rowMeans(twistyM$sims.list$sigma) - sigma)/sigma)) * 100)
    bias[5, 4, subiter, iter] <- sigma
    bias[5, 5, subiter, iter] <- mean(altM$mean$sigma)
    bias[5, 6, subiter, iter] <- (abs(mean((rowMeans(altM$sims.list$sigma) - sigma)/sigma)) * 100)
    
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

save(bias, file ="sim.save.R")

summary <- array(NA, dim = c(5,2))

summary <- apply(bias[,c(3,6),,], MARGIN = c(1,2), FUN = mean)

colnames(summary) <- c("Winding", "Straight")

save(summary, file = "summary.R")
