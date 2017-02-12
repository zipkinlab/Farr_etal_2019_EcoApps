# RealvsAlt
Matthew Farr  
February 12th, 2017  

I thought I would address the sampling design in the methods, but I plan on putting the majority of this in an appendix. Below isn’t necessarily what I would like to put in the manuscript, but just my thoughts on the topic. Let me know what you think.

The purpose of the simulation study is to compare standard distance sampling survey design to our study’s twisty distance sampling survey design. Our study implemented twisty distance sampling survey design due to sampling constraints. This sampling design may create biases by overinflating the estimation of detection. Hiby and Krishna (2001) showed that twisty sampling design had a minimum impact on estimating detection, and we will confirm this notion by comparing the relative bias of twisty sampling to standard sampling. 

I have a simulation with multiple iterations up and running. I plan to have multiple iterations with sub-iterations for resampling. For each iteration, observations (groups not individuals) are distributed uniformly across the sampling area. For each sub-iteration, the iteration of uniformly distributed observations is resampled. Below is the annotated code for the simulation study. This Rmarkdown file only shows a single iteration / sub-iteration of the simulation. 

Rmarkdown file last edited by Matthew Farr (2/12/17).

Set seed

```r
set.seed(1985)
```

Set working directory

```r
setwd("C:/Users/farrm/Documents/GitHub/CDSM/Simulations")
```

Load R packages

```r
library(rgdal)
library(sp)
library(dplyr)
library(tidyr)
library(jagsUI)
```

Create sampling boundary where individuals will be simulated

```r
#Easting UTM
xlim <- c(715304, 752393)

#Northing UTM
ylim <- c(9831970, 9857296)
```

Simulate parameter values

```r
#Number of groups
N <- 1000

#Simulate UTM coordinates of groups
u1 <- runif(N, xlim[1], xlim[2]) 
u2 <- runif(N, ylim[1], ylim[2])

#Group size
lambda.group <- 2
cs <- rpois(N, lambda.group) + 1

#Abundance
Ntotal <- sum(cs)

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
```

Begin real sampling

Import transect shapefiles

```r
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
```

![](https://github.com/farrmt/CDSM/blob/master/Simulations/Image1.png)<!-- -->

Sample coordintaes from transect. Used to calculate distances of observed groups.

```r
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
```

Combine site coordinates

```r
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
```

![](https://github.com/farrmt/CDSM/blob/master/Simulations/Image2.png)<!-- -->

Initialize values

```r
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
```

Simulate distances and site of groups

```r
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
```

Harvest data

```r
#Dataframe that includes information on all groups
Dtot <- cbind(y, index, u1, u2, site, cs)

#Dataframe containing only groups within 650 meters to transect
Din <- Dtot[complete.cases(Dtot),]

#Number of groups within 650 meters
Nin <- length(Din[,1])

#Abundance within 650 meters
Nintotal <- sum(Din[,6])
```

Initialize data

```r
#Remove groups not within 650 meters
index <- index[y==1]
index <- index[!is.na(index)]

#Detection Probability
p <- NULL

#Number of captured ("detected") groups
ncap <- rep(NA, Nin)

#Distance Class
dclass <- rep(NA, Nin)
```

Simulate detection of groups less than 650 meters

```r
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
```

Harvest data

```r
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
```

Calculate the number of groups that are within more than 1 site (overlap)

```r
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
```

Create offset for sites with longer transects and sampling area

```r
#Search area (meters squared) of each site
A.site <- as.vector(c(11.6542, 11.9619, 12.4702, 12.5182, 10.7843, 10.2384, 10.7495, 
                      12.0545, 9.0114, 11.2589, 10.4075, 9.7834, 11.8226, 10.5295,
                      11.5376, 14.8511, 14.0352))
```

BUGS Model

```r
cat("
    model{ 
    
    ##Priors
    
    for(j in 1:nsites){
    
    #Abundance prior
    alpha[j] ~ dnorm(0, 0.01)
    
    #Detection prior
    sigma[j] ~ dunif(0, 500)
    
    }#End j loop
    
    #Group size prior
    beta ~ dunif(0, 50)
    
    ##Likelihood
    
    #Multinomial detection component
    for(i in 1:nobs){
    
    dclass[i] ~ dcat(fc[1:nG, site[i]])
    
    }#End i loop
    
    for(j in 1:nsites){
    
    #Construct cell probabilities for nG cells
    for(k in 1:nG){  
    
    #Half normal detection function at midpt (length of rectangle)
    p[k,j] <- exp(- midpt[k] * midpt[k] / (2 * sigma[j] * sigma[j])) 
    
    #Probability of x in each interval (width of rectangle)
    pi[k,j] <- v/B 
    
    #Detection probability for each interval (area of each rectangle)
    f[k,j] <- p[k,j] * pi[k,j] 
    
    #Conditional detection probability (scale to 1)
    fc[k,j] <- f[k,j] / pcap[j] 
    
    }#End k loop
    
    #Detection probability at each site (sum of rectangles)
    pcap[j] <- sum(f[1:nG,j])               
    
    #Observation process
    y[j] ~ dbin(pcap[j], N[j])
    
    #Description of latent number of groups
    N[j] ~ dpois(lambda[j])
    
    #Linear model for number of groups
    lambda[j] <- exp(alpha[j] + log(offset[j]))

    #Group size
    gs.lam[j] <- exp(beta)
    
    }#End j loop
    
    for(i in 1:nobs){
    
    gs[i] ~ dpois(gs.lam[site[i]]) T(1,)
    
    }#End i loop
    
    ##Derived quantities
    
    #Number of groups within
    Nin <- sum(N[1:nsites])
    
    for(j in 1:nsites){
    
    #Abundance within
    Ntotal[j] <- N[j] * gs.lam[j]
    
    } #End j loop
    
    Nintotal <- sum(Ntotal[])
    
    #Density
    D <- (939.316/164.4837)
    
    #Number of total groups in sampling boundary
    Nreal <- Nin * D
    
    #Abundance in sampling boundary
    Nrealtotal <- Nintotal * D

    }",fill=TRUE, file="ssds.txt")
```

Compile BUGS data

```r
#Input data
str(realD <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                 nobs = nobs, dclass = dclass, nsites = nsites, gs = gs, offset = A.site))
```

```
## List of 11
##  $ nG    : int 26
##  $ v     : num 25
##  $ site  : num [1:90] 2 6 12 1 14 13 9 7 8 12 ...
##  $ y     : num [1:17] 4 7 4 5 1 2 5 5 6 2 ...
##  $ B     : num 650
##  $ midpt : num [1:26] 12.5 37.5 62.5 87.5 112.5 ...
##  $ nobs  : num 90
##  $ dclass: num [1:90] 7 9 14 13 5 7 6 1 5 7 ...
##  $ nsites: num 17
##  $ gs    : num [1:90] 2 1 5 2 3 2 2 4 2 3 ...
##  $ offset: num [1:17] 11.7 12 12.5 12.5 10.8 ...
```

```r
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
```


```r
realM <- jags(data = realD, inits = inits, parameters.to.save = params, model.file = "ssds.txt", 
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 197
##    Unobserved stochastic nodes: 52
##    Total graph size: 3398
## 
## Initializing model
## 
## Adaptive phase, 100 iterations x 3 chains 
## If no progress bar appears JAGS has decided not to adapt 
##  
## 
##  Burn-in phase, 2000 iterations x 3 chains 
##  
## 
## Sampling from joint posterior, 10000 iterations x 3 chains 
##  
## 
## Calculating statistics....... 
## 
## Done.
```

Save and remove data for next sampling

```r
realVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal)

rm(X, Y, nsites, J, si, di, dclass, dst, q, 
   site, d, y, index, Dtot, Din, Nin, Nintotal,
   p, ncap, Dcap, y.new, miss, yobs, nobs, gs,
   N.in, inits)
```

Begin alternative sampling

Create alternative sampling transect. There are 10 transects that run north to south.

```r
#Sampling area middle UTM coordinate
mdE <- 733848.5
mdN <- 9844633

#Sampling area left corner UTM coordinate
Et <- mdE - (13000/2)
Nt <- mdN + (12650/2)

#Sample points from alternative transects
Ep <- seq((Et + 650), (Et + (13000 - 650)), 1300)
Np <- seq(Nt, (Nt - 12650), -253)
X <- rep(Ep, rep(length(Np), length(Ep)))
Y <- rep(Np, length(Ep))
```

![](https://github.com/farrmt/CDSM/blob/master/Simulations/Image3.png)<!-- -->

inititalize values

```r
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
```

Simulate data for distances and site for groups

```r
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
```

Harvest data

```r
#Dataframe that includes information on all groups
Dtot <- cbind(y, index, u1, u2, site, cs)

#Dataframe containing only groups within 650 meters to transect
Din <- Dtot[complete.cases(Dtot),]

#Number of groups within 650 meters
Nin <- length(Din[,1])

#Abundance within 650 meters
Nintotal <- sum(Din[,6])
```

Initialize data

```r
#Remove groups not within 650 meters
index <- index[y==1]
index <- index[!is.na(index)]

#Detection Probability
p <- NULL

#Number of captured ("detected") groups
ncap <- rep(NA, Nin)

#Distance Class
dclass <- rep(NA, Nin)
```

Simulate detection of groups less than 650 meters

```r
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
```

Harvest data

```r
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
```

Compile BUGS data. Reuse BUGS model, parameters to save, and MCMC settings.

```r
#Input data
str(altD <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                 nobs = nobs, dclass = dclass, nsites = nsites, gs = gs, offset = rep(1, nsites)))
```

```
## List of 11
##  $ nG    : int 26
##  $ v     : num 25
##  $ site  : num [1:88] 2 7 8 10 2 1 2 2 3 6 ...
##  $ y     : num [1:10] 13 17 10 7 8 6 4 8 5 10
##  $ B     : num 650
##  $ midpt : num [1:26] 12.5 37.5 62.5 87.5 112.5 ...
##  $ nobs  : num 88
##  $ dclass: num [1:88] 5 15 1 7 12 19 7 1 11 18 ...
##  $ nsites: num 10
##  $ gs    : num [1:88] 4 4 3 2 4 2 1 3 8 3 ...
##  $ offset: num [1:10] 1 1 1 1 1 1 1 1 1 1
```

```r
#Initial values
N.in <- yobs + 1

inits <- function(){list(N = N.in, sigma = runif(10, 50, 350))}
```

Run BUGS model.

```r
altM <- jags(data = altD, inits = inits, parameters.to.save = params, model.file = "ssds.txt", 
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, store.data = TRUE)
```

```
## 
## Processing function input....... 
## 
## Done. 
##  
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 186
##    Unobserved stochastic nodes: 31
##    Total graph size: 2170
## 
## Initializing model
## 
## Adaptive phase, 100 iterations x 3 chains 
## If no progress bar appears JAGS has decided not to adapt 
##  
## 
##  Burn-in phase, 2000 iterations x 3 chains 
##  
## 
## Sampling from joint posterior, 10000 iterations x 3 chains 
##  
## 
## Calculating statistics....... 
## 
## Done.
```

Save and remove data

```r
altVals <- list(cbind(Din[,2], Din[,3]), cbind(Dcap[,2], Dcap[,3]), Nin, Nintotal,
                cbind(X, Y))

rm(X, Y, nsites, J, si, di, dclass, dst, q, 
   site, d, y, index, Dtot, Din, Nin, Nintotal,
   p, ncap, Dcap, y.new, miss, yobs, nobs, gs,
   N.in, inits)
```

Absoulte relative bias estimates

```r
bias <- t(matrix(data = c(

#Number of groups in search area

#Real True
realVals[[3]],

#Real Estimate
realM$mean$Nin,

#Real Bias
(abs(mean((realM$sims.list$Nin - realVals[[3]])/realVals[[3]])) * 100),

#Alternative True
altVals[[3]],

#Alternative Estimate
altM$mean$Nin,

#Alternative Bias
(abs(mean((altM$sims.list$Nin - altVals[[3]])/altVals[[3]])) * 100),

#Abundance in search area

#Real True
realVals[[4]],

#Real Estimate
realM$mean$Nintotal,

#Real Bias
(abs(mean((realM$sims.list$Nintotal - realVals[[4]])/realVals[[4]])) * 100),

#Alternative True
altVals[[4]],

#Alternative Estimate
altM$mean$Nintotal,

#Alternative Bias
(abs(mean((altM$sims.list$Nintotal - altVals[[4]])/altVals[[4]])) * 100),

#Number of groups in survey boundary

#Real True
N,

#Real Estimate
realM$mean$Nreal,

#Real Bias
(abs(mean((realM$sims.list$Nreal - N)/N)) * 100),

#Alternative True
N,

#Alternative Estimate
altM$mean$Nreal,

#Alternative Bias
(abs(mean((altM$sims.list$Nreal - N)/N)) * 100),

#Abundance in survey boundary

#Real True
Ntotal,

#Real Estimate
realM$mean$Nrealtotal,

#Real Bias
(abs(mean((realM$sims.list$Nrealtotal - Ntotal)/Ntotal)) * 100),

#Alternative True
Ntotal,

#Alternative Estimate
altM$mean$Nrealtotal,

#Alternative Bias
(abs(mean((altM$sims.list$Nrealtotal - Ntotal)/Ntotal)) * 100),

#Detection parameter

#Real True
sigma,

#Real Estimate
mean(realM$mean$sigma),

#Real Bias
(abs(mean((rowMeans(realM$sims.list$sigma) - sigma)/sigma)) * 100),

#Alternative True
sigma,

#Alternative Estimate
mean(altM$mean$sigma),

#Alternative Bias
(abs(mean((rowMeans(altM$sims.list$sigma) - sigma)/sigma)) * 100)),

nrow = 6, ncol = 5))

colnames(bias) <- c("Real True", "Real Est", "Real Bias", "Alt True", "Alt Est", "Alt Bias" )
rownames(bias) <- c("#GroupsWithin", "AbundanceWithin", "#Groups", "Abundance", "Sigma")
```

Results

Real sampling output

```
## JAGS output for model 'ssds.txt', generated by jagsUI.
## Estimates based on 3 chains of 12000 iterations,
## burn-in = 2000 iterations and thin rate = 4,
## yielding 7500 total samples from the joint posterior. 
## MCMC ran for 0.516 minutes at time 2017-02-12 13:35:19.
## 
##                mean      sd     2.5%      50%    97.5% overlap0 f  Rhat
## sigma[1]    389.076  73.481  233.559  398.159  495.753    FALSE 1 1.000
## sigma[2]    416.385  59.462  286.545  425.777  496.578    FALSE 1 1.001
## sigma[3]    325.820  93.916  164.895  319.885  490.519    FALSE 1 1.001
## sigma[4]    270.756  95.259  133.835  251.800  478.308    FALSE 1 1.002
## sigma[5]    284.752 116.283   95.162  277.577  487.760    FALSE 1 1.000
## sigma[6]    290.353 107.980  117.772  280.432  486.910    FALSE 1 1.000
## sigma[7]    313.222  92.816  163.957  303.380  486.322    FALSE 1 1.000
## sigma[8]    244.615  91.897  120.105  222.131  464.530    FALSE 1 1.004
## sigma[9]    251.395  93.129  126.296  227.629  472.552    FALSE 1 1.001
## sigma[10]   247.970 112.752   88.222  227.274  480.231    FALSE 1 1.000
## sigma[11]   313.270  95.864  155.747  304.402  489.057    FALSE 1 1.000
## sigma[12]   241.948  77.432  139.772  222.467  442.054    FALSE 1 1.002
## sigma[13]   258.335  78.388  152.681  239.257  459.375    FALSE 1 1.000
## sigma[14]   319.019 101.923  141.369  317.755  490.973    FALSE 1 1.000
## sigma[15]   370.050  74.956  230.164  371.478  492.881    FALSE 1 1.000
## sigma[16]   286.559  99.615  133.445  270.339  484.758    FALSE 1 1.001
## sigma[17]   253.567  76.683  149.746  235.691  451.414    FALSE 1 1.007
## Nin         184.935  20.562  149.000  184.000  230.000    FALSE 1 1.001
## Nintotal    480.679  63.177  371.372  475.968  617.387    FALSE 1 1.001
## Nreal      1056.105 117.425  850.893 1050.768 1313.460    FALSE 1 1.001
## Nrealtotal 2745.008 360.783 2120.794 2718.108 3525.706    FALSE 1 1.001
## deviance    892.781   6.434  881.542  892.332  907.060    FALSE 1 1.001
##            n.eff
## sigma[1]    5926
## sigma[2]    1974
## sigma[3]    1764
## sigma[4]     989
## sigma[5]    6492
## sigma[6]    3946
## sigma[7]    7500
## sigma[8]     557
## sigma[9]    1752
## sigma[10]   3430
## sigma[11]   7500
## sigma[12]   1393
## sigma[13]   7500
## sigma[14]   6000
## sigma[15]   7500
## sigma[16]   3122
## sigma[17]    395
## Nin         1345
## Nintotal    1847
## Nreal       1345
## Nrealtotal  1847
## deviance    1932
## 
## Successful convergence based on Rhat values (all < 1.1). 
## Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
## For each parameter, n.eff is a crude measure of effective sample size. 
## 
## overlap0 checks if 0 falls in the parameter's 95% credible interval.
## f is the proportion of the posterior with the same sign as the mean;
## i.e., our confidence that the parameter is positive or negative.
## 
## DIC info: (pD = var(deviance)/2) 
## pD = 20.7 and DIC = 913.466 
## DIC is an estimate of expected predictive error (lower is better).
```

![](https://github.com/farrmt/CDSM/blob/master/Simulations/Image4.png)<!-- -->

Alternative sampling output

```r
altM
```

```
## JAGS output for model 'ssds.txt', generated by jagsUI.
## Estimates based on 3 chains of 12000 iterations,
## burn-in = 2000 iterations and thin rate = 4,
## yielding 7500 total samples from the joint posterior. 
## MCMC ran for 0.392 minutes at time 2017-02-12 13:35:52.
## 
##                mean      sd     2.5%      50%    97.5% overlap0 f  Rhat
## sigma[1]    351.788  75.634  221.015  346.438  490.079    FALSE 1 1.001
## sigma[2]    376.664  69.642  248.988  376.323  493.149    FALSE 1 1.001
## sigma[3]    378.670  73.096  239.335  382.402  493.885    FALSE 1 1.002
## sigma[4]    235.224  85.667  123.081  213.355  455.628    FALSE 1 1.001
## sigma[5]    377.701  74.010  234.794  382.221  494.335    FALSE 1 1.001
## sigma[6]    384.479  73.506  236.607  391.456  494.424    FALSE 1 1.000
## sigma[7]    390.490  72.861  237.450  399.899  495.055    FALSE 1 1.000
## sigma[8]    295.936  87.332  164.010  280.666  482.403    FALSE 1 1.001
## sigma[9]    340.675  87.958  185.739  338.314  491.342    FALSE 1 1.000
## sigma[10]   192.056  60.652  115.981  178.638  362.005    FALSE 1 1.003
## Nin         162.500  17.711  132.000  161.000  201.000    FALSE 1 1.001
## Nintotal    447.154  57.711  347.148  442.144  572.064    FALSE 1 1.000
## Nreal       927.989 101.144  753.812  919.422 1147.849    FALSE 1 1.001
## Nrealtotal 2553.561 329.568 1982.457 2524.951 3266.883    FALSE 1 1.000
## deviance    878.253   5.441  869.378  877.682  890.427    FALSE 1 1.000
##            n.eff
## sigma[1]    2303
## sigma[2]    2742
## sigma[3]    1058
## sigma[4]    7500
## sigma[5]    3085
## sigma[6]    3119
## sigma[7]    7297
## sigma[8]    1633
## sigma[9]    7500
## sigma[10]   3865
## Nin         2406
## Nintotal    3918
## Nreal       2406
## Nrealtotal  3918
## deviance    7500
## 
## Successful convergence based on Rhat values (all < 1.1). 
## Rhat is the potential scale reduction factor (at convergence, Rhat=1). 
## For each parameter, n.eff is a crude measure of effective sample size. 
## 
## overlap0 checks if 0 falls in the parameter's 95% credible interval.
## f is the proportion of the posterior with the same sign as the mean;
## i.e., our confidence that the parameter is positive or negative.
## 
## DIC info: (pD = var(deviance)/2) 
## pD = 14.8 and DIC = 893.056 
## DIC is an estimate of expected predictive error (lower is better).
```

```r
plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, 
     yaxt = "n", xaxt = "n", ylab = "", xlab = "")
points(cbind(u1, u2), col = "black", pch = 20)
points(altVals[[5]], col = "blue")
points(altVals[[1]], col = "green")
points(altVals[[2]], col = "red", pch = 20)
```

![](https://github.com/farrmt/CDSM/blob/master/Simulations/Image5.png)<!-- -->

Bias output

```r
bias
```

```
##                 Real True  Real Est Real Bias Alt True   Alt Est  Alt Bias
## #GroupsWithin         170  184.9347 8.7850980      184  162.5003 11.684638
## AbundanceWithin       484  480.6786 0.6862402      534  447.1542 16.263254
## #Groups              1000 1056.1052 5.6105203     1000  927.9892  7.201078
## Abundance            2963 2745.0082 7.3571329     2963 2553.5607 13.818403
## Sigma                 300  298.6524 0.4491995      300  332.3682 10.789402
```

Overlap correlation with detection (sigma)

```r
cor(OVA[,2], (abs(colMeans((realM$sims.list$sigma - sigma)/sigma)) * 100))
```

```
##                 [,1]
## overlaps.x -0.160825
```
