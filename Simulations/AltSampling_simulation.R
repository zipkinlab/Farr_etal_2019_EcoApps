#sampling area
mdE <- 733848.5
mdN <- 9844633

Et <- mdE - (13000/2)
Nt <- mdN + (12650/2)

Ep <- seq((Et + 650), (Et + (13000 - 650)), 1300)
Np <- seq(Nt, (Nt - 12650), -253)

X <- rep(Ep, rep(length(Np), length(Ep)))
Y <- rep(Np, length(Ep))

#--------------------------#
#-Create Sampling Boundary-#
#--------------------------#
  
#Easting
xlim <- c(715304, 752393)

#Northing
ylim <- c(9831970, 9857296)

#---------------------------#
#-Simulate Parameter Values-#
#---------------------------#

#Number of groups
N <- 1000

#Simulate coordinates of groups
u1 <- runif(N, xlim[1], xlim[2]) 
u2 <- runif(N, ylim[1], ylim[2])

#Group size
lambda.group <- 2
gs <- rpois(N, lambda.group) + 1

#Abundance
Ntotal <- sum(gs)

#Half normal detection parameter
sigma <- 300

#Mid point of each distance class
midpt <- seq(25, 650, 25)

#Index for distance class
nG <- length(midpt)

#Width of distance class
v <- 25

#Transect half width
B <- 650

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
Dtot <- cbind(y, index, u1, u2, site, gs)

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

#---------------#
#-Model Offsets-#
#---------------#

#Search area of each site
A.site <- as.vector(c(11.6542, 11.9619, 12.4702, 12.5182, 10.7843, 10.2384, 10.7495, 
                      12.0545, 9.0114, 11.2589, 10.4075, 9.7834, 11.8226, 10.5295,
                      11.5376, 14.8511, 14.0352))

#------------#
#-BUGS Model-#
#------------#

sink("ssds.txt")
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
    
    #Fit statistic for detection
    dclassnew[i] ~ dcat(fc[1:nG, site[i]])
    Tobs[i] <- pow(1 - sqrt(fc[dclass[i], site[i]]), 2)
    Tobsnew[i] <- pow(1 - sqrt(fc[dclassnew[i], site[i]]), 2)
    
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
    lambda[j] <- exp(alpha[j])
    
    #Fit statistic for number of groups
    Nnew[j] ~ dpois(lambda[j])
    Tn[j] <- pow(sqrt(N[j]) - sqrt(lambda[j]), 2)
    Tnnew[j] <- pow(sqrt(Nnew[j]) - sqrt(lambda[j]), 2)
    
    #Chi squared fit statistic
    eval[j] <- pcap[j] * N[j]
    E[j] <- pow((y[j] - eval[j]),2)/(eval[j] + 0.5)
    y.new[j] ~ dbin(pcap[j], N[j])
    E.new[j] <- pow((y.new[j] - eval[j]),2)/(eval[j] + 0.5)
    
    #Group size
    gs.lam[j] <- exp(beta)
    
    }#End j loop
    
    for(i in 1:nobs){
    
    gs[i] ~ dpois(gs.lam[site[i]]) T(1,)
    
    #Group fit statistic
    gsnew[i] ~ dpois(gs.lam[site[i]]) T(1,)
    Tg[i] <- pow(gs[i] - gs.lam[site[i]], 2) / gs.lam[site[i]]
    Tgnew[i] <- pow(gsnew[i] - gs.lam[site[i]], 2) / gs.lam[site[i]]
    
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
    
    ##Fit statistics
    
    #Detection
    
    fit.obs <- sum(Tobs[])
    fit.obs.new <- sum(Tobsnew[])
    
    #Abundance
    fit.ab <- sum(Tn[])
    fit.ab.new <- sum(Tnnew[])
    
    #Group size
    fit.gs <- sum(Tg[])
    fit.gs.new <- sum(Tgnew[])
    
    
    fit <- sum(E[])
    fit.new <- sum(E.new[])
    }
    
    
    ",fill=TRUE)
sink()


#-------------------#
#-Compile BUGS data-#
#-------------------#

#Input data
str(data <- list(nG = nG, v = v, site = site, y = yobs, B = B, midpt = midpt,
                 nobs = nobs, dclass = dclass, nsites = nsites, gs = gs))

#Initial values
N.in <- yobs + 1

inits <- function(){list(N = N.in, sigma = runif(10, 50, 350))} 

#Parameters to monitor
params<-c('gs.lam', 'sigma', 'Nin', 'Nintotal', 'Nreal', 'Nrealtotal', 
          'fit', 'fit.new', 'fit.obs', 'fit.obs.new', 'fit.ab', 
          'fit.ab.new', 'fit.gs', 'fit.gs.new')

#MCMC settings

nc <- 3
ni <- 25000
nb <- 5000
nt <- 5

#----------------#
#-Run BUGS Model-#
#----------------#

ssds <- jags(data = data, inits = inits, parameters.to.save = params, model.file = "ssds.txt", 
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, store.data = TRUE)

#Visualize
plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, 
     yaxt = "n", xaxt = "n", ylab = "", xlab = "")
points(cbind(u1, u2), col = "black", pch = 20)
points(cbind(X,Y), col = "blue")
points(cbind(Din[,2], Din[,3]), col = "green")
points(cbind(Dcap[,2], Dcap[,3]), col = "red", pch = 20)
