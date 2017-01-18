#-----------------------------------------------------------#
#----Simulation of distance sampling data to check model.---#
#----Both models are similar to the model described in------#
#----AHM (Kery and Royle 2015) 8.5.3------------------------#
#----Script last edited by Matthew Farr (1/18/17)-----------#
#-----------------------------------------------------------#

#----------------#
#-Set Parameters-#
#----------------#

#Number of sites
nsites <- 100

#Mean group size (rate parameter)
lambda <- 5

#Detection parameter
sigma <- 225

#Transect half-width
B <- 650

#------------------#
#-Number of Groups-#
#------------------#

N <- rpois(nsites, lambda)
N.true <- N

#------------#
#-Group Size-#
#------------#

lambda.group <- 2

#---------------#
#-Simulate Data-#
#---------------#

data <- NULL
for(i in 1:nsites){
  if (N[i] == 0) {
    data <- rbind(data, c(i, NA, NA, NA, NA))
    next
  }
  d <- runif(N[i], 0, B)
  p <- exp(-d * d / (2 * sigma * sigma))
  gs <- rpois(N[i], lambda.group) + 1
  y <- rbinom(N[i], 1, p)
  d <- d[y == 1]
  ncap <- y[y == 1]
  gs <- gs[y == 1]
  if(sum(y) > 0)
    data <- rbind(data, cbind(rep(i, sum(ncap)), ncap, d, gs))
  else data <- rbind(data, c (i, NA, NA, NA))
}

ncap <- table(data[,1])            # ncap = 1 if no individuals captured
sites0 <- data[is.na(data[,2]),][,1] # sites where nothing was seen
ncap[as.character(sites0)] <- 0    # Fill in 0 for sites with no detections
ncap <- as.vector(ncap)            # Number of individuals detected per site

# Other data
site <- data[!is.na(data[,2]),1]   # Site ID of each observation
delta <- 25                       # Distance bin width for rect. approx.
midpt <- seq(delta/2, B, delta)    # Make mid-points and chop up data
dclass <- data[,3] %/% delta + 1   # Convert distance to distance category
nD <- length(midpt)                # Number of distance intervals
dclass <- dclass[!is.na(data[,2])] # Observed categorical observations
nind <- length(dclass)             # Total number of individuals detected

library(jagsUI)

sink("ssds_mdpt.txt")
cat("
    model{ 
    # Priors
    
    for(j in 1:nsites){
    alpha[j] ~ dnorm(0, 0.01)
    sigma[j] ~ dunif(0, 500)
    }
    
    # Multinomial Component
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[1:nG, site[i]]) # Part 1 of HM
    }
    
    for(j in 1:nsites){
    
    # construct cell probabilities for nG cells
    for(k in 1:nG){  
    
    # half normal detection function at midpt (length of rectangle)
    p[k,j] <- exp(- midpt[k] * midpt[k] / (2 * sigma[j] * sigma[j])) 
    
    # probability of x in each interval (width of rectangle)
    pi[k,j] <- v/B 
    
    # detection probability for each interval (area of each rectangle)
    f[k,j] <- p[k,j] * pi[k,j] 
    
    # conditional detection probability (scale to 1)
    fc[k,j] <- f[k,j] / pcap[j] 
    }
    
    # detection probability at each site (sum of rectangles)
    pcap[j] <- sum(f[1:nG,j])               
    
    
    y[j] ~ dbin(pcap[j], N[j])   # Part 2 of HM
    N[j] ~ dpois(lambda[j]) # Part 3 of HM
    lambda[j] <- exp(alpha[j])    # linear model for abundance
    
    #Chi squared fit statistic
    eval[j] <- pcap[j] * N[j]
    E[j] <- pow((y[j] - eval[j]),2)/(eval[j] + 0.5)
    y.new[j] ~ dbin(pcap[j], N[j])
    E.new[j] <- pow((y.new[j] - eval[j]),2)/(eval[j] + 0.5)
    }
    
    # Derived params
    Ntotal <- sum(N[1:nsites])
    D <- (Ntotal/164.4837)
    Nreal <- D * 939.316
    
    #Chi squared fit statistic
    fit <- sum(E[])
    fit.new <- sum(E.new[])
    }
    
    
    ",fill=TRUE)
sink()



### compile data for JAGS model

str(data1 <- list(nG = nD, v = delta, site = site, y = ncap, B = B, midpt = midpt,
                 nind = nind, dclass = dclass, nsites = nsites))

### create initial values
N.in <- ncap + 1

inits1 <- function(){list(N = N.in, sigma = runif(nsites, 150, 300))} 


### set parameters to monitor
params1<-c('Ntotal', 'N')

### mcmc settings

nc <- 3
ni <- 12000
nb <- 2000
nt <- 1

### run model

ssds <- jags(data = data1, inits = inits1, parameters.to.save = params1, model.file = "ssds_mdpt.txt", 
             n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, store.data = TRUE)

cat("
    model{ 
    # Priors
    
    for(j in 1:nsites){
    alpha0[j] ~ dunif(-10,10)
    sigma[j] ~ dunif(0, 500)
    }
    
    for(i in 1:nind){
    dclass[i] ~ dcat(fc[site[i],]) # Part 1 of HM
    }
    for(s in 1:nsites){
    # Construct cell probabilities for nD distance bands
    for(g in 1:nD){                # midpt = mid-point of each band
    log(p[s,g]) <- -midpt[g] * midpt[g] / (2 * sigma[s] * sigma[s])
    pi[s,g] <- ((2 * midpt[g] ) / (B * B)) * delta # prob. per interval
    f[s,g] <- p[s,g] * pi[s,g]
    fc[s,g] <- f[s,g] / pcap[s]
    }
    pcap[s] <- sum(f[s,])           # Pr(capture): sum of rectangular areas
    
    ncap[s] ~ dbin(pcap[s], N[s])   # Part 2 of HM
    N[s] ~ dpois(lambda[s])         # Part 3 of HM
    log(lambda[s]) <- alpha0[s] # linear model abundance
    }
    
    # Derived parameters
    Ntotal <- sum(N[])
    }
    ",fill=TRUE, file="model4.txt")

str(data2 <- list(nD = nD, delta = delta, site = site, ncap = ncap, B = B, midpt = midpt,
                  nind = nind, dclass = dclass, nsites = nsites))

ssds2 <- jags(data = data2, inits = inits1, parameters.to.save = params1, model.file = "model4.txt", 
              n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, store.data = TRUE)
