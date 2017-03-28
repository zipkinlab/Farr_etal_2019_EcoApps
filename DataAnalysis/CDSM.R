#---------------------------------------------------------#
#----Community distance sampling model -------------------#
#----Script lasted edited by Matthew Farr (3/03/17)-------#
#---------------------------------------------------------#

#-----------------------#
#-Set Working Directory-#
#-----------------------#

setwd("C:/Users/farrm/Documents/GitHub/CDSM/DataAnalysis")

#---------------#
#-Load Libaries-#
#---------------#

library(jagsUI)

#-----------#
#-Load Data-#
#-----------#

DSdata <- dget("C:/Users/farrm/Documents/GitHub/CDSM/DataFormat/DSdata")

#---------------#
#-Attach DSdata-#
#---------------#

attach(DSdata)

#------------#
#-BUGS Model-#
#------------#

sink("CDSM.txt")
cat("
    model{
    
    #Priors
    
    #Species-specific Parameters
    
    for(s in 1:nspec){
    
    #Abundance parameters
    alpha0[s] ~ dnorm(mu_a0, tau_a0) #intercept of lambda linear predictor
    
    alpha1[s] ~ dnorm(mu_a1, tau_a1) #coefficient for region covariate
    
    #Detection parameter
    asig[s] ~ dnorm(mu_s, tau_s) #intercept of sigma linear predictor
    
    #Group size parameter
    beta0[s] ~ dunif(0,50) #intercept of lam.gs linear predictor
    
    beta1[s] ~ dnorm(0, 0.01) #coefficient for region covariate
    
    }
    
    #Hyperparameters
    
    #Sigma
    mu_s ~ dnorm(0, 0.01)
    tau_s <- 1/(sig_s * sig_s)
    sig_s ~ dunif(0, 500)
    
    #Alpha0
    mu_a0 ~ dnorm(0, 0.01)
    sig_a0 <- 1/sqrt(tau_a0)
    tau_a0 ~ dgamma(0.1, 0.1)
    
    #Alpha1
    mu_a1 ~ dnorm(0, 0.01)
    sig_a1 <- 1/sqrt(tau_a1)
    tau_a1 ~ dgamma(0.1, 0.1)
    
    #Group Size Gamma Parameter
    
    r.N ~ dunif(0,100)
    
    #Likelihood
    
    ## Part 1
    
    for(i in 1:nobs){
    
    dclass[i] ~ dcat(fc[1:nD, rep[i], site[i], spec[i]])
    
    }#end i loop
    
    for(s in 1:nspec){
    
    for(j in 1:nsites){
    
    for(t in 1:nreps[j]){
    
    # construct cell probabilities for nG cells using numerical integration
    # sum of the area (rectangles) under the detection function
    for(k in 1:nD){
    
    # half normal detection function at midpt (length of rectangle)
    g[k,t,j,s] <- exp(-mdpt[k]*mdpt[k]/(2*sigma[s]*sigma[s]))
    
    # probability of x in each interval (width of rectangle) for both sides of the transect
    pi[k,t,j,s] <- v/B
    
    # detection probability for each interval (area of each rectangle)
    f[k,t,j,s] <- g[k,t,j,s] * pi[k,t,j,s]
    
    # conditional detection probability (scale to 1)
    fc[k,t,j,s] <- f[k,t,j,s]/pcap[t,j,s]
    
    }#end k loop
    
    # detection probability at each site (sum of rectangles)
    pcap[t,j,s] <- sum(f[1:nD,t,j,s])
    
    ## Part 2
    
    #Observed population @ each t,j,s
    y[t,j,s] ~ dbin(pcap[t,j,s], N[t,j,s])
    
    ## Part 3
    
    #Population @ each t,j,s
    N[t,j,s] ~ dpois(lambda[t,j,s])
    
    #Linear Predictor for Lambda
    lambda[t,j,s] <- exp(alpha0[s] + alpha1[s] * region[j])
    
    }#end t loop
    
    }#end j loop
    
    #Linear Predictor for Sigma
    sigma[s] <- exp(asig[s])
    
    }#end s loop
    
    ## Group size
    
    for(s in 1:nspec){
    
    for(j in 1:nsites){

    for(t in 1:nreps[j]){
    
    gs.lam[t,j,s] <- exp(beta0[s] + beta1[s] * region[j] ) 
    gs.lam.star[t,j,s] <- gs.lam[t,j,s] * gs.rho[t,j,s]
    gs.rho[t,j,s] ~ dgamma(r.N, r.N)
    
    }#end t loop

    }#end j loop
    
    }#end s loop
    
    for(i in 1:nobs){
    
    gs[i] ~ dpois(gs.lam.star[rep[i], site[i], spec[i]]) T(1,)
    
    }#end i loop
    
    ## Derived quantities
    
    for(s in denspec){
    
    for(j in 1:nsites){
    
    for(t in 1:nreps[j]){

    GSrep[t,j,s] <- N[t,j,s] * gs.lam.star[t,j,s]

    }#end t loop

    GSsite[j,s] <- mean(GSrep[1:nreps[j], j, s])
    DenGSsite[j,s] <- GSsite[j,s] / area[j]
    
    }#end j loop
    
    GS[s] <- sum(GSsite[1:nsites, s])
    
    DenGS[s] <- mean(DenGSsite[1:nsites, s]) * 100 #per 100 km2

    RegDen[s,1] <- mean(DenGSsite[1:13, s]) * 100 #per 100 km2
    RegDen[s,2] <- mean(DenGSsite[14:17, s]) * 100 #per 100 km2
    
    }#end s loop
    
    }
    ",fill=TRUE)
sink()

#-------------------#
#-Compile BUGS data-#
#-------------------#

data <- list(nspec = nspec, nD = nD, v = v, area = area, site = site, rep = rep, spec = spec, 
             y = y, B = B, mdpt = mdpt, nobs = nobs, dclass = dclass, nsites = nsites, 
             nreps = nreps, gs = gs, region = region, denspec = denspec)

#---------------#
#-Inital values-#
#---------------#

Nst <- y + 1

inits <- function(){list(N = Nst, mu_a0 = runif(1, 0, 1), tau_a0 = runif(1, 0, 1),
                          mu_s = runif(1, 5, 6), sig_s=runif(1),asig = runif(nspec, 5, 6), 
                          mu_a1 = runif(1, 0, 1), tau_a1 = runif(1, 0, 1))}

#--------------------#
#-Parameters to save-#
#--------------------#

params<-c('mu_a0', 'mu_a1', 'mu_s', 'tau_a0', 'tau_a1', 'tau_s', 
          'alpha0', 'alpha1', 'beta0', 'beta1', 'GS', 'DenGS', 'RegDen')

#---------------#
#-MCMC settings-#
#---------------#

nc <- 3
ni <- 250000
nb <- 50000
nt <- 20

CDSM <- jags(data = data, inits = inits, parameters.to.save = params, model.file = "CDSM.txt", 
                         n.chains = nc, n.iter = ni, n.burnin = nb, n.thin = nt, store.data = FALSE, parallel = TRUE)

