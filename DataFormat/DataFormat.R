#----------------------------------------------#
#----Formatting CDSM data from raw csv file----#
#----Last edited by Matthew Farr 2/17/17-------#
#----------------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#
setwd("C:/Users/farrm/Documents/GitHub/CDSM/DataFormat")

#----------------#
#-Load libraries-#
#----------------#

library(dplyr)
library(tidyr)

#------------#
#-Import CSV-#
#------------#

raw <- read.csv("C:/Users/farrm/Documents/GitHub/CDSM/RawData/CDSM_rawdata.csv", header=TRUE)
raw <- tbl_df(raw)

#Sort by species
raw <- arrange(raw, Animal)

#-------------------------------------------------#
#-Create observation array (rep x site x species)-#
#-------------------------------------------------#

#Filter out unnecessary data
D <- data.frame(raw$Animal, raw$Site_ID, raw$Sample_ID)
colnames(D) <- c("Animal", "Site_ID", "Sample_ID")

#Vector of species
name <- c("BandedMongoose", "BatEaredFox", "BlackBackedJackal", 
          "Caracal", "Cheetah", "Hyena", "Leopard", "Lion", "Serval",
          "SideStripedJackal", "SlenderMongoose")

#Initialize observation array (rep x site x species)
y <- array(NA, dim = c(16,17,11))

#Generate observation array
for(s in 1:11){
  J <- (filter(D, Animal == name[s]))
  J <- group_by(J, Site_ID, Sample_ID, Animal)%>%summarize(n())
  W <- data.frame(rep(1:17, rep(16, 17)), rep(1:16, 17))
  colnames(W) <- c("Site_ID", "Sample_ID")
  Y <- full_join(W, J, by.x = c("Site_ID", "Sample_ID"), by.y = c("Site_ID", "Sample_ID"))
  Y$`n()`[is.na(Y$`n()`)] = 0
  X <- split(Y$`n()`, f = Y$Site_ID)
  X <- do.call(cbind, X)
  for(j in 14:17){
    for(t in 14:16){
      X[t,j] = NA
    }
  }
  y[,,s] <- X
}

#-------------------------#
#-Create distance classes-#
#-------------------------#

#Distance class were created in ArcGIS using minimum distance
Dst <- raw$DSclass

#Width of distance classes
v <- 50 #meters

#Transect half-width
B <- 650 #meters

#------------------------#
#-Create IDs and indices-#
#------------------------#

#Site ID
site <- raw$Site_ID

#Replicate ID
rep <- raw$Sample_ID

#Species ID
spec <- as.integer(raw$Animal)

#Distance class midpoint ID
mdpt <- seq(25, 625, 50)

#Number of sites
nsites <- as.integer(max(raw$Site_ID))

#Number of reps
nreps <- as.vector(c(rep(16, 13), rep(13, 4)))

#Number of speceis
nspec <- 11

#Number of observations
nobs <- sum(y, na.rm = TRUE)

#Number of distance classes
nD <- length(mdpt)


#-----------------------------------#
#-Import group size of observations-#
#-----------------------------------#

gs <- raw$count

#----------------------------#
#-Sampling area of each site-#
#----------------------------#

area <- as.vector(c(11.6542, 11.9619, 12.4702, 12.5182, 10.7843, 10.2384, 10.7495, 
                    12.0545, 9.0114, 11.2589, 10.4075, 9.7834, 11.8226, 10.5295,
                    11.5376, 14.8511, 14.0352))

#-------------------------#
#-Create Region Covariate-#
#-------------------------#

region <- c(rep(0, 13), rep(1, 4))

#--------------#
#-Combine data-#
#--------------#

DSdata <- list(y, Dst, v, B, site, rep, spec, mdpt, nsites, nreps, nspec, nobs, 
               nD, gs, area, region)
heads <- c("y", "Dst", "v", "B", "site", "rep", "spec", "mdpt", "nsites", "nreps", 
                          "nspec", "nobs", "nD", "gs", "area", "region")
DSdata <- setNames(DSdata, nm = heads)

#-------------#
#-Export data-#
#-------------#

dput(DSdata, file = "C:/Users/farrm/Documents/GitHub/CDSM/DataFormat/DSdata")