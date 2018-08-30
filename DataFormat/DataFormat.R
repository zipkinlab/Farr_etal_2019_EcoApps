#-------------------------------------------#
#----Formatting DSdata from raw csv file----#
#----Created by Matthew Farr----------------#
#-------------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("./DataFormat")

#----------------#
#-Load libraries-#
#----------------#

library(dplyr)
library(tidyr)

#------------#
#-Import CSV-#
#------------#

raw <- read.csv("../RawData/HMSDS_rawdata.csv", header=TRUE)
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
dclass <- raw$DSclass

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

#Index for nine species
aspec <- c(1,2,3,6,8)

social <- c(1,2,3,5,6,8,11)

#-----------------------------------#
#-Import group size of observations-#
#-----------------------------------#

#Filter out solitary species
G <- filter(raw, raw$Animal != "Caracal", raw$Animal != "Leopard",
            raw$Animal != "Serval", raw$Animal != "SideStripedJackal") 
#Group size
gs <- G$Count

gs1 <- raw$Count

#Site ID for social species
s.site <- G$Site_ID

#Replicate ID
s.rep <- G$Sample_ID

#Species ID
s.spec <- as.integer(G$Animal)

#Number of social species observations
nsoc <- length(gs)

#--------------------#
#-Average group size-#
#--------------------#

#NEITHER USED IN FINAL DETECTION MODEL

#Observed from data
ags <- rep(1, 11)
for(i in social){ags[i] <- mean(raw$Count[spec==i])}
ags <- scale(ags)

#Reported from Gittleman 1989 pg 189 - 191
ags2 <- c(17, 2, 2, 1, 1, 6, 1, 9, 1, 2, 2)
ags2 <- scale(ags2)

#-------------------#
#-Average body size-#
#-------------------#

#Reported from Gittleman 1989 pg 189 - 191
abs <- c(1.26, 3.94, 7.69, 11.59, 58.56, 51.94, 52.46, 166.02, 11.70, 11.25, 0.49)
abs <- scale(abs)

#----------------------------#
#-Offset for transect length-#
#----------------------------#

offset <- as.vector(c(1, 1, 1, 1, 1, 1, 1, 1.080,
                      0.878, 1, 1, 1, 1, 1, 1.100,
                      1.300, 1.237))

#-------------------------#
#-Create Region Covariate-#
#-------------------------#

region <- c(rep(0, 13), rep(1, 4))

#-------------------#
#-Vehicle Covariate-#
#-------------------#

#NOT USED IN FINAL DETECTION MODEL

car <- array(1, dim = c(16,17))
car[1, 1:8] <- 3
car[5, 1:13] <- 1
car[3, 9:13] <- 2
car[10, 9:13]  <- 3
car[c(1, 7:9), 14:17] <- 2
car[2, 14:17] <- 2
car[c(3,4,6,11:13), 14:17] <- 3
car[10, 14:17] <- 3

#--------------#
#-Combine data-#
#--------------#

DSdata <- list(y, dclass, v, B, site, rep, spec, mdpt, nsites, nreps, nspec, nobs, 
               nD, aspec, gs, gs1, s.site, s.rep, s.spec, nsoc, region, offset, car, social, ags, ags2, abs)
heads <- c("y", "dclass", "v", "B", "site", "rep", "spec", "mdpt", "nsites", "nreps", 
           "nspec", "nobs", "nD", "aspec", "gs", "gs1", "s.site", "s.rep", "s.spec", "nsoc", "region",
           "offset", "car", "social", "ags", "ags2", "abs")
DSdata <- setNames(DSdata, nm = heads)

#-------------#
#-Export data-#
#-------------#

dput(DSdata, file = "DSdata")
