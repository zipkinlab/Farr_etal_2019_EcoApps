#------------------------------------------------------#
#----Hierarchical community distance sampling model----#
#----Post analysis of HMSDS model and generation of----#
#----figures.------------------------------------------#
#----Created by Matthew Farr---------------------------#
#------------------------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("./PostAnalysis")

#----------------#
#-Load libraries-#
#----------------#

library(jagsUI)
library(ggplot2)
library(ggthemes)
library(grid)
library(gridExtra)
library(extrafont)
loadfonts()

#------------#
#-Load Files-#
#------------#

#Formatted Data
DSdata <- dget("../DataFormat/DSdata")
attach(DSdata)
#Model Output
load("../DataAnalysis/HMSDS.Rdata")

#---------#
#-Table 1-#
#---------#

table1 <- data.frame(cbind(spec, gs1, ifelse(site > 13, 1, 0)))
tr <- cbind(aggregate(gs1[V3 == 1] ~ spec[V3 == 1], table1, length)[,2],
            aggregate(gs1[V3 == 1] ~ spec[V3 == 1], table1, sum)[,2])
table1 <- data.frame(
  cbind(rbind(
    tr[1:3,], rep(0, 2), tr[4:5,], rep(0, 2), tr[6:9,]
  ),
  aggregate(gs1[V3 == 0] ~ spec[V3 == 0], table1, length)[,2],
  aggregate(gs1[V3 == 0] ~ spec[V3 == 0], table1, sum)[,2],
  aggregate(gs1 ~ spec, table1, length)[,2],
  aggregate(gs1 ~ spec, table1, sum)[,2]))

rm(tr)

#Rearrange species order
table1 <- table1[c(8,1:5,7,9:11,6),]

#NAs to non-grouping species
table1[c(5,7,8,9),c(1,3,5)] <- NA

#Total
table1[12,] <- apply(table1, 2, sum, na.rm = TRUE)
colnames(table1) <- c("Talek Groups", "Talek Ind", "MT Group", "MT Ind", "Group", "Ind")
rownames(table1) <- c("African Lion", "Banded Mongoose", "Bat-eared Fox", "Black-backed Jackal", "Caracal",
                   "Cheetah", "Leopard", "Serval", "Side-striped Jackal", 
                   "Slender Mongoose", "Spotted Hyena", "Total")

write.csv(table1, file = "Table1.csv")

#-----------------------#
#-Appendix S3: Table S1-#
#-----------------------#

tableS1 <- HMSDS$summary[,-c(4:6,8:11)]
tableS1 <- round(tableS1, digits = 2)
write.csv(tableS1, file = "AppendixS3TableS1.csv")

#---------#
#-Results-#
#---------#

#Listed in order as published

#Observation summaries
table1[12,6] #Number of individuals
length(unique(spec)) #Number of species
nobs #Number of observations = number of groups obs + number of ind obs

#Community parameters
tableA['mu_a1',c(1,3,4)] #Community-level effect of management regime on the expected number of groups (log scale)
tableA['mu_b1',c(1,3,4)] #Community-level effect of management regime on expected group size (log scale)

#Number of groups (species in groups) OR Abundance (solitary species)
mean((HMSDS$sims.list$alpha1[,8]) < 0) #Probability of expected number of lion groups per transect being higher in MT than TR
tableA['alpha1[8]',c(1,3,4)] #Expected number of lion groups per transect (log scale)
tableA['alpha1[7]',c(1,3,4)] #Expected leopard abundance per transect (log scale)
tableA['alpha1[9]',c(1,3,4)] #Expected serval abundance per transect (log scale)
mean((HMSDS$sims.list$alpha1[,7]) < 0) #Probability of expected abundance of leopards per transect being higher in MT than TR
mean((HMSDS$sims.list$alpha1[,9]) < 0) #Probability of expected abundance of servals per transect being higher in MT than TR
tableA['alpha1[3]',c(1,3,4)] #Expected number of black-backed jackal groups per transect (log scale)
tableA['alpha1[6]',c(1,3,4)] #Expected number of spotted hyena groups per transect (log scale)
mean((HMSDS$sims.list$alpha1[,3]) > 0) #Probability of expected number of black-backed jackal groups per transect being higher in TR than MT
mean((HMSDS$sims.list$alpha1[,6]) > 0) #Probability of expected number of spotted hyena groups per transect being higher in TR than MT

#Group size
tableA['beta1[6]',c(1,3,4)] #Expected group size of spotted hyena per transect (log scale)

#Density (only species with 20 or more observations)
mean((HMSDS$sims.list$RegGS[,8,1] - HMSDS$sims.list$RegGS[,8,2]) > 0) #Probability of lion density being higher in MT than TR
tableA['RegGS[8,2]', c(1,3,4)] #Lion density in TR
tableA['RegGS[8,1]', c(1,3,4)] #Lion density in MT
mean((HMSDS$sims.list$RegGS[,2,1] - HMSDS$sims.list$RegGS[,2,2]) > 0) #Probability of bat-eared fox density being higher in MT than TR
tableA['RegGS[2,2]', c(1,3,4)] #Bat-eared fox density in TR
tableA['RegGS[2,1]', c(1,3,4)] #Bat-eared fox density in MT
mean((HMSDS$sims.list$RegGS[,1,1] - HMSDS$sims.list$RegGS[,1,2]) > 0) #Probability of banded Mongoose density being higher in MT than TR
tableA['RegGS[1,2]', c(1,3,4)] #Bat-eared fox density in TR
tableA['RegGS[1,1]', c(1,3,4)] #Bat-eared fox density in MT
mean((HMSDS$sims.list$RegGS[,3,1] - HMSDS$sims.list$RegGS[,3,2]) < 0) #Probability of black-backed jackal density being higher in TR than MT
tableA['RegGS[3,2]', c(1,3,4)] #Black-backed jackal density in TR
tableA['RegGS[3,1]', c(1,3,4)] #Black-backed jackal density in MT
tableA['RegGS[6,2]', c(1,3,4)] #Spotted hyena density in TR
tableA['RegGS[6,1]', c(1,3,4)] #Spotted hyena density in MT
mean((HMSDS$sims.list$RegGS[,6,1] - HMSDS$sims.list$RegGS[,6,2]) < 0) #Probability of spotted hyena density being higher in TR than MT

#Detection probability
c(exp(mean(HMSDS$sims.list$mu_s)),  
exp(quantile(HMSDS$sims.list$mu_s, probs = 0.025)),
exp(quantile(HMSDS$sims.list$mu_s, probs = 0.975))) #Community-level intercept parameter for detection (inverse log scale)
sqrt(-2*exp(mean(HMSDS$sims.list$mu_s))^2*log(0.5)) #Average species' distance at 50% detection probability 
sqrt(-2*exp(mean(HMSDS$sims.list$mu_s))^2*log(0.1)) #Average species' distance at 10% detection probability
tableA['bsig2', c(1,3,4)] #Effect of average body size on detection
tableA['bsig1', c(1,3,4)] #Effect of management regime on detection
c(exp(mean(HMSDS$sims.list$asig[,11] + HMSDS$sims.list$bsig1*DSdata$abs[11])),
exp(quantile(HMSDS$sims.list$asig[,11] + HMSDS$sims.list$bsig1*DSdata$abs[11], probs = 0.025)),
exp(quantile(HMSDS$sims.list$asig[,11] + HMSDS$sims.list$bsig1*DSdata$abs[11], probs = 0.975))) #Scale parameter of slender mongoose
sqrt(-2*exp(mean(HMSDS$sims.list$asig[,11] + HMSDS$sims.list$bsig1*DSdata$abs[11]))^2*log(0.5)) #Slender mongoose distance at 50% detection probability
sqrt(-2*exp(mean(HMSDS$sims.list$asig[,11] + HMSDS$sims.list$bsig1*DSdata$abs[11]))^2*log(0.1)) #Slender mongoose distance at 10% detection probability
c(exp(mean(HMSDS$sims.list$asig[,8] + HMSDS$sims.list$bsig1*DSdata$abs[8])),
  exp(quantile(HMSDS$sims.list$asig[,8] + HMSDS$sims.list$bsig1*DSdata$abs[8], probs = 0.025)),
  exp(quantile(HMSDS$sims.list$asig[,8] + HMSDS$sims.list$bsig1*DSdata$abs[8], probs = 0.975))) #Scale parameter of lion
sqrt(-2*exp(mean(HMSDS$sims.list$asig[,8] + HMSDS$sims.list$bsig1*DSdata$abs[8]))^2*log(0.5)) #Lion distance at 50% detection probability
sqrt(-2*exp(mean(HMSDS$sims.list$asig[,8] + HMSDS$sims.list$bsig1*DSdata$abs[8]))^2*log(0.1)) #Lion distance at 10% detection probability

#----------#
#-Figure 1-#
#----------#

#Created in ArcGIS

#----------#
#-Figure 2-#
#----------#

alpha1.val <- cbind(c(HMSDS$q2.5$alpha1[c(8,1:5,7,9:11,6)], HMSDS$q2.5$mu_a1),
                    c(HMSDS$q25$alpha1[c(8,1:5,7,9:11,6)], HMSDS$q25$mu_a1),
                    c(HMSDS$mean$alpha1[c(8,1:5,7,9:11,6)], HMSDS$mean$mu_a1),
                    c(HMSDS$q75$alpha1[c(8,1:5,7,9:11,6)], HMSDS$q75$mu_a1),
                    c(HMSDS$q97.5$alpha1[c(8,1:5,7,9:11,6)], HMSDS$q97.5$mu_a1),
                    c(HMSDS$overlap0$alpha1[c(8,1:5,7,9:11,6)],HMSDS$overlap0$mu_a1))

sppnames <- c("African Lion", "Banded Mongoose", "Bat-eared Fox", "Black-backed Jackal", "Caracal",
              "Cheetah", "Leopard", "Serval", "Side-striped Jackal", 
              "Slender Mongoose", "Spotted Hyena", "The Community")

values <- data.frame(sppnames, alpha1.val)
colnames(values) <- c("species", "lower.alpha", "l25.alpha", "mean.alpha", "u75.alpha", "upper.alpha", "alpha.sig")

values$species <- factor(values$species, levels = values$species)

values$nodiff <- values$l25.alpha < 0 & values$u75.alpha > 0
values$H95 <- values$lower.alpha > 0
values$H50 <- values$lower.alpha < 0 & values$l25.alpha > 0
values$L95 <- values$upper.alpha < 0
values$L50 <- values$upper.alpha > 0 & values$u75.alpha < 0

Figure2 <- ggplot() + 
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = mean.alpha, ymax = mean.alpha, color = "black"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = lower.alpha, ymax = upper.alpha, color = "black"), 
                width = 0, size = 1.25) +
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = l25.alpha, ymax = u75.alpha, color = "black"), 
                width = 0, size = 3.5) +
  geom_errorbar(data = subset(values, H50 == TRUE), aes(x = species, ymin = mean.alpha, ymax = mean.alpha, color = "orangegrey"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, H50 == TRUE), aes(x = species, ymin = lower.alpha, ymax = upper.alpha, color = "orangegrey"), 
                size = 1.25, width = 0) +
  geom_errorbar(data = subset(values, H50 == TRUE), aes(x = species, ymin = l25.alpha, ymax = u75.alpha, color = "orangegrey"), 
                size = 3.5, width = 0) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = mean.alpha, ymax = mean.alpha, color = "greengrey"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = lower.alpha, ymax = upper.alpha, color = "greengrey"), 
                size = 1.25, width = 0) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = l25.alpha, ymax = u75.alpha, color = "greengrey"), 
                size = 3.5, width = 0) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = mean.alpha, ymax = mean.alpha, color = "green"),
                width = 0.25) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = lower.alpha, ymax = upper.alpha, color = "green"),
                size = 1.25, width = 0) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = l25.alpha, ymax = u75.alpha, color = "green"),
                size = 3.5, width = 0) +
  annotate("text", x = 0.5, y = 3.375, hjust = 0, size = 7, family = "Times New Roman", 
           label = "Greater under passive enforcement (Talek region)") + 
  annotate("text", x = 0.5, y = -3.375, hjust = 0, size = 7, family = "Times New Roman", 
           label = "Greater under active enforcement (Mara Triangle)") + 
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  geom_vline(xintercept = (which(values$species == "The Community") - 0.5), linetype = "dotted") +
  scale_color_manual(name = "", values = c("black" = "black", "orange" = "#ff8000", "orangegrey" = "#CC6600", "green" = "#33ff77", "greengrey" = "#009933"),
                     labels = c("No difference", "Higher\nin Talek region (95%)", "Higher\nin Talek region (50%)", "Lower\nin Talek region (95%)", "Lower\nin Talek region (50%)")) +
  theme_few() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(family = "Times New Roman", size = 24),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(labels = c("African Lion" = "African\nlion", "Banded Mongoose" = "Banded\nmongoose", "Bat-eared Fox" = "Bat-eared\nfox", "Black-backed Jackal" = "Black-backed\njackal",
                              "Side-striped Jackal" = "Side-striped\njackal", "Slender Mongoose" = "Slender\nmongoose", "Spotted Hyena" = "Spotted\nhyena",
                              "The Community" = "The\nCommunity")) +
  labs(y = expression("Management effect"), x = expression())

ggsave(file = "Figure2.png", bg = "transparent", width = 12, height = 6)

#----------#
#-Figure 3-#
#----------#

beta1.val <- cbind(c(HMSDS$q2.5$beta1[c(8,1:3,5,11,6)], HMSDS$q2.5$mu_b1),
                    c(HMSDS$q25$beta1[c(8,1:3,5,11,6)], HMSDS$q25$mu_b1),
                    c(HMSDS$mean$beta1[c(8,1:3,5,11,6)], HMSDS$mean$mu_b1),
                    c(HMSDS$q75$beta1[c(8,1:3,5,11,6)], HMSDS$q75$mu_b1),
                    c(HMSDS$q97.5$beta1[c(8,1:3,5,11,6)], HMSDS$q97.5$mu_b1),
                    c(HMSDS$overlap0$beta1[c(8,1:3,5,11,6)],HMSDS$overlap0$mu_b1))

sppnames <- c("African Lion", "Banded Mongoose", "Bat-eared Fox", "Black-backed Jackal",
              "Cheetah", "Slender Mongoose", "Spotted Hyena", "The Community")

values <- data.frame(sppnames, beta1.val)
colnames(values) <- c("species", "lower.beta", "l25.beta", "mean.beta", "u75.beta", "upper.beta", "beta.sig")

values$species <- factor(values$species, levels = values$species)

values$nodiff <- values$l25.beta < 0 & values$u75.beta > 0
values$H95 <- values$lower.beta > 0
values$H50 <- values$lower.beta < 0 & values$l25.beta > 0
values$L95 <- values$upper.beta < 0
values$L50 <- values$upper.beta > 0 & values$u75.beta < 0

Figure3 <- ggplot() + 
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = mean.beta, ymax = mean.beta, color = "black"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = lower.beta, ymax = upper.beta, color = "black"), 
                width = 0, size = 1.25) +
  geom_errorbar(data = subset(values, nodiff == TRUE), aes(x = species, ymin = l25.beta, ymax = u75.beta, color = "black"), 
                width = 0, size = 3.5) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = mean.beta, ymax = mean.beta, color = "green"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = lower.beta, ymax = upper.beta, color = "green"), 
                size = 1.25, width = 0) +
  geom_errorbar(data = subset(values, L95 == TRUE), aes(x = species, ymin = l25.beta, ymax = u75.beta, color = "green"),
                size = 3.5, width = 0) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = mean.beta, ymax = mean.beta, color = "greengrey"), 
                width = 0.25) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = lower.beta, ymax = upper.beta, color = "greengrey"),
                size = 1.25, width = 0) +
  geom_errorbar(data = subset(values, L50 == TRUE), aes(x = species, ymin = l25.beta, ymax = u75.beta, color = "greengrey"),
                size = 3.5, width = 0) +
  annotate("text", x = 0.5, y = 3.125, hjust = 0, size = 7, family = "Times New Roman", 
           label = "Greater under passive enforcement (Talek region)") + 
  annotate("text", x = 0.5, y = -3.125, hjust = 0, size = 7, family = "Times New Roman", 
           label = "Greater under active enforcement (Mara Triangle)") + 
  coord_cartesian(ylim = c(-3.25, 3.25)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  geom_vline(xintercept = (which(values$species == "The Community") - 0.5), linetype = "dotted") +
  scale_color_manual(name = "", values = c("black" = "black", "orange" = "#ff8000", "orangegrey" = "#cc6600", "green" = "#33ff77", "greengrey" = "#009933"),
                     labels = c("No difference", "Higher\nin Talek region (95%)", "Higher\nin Talek region (50%)", "Lower\nin Talek region (95%)", "Lower\nin Talek region (50%)")) +
  theme_few() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(family = "Times New Roman", size = 24),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(labels = c("African Lion" = "African\nlion", "Banded Mongoose" = "Banded\nmongoose", "Bat-eared Fox" = "Bat-eared\nfox", "Black-backed Jackal" = "Black-backed\njackal", 
                              "Side-striped Jackal" = "Side-striped\njackal", "Slender Mongoose" = "Slender\nmongoose", "Spotted Hyena" = "Spotted\nhyena", 
                              "The Community" = "The\nCommunity")) + 
  labs(y = expression("Management effect"), x = expression())

ggsave(file = "Figure3.png", bg = "transparent", width = 12, height = 6)

#----------#
#-Figure 4-#
#----------#

data <- NULL
Denplot <- NULL
regnames <- c("Talek region", "Mara Triangle")
plotter <- c(8,1:3,6)
sppnames <- c("African\nlion", "Banded\nmongoose", "Bat-eared\nfox", "Black-backed\njackal", "Spotted\nhyena")
for(i in 1:5){
  den.val <- cbind(c(HMSDS$mean$RegGS[plotter[i], 2],  HMSDS$mean$RegGS[plotter[i], 1]), 
                             c(HMSDS$q2.5$RegGS[plotter[i], 2],  HMSDS$q2.5$RegGS[plotter[i], 1]),
                             c(HMSDS$q25$RegGS[plotter[i], 2],  HMSDS$q25$RegGS[plotter[i], 1]),
                             c(HMSDS$q75$RegGS[plotter[i], 2],  HMSDS$q75$RegGS[plotter[i], 1]),
                             c(HMSDS$q97.5$RegGS[plotter[i], 2],  HMSDS$q97.5$RegGS[plotter[i], 1]))
  colornames <- c()
  data <- data.frame(den.val, regnames)
  colnames(data) <- c("mean", "l2.5", "l25", "u75", "u97.5", "Region")
  Denplot[[i]] <- ggplot(data, aes(x = Region, y = mean, color = Region)) +
                         geom_point(aes(color = Region), size = 0.0001) +
                         geom_errorbar(aes(x = Region, ymax = mean, ymin = mean), 
                                      width = 0.25) +
                         geom_errorbar(aes(x = Region, ymax = u97.5, ymin = l2.5), 
                                       width = 0, size = 0.5) +
                         geom_errorbar(aes(x = Region, ymax = u75, ymin = l25), 
                                       width = 0, size = 1.5) +
                         scale_color_manual(values = c("#009933", "#CC6600")) +
                         scale_y_continuous(breaks = scales::pretty_breaks(4)) +
                         theme_few() +
                         theme(text = element_text(family = "Times New Roman", size = 12),
                               panel.background = element_rect(fill = "transparent", color = NA),
                               plot.background = element_rect(fill = "transparent", color = NA),
                               panel.border = element_blank(),
                               axis.text.x = element_blank(),
                               axis.ticks.x=element_blank(),
                               legend.position = "none", 
                               plot.margin = unit(c(0.25,0,0.125,0.25), "in"),
                               axis.line = element_line(color = "black")) + 
                         labs(y = NULL, x = sppnames[i], color = "", title = NULL) +
                         coord_cartesian(ylim = c(0, max(data$u97.5) + 1))
  rm(data)
  Denplot[[i]]
}

Figure4 <- grid.arrange(textGrob(expression(Density~at~13~km^2), gp=gpar(fontfamily = "Times New Roman", fontsize = 14), rot = 90),
                        arrangeGrob(grobs = Denplot,  nrow = 1), ncol = 2, widths = c(0.05, 0.95))
ggsave(file = "Figure4.png", plot = Figure4, width = 6.5, height = 3)

#----------#
#-Figure 5-#
#----------#

#Simulate distance across transect half-width
dist.sim <- seq(0, 650, 0.1)

#Harvest mean sigma value
sigma <- c(exp(HMSDS$mean$asig + as.numeric(DSdata$abs %*% HMSDS$mean$bsig1)), exp(HMSDS$mean$mu_s))

#Calculate detection probability across distances
distfunc <- matrix(NA, nrow = 12, ncol = length(dist.sim))
for(i in 1:12){
  for(j in 1:length(dist.sim)){
    distfunc[i,j] <- exp(-dist.sim[j]*dist.sim[j]/(2*sigma[i]*sigma[i]))
  }
}

#Plot detection probability vs observed distance
Figure5 <- ggplot() +  
  geom_line(aes(x = dist.sim, y = distfunc[1,], color = "Banded mongoose"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[2,], color = "Bat-eared fox"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[3,], color = "Black-backed jackal"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[4,], color = "Caracal"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[5,], color = "Cheetah"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[6,], color = "Spotted hyena"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[7,], color = "Leopard"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[8,], color = "African lion"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[9,], color = "Serval"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[10,], color = "Side-striped jackal"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[11,], color = "Slender mongoose"), size = 1.5) +
  geom_line(aes(x = dist.sim, y = distfunc[12,], color = "The Community"), size = 1.5) + 
  scale_color_manual(name = "", values = c("Banded mongoose" = "#d9d9d9", 
                                           "Bat-eared fox" = "#fccde5", 
                                           "Black-backed jackal" = "#8dd3c7", 
                                           "Caracal" = "#fdb462", 
                                           "Cheetah" = "#fb8072", 
                                           "Spotted hyena" = "#bebada",
                                           "Leopard" = "#bc80bd", 
                                           "African lion" = "#80b1d3", 
                                           "Serval" = "#b3de69", 
                                           "Side-striped jackal" = "#ffffb3", 
                                           "Slender mongoose" = "#ccebc5",
                                           "The Community" = "black")) +
  labs(
    title = "",
    x = "Distance (meters)",
    y = "Detection probability\n") + 
  theme_few() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "in"),
        text = element_text(family = "Times New Roman", size = 24))

print(Figure5)
ggsave(file = "Figure5.png", plot = Figure5, width = 12, height = 6)