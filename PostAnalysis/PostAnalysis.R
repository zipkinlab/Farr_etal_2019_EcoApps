g <- array(NA, dim = c(nD, 30000))
g2 <- NULL

  for(i in 1:30000){
  for(k in 1:nD){
  g[k,i] <- exp(-mdpt[k]*mdpt[k]/(2*exp(CDSM$sims.list$mu_s[i])*exp(CDSM$sims.list$mu_s[i]))) * v/B
  }
  g2[i] <- sum(g[,i])
  }

coef.vect_g2 <- mean(g2)
lower95_g2 <- quantile(g2, probs = c(0.025))
upper95_g2 <- quantile(g2, probs = c(0.975))
##########################################################################################################
### plot parameter values and CI
alpha1.val <- cbind(CDSM$q2.5$alpha1, CDSM$q25$alpha1, CDSM$mean$alpha1, 
                    CDSM$q75$alpha1, CDSM$q97.5$alpha1, CDSM$overlap0$alpha1)

beta1.val <- cbind(CDSM$q2.5$beta1, CDSM$q25$beta1, CDSM$mean$beta1, 
                   CDSM$q75$beta1, CDSM$q97.5$beta1, CDSM$overlap0$beta1)

sppnames <- c("Banded Mongoose", "Bat-eared Fox", "Black-backed Jackal", "Caracal",
              "Cheetah", "Hyena", "Leopard", "Lion", "Serval", "Side-striped Jackal", 
              "Slender Mongoose")

values <- data.frame(sppnames, alpha1.val, beta1.val)
colnames(values) <- c("species", "lower.alpha", "l25.alpha", "mean.alpha", "u75.alpha", "upper.alpha", "alpha.sig",
                      "lower.beta", "l25.beta", "mean.beta", "u75.beta", "upper.beta", "beta.sig")

values$species <- factor(values$species, levels = values$species)


library(ggplot2)
library(ggthemes)


###Alpha 1

Alpha1plot <- ggplot() + geom_point(data = subset(values, alpha.sig == 1), aes(x = mean.alpha, y = species), size = 3, fill = "white")+ 
  geom_errorbarh(data = subset(values, alpha.sig == 1), aes(x = mean.alpha, y = species, xmin = lower.alpha, xmax = upper.alpha), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = subset(values, alpha.sig == 1), aes(x = mean.alpha, y = species, xmin = l25.alpha, xmax = u75.alpha), 
                 size = 1.5, height = 0, alpha = 0.8) +
  geom_point(data = subset(values, alpha.sig == 0), aes(y = species, x = mean.alpha), color = "red", size = 3, fill = "white") +
  geom_errorbarh(data = subset(values, alpha.sig == 0), aes(x = mean.alpha, y = species, xmin = lower.alpha, xmax = upper.alpha), 
                 color = "red", size = 0.3, height = 0) +
  geom_errorbarh(data = subset(values, alpha.sig == 0), aes(x = mean.alpha, y = species, xmin = l25.alpha, xmax = u75.alpha), 
                 color = "red", size = 1.5, height = 0) +
  coord_cartesian(xlim = c(-4, 4)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  theme_few() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = "Mara Triangle                                 Talek Region")

###Beta 1

Beta1plot <- ggplot() + geom_point(data = subset(values, beta.sig == 1), aes(x = mean.beta, y = species), size = 3, fill = "white")+ 
  geom_errorbarh(data = subset(values, beta.sig == 1), aes(x = mean.beta, y = species, xmin = lower.beta, xmax = upper.beta), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = subset(values, beta.sig == 1), aes(x = mean.beta, y = species, xmin = l25.beta, xmax = u75.beta), 
                 size = 1.5, height = 0, alpha = 0.8) +
  geom_point(data = subset(values, beta.sig == 0), aes(y = species, x = mean.beta), color = "red", size = 3, fill = "white") +
  #geom_errorbarh(data = subset(values, beta.sig == 0), aes(x = mean.beta, y = species, xmin = lower.beta, xmax = upper.beta), 
  # color = "red", size = 0.3, height = 0) +
  #geom_errorbarh(data = subset(values, beta.sig == 0), aes(x = mean.beta, y = species, xmin = l25.beta, xmax = u75.beta), 
  #color = "red", size = 1.5, height = 0) +
  coord_cartesian(xlim = c(-4, 4)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  theme_few() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = "Mara Triangle                                 Talek Region")

###Density difference

den.val <- cbind(CDSM$q2.5$Ddiff, CDSM$q25$Ddiff, CDSM$mean$Ddiff, 
                 CDSM$q75$Ddiff, CDSM$q97.5$Ddiff, CDSM$overlap0$Ddiff)

sppnames <- c("Banded Mongoose", "Bat-eared Fox", "Black-backed Jackal", 
              "Cheetah", "Hyena", "Lion", "Serval", "Side-striped Jackal", 
              "Slender Mongoose")

Dvalues <- data.frame(sppnames, den.val)
colnames(Dvalues) <- c("species", "lower.den", "l25.den", "mean.den", "u75.den", "upper.den", "den.sig")

Dvalues$species <- factor(Dvalues$species, levels = Dvalues$species)

Ddiffplot <- ggplot() + geom_point(data = subset(Dvalues, den.sig == 1), aes(x = mean.den, y = species), size = 3, fill = "white")+ 
  geom_errorbarh(data = subset(Dvalues, den.sig == 1), aes(x = mean.den, y = species, xmin = lower.den, xmax = upper.den), 
                 size = 0.3, height = 0, alpha = 0.8) +
  geom_errorbarh(data = subset(Dvalues, den.sig == 1), aes(x = mean.den, y = species, xmin = l25.den, xmax = u75.den), 
                 size = 1.5, height = 0, alpha = 0.8) +
  geom_point(data = subset(Dvalues, den.sig == 0), aes(y = species, x = mean.den), color = "red", size = 3, fill = "white") +
  geom_errorbarh(data = subset(Dvalues, den.sig == 0), aes(x = mean.den, y = species, xmin = lower.den, xmax = upper.den), 
                 color = "red", size = 0.3, height = 0) +
  geom_errorbarh(data = subset(Dvalues, den.sig == 0), aes(x = mean.den, y = species, xmin = l25.den, xmax = u75.den), 
                 color = "red", size = 1.5, height = 0) +
  coord_cartesian(xlim = c(-1.25, 1.25)) +
  geom_vline(xintercept = 0, alpha = 0.3) +
  theme_few() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  labs(x = "Mara Triangle                                 Talek Region")
