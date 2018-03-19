# Marina Golivets
# March 2018
# Code to create figures for the competition-relatedness hypothesis test


# load packages
library(metafor)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)
library(bayesplot)
library(grid)

par(mfrow = c(2, 2), oma=c(3.5,0,0,0), mar = c(3, 3, 1, 2), mgp = c(2, .4, .4), las = 1, tck = -0.03)


load("compet_related_hypoth/additive.PairMean.CompRelated.RData")
load("compet_related_hypoth/additive.PairDifs.CompRelated.RData")


load("compet_related_hypoth/replacement.PairMean.CompRelated.RData")
load("compet_related_hypoth/replacement.PairDifs.CompRelated.RData")

dat <- read.csv("compet_related_hypoth/dat.ad.PairMean.CompRelated.csv")
dat <- read.csv("compet_related_hypoth/dat.ad.PairDifs.CompRelated.csv")

dat <- read.csv("compet_related_hypoth/dat.rp.PairMean.CompRelated.csv")
dat <- read.csv("compet_related_hypoth/dat.rp.PairDifs.CompRelated.csv")


X0 <- c(seq(0,2, 0.001), seq(2, 400, 1))
Y1 <- Y2 <- matrix(NA, ncol = length(rstan::extract(m2)$a[,1]), nrow = length(X0))
a1 <- rstan::extract(m2)$a[ ,1]
b1 <- rstan::extract(m2)$beta1[ ,1]
a2 <- rstan::extract(m2)$a[ ,2]
b2 <- rstan::extract(m2)$beta1[ ,2]
for(i in 1:length(X0)){
    Y1[i, ] <- a1 * (1 - exp(- b1* X0[i]))
    Y2[i, ] <- a2 * (1 - exp(-b2 * X0[i]))
}

Y1_CI95 <- apply(Y1, 1, quantile, c(0.025, 0.975))
Y2_CI95 <- apply(Y2, 1, quantile, c(0.025, 0.975))

d_pred <- apply(rstan::extract(m2)$d_pred, 2, mean)



my_colors <- rep(NA, nrow(dat))
for (i in 1:nrow(dat)) {my_colors[i] <- c("#a6611a", "#018571")[as.numeric(as.factor(dat$monocot1))[i]] }


plot(apply(Y2, 1, mean) ~ X0, type = "l", col = "#a6611a", lwd = 2,  ylim = c(0, 15),
     frame.plot = FALSE, xlab = "Phylogenetic distance", ylab = "Mean ANNE", cex.axis = .8, 
     cex.label = .9)
lines(apply(Y1, 1, mean) ~ X0, col = "#018571", lwd = 2)
lines(Y2_CI95[1, ] ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(Y2_CI95[2, ] ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(Y1_CI95[1, ] ~ X0, col = "#018571", lty = 2, lwd = 1)
lines(Y1_CI95[2, ] ~ X0, col = "#018571", lty = 2, lwd = 1)
points((dat$yi - min(dat$yi)) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points((dat$yi - min(dat$yi)) ~ dat$pair.dist, col = my_colors)
points(d_pred ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(d_pred ~ dat$pair.dist, col = my_colors, pch = 16)
mtext("(a)", side = 1, line = -9.5, cex = .75, at = -80, font = 2)


plot(apply(Y2, 1, mean) ~ X0, type = "l", col = "#a6611a", lwd = 2,  ylim = c(0, 30),
     frame.plot = FALSE, xlab = "Phylogenetic distance", ylab = "Mean RNNE", cex.axis = .8, 
     cex.label = .9)
lines(apply(Y1, 1, mean) ~ X0, col = "#018571", lwd = 2)
lines(Y2_CI95[1, ] ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(Y2_CI95[2, ] ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(Y1_CI95[1, ] ~ X0, col = "#018571", lty = 2, lwd = 1)
lines(Y1_CI95[2, ] ~ X0, col = "#018571", lty = 2, lwd = 1)
points((dat$yi - min(dat$yi)) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points((dat$yi - min(dat$yi)) ~ dat$pair.dist, col = my_colors)
points(d_pred ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(d_pred ~ dat$pair.dist, col = my_colors, pch = 16)
mtext("(b)", side = 1, line = -9.5, cex = .75, at = -80, font = 2)


plot(log(apply(Y2, 1, mean)) ~ X0, type = "l", col = "#a6611a", ylim = c(-7, 5),
     frame.plot = FALSE, cex.axis = .8, 
     cex.label = .9,
     lwd = 2, ylab = "log (ANNE asymmetry)", xlab = "Phylogenetic distance" ) 
lines(log(apply(Y1, 1, mean)) ~ X0, col = "#018571", lwd = 2)
lines(log(Y2_CI95[1, ]) ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(log(Y2_CI95[2, ]) ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(log(Y1_CI95[1, ]) ~ X0, col = "#018571", lty = 2, lwd = 1)
lines(log(Y1_CI95[2, ]) ~ X0, col = "#018571", lty = 2, lwd = 1)
points(log(dat$yi) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(log(dat$yi) ~ dat$pair.dist, col = my_colors)
points(log(d_pred) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(log(d_pred) ~ dat$pair.dist, col = my_colors, pch = 16)
mtext("(c)", side = 1, line = -9.5, cex = .75, at = -80, font = 2)


plot(log(apply(Y2, 1, mean)) ~ X0, type = "l", col = "#a6611a", ylim = c(-7, 5),
     frame.plot = FALSE, cex.axis = .8, 
     cex.label = .9,
     lwd = 2, ylab = "log (RNNE asymmetry)", xlab = "Phylogenetic distance" ) 
lines(log(apply(Y1, 1, mean)) ~ X0, col = "#018571", lwd = 2)
lines(log(Y2_CI95[1, ]) ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(log(Y2_CI95[2, ]) ~ X0, col = "#a6611a", lty = 2, lwd = 1)
lines(log(Y1_CI95[1, ]) ~ X0, col = "#018571", lty = 2, lwd = 1)
lines(log(Y1_CI95[2, ]) ~ X0, col = "#018571", lty = 2, lwd = 1)
points(log(dat$yi) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(log(dat$yi) ~ dat$pair.dist, col = my_colors)
points(log(d_pred) ~ dat$pair.dist, col = "white", pch = 16, cex = 1.2)
points(log(d_pred) ~ dat$pair.dist, col = my_colors, pch = 16)
mtext("(d)", side = 1, line = -9.5, cex = .75, at = -60, fon = 2)

legend(x=-600, y=-12,legend=c("eudicots", "monocots"), col = c("#a6611a", "#018571"), pch = 15, cex = .8, bty = "n", xpd = NA, ncol = 2)
legend(x=-230, y=-12,legend=c("observed value", "predicted value"), ncol=2, pch = c(1, 16), cex = .8, bty = "n", xpd = NA)
       
legend(x=-500, y=-14,legend=c("fitted non-linear curve"), cex = 0.8,
       lwd = 2, lty = 1, bty = "n", xpd = NA)
legend(x=-110, y=-14,legend=c("95% credible interval"), cex = 0.8,
       lty = 2, lwd = 1, bty = "n", xpd = NA)


