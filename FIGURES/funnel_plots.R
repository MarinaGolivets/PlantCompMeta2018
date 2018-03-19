# Marina Golivets
# March 2018
# Code to create funnel plots


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


load("additive.subs20.Final.20000iter.RData")
Fin_ad <- m2

load("replacement.subs20.Final.run4.20000iter.RData")
Fin_rp <- m2

load("additive.subs20.Pairs.Final.20000iter.RData")
Fin_ad_pair <- m2

load("replacement.subs20.Pairs.Final.run4.20000iter.RData")
Fin_rp_pair <- m2

# load data
dat <- read.csv(file = "MGolivets_plant_biomass_additive_subs20_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_ad <- dat

dat <- read.csv(file = "MGolivets_plant_biomass_replacement_subs20_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_rp <- dat


dat <- read.csv(file = "MGolivets_plant_biomass_additive_subs20_Pairs_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_ad_pair <- dat

dat <- read.csv(file = "DATA/MGolivets_plant_biomass_replacement_subs20_Pairs_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_rp_pair <- dat



# publication bias

res_ad <- apply(rstan::extract(Fin_ad)$res, 2, mean)
res_ad <- res_ad / sd(res_ad)
shapiro.test(res_ad)
qqnorm(res_ad); qqline(res_ad, col = 2)
s_ad <- 1/sqrt(dat_ad$vi)
y_ad <- res_ad*s_ad
summary(lm(y_ad ~ s_ad))
col_ad <- as.numeric(as.factor(dat_ad$group))

table(res_ad <= -1.3)
incl <- which(res_ad > -1.3)
summary(lm(y_ad[incl] ~ s_ad[incl]))

res_rp <- apply(rstan::extract(Fin_rp)$res, 2, mean)
res_rp <- res_rp / sd(res_rp)
shapiro.test(res_rp)
s_rp <- 1/sqrt(dat_rp$vi)
y_rp <- res_rp*s_rp
summary(lm(y_rp ~ s_rp))
col_rp <- as.numeric(as.factor(dat_rp$group))

table(res_rp <= -2.2)
incl <- which(res_rp > -2.2)
summary(lm(y_rp[incl] ~ s_rp[incl]))

res_ad_pair <- apply(rstan::extract(Fin_ad_pair)$res, 2, mean)
res_ad_pair <- res_ad_pair/ sd(res_ad_pair)
shapiro.test(res_ad_pair)
s_ad_pair <- 1/sqrt(dat_ad_pair$vi)
y_ad_pair <- res_ad_pair*s_ad_pair
summary(lm(y_ad_pair ~ s_ad_pair))
col_ad_apir <- as.numeric(as.factor(dat_ad_pair$group))

res_rp_pair <- apply(rstan::extract(Fin_rp_pair)$res, 2, mean)
res_rp_pair <- res_rp_pair / sd(res_rp_pair)
shapiro.test(res_rp_pair)
s_rp_pair <- 1/sqrt(dat_rp_pair$vi)
y_rp_pair <- res_rp_pair*s_rp_pair
summary(lm(y_rp_pair ~ s_rp_pair))
col_rp_pair <- as.numeric(as.factor(dat_rp_pair$group))


# ggplot(as.data.frame(cbind(y_egger, error_egger)), aes(error_egger, y_egger)) + 
#     geom_point() + 
#     geom_abline(aes(slope = 1, intercept = 0)) + 
#     theme_bw()

library(meta)

par(mfrow = c(2,2), mar = c(3.2, 5, 3.2, 0), mgp = c(1.5, .3, .3), las = 1, tck = -0.02)

funnel.plot.function <- function (x = res_ad, y = dat_ad$vi, plot.num = "(a)", col = col_ad,
                                  at.plot.num = -12, b0 = -0.50, at.b0 = -7.3, xlim = c(-10, 4))
                                   {
    funnel(x = x, y = y, 
           contour.levels = c(.05, .5, .95),
           col.contour = c("#80cdc1", "darkgrey", "white"),
           comb.random = T,
           yaxis = "invse", sm = "SMD", 
           pch = ".", cex = 2, col = "#a6611a",
           cex.axis = .8, xlim = xlim, ylim = c(0, 50),
           xlab = "", ylab = "", font.lab = 2, cex.lab = .9,
           frame.plot = FALSE)
    mtext(plot.num, side = 1, line = -9, cex = .8, at = at.plot.num, font = 2, adj = 1)
    text(bquote(paste(italic(b[0]), " = ", .(b0), ", ", italic(P), " < 0.001")), 
         x = at.b0, y = 41, cex = 0.85)
}

funnel.plot.function()
funnel.plot.function(x = res_rp, y = dat_rp$vi, plot.num = "(b)", 
                     at.plot.num = -13, b0 = -0.33, at.b0 = -6.1, xlim = c(-10, 10))
funnel.plot.function(x = res_ad_pair, y = dat_ad_pair$vi, plot.num = "(c)", 
                     b0 = -0.63) 
funnel.plot.function(x = res_rp_pair, y = dat_rp_pair$vi, plot.num = "(d)", 
                     at.plot.num = -13, b0 = -0.33, at.b0 = -6.1, xlim = c(-10, 10))
mtext("Meta-analytic residuals", side = 3, line = -10.5, cex = .8, at = -15, font = 2)
mtext("Reciprocal effects", side = 1, line = -11.5, cex = .8, at = -16, font = 2, padj = 1)
mtext("All effects", side = 1, line = -25.5, cex = .8, at = -17, font = 2, padj = 1)
mtext("Precision", side = 1, line = -14, cex = .8, at = -44, font = 2, padj = 1, las = 2)

cor(res_rp, dat_rp$yi)
hist(res_rp)

