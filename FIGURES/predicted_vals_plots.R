# Marina Golivets
# March 2018
# Code to plots predicted vs observed effect sizes


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

# load Stan outputs
load("additive.subs20.Final.20000iter.RData")
Fin_ad <- m2

load("replacement.subs20.Final.run4.20000iter.RData")
Fin_rp <- m2

load("additive.subs20.Pairs.Final.20000iter.RData")
Fin_ad_pair <- m2

load("replacement.subs20.Pairs.Final.run4.20000iter.RData")
Fin_rp_pair <- m2


# load raw data
dat_ad <- read.csv(file = "MGolivets_plant_biomass_additive_subs20_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)

dat_rp <- read.csv(file = "MGolivets_plant_biomass_replacement_subs20_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)

dat_ad_pair <- read.csv(file = "MGolivets_plant_biomass_additive_subs20_Pairs_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)

dat_rp_pair <- read.csv(file = "MGolivets_plant_biomass_replacement_subs20_Pairs_12152017.csv",
                header = TRUE, stringsAsFactors = FALSE)


# create a data frames with observed and predicted effect size values
df_pred_ad <- data.frame(x = 1:nrow(dat_ad),
                         d = dat_ad$yi, 
                         d_pred = apply(rstan::extract(Fin_ad)$d_pred, 2, mean), 
                         d_pred_CrI25 = apply(rstan::extract(Fin_ad)$d_pred, 2, quantile, .25),
                         d_pred_CrI75 = apply(rstan::extract(Fin_ad)$d_pred, 2, quantile, .75))


df_pred_rp <- data.frame(x = 1:nrow(dat_rp),
                         d = dat_rp$yi, 
                         d_pred = apply(rstan::extract(Fin_rp)$d_pred, 2, mean), 
                         d_pred_CrI25 = apply(rstan::extract(Fin_rp)$d_pred, 2, quantile, .25),
                         d_pred_CrI75 = apply(rstan::extract(Fin_rp)$d_pred, 2, quantile, .75))


df_pred_ad_pair <- data.frame(x = 1:nrow(dat_ad_pair),
                              d = dat_ad_pair$yi, 
                              d_pred = apply(rstan::extract(Fin_ad_pair)$d_pred, 2, mean), 
                              d_pred_CrI25 = apply(rstan::extract(Fin_ad_pair)$d_pred, 2, quantile, .25),
                              d_pred_CrI75 = apply(rstan::extract(Fin_ad_pair)$d_pred, 2, quantile, .75))


df_pred_rp_pair <- data.frame(x = 1:nrow(dat_rp_pair),
                              d = dat_rp_pair$yi, 
                              d_pred = apply(rstan::extract(Fin_rp_pair)$d_pred, 2, mean), 
                              d_pred_CrI25 = apply(rstan::extract(Fin_rp_pair)$d_pred, 2, quantile, .25),
                              d_pred_CrI75 = apply(rstan::extract(Fin_rp_pair)$d_pred, 2, quantile, .75))


# df_pred <- cbind(rbind(df_pred_ad, df_pred_rp), dataset = c(rep("ad", nrow(dat_ad)), rep("rp", nrow(dat_rp))))

# calculate proportions of times when observed values are within 50% posterior predictive intervals

count.pred.function <- function(data){
    counts <- rep(NA, nrow(data))
    for(i in 1:nrow(data)){
        counts[i] <- ifelse(data$d_pred_CrI25[i] <= data$d[i] & 
                                  data$d_pred_CrI75[i] >= data$d[i], 1, 0)
    } 
    counts <- round(sum(counts)/nrow(data), 2) 
    return(counts)
}

count.pred.function(df_pred_ad)
count.pred.function(df_pred_ad_pair)

count.pred.function(df_pred_rp)
count.pred.function(df_pred_rp_pair)


# plot observed data along with 50% posterior predictive intervals

library(grid)

plot.prediction.function <- function (data = df_pred_ad, y_lab = "ANNE", plot.num = "(a)") {
    p <-  ggplot(data, aes(x = x, y = d)) +
        geom_errorbar(aes(ymin = d_pred_CrI25,
                          ymax = d_pred_CrI75, color = "col1", size = "size1")) +
        geom_point(aes(color = "col2", shape = "shape1"), size = .7) +
        scale_color_manual(name= "", values = c(col1 = "#80cdc1", col2 = "#a6611a")) +
        scale_shape_manual(name= "", values = c(shape1=3), labels = c("Observed value")) + 
        scale_fill_manual(name= "", values = c(col2 = "#a6611a"), labels = c("Observed value")) + 
        scale_size_manual(name= "", values = c(size1 = .5), labels = c("50% PPI")) +
        
        guides(color = FALSE,  
               shape = guide_legend(order = 1, override.aes = list(color = c("#a6611a"))),
               fill = guide_legend(order = 1),  
               size = guide_legend(order = 2, override.aes = list(color = c("#80cdc1")))) +
    
        theme(strip.text.x = element_blank(), 
              axis.title.x = element_text(family = "Calibri", size = 9, face = "plain"),
              axis.title.y = element_text(family = "Calibri", size = 9, face = "plain"),
              axis.text.x = element_text(family = "Calibri", size = 8, face = "plain", color = "black"),
              axis.text.y = element_text(family = "Calibri", size = 8, face = "plain", color = "black"),
              legend.title = element_text(family = "Calibri", color = "black", size = 9, face="bold"),
              legend.text =  element_text(family = "Calibri", color = "black", size = 9, face="plain"),                           
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              legend.position="bottom",
              legend.box = "horizontal",
              axis.ticks.length = unit(.1, "cm"),
              plot.margin = unit(c(2, 1, 1, 1), "lines")) +
        
        labs(y = y_lab, x = "Observation no.")

    return(p)
}
plot.prediction.function()

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend) }

mylegend <- g_legend(plot.prediction.function())

grid.newpage()
p3 <- grid.arrange(arrangeGrob(plot.prediction.function() + theme(legend.position="none"),
                               plot.prediction.function(data = df_pred_rp, y_lab = "RNNE") + theme(legend.position="none"),
                               plot.prediction.function(data = df_pred_ad_pair) + theme(legend.position="none"),
                               plot.prediction.function(data = df_pred_rp_pair, y_lab = "RNNE") + theme(legend.position="none"),
                               nrow=2),
                   mylegend, nrow=3, heights=c(10, 1, 0))




grid.text(x=0.05, y=.93, label= "(a)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.55, y=.93, label= "(b)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.05, y=.48, label= "(c)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.55, y=.48, label= "(d)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.99, label= "All effects", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.54, label= "Reciprocal effects", gp = gpar(cex = .8, fontface="bold"))
