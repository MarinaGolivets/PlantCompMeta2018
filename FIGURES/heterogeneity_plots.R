# Marina Golivets
# March 2018
# Code to create a plot showing heterogeneity partitioning



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


mean.sd.fn <- function (ad = MEAN_ad, rp = MEAN_rp, dataset.name = "all data") {
    
    # additive
    I2_sdep <- mean(rstan::extract(ad)$I2_sdep)
    I2_study <- mean(rstan::extract(ad)$I2_study)
    I2_nb <- mean(rstan::extract(ad)$I2_nb)
    I2_tg <- mean(rstan::extract(ad)$I2_tg)
    I2_obs <- mean(rstan::extract(ad)$I2_obs)
    
    sd_I2_sdep <- sd(rstan::extract(ad)$I2_sdep)
    sd_I2_study <- sd(rstan::extract(ad)$I2_study)
    sd_I2_nb <- sd(rstan::extract(ad)$I2_nb)
    sd_I2_tg <- sd(rstan::extract(ad)$I2_tg)
    sd_I2_obs <- sd(rstan::extract(ad)$I2_obs)
    
    
    I2.ad <- cbind(c(I2_study, I2_sdep, I2_tg, I2_nb, I2_obs),
                   c(sd_I2_study, sd_I2_sdep, sd_I2_tg, sd_I2_nb, sd_I2_obs))
    
    
    # replacement
    I2_sdep <- mean(rstan::extract(rp)$I2_sdep)
    I2_study <- mean(rstan::extract(rp)$I2_study)
    I2_nb <- mean(rstan::extract(rp)$I2_nb)
    I2_tg <- mean(rstan::extract(rp)$I2_tg)
    I2_obs <- mean(rstan::extract(rp)$I2_obs)
    
    
    sd_I2_sdep <- sd(rstan::extract(rp)$I2_sdep)
    sd_I2_study <- sd(rstan::extract(rp)$I2_study)
    sd_I2_nb <- sd(rstan::extract(rp)$I2_nb)
    sd_I2_tg <- sd(rstan::extract(rp)$I2_tg)
    sd_I2_obs <- sd(rstan::extract(rp)$I2_obs)
    
    
    I2.rp <- cbind(c(I2_study, I2_sdep, I2_tg, I2_nb, I2_obs),
                   c(sd_I2_study, sd_I2_sdep, sd_I2_tg, sd_I2_nb, sd_I2_obs))
    
    
    df <- data.frame(Model = rep(c("ANNEs", "RNNEs"), each = 5), 
                     effect_name = rep(c("study", "sampling group", "target species", "neighbor species", "observation"), 2), 
                     effect_ratio = c(I2.ad[ ,1], I2.rp[ ,1]), 
                     effect_error = c(I2.ad[ ,2], I2.rp[ ,2]), dataset = rep(dataset.name, 10))
    
    return (df)
}    

mean.sd.phylo.fn <- function (ad = MEAN_ad_phylo, rp = MEAN_rp_phylo, dataset.name = "all data") {
    
    # additive
    I2_sdep <- mean(rstan::extract(ad)$I2_sdep)
    I2_study <- mean(rstan::extract(ad)$I2_study)
    I2_nb <- mean(rstan::extract(ad)$I2_nb)
    I2_tg <- mean(rstan::extract(ad)$I2_tg)
    I2_obs <- mean(rstan::extract(ad)$I2_obs)
    I2_nb_phylo <- mean(rstan::extract(ad)$I2_nb_phylo)
    I2_tg_phylo <- mean(rstan::extract(ad)$I2_tg_phylo)
    
    sd_I2_sdep <- sd(rstan::extract(ad)$I2_sdep)
    sd_I2_study <- sd(rstan::extract(ad)$I2_study)
    sd_I2_nb <- sd(rstan::extract(ad)$I2_nb)
    sd_I2_tg <- sd(rstan::extract(ad)$I2_tg)
    sd_I2_obs <- sd(rstan::extract(ad)$I2_obs)
    sd_I2_nb_phylo <- sd(rstan::extract(ad)$I2_nb_phylo)
    sd_I2_tg_phylo <- sd(rstan::extract(ad)$I2_tg_phylo)
    
    
    I2.ad <- cbind(c(I2_study, I2_sdep, I2_tg_phylo, I2_nb_phylo, I2_tg, I2_nb, I2_obs),
                   c(sd_I2_study, sd_I2_sdep, sd_I2_tg_phylo, sd_I2_nb_phylo, sd_I2_tg, sd_I2_nb, sd_I2_obs))
    
    
    # replacement
    I2_sdep <- mean(rstan::extract(rp)$I2_sdep)
    I2_study <- mean(rstan::extract(rp)$I2_study)
    I2_nb <- mean(rstan::extract(rp)$I2_nb)
    I2_tg <- mean(rstan::extract(rp)$I2_tg)
    I2_obs <- mean(rstan::extract(rp)$I2_obs)
    I2_nb_phylo <- mean(rstan::extract(rp)$I2_nb_phylo)
    I2_tg_phylo <- mean(rstan::extract(rp)$I2_tg_phylo)
    
    sd_I2_sdep <- sd(rstan::extract(rp)$I2_sdep)
    sd_I2_study <- sd(rstan::extract(rp)$I2_study)
    sd_I2_nb <- sd(rstan::extract(rp)$I2_nb)
    sd_I2_tg <- sd(rstan::extract(rp)$I2_tg)
    sd_I2_obs <- sd(rstan::extract(rp)$I2_obs)
    sd_I2_nb_phylo <- sd(rstan::extract(rp)$I2_nb_phylo)
    sd_I2_tg_phylo <- sd(rstan::extract(rp)$I2_tg_phylo)
    
    I2.rp <- cbind(c(I2_study, I2_sdep, I2_tg_phylo, I2_nb_phylo, I2_tg, I2_nb, I2_obs),
                   c(sd_I2_study, sd_I2_sdep, sd_I2_tg_phylo, sd_I2_nb_phylo, sd_I2_tg, sd_I2_nb, sd_I2_obs))
    
    df <- data.frame(Model = rep(c("ANNEs", "RNNEs"), each = 7), 
                     effect_name = rep(c("study", "sampling group", "target phylogeny", "neighbor phylogeny", "target species", "neighbor species", "observation"), 2), 
                     effect_ratio = c(I2.ad[ ,1], I2.rp[ ,1]), 
                     effect_error = c(I2.ad[ ,2], I2.rp[ ,2]), dataset = rep(dataset.name, 14))
    return (df)
}

# load data (one at a time)
# all data; target ratio
load("all_data/additive.subs20.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("all_data/m2.mrem.replacement.subs20.TargetRatio.Mean.20000iter.RData")
MEAN_rp <- m2

df1 <- mean.sd.fn()
df1$effect_ratio <- df1$effect_ratio*100
df1$effect_error <- df1$effect_error*100

# pairs; target ratio
load("reciprocal_effects/additive.subs20.Pairs.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("reciprocal_effects/m2.mrem.replacement.subs20.Pairs.TargetRatio.Mean.20000iter.RData")
MEAN_rp <- m2

df2 <- mean.sd.fn(dataset.name = "pairs")
df2$effect_ratio <- df2$effect_ratio*100
df2$effect_error <- df2$effect_error*100


df_all <- rbind(df1, df2)

df_all$effect_name <- factor(df_all$effect_name, levels = c("study", "sampling group", "target species", "neighbor species", "observation"))


p1 <- ggplot(df_all, aes(x = effect_name, y = effect_ratio, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = effect_ratio - effect_error, ymax = effect_ratio + effect_error), width = .2,
                  position = position_dodge(.9)) + 
    facet_grid(. ~ dataset) + #annotate("text", label = c("(a)", "(b)", "(c)"), size = 3, x = 5.5, y = .5) + 
    coord_flip() + 
    scale_fill_manual(values = c("#dfc27d", "#a6611a")) + 
    theme_classic() + 
    # guides(fill = FALSE) +
    guides(fill = guide_legend(title = "Model without phylogenetic correction")) + 
    # labs(y = "% Total heterogeneity", x = "Source of heterogeneity") + 
    theme(strip.text.x = element_blank(), 
          # axis.title = element_text(size = 10),
          # axis.title.x = element_text(family = "Calibri", size = 9, face = "bold"),
          # axis.title.y = element_text(family = "Calibri", size = 9, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = "Calibri", size = 8, face = "plain", colour = "black"),
          axis.text.y = element_text(family = "Calibri", size = 8, face = "plain", colour = "black"),
          plot.margin = unit(c(1, .5, .5, 2.1), "lines"),
          panel.spacing = unit(2, "lines"),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(family = "Calibri", color = "black", size = 8, face="bold"),
          legend.text =  element_text(family = "Calibri", color = "black", size = 8, face="plain"))
p1


# phylo
load("all_data_phylo/additive.phylo2.sdep2.subs20.TargetRatio.Mean.RData")
MEAN_ad_phylo <- m2

load("all_data_phylo/m2.phylo2.mrem.sdep2.replacement.subs20.TargRatio.Mean.20000iter.RData")
MEAN_rp_phylo <- m2

df3 <- mean.sd.phylo.fn()
df3$effect_ratio <- df3$effect_ratio*100
df3$effect_error <- df3$effect_error*100

load("reciprocal_effects_phylo/additive.phylo2.sdep2.subs20.Pairs.TargetRatio.Mean.RData")
MEAN_ad_phylo <- m2

load("reciprocal_effects_phylo/m2.phylo2.mrem.sdep2.Pairs.replacement.subs20.TargRatio.Mean.20000iter.RData")
MEAN_rp_phylo <- m2

df4 <- mean.sd.phylo.fn(dataset.name = "pairs")
df4$effect_ratio <- df4$effect_ratio*100
df4$effect_error <- df4$effect_error*100

df_all <- rbind(df3, df4)
df_all$effect_name <- factor(df_all$effect_name, 
                             levels = c("study", "sampling group", "target phylogeny", "neighbor phylogeny", 
                                        "target species", "neighbor species", "observation"))


p2 <- ggplot(df_all, aes(x = effect_name, y = effect_ratio, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = effect_ratio - effect_error, ymax = effect_ratio + effect_error), width = .2,
                  position = position_dodge(.9)) + 
    facet_grid(. ~ dataset) + #annotate("text", label = c("(a)", "(b)", "(c)"), size = 3, x = 5.5, y = .5) + 
    coord_flip() + 
    scale_fill_manual(values = c("#80cdc1", "#018571")) + 
    theme_classic() + 
    # guides(fill = FALSE) +
    guides(fill = guide_legend(title = "Model with phylogenetic correction      ")) + 
    # labs(y = "% Total heterogeneity", x = "Source of heterogeneity") + 
    theme(strip.text.x = element_blank(), 
          # axis.title = element_text(size = 10),
          # axis.title.x = element_text(family = "Calibri", size = 9, face = "bold"),
          # axis.title.y = element_text(family = "Calibri", size = 9, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(family = "Calibri", size = 8, face = "plain", colour = "black"),
          axis.text.y = element_text(family = "Calibri", size = 8, face = "plain", colour = "black"),
          plot.margin = unit(c(1, .5, 2, 1.5), "lines"),
          panel.spacing = unit(2, "lines"),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_text(family = "Calibri", color = "black", size = 8, face = "bold"),
          legend.text =  element_text(family = "Calibri", color = "black", size = 8, face = "plain"))
p2


grid.newpage()
# pushViewport(viewport(layout=grid.layout(2, 2)))#, widths = c(0.5, 0.5), heights = c(.2, .2, .3, .3))))
# print(plot5, vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
# print(plot6, vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
# print(plot7, vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
# print(plot8, vp=viewport(layout.pos.row = 2,layout.pos.col = 2))

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend) }

mylegend1 <- g_legend(p1)
mylegend2 <- g_legend(p2)

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               nrow=1),
                   arrangeGrob(p2 + theme(legend.position="none"),
                               nrow=1),
                   mylegend1, mylegend2, nrow=4, heights=c(12, 16, 2, 2))

grid.text(x=0.1, y=.98, label= "(a)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.1, y=.6, label= "(c)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.62, y=.98, label= "(b)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.62, y=.6, label= "(d)", gp = gpar(cex = .8, fontface="bold"))

grid.text(x=0.55, y=.18, label= "% Total heterogeneity", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.03, y=.6, label= "Source of heterogeneity", rot = 90, gp = gpar(cex = .8, fontface="bold"))

