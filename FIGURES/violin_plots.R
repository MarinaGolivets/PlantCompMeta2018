# Marina Golivets
# March 2018
# Code to create violin plots



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
library(brms)


# function for summarizing model outputs (i.e. calculate the mean and 95% CrI)
sampSummary <- function(stan.obj) {
    
    # load stan model output (only betas in this case):
    mcmcSamp <- rstan::extract(stan.obj)$beta
    
    # calculate the mean and credible intervals:
    summar <- rbind(apply(mcmcSamp, 2, mean), apply(mcmcSamp, 2, quantile, c(.025, .975)))
    
    # return the latter:
    return(summar)
    
}


group_means_fn <- function(stan_data, n) {
    group_means <- as.data.frame(matrix(NA, nrow = 1500*4, ncol = 2), stringsAsFactors = F)
    
    colnames(group_means) <- c("group_value", "group_name")
    
    k <- 0
    group_names <- c(paste("non-native\non non-native\nn = ", n[1]),
                    paste("native on\nnon-native\nn = ", n[2]),
                    paste("non-native\non native\nn = ", n[3]),
                    paste("native on\nnative\nn = ", n[4]))
               
    for(i in 1:4) {
        group_means[(k+1):(k+1500), ] <- cbind(as.vector((as.array(stan_data)[,,i])), 
                                               rep(group_names[i], 1500))
        k <- 1500 * i
    }    
    
    group_means$group_value <- as.numeric(group_means$group_value)
    
    return(group_means)
}

pairdif_means_fn <- function(stan_data, n) {
    group_means <- as.data.frame(matrix(NA, nrow = 1500*3, ncol = 2), stringsAsFactors = F)
    
    colnames(group_means) <- c("group_value", "group_name")
    
    k <- 0
    group_names <- c(paste("non-native\nvs. non-native\nn = ", n[1]),
                     paste("native vs.\nnon-native\nn = ", n[2]),
                     paste("native vs.\nnative\nn = ", n[3]))
    
    for(i in 1:3) {
        group_means[(k+1):(k+1500), ] <- cbind(as.vector((as.array(stan_data)[,,i])), 
                                               rep(group_names[i], 1500))
        k <- 1500 * i
    }    
    
    group_means$group_value <- as.numeric(group_means$group_value)
    
    return(group_means)
}

violin_plot1 <- function(df = group_means_ad, 
                         y_lab = "ANNE",
                         difs = c("", "1,2", "1", "2"),
                         overall_mean = Mean_ad.summary[1, ],
                         y_add = rep(.1, 2)) {
    
    p <-  ggplot(df, aes(x = group_name, y = group_value) ) + 
        
        stat_ydensity(geom = "violin", trim = T, draw_quantiles = c(0.025, 0.5, 0.975), 
                      color = colors()[205], fill = "#dfc27d") +
        
        stat_summary(fun.y = mean, geom = "point", size = 2, color = "#a6611a") + 
        
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                     geom = "pointrange", color = "#a6611a") + 
        
        theme(strip.text.x = element_blank(), 
              axis.title.x = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.title.y = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.text.x = element_text(family = "Calibri", size = 8.5, face = "plain", colour = "black"),
              axis.text.y = element_text(family = "Calibri", size = 8.5, face = "plain", colour = "black"),
              axis.title = element_text(size = 10.5),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              plot.margin = unit(c(2.5, .5, 0, .5), "lines")) +
        
        labs(y = y_lab, x = "") +
        
        annotate("text", x = c(1, 2, 3, 4), 
                 y =  aggregate(df$group_value ~ as.factor(df$group_name), FUN = max)[ ,2] + y_add, 
                 label = difs, size = 2) +
        
        geom_hline(yintercept = overall_mean, color = "#a6611a", linetype = "solid", size = .75) 
        
}  

violin_plot11 <- function(df = pairdif_means_ad, 
                          y_lab = "ANNE asymmetry",
                          difs = c("1,2", "1", "2"),
                          overall_mean = Mean_ad.summary[1, ],
                          y_add = rep(.1, 3)) {
    
    p <-  ggplot(df, aes(x = group_name, y = group_value) ) + 
        
        stat_ydensity(geom = "violin", trim = T, draw_quantiles = c(0.025, 0.5, 0.975), 
                      color = colors()[205], fill = "#dfc27d") +
        
        stat_summary(fun.y = mean, geom = "point", size = 2, color = "#a6611a") + 
        
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                     geom = "pointrange", color = "#a6611a") + 
        
        theme(strip.text.x = element_blank(), 
              axis.title.x = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.title.y = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.text.x = element_text(family = "Calibri", size = 8.5, face = "plain", colour = "black"),
              axis.text.y = element_text(family = "Calibri", size = 8.5, face = "plain", colour = "black"),
              axis.title = element_text(size = 10.5),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              plot.margin = unit(c(1, .5, 0, .5), "lines")) +
        
        labs(y = y_lab, x = "") +
        
        annotate("text", x = c(1, 2, 3), 
                 y =  aggregate(df$group_value ~ as.factor(df$group_name), FUN=max)[ ,2] + y_add, 
                 label = difs, size = 2) +
        
        geom_hline(yintercept = overall_mean, color = "#a6611a", linetype = "solid", size = .75) 
    
} 


violin_plot2 <- function(df = group_means_ad, 
                         y_lab = "ANNE",
                         difs = c("1", "3", "1,2", "3,4", "2", "4", "", ""),
                         overall_mean = Mean_ad.summary[1, ],
                         overall_mean_phylo = Mean_ad_phylo.summary[1,],
                         y_add = rep(.2, 8)) {

    df$new_groups <- factor(interaction(df$group_name, df$Model), 
                            levels = sort(unique(as.character(interaction(df$group_name, df$Model))), 
                                          decreasing = T))
    
    p <-  ggplot(df, aes(x = group_name, y = group_value, fill = Model) ) +
        
        stat_ydensity(geom = "violin", trim = T, draw_quantiles = c(0.025, 0.5, 0.975), 
                      color = colors()[205], 
                      position = position_dodge(1), adjust = 1) +
        
        scale_fill_manual(values = c("#dfc27d", "#80cdc1")) +
        
        # guides(fill = FALSE, color = FALSE) +
        
        stat_summary(fun.y = mean, geom = "point", size = 2, 
                     aes(color = Model), position = position_dodge(1)) + 
        
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                     geom = "pointrange", aes(color = Model), 
                     position = position_dodge(1)) + 
        
        scale_color_manual(values = c("#a6611a", "#018571")) +
        
        theme(strip.text.x = element_blank(), 
              axis.title.x = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.title.y = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.text.x = element_text(family = "Calibri", size = 8.5, face = "plain", color = "black"),
              axis.text.y = element_text(family = "Calibri", size = 8.5, face = "plain", color = "black"),
              legend.title = element_text(family = "Calibri", color = "black", size = 9.3, face="bold"),
              legend.text =  element_text(family = "Calibri", color = "black", size = 9.3, face="plain"),                           
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              plot.margin = unit(c(2.5, .5, 0, .5), "lines"),
              legend.position="bottom",
              legend.box = "horizontal") +
        
        labs(y = y_lab, x = "") +
        
        annotate("text", x = c(.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25), 
                 y =  aggregate(df$group_value ~ df$new_groups, FUN = max)[ ,2][c(7,8,5,6,3,4,1,2)] + y_add,
                 label = difs, size = 2) +
        
        geom_hline(yintercept = overall_mean_phylo, color = "#018571", linetype = "solid", size = .75) + 
        
        geom_hline(yintercept = overall_mean, color =  "#a6611a", linetype = "solid", size = .75) 

} 

violin_plot22 <- function(df = pairdif_means_ad, 
                          y_lab = "ANNE asymmetry",
                          difs = c("1", "2", "", "", "1", "2"),
                          overall_mean = Mean_ad.summary[1, ],
                          overall_mean_phylo = Mean_ad_phylo.summary[1, ],
                          y_add = rep(.2, 6)) {
    
    df$new_groups <- factor(interaction(df$group_name, df$Model), 
                            levels = sort(unique(as.character(interaction(df$group_name, df$Model))), 
                                          decreasing = T))
    
    p <-  ggplot(df, aes(x = group_name, y = group_value, fill = Model) ) +
        
        stat_ydensity(geom = "violin", trim = T, draw_quantiles = c(0.025, 0.5, 0.975), 
                      color = colors()[205], 
                      position = position_dodge(1), adjust = 1) +
        
        scale_fill_manual(values = c("#dfc27d", "#80cdc1")) +
        
        # guides(fill = FALSE, color = FALSE) +
        
        stat_summary(fun.y = mean, geom = "point", size = 2, 
                     aes(color = Model), position = position_dodge(1)) + 
        
        stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                     geom = "pointrange", aes(color = Model), 
                     position = position_dodge(1)) + 
        
        scale_color_manual(values = c("#a6611a", "#018571")) +
        
        theme(strip.text.x = element_blank(), 
              axis.title.x = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.title.y = element_text(family = "Calibri", size = 9.5, face = "plain"),
              axis.text.x = element_text(family = "Calibri", size = 8.5, face = "plain", color = "black"),
              axis.text.y = element_text(family = "Calibri", size = 8.5, face = "plain", color = "black"),
              legend.title = element_text(family = "Calibri", color = "black", size = 9.3, face="bold"),
              legend.text =  element_text(family = "Calibri", color = "black", size = 9.3, face="plain"),                           
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), 
              panel.background = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"),
              plot.margin = unit(c(1, .5, 0, .5), "lines"),
              legend.position="bottom",
              legend.box = "horizontal") +
        
        labs(y = y_lab, x = "") +
        
        annotate("text", x = c(.75, 1.25, 1.75, 2.25, 2.75, 3.25), 
                 y =  aggregate(df$group_value ~ df$new_groups, FUN = max)[ ,2][c(5,6,3,4,1,2)] + y_add,
                 label = difs, size = 2) +
        
        geom_hline(yintercept = overall_mean_phylo, color = "#018571", linetype = "solid", size = .75) + 
        
        geom_hline(yintercept = overall_mean, color =  "#a6611a", linetype = "solid", size = .75) 
    
} 

################     ALL DATA    ###################

# load data
dat <- read.csv(file = "MGolivets_plant_biomass_additive.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat <- subset(dat, dat$yi > -20 & dat$yi < 20)
dat_ad <- dat

dat <- read.csv(file = "MGolivets_plant_biomass_replacement.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat <- subset(dat, dat$yi > -20 & dat$yi < 20)
dat_rp <- dat


# load workspaces (from portable hard drive)
load("additive.subs20.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.TargetRatio.BOG.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("m2.mrem.replacement.subs20.TargetRatio.Mean.20000iter.RData")
MEAN_rp <- m2


# extract means
BOG_ad.summary <- sampSummary(BOG_ad)[, 1:4]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:4]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 


# BOG_ad.summary[1,] - BOG_rp.summary[1,]
# 
# # extract p-values
# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad)$p_g34)
# 
# hypothesis(BOG_ad, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                  "beta[2] = beta[3]", "beta[4] = beta[2]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp)$p_g34)
# 
# hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# hypothesis(BOG_rp, c("beta[1] = 0", "beta[2] = 0", "beta[3] = 0",
#                      "beta[4] = 0"), alpha = 0.05, seed = 500)
           


# violin plots
group_means_ad <- group_means_fn(stan_data = BOG_ad, n = as.numeric(table(dat_ad$group)))

group_means_rp <- group_means_fn(stan_data = BOG_rp, n = as.numeric(table(dat_rp$group)))



plot1 <- violin_plot1()
plot2 <- violin_plot1(df = group_means_rp, 
                   y_lab = "RNNE",
                   difs = c("1", "1,2,3", "2", "3"),
                   overall_mean = Mean_rp.summary[1, ],
                   y_add = rep(.2, 4))
plot2 <- plot2 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")


##################     DATA ON RECIPROCAL EFFECTS     ####################

# load data
dat_ad <- arrange(dat_ad, pb.id, pair)
dat_ad <- dat_ad[!is.na(dat_ad$pair), ]
names <- names(which(table(dat_ad$pair.id) != 2))
dat_ad <- dat_ad[-which(dat_ad$pair.id %in% names), ] 

dat_rp <- arrange(dat_rp, pb.id, pair)
dat_rp <- dat_rp[!is.na(dat_rp$pair), ]
names <- names(which(table(dat_rp$pair.id) != 2))
dat_rp <- dat_rp[-which(dat_rp$pair.id %in% names), ] 


# subs 20; pairs; target ratio
load("additive.subs20.Pairs.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.Pairs.TargetRatio.BOG.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.Pairs.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("m2.mrem.replacement.subs20.Pairs.TargetRatio.Mean.20000iter.RData")
MEAN_rp <- m2


# extract means
BOG_ad.summary <- sampSummary(BOG_ad)[, 1:4]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:4]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 

# BOG_ad.summary[1,] - BOG_rp.summary[1,]
# 
# # # extract p-values
# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad)$p_g34)
# 
# hypothesis(BOG_ad, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# # 
# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp)$p_g34)
# 
# hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)


group_means_ad <- group_means_fn(stan_data = BOG_ad, n = as.numeric(table(dat_ad$group)))

group_means_rp <- group_means_fn(stan_data = BOG_rp, n = as.numeric(table(dat_rp$group)))



plot3 <- violin_plot1(difs = c("", "1,2", "1", "2"))
plot4 <- violin_plot1(df = group_means_rp, 
                   y_lab = "RNNE",
                   difs = c("1", "1,2", "2", ""),
                   overall_mean = Mean_rp.summary[1,],
                   y_add = rep(.2, 4))
plot4 <- plot4 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")

grid.newpage()
pushViewport(viewport(layout=grid.layout(2, 2)))#, widths = c(0.5, 0.5), heights = c(.2, .2, .3, .3))))
print(plot1, vp=viewport(layout.pos.row = 1,layout.pos.col = 1))
print(plot2, vp=viewport(layout.pos.row = 1,layout.pos.col = 2))
print(plot3, vp=viewport(layout.pos.row = 2,layout.pos.col = 1))
print(plot4, vp=viewport(layout.pos.row = 2,layout.pos.col = 2))

grid.text(x=0.02, y=.93, label= "(a)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.02, y=.43, label= "(c)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.93, label= "(b)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.43, label= "(d)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.97, label= "All effects", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.47, label= "Reciprocal effects", gp = gpar(cex = .8, fontface="bold"))


####################    ALL DATA PHYLO      #################

# load data
dat <- read.csv(file = "MGolivets.dat.phylo2.ad.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat$group1 <- dat$group
dat[dat$group == "nntv-ntv", ]$group1 <- "ntv-nntv"
table(dat$group1)
dat_ad <- dat

dat <- read.csv(file = "MGolivets.dat.phylo2.rp.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat$group1 <- dat$group
dat[dat$group == "nntv-ntv", ]$group1 <- "ntv-nntv"
table(dat$group1)
dat_rp <- dat


# load workspace
# subs 20; target ratio
load("additive.phylo2.sdep2.subs20.TargetRatio.BOG.RData")
BOG_ad_phylo <- m2

load("replacement.subs20.TargRatio.BOG.20000iter.RData")
BOG_rp_phylo <- m2

load("additive.phylo2.sdep2.subs20.TargetRatio.Mean.RData")
MEAN_ad_phylo <- m2

load("replacement.phylo2.subs20.Mean.RData")
MEAN_rp_phylo <- m2


# subs 20; subsPhylo2; target ratio
load("additive.subs20.subsPhylo.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.subsPhylo2.TargRatio.BOG.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.subsPhylo.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("replacement.subs20.subsPhylo2.TargRatio.Mean20000iter.RData")
MEAN_rp <- m2


# extract means
BOG_ad_phylo.summary <- sampSummary(BOG_ad_phylo)[, 1:4]
Mean_ad_phylo.summary <- as.matrix(sampSummary(MEAN_ad_phylo)[, 1]) 

BOG_rp_phylo.summary <- sampSummary(BOG_rp_phylo)[, 1:4]
Mean_rp_phylo.summary <- as.matrix(sampSummary(MEAN_rp_phylo)[, 1]) 

BOG_ad.summary <- sampSummary(BOG_ad)[, 1:4]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:4]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 


# # extract p-values
# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad)$p_g34)
# 
# p_ad_phyloSubs <- c(p1, p2, p3, p4, p5, p6)
# p_ad_phyloSubs
# 
# hypothesis(BOG_ad, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# 
# p1 <- mean(rstan::extract(BOG_ad_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad_phylo)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad_phylo)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad_phylo)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad_phylo)$p_g34)
# 
# p_ad_phylo <- c(p1, p2, p3, p4, p5, p6)
# p_ad_phylo
# 
# hypothesis(BOG_ad_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# 
# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp)$p_g34)
# 
# p_rp_phyloSubs <- c(p1, p2, p3, p4, p5, p6)
# p_rp_phyloSubs
# 
# hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                            "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# p1 <- mean(rstan::extract(BOG_rp_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp_phylo)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp_phylo)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp_phylo)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp_phylo)$p_g34)
# 
# p_rp_phylo <- c(p1, p2, p3, p4, p5, p6)
# p_rp_phylo
# 
# hypothesis(BOG_rp_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)


group_means_ad <- cbind(rbind(group_means_fn(stan_data = BOG_ad_phylo, n = as.numeric(table(dat_ad$group))), 
                              group_means_fn(stan_data = BOG_ad, n = as.numeric(table(dat_ad$group)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*4), 
                                       levels = c("without phylogenetic correction", "with phylogenetic correction")))

group_means_rp <- cbind(rbind(group_means_fn(stan_data = BOG_rp_phylo, n = as.numeric(table(dat_rp$group))), 
                              group_means_fn(stan_data = BOG_rp, n = as.numeric(table(dat_rp$group)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*4), 
                                       levels = c("without phylogenetic correction", "with phylogenetic correction")))


plot5 <- violin_plot2(difs = c("", "", "1", "2,3", "", "2", "1", "3"),
                      y_add = rep(.15, 8))
plot5 <- plot5 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")

plot6 <- violin_plot2(df = group_means_rp, 
                   y_lab = "RNNE",
                   difs = c("1", "3", "1,2", "3,4", "2", "4", "", ""),
                   overall_mean_phylo = Mean_rp_phylo.summary[1,],
                   overall_mean = Mean_rp.summary[1,],
                   y_add = rep(.3, 8))
plot6 <- plot6 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")



############    RECIPROCAL EFFECTS PHYLO    ##############

# load data
dat <- read.csv(file = "MGolivets.dat.phylo2.ad.Pairs.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_ad <- dat

dat <- read.csv(file = "MGolivets.dat.phylo2.rp.Pairs.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat_rp <- dat


# load workspace
# subs 20; target ratio
load("additive.phylo2.sdep2.subs20.Pairs.TargetRatio.BOG.RData")
BOG_ad_phylo <- m2

load("replacement.phylo2.subs20.Pairs.BOG.RData")
BOG_rp_phylo <- m2

load("additive.phylo2.sdep2.subs20.Pairs.TargetRatio.Mean.RData")
MEAN_ad_phylo <- m2

load("replacement.subs20.TargRatio.Pairs.Mean.20000iter.RData")
MEAN_rp_phylo <- m2


# subs 20; subsPhylo2; target ratio
load("additive.subs20.subsPhylo.Pairs.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.subsPhylo2.TargRatio.Pairs.BOG.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.subsPhylo.Pairs.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("replacement.subs20.subsPhylo2.TargRatio.Pairs.Mean.20000iter.RData")
MEAN_rp <- m2


# extract means
BOG_ad_phylo.summary <- sampSummary(BOG_ad_phylo)[, 1:4]
Mean_ad_phylo.summary <- as.matrix(sampSummary(MEAN_ad_phylo)[, 1]) 

BOG_rp_phylo.summary <- sampSummary(BOG_rp_phylo)[, 1:4]
Mean_rp_phylo.summary <- as.matrix(sampSummary(MEAN_rp_phylo)[, 1]) 

BOG_ad.summary <- sampSummary(BOG_ad)[, 1:4]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:4]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 

# # extract p-values
# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad)$p_g34)
# 
# p_ad_phyloSubs <- c(p1, p2, p3, p4, p5, p6)
# p_ad_phyloSubs
# 
# hypothesis(BOG_ad, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# 
# p1 <- mean(rstan::extract(BOG_ad_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad_phylo)$p_g14)
# p4 <- mean(rstan::extract(BOG_ad_phylo)$p_g23)
# p5 <- mean(rstan::extract(BOG_ad_phylo)$p_g24)
# p6 <- mean(rstan::extract(BOG_ad_phylo)$p_g34)
# 
# p_ad_phylo <- c(p1, p2, p3, p4, p5, p6)
# p_ad_phylo
# 
# hypothesis(BOG_ad_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                            "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# 
# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp)$p_g34)
# 
# p_rp_phyloSubs <- c(p1, p2, p3, p4, p5, p6)
# p_rp_phyloSubs
# 
# hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                      "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)
# 
# p1 <- mean(rstan::extract(BOG_rp_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp_phylo)$p_g14)
# p4 <- mean(rstan::extract(BOG_rp_phylo)$p_g23)
# p5 <- mean(rstan::extract(BOG_rp_phylo)$p_g24)
# p6 <- mean(rstan::extract(BOG_rp_phylo)$p_g34)
# 
# p_rp_phylo <- c(p1, p2, p3, p4, p5, p6)
# p_rp_phylo
# 
# hypothesis(BOG_rp_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]", "beta[1] = beta[4]",
#                            "beta[2] = beta[3]", "beta[2] = beta[4]", "beta[3] = beta[4]"), 
#            alpha = 0.05, seed = 500)


group_means_ad <- cbind(rbind(group_means_fn(stan_data = BOG_ad_phylo, n = as.numeric(table(dat_ad$group))), 
                              group_means_fn(stan_data = BOG_ad, n = as.numeric(table(dat_ad$group)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*4), 
                        levels = c("without phylogenetic correction", "with phylogenetic correction")))

group_means_rp <- cbind(rbind(group_means_fn(stan_data = BOG_rp_phylo, n = as.numeric(table(dat_rp$group))), 
                              group_means_fn(stan_data = BOG_rp, n = as.numeric(table(dat_rp$group)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*4), 
                                          levels = c("without phylogenetic correction", "with phylogenetic correction")))

plot7 <- violin_plot2(difs = c("", "", "1,2", "3,4", "1", "3", "2", "4"),
                      y_add = c(.15, .15, .15, .15, .22, .15, .15, .15))
plot7 <- plot7 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")
print(plot7)
plot8 <- violin_plot2(df = group_means_rp, 
                      y_lab = "RNNE",
                      difs = c("1", "3", "1,2", "3,4", "2", "4", "", ""),
                      overall_mean_phylo = Mean_rp_phylo.summary[1,],
                      overall_mean = Mean_rp.summary[1,],
                      y_add = rep(.3, 8))
plot8 <- plot8 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")


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

mylegend <- g_legend(plot8)

p3 <- grid.arrange(arrangeGrob(plot5 + theme(legend.position="none"),
                               plot6 + theme(legend.position="none"),
                               plot7 + theme(legend.position="none"),
                               plot8 + theme(legend.position="none"),
                               nrow=2),
                   mylegend, nrow=3, heights=c(10, 1, 0))
                    
grid.text(x=0.02, y=.93, label= "(a)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.02, y=.48, label= "(c)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.93, label= "(b)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.48, label= "(d)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.97, label= "All effects", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.5, y=.52, label= "Reciprocal effects", gp = gpar(cex = .8, fontface="bold"))


####################    COMPETITIVE INEQUALITY AND RELEASE    ###############
# load data
dat <- read.csv(file = "DATA/MGolivets_plant_biomass_additive_subs20_PairDifs.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat$group1 <- dat$group
dat[dat$group == "nntv-ntv", ]$group1 <- "ntv-nntv"
table(dat$group1)
dat_ad <- dat

dat <- read.csv(file = "MGolivets_plant_biomass_replacement_subs20_PairDifs.csv",
                header = TRUE, stringsAsFactors = FALSE)
dat$group1 <- dat$group
dat[dat$group == "nntv-ntv", ]$group1 <- "ntv-nntv"
table(dat$group1)
dat_rp <- dat


# subs 20; pair diffs; target ratio
load("additive.subs20.PairDifs.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.TargetRatio.RYT.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.PairDifs.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("replacement.subs20.TargetRatio.RYT.Mean.20000iter.RData")
MEAN_rp <- m2

# extract means
BOG_ad.summary <- sampSummary(BOG_ad)[, 1:3]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:3]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 


# # extract p-values
# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g23)
# p <- c(p1, p2, p3)
# p
# hypothesis(BOG_ad, c("beta[2] = beta[1]", "beta[1] = beta[3]",
#                      "beta[2] = beta[3]"), 
#            alpha = 0.05, seed = 500)
# 
# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g23)
# p_rp <- c(p1, p2, p3)
# p_rp
# hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]",
#                      "beta[2] = beta[3]"), 
#            alpha = 0.05, seed = 500)


pairdif_means_ad <- pairdif_means_fn(BOG_ad, n=as.numeric(table(dat_ad$group1)))
pairdif_means_rp <- pairdif_means_fn(BOG_rp, n=as.numeric(table(dat_rp$group1)))



plot9 <- violin_plot11(difs = c("1", "", "1"), y_lab = "Competitive inequality")
plot10 <- violin_plot11(df = pairdif_means_rp, 
                        y_lab = "Competitive release",
                        difs = c("1", "1", ""),
                        overall_mean = Mean_rp.summary[1, ],
                        y_add = rep(.2, 3))
plot10 <- plot10 + geom_hline(yintercept = 0, color = "black", linetype = "dotted")
plot10


####################    COMPETITIVE INEQUALITY AND RELEASE PHYLO   ###############
# load data
dat_ad <- read.csv(file = "MGolivets_plant_biomass_additive_phylo_subs20_PairDifs.csv", header = TRUE, stringsAsFactors = FALSE)
dat_rp <- read.csv(file = "OneDrive/META-ANALYSIS_PLANT_COMP/Meta_Pairwise_Plant_Competition/DATA/MGolivets_plant_biomass_replacement_phylo_subs20_PairDifs.csv", header = TRUE, stringsAsFactors = FALSE)

# subs 20; targ ratio
load("additive.phylo2.sdep2.subs20.PairDifs.TargetRatio.BOG.RData")
BOG_ad_phylo <- m2

load("replacement.subs20.TargetRatio.RYT.Phylo.BOG.20000iter.RData")
BOG_rp_phylo <- stan.output

load("additive.phylo2.sdep2.subs20.PairDifs.TargetRatio.Mean.RData")
MEAN_ad_phylo <- m2

load("replacement.subs20.TargetRatio.RYT.Phylo.Mean.20000iter.RData")
MEAN_rp_phylo <- stan.output


# subs 20; subsPhylo2; target ratio
load("additive.subs20.subsPhylo.PairDifs.TargetRatio.BOG.20000iter.RData")
BOG_ad <- m2

load("replacement.subs20.TargetRatio.PhyloSubs.RYT.BOG.20000iter.RData")
BOG_rp <- m2

load("additive.subs20.subsPhylo.PairDifs.TargetRatio.Mean.20000iter.RData")
MEAN_ad <- m2

load("replacement.subs20.TargetRatio.PhyloSubs.RYT.Mean.20000iter.RData")
MEAN_rp <- m2


# extract means
BOG_ad_phylo.summary <- sampSummary(BOG_ad_phylo)[, 1:3]
Mean_ad_phylo.summary <- as.matrix(sampSummary(MEAN_ad_phylo)[, 1]) 

BOG_rp_phylo.summary <- sampSummary(BOG_rp_phylo)[, 1:3]
Mean_rp_phylo.summary <- as.matrix(sampSummary(MEAN_rp_phylo)[, 1]) 

BOG_ad.summary <- sampSummary(BOG_ad)[, 1:3]
Mean_ad.summary <- as.matrix(sampSummary(MEAN_ad)[, 1]) 

BOG_rp.summary <- sampSummary(BOG_rp)[, 1:3]
Mean_rp.summary <- as.matrix(sampSummary(MEAN_rp)[, 1]) 

pairdif_means_ad <- cbind(rbind(pairdif_means_fn(stan_data = BOG_ad_phylo, n = as.numeric(table(dat_ad$group1))), 
                              pairdif_means_fn(stan_data = BOG_ad, n = as.numeric(table(dat_ad$group1)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*3), 
                                       levels = c("without phylogenetic correction", "with phylogenetic correction")))

pairdif_means_rp <- cbind(rbind(pairdif_means_fn(stan_data = BOG_rp_phylo, n = as.numeric(table(dat_rp$group1))), 
                              pairdif_means_fn(stan_data = BOG_rp, n = as.numeric(table(dat_rp$group1)))), 
                        Model = factor(rep(c("with phylogenetic correction", "without phylogenetic correction"), each = 1500*3), 
                                       levels = c("without phylogenetic correction", "with phylogenetic correction")))

# p1 <- mean(rstan::extract(BOG_ad)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad)$p_g23)
# p_ad <- c(p1, p2, p3)
# p_ad
hypothesis(BOG_ad, c("beta[1] = beta[2]", "beta[1] = beta[3]",
                     "beta[2] = beta[3]"), 
           alpha = 0.05, seed = 500)

# p1 <- mean(rstan::extract(BOG_rp)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp)$p_g23)
# p_rp <- c(p1, p2, p3)
# p_rp
hypothesis(BOG_rp, c("beta[1] = beta[2]", "beta[1] = beta[3]",
                     "beta[2] = beta[3]"), 
           alpha = 0.05, seed = 500)

# p1 <- mean(rstan::extract(BOG_ad_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_ad_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_ad_phylo)$p_g23)
# p_ad <- c(p1, p2, p3)
# p_ad
hypothesis(BOG_ad_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]",
                     "beta[2] = beta[3]"), 
           alpha = 0.05, seed = 500)

# p1 <- mean(rstan::extract(BOG_rp_phylo)$p_g12)
# p2 <- mean(rstan::extract(BOG_rp_phylo)$p_g13)
# p3 <- mean(rstan::extract(BOG_rp_phylo)$p_g23)
# p_rp <- c(p1, p2, p3)
# p_rp
hypothesis(BOG_rp_phylo, c("beta[1] = beta[2]", "beta[1] = beta[3]",
                     "beta[2] = beta[3]"), 
           alpha = 0.05, seed = 500)



plot11 <- violin_plot22(df = subset(pairdif_means_ad, group_value >=0), y_lab = "Competitive inequality",
                        difs = c("", "", "", "", "", "")) + 
  geom_hline(yintercept = 0, color = "black", linetype = "dotted")
plot12 <- violin_plot22(df = pairdif_means_rp,  
                        y_lab = "Competitive release",
                        difs = rep("", 6),
                        overall_mean = Mean_rp.summary[1, ],
                        overall_mean_phylo = Mean_rp_phylo.summary[1, ],
                        y_add = rep(.2, 6)) + geom_hline(yintercept = 0, color = "black", linetype = "dotted")

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

mylegend <- g_legend(plot12)

p3 <- grid.arrange(arrangeGrob(plot9 + theme(legend.position="none"),
                               plot10 + theme(legend.position="none"),
                               nrow=1),
                   arrangeGrob(plot11 + theme(legend.position="none"),
                               plot12 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=3, heights=c(13, 15, 2))

grid.text(x=0.02, y=.99, label= "(a)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.02, y=.54, label= "(c)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.99, label= "(b)", gp = gpar(cex = .8, fontface="bold"))
grid.text(x=0.52, y=.54, label= "(d)", gp = gpar(cex = .8, fontface="bold"))

         
