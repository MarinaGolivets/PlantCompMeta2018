# Marina Golivets
# 6/14/2017
# Function to calculate competitive asymmetry and subset data 


dat.comp.asym <- function (dat = dat) {
    
    # arrange data
    dat <- arrange(dat, pb.id, pair, yi)
    
    # create a sequance {1, 3, 5, ...}
    sequence <- seq(1, nrow(dat), by = 2)

    # create a data frame to store effect sizes
    pair.es <- as.data.frame(cbind(rep(NA, nrow(dat)), rep(NA, nrow(dat))))

    # calculate effect sizes
    for (i in sequence) {
    
        # delta Hedges d
        pair.es[i, 1] <- dat[i, ]$yi - dat[i+1, ]$yi
        pair.es[i+1, 1] <- dat[i+1, ]$yi - dat[i, ]$yi
    
        # its variance
        pair.es[i, 2] <- dat[i, ]$vi + dat[i+1, ]$vi
        pair.es[i+1, 2] <- dat[i, ]$vi + dat[i+1, ]$vi
    
    }

dat$yi <- pair.es[, 1]
dat$vi <- pair.es[, 2]

dat <- subset(dat, yi >= 0) # size should reduce to n/2

return(dat)

}