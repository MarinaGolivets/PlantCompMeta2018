# Marina Golivets
# 6/14/2017
# Function to calculate biomass ratio in competition treatment


biomass.ratio.fn <- function (dat = dat) {
    
    # create an empty object to store data
    biomass.ratio <- NULL
    
    # run a for loop to calculate biomass ratios
    for (i in 1:length(unique(dat$pair.id))) {

            
        sp2.bms <- dat[dat$pair.id == unique(dat$pair.id)[i], ]$trt.mean[2] * 
                dat[dat$pair.id == unique(dat$pair.id)[i], ]$target.ratio[2] 
            
        sp1.bms <- dat[dat$pair.id == unique(dat$pair.id)[i], ]$trt.mean[1] * 
                dat[dat$pair.id == unique(dat$pair.id)[i], ]$target.ratio[1]
        
        tot.bms <- sp1.bms + sp2.bms
        
        biomass.ratio[dat$pair.id == unique(dat$pair.id)[i]] <- 
            c(sp2.bms/tot.bms, sp1.bms/tot.bms)
        
    } 
    
    return(biomass.ratio)
}