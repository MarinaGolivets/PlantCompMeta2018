# Marina Golivets
# 6/19/2017
# Function to standardize covariates

std.covs.fn <- function (dat) {
    
    dat$neighbor.life.span[dat$neighbor.life.span == "an, bn, pr"] <- "an"
    dat$neighbor.life.span[dat$neighbor.life.span == "an, bn"] <- "an"
    dat$neighbor.life.span[dat$neighbor.life.span == "bn, pr"] <- "bn"
    dat$neighbor.life.span[dat$neighbor.life.span == "an, pr"] <- "pr"

    dat$target.life.span[dat$target.life.span == "an, bn, pr"] <- "an"
    dat$target.life.span[dat$target.life.span == "an, bn"] <- "an"
    dat$target.life.span[dat$target.life.span == "bn, pr"] <- "bn"
    dat$target.life.span[dat$target.life.span == "an, pr"] <- "pr"

    dat$neighbor.life.habit[dat$neighbor.life.habit == "ssb, fb"] <- "fb"
    dat$neighbor.life.habit[dat$neighbor.life.habit == "tr, sb"] <- "sb"
    dat$neighbor.life.habit[dat$neighbor.life.habit == "sb, fb"] <- "fb"

    dat$target.life.habit[dat$target.life.habit == "ssb, fb"] <- "fb"
    dat$target.life.habit[dat$target.life.habit == "tr, sb"] <- "sb"
    dat$target.life.habit[dat$target.life.habit == "sb, fb"] <- "fb"
    dat$target.life.habit[dat$target.life.habit == "ssb, sb, fb"] <- "fb"

    dat$climate.group <- rep(NA, nrow(dat))
    dat$climate.group[dat$KG.climate == "Dfb" | dat$KG.climate == "Dfa" | dat$KG.climate == "Dfc" | dat$KG.climate == "Dwa" | dat$KG.climate == "Dsb" | dat$KG.climate == "Dsc"] <- "Continental"
    dat$climate.group[dat$KG.climate == "Cfa" | dat$KG.climate == "Csa" | dat$KG.climate == "Cfb" | dat$KG.climate == "Cfc" | dat$KG.climate == "Csb" | dat$KG.climate == "Cwa" | dat$KG.climate == "Cwb"] <- "Temperate"
    dat$climate.group[dat$KG.climate == "Bsk" |  dat$KG.climate == "BSk"  | dat$KG.climate == "BSh" | dat$KG.climate == "BWh" | dat$KG.climate == "Bwh"] <- "Dry"  
    dat$climate.group[dat$KG.climate == "Am" |  dat$KG.climate == "Aw" | dat$KG.climate == "Af" ] <- "Tropical" 


    dat$hab.type <- dat$hab.type1
    dat$hab.type[dat$hab.type1 == "salt marsh" |  dat$hab.type1 == "mangrove forest"  | dat$hab.type1 == "bog"] <- "marsh"
    dat$hab.type[dat$hab.type1 == "serpentine seep"] <- "grassland"
    dat$hab.type[dat$hab.type == "marsh"] <- "wetland"
    dat$hab.type[dat$hab.type == "old field"] <- "ruderal/disturbed"
    # dat$hab.type <- factor(dat$hab.type, levels = unique(dat$hab.type)[c(4,3,2,5,7,6,1,8,9)])

    dat$group1 <- dat$group
    dat[dat$group == "nntv-ntv", ]$group1 <- "ntv-nntv"

    
    dat[dat$target.glob.invas.status == "naturalized" |
        dat$target.glob.invas.status == "casual" |
        dat$target.glob.invas.status == "none", ]$target.glob.invas.status <- "non-invasive"

    dat[dat$neighbor.glob.invas.status == "naturalized" |
        dat$neighbor.glob.invas.status == "casual" |
        dat$neighbor.glob.invas.status == "none", ]$neighbor.glob.invas.status <- "non-invasive"

    dat$target.glob.invas.status <- factor(dat$target.glob.invas.status, levels = c("non-invasive", "invasive"))
    dat$neighbor.glob.invas.status <- factor(dat$neighbor.glob.invas.status, levels = c("non-invasive", "invasive"))

    
table(dat$neighbor.monocot)
dat[dat$neighbor.monocot == "none", ]$neighbor.monocot <- "n"
dat[dat$target.monocot == "none", ]$target.monocot <- "n"

return(dat)

}