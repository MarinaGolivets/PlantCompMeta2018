# Marina Golivets
# Last modified: 10/24/2017

# The function calculates:
# a) non-bias-corrected Hedges d for observations that share control treatments following Gleser & Olkin (2009)

# b) bias-corrected Hedges d for observations that share control treatments following Stevens & Taylor (2008)

# Gleser, L. and Olkin, I. (2009). Stochastically dependent effect sizes. The handbook of research 
# synthesis and meta-analysis, 357.

# Stevens, J.R. & Taylor, A.M. (2008) Hierarchical dependence in meta-analysis. 
# J. Educ. Behav. Stat. 


HedgesD.fn <- function (Data, biasCorrection = FALSE) {
  
  # arrange data in alphabethical order
  Data <- arrange(Data, s.dep)
  
  # calculate total sample sizes for blocks (groups) of correlated studies
  Data$n.tot  <- unlist(lapply(split(Data, Data$s.dep), 
                              function(x) rep(sum(x$trt.s.size) + unique(x$ctrl.s.size), each = nrow(x))))
  
  # calculate deegres of freedom for blocks (groups) of correlated studies
  Data$d.f  <- unlist(lapply(split(Data, Data$s.dep), 
                               function(x) rep(sum(x$trt.s.size) + unique(x$ctrl.s.size) - nrow(x) - 1, 
                                               each = nrow(x))))
  
  # calculate pooled variance
  Data$Sp <- unlist(lapply(split(Data, Data$s.dep), 
                          function(x) sqrt((sum(x$trt.sd^2 * (x$trt.s.size - 1)) + 
                                                  unique(x$ctrl.sd)^2 * (unique(x$ctrl.s.size) - 1)) / x$d.f)) )
                                                  
  if (biasCorrection == TRUE) {
    
    # correction for small sample size
    Data$J <- with(Data, 1 - 3 / (4 * (d.f - 1)) )
    
    # Hedges' d
    Data$yi <- with(Data, (trt.mean - ctrl.mean) / Sp * J )
    
    # p (Steven & Taylor 2008)
    Data$p <- with(Data, J^2 * d.f / (d.f - 2) )
      
    # Hedges' d variance
    Data$vi <- with(Data, p * (1/ctrl.s.size + 1/trt.s.size) + (p - 1) * yi^2 )
    
  } else {
    
    # Hedges' d (without bias correction)
    Data$yi <- with(Data, (trt.mean - ctrl.mean) / Sp )
    
    # Hedges' d variance (without bias correction)
    Data$vi <- with(Data, 1/ctrl.s.size + 1/trt.s.size + yi^2 / (2 * n.tot) )
    
  }
  

  return(Data)
  
}