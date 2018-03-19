# Marina Golivets
# Last modified: 10/24/2017

# The function calculates a variance-covariance matrix for standardized mean differences following 
# a) Gleser & Olkin (2009)
# b) Stevens & Taylor (2008) (bias-corrected method)

# Gleser, L. and Olkin, I. (2009). Stochastically dependent effect sizes. The handbook of research 
# synthesis and meta-analysis, 357.

# Stevens, J.R. & Taylor, A.M. (2008) Hierarchical dependence in meta-analysis. 
# J. Educ. Behav. Stat. 

# code modified from: http://www.metafor-project.org/doku.php/analyses:gleser2009


calcV.fn <- function (Data, biasCorrection = FALSE, returnArray = FALSE) {
    
    # arrange data
    Data <- arrange(Data, s.dep)
    
    # calculate a variance-covariance matrix for a group of ESs with shared control:
    
    calc.v1 <- function (x) {
    
      v <- matrix((x$p[1] - 1) * outer(x$yi, x$yi, "*"), nrow = nrow(x), ncol = nrow(x))
        
      diag(v) <- x$vi
      
      return(v)
    
    }
    
    
    calc.v2 <- function (x) {
        
      v <- matrix(1/x$ctrl.s.size[1] + outer(x$yi, x$yi, "*") / (2*x$n.tot[1]), nrow = nrow(x), ncol = nrow(x))
        
      diag(v) <- x$vi
      
      return(v)
      
    }
    
    # calculate the entire variance-covariance matrix
    
    if (biasCorrection == TRUE) {
      
      V <- lapply(split(Data, Data$s.dep), calc.v1)
      
    } else {
      
      V <- lapply(split(Data, Data$s.dep), calc.v2)
      
    }
    
    # calculate a correlation matrix
    # R <- solve(sqrt(diag(Data$vi))) %*% V %*% solve(sqrt(diag(Data$vi)))
    
    # create an array where each element is of size == max size of sampling group
    if (returnArray == TRUE) {
      
       V.array <- array(data = 0, dim = c(length(unique(Data$s.dep)), 
                                       max(table(Data$s.dep)), max(table(Data$s.dep))))
    
      for (i in 1:length(unique(Data$s.dep))) {
        
       k <- length(which(Data$s.dep == unique(Data$s.dep)[i]))
       
       V.array[i, 1:k, 1:k] <- V[[i]]
       
      }
     
     return(V.array)
     
    } else {
      
      return(bldiag(V))
      
    }
    
}

#### END OF CODE ####    