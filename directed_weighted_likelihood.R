#Function to calculate likelihood for directed weighted nets with constrOptim using find_latent


directed_weighted_lkl <- function(par, data){
  #Test:
  #test_data <- matrix(rpois(100,8),nrow=10,ncol=10)
  #data <- c(apply(test_data,1,sum),apply(test_data,2,sum))
  #par <- rep(0.5,length(data))
  #Here we want the data to be in the form:
  #N is the whole length, the double of the number 
  #of nodes. Such that each unique indeg is in the
  #first half and the unique outdeg is in the second
  
  
  #We walso want the the parameters to be the doubled
  #size.
  

  N <- length(par)
  M <- N/2
  
  x <- par[1:M]
  y <- par[(M+1):N]
  
  deg_out <- data[1:M]
  deg_in <- data[(M+1):N]
  
  lkl <- c()
  
  #Make this vectorized!
  for(i in 1:M){
    lkl[i] <-  deg_out[i]*log(x[i]) + deg_in[i]*log(y[i])
    for(j in as.vector(1:M)[-i]){
      lkl[i] <- lkl[i] + log(1-x[i]*y[j])
    }
  }
  
  lkl <- -sum(lkl)
  
  return(lkl)
}
