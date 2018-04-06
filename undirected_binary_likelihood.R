#Function to calculate the likelihood, to be used in optim with find_latent.

undirected_binary_lkl  <- function(par, data){
  #Test:
  #par=runif(n = length(k))
  #data = k
  
  data <- as.vector(data)
  
  #Reparametrize to get to restrict to non-neg:
  #Goes together with returning a log
  #and evaluating as exp
  par <- exp(par) + 10^-10
  
  N <- length(par)
  lkl <- c()
  
  for(i in 1:(N-1)){
    lkl[i] <-  log(par[i]^data[i])
    for(j in (i+1):N){
      lkl[i] <- lkl[i] - log(1+par[i]*par[j])
    }
  }
  
  lkl[N] <-  log(par[N]^data[N])
  lkl <- -sum(lkl)
  
  return(log(lkl))
}