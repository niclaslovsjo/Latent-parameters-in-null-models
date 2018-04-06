#This takes a vector of latent variables and returns expected values

expected_weighted_directed <- function(param){
  
  n <- length(param)
  indeg <- param[1:(n/2)]
  outdeg <- param[(n/2+1):n]
  N <- length(indeg)
  
  s_out <- c()
  s_in <- c()
  for(i in 1:N){
    temp_in <- c()
    temp_out <- c()
    for(j in as.vector(1:N)[-i]){
      temp_in <- c(temp_in,indeg[i]*outdeg[j]/(1-indeg[i]*outdeg[j])) 
      temp_out <- c(temp_out,outdeg[i]*indeg[j]/(1-outdeg[i]*indeg[j])) 
    }
    s_out[i] <- sum(temp_in)
    s_in[i] <- sum(temp_out)
  }
  s <- data.table(expected_in_degree = s_in,expected_out_degree=s_out)
  return(s)
}