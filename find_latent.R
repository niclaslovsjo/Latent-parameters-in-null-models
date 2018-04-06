library(data.table)

#----------------------------------------------------------------
#Test-datasets
# test_data <- matrix(rpois(100,8),nrow=10,ncol=10)
# diag(test_data)<-0

# test_data <- fread("C:/Users/niclo06/Dropbox/valued null networks/data2/7cont_nodiag.txt")
# test_data <- as.matrix(test_data[,-1])

library(igraphdata)
library(igraph)
data("karate")
test_data <- as_adjacency_matrix(karate)
test_data <- as.matrix(test_data)

#-----------------------------------------------------------------


find_latent <- function(data,
                        type="weighted_directed"){
  

  
  #---------------------------------------------------------------
  if(type == "weighted_directed"){
    
    
    if(is.matrix(data)){
      data <- as.vector(c(apply(data,1,sum),apply(data,2,sum)))
    } else {cat("Assuming data is a vector of in-degrees and out-degrees for each node.")}
    
    
    init_values = rep(0.5,length(data))
    
    N <- length(data)
    n <- N/2
    
    
    #Set up constraints for constrOptim
    #This is simply a design matrix that follows
    #ui %*%theta >= ci
    
    ui1 <- matrix(0,ncol=N,nrow=N)
    diag(ui1)<-1
    ui2 <- matrix(0,ncol=N,nrow=N)
    diag(ui2)<- -1
    ui <- rbind(ui1,ui2)
    ci <- c(rep(0,N),rep(-1,N))
    
    #----------------------------------------------------------
    #Gradient function for directed weighted
    #Space for creating a gradient function
    
    
    #---------------------------------------------------------------
    #Optimization 
    #For now only by numerically approximated gradients
    
    result <- constrOptim(theta = init_values,
                          data = data,
                          f = directed_weighted_lkl,
                          grad = NULL,
                          ui = ui,
                          ci = ci,
                          method="Nelder-Mead",
                          control=list(maxit=100000))
    
    result <- data.table(indeg=result$par[1:n],outdeg=result$par[(n+1):N])
    }
  else if(type == "binary_undirected"){
    
    if(is.matrix(data)){
      data <- as.vector(c(apply(data,1,sum)))
    } else {cat("Assuming data is a vector of in-degrees and out-degrees for each node.")}
    
    
    init_values = rep(0.5,length(data))
    
    N <- length(data)
    
    #Unconstrained optim, instead taking exponentials to force latent>0
    result <- optim(par = init_values,
                    data = data,
                    fn=undirected_binary_lkl,
                    method="CG",
                    control=list(maxit=1000))
    result <- data.table(degree = exp(result$par))
  }
  
  return(result)
}

#first <- find_latent(data = test_data,type = "weighted_directed")
