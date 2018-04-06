#Constraints for constrOptim
#On the form: ui%*%theta >= ci
# test_data <- matrix(rpois(100,8),nrow=10,ncol=10)
# diag(test_data)<-0
# 
# test_data <- fread("C:/Users/niclo06/Dropbox/valued null networks/data2/7cont_nodiag.txt")
# test_data <- as.matrix(test_data[,-1])

find_latent <- function(data, 
                        optim_method="Nelder-Mead",
                        init_values = rep(0.5,length(data)),
                        type="weighted_directed"
                        ){
  
  if(is.matrix(data)){
    data <- as.vector(c(apply(data,1,sum),apply(data,2,sum)))
  } else {cat("Assuming data is a vector of in-degrees and out-degrees for each node.")}
  
  
  N <- length(data)
  
  #---------------------------------------------------------------
  if(type=="weighted_directed"){
    
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
                          f = dir_weighted_lkl,
                          grad = NULL,
                          ui = ui,
                          ci = ci,
                          method=optim_method,
                          control=list(maxit=100000))
    
    
  }else if(type="binary_undirected"){
    
    #Unconstrained optim, instead taking exponentials to force latent>0
    result <- optim(par = init_values,
                    data = data,
                    fn=undir_bin_lkl,
                    method="CG",
                    control=list(maxit=1000))
  }

  
  if(result$convergence==0){
    cat("Converged!\n")
    cat("Used", result$outer.iterations, "iterations.\n")
  }
  
  
  return(result$par)
}


