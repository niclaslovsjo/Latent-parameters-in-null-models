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




find_latent(data = test_data)->mm









res <- data.table(observed_in=data[1:(length(data)/2)],
                  estimated_in=result$par[1:(length(data)/2)],
                  observed_out = data[(length(data)/2+1):length(data)],
                  estimated_out = result$par[(length(data)/2+1):length(data)])
res <- res[order(estimated,decreasing = T)]



exp_out(indeg=result$par[1:(N/2)],outdeg=result$par[(N/2+1):N])
data[1:(N/2)]

exp_in(indeg=result$par[1:(N/2)],outdeg=result$par[(N/2+1):N])
data[(N/2+1):N]

#-----------------------------------------------------------
#----------------------------------------------------------------------
#Generate new nets
n <- length(result$par)
probs <- unlist(lapply(result$par[1:(n/2)],function(el1){lapply(result$par[(n/2+1):n],function(el2)el1*el2)}))

probs <- unlist(lapply(result$par[1:(n/2)],function(el1){lapply(result$par[(n/2+1):n],function(el2)el1*el2/(1-el1*el2))}))
matrix(probs,nrow=n/2,ncol=n/2,byrow=T)->carl

for(i in 1:(n/2)){
  carl[i,] <- carl[i,]*data[1:7]
}

gen_net <-lapply(1:1000,function(i){
  out_vec <- rbinom(length(probs),1,probs)
  data.table(sim_vec = out_vec,gen_nr=i,element_id=1:length(out_vec))
  #out_mat <- matrix(out_vec, nrow=N,ncol=N,byrow=TRUE)
})

gen_net <- rbindlist(gen_net)
gen_net[,means:=mean(sim_vec),by="element_id"]
gen_net[,st_dev:=sd(sim_vec),by="element_id"]

#Test if equal
mean(gen_net$means)
edge_density(graph_from_adjacency_matrix(A))



#------------------------------------------

exp_out <- function(indeg,outdeg){
  N <- length(indeg)
  s_out <- c()
  for(i in 1:N){
    temp_in <- c()
    for(j in as.vector(1:N)[-i]){
       temp_in <- c(temp_in,indeg[i]*outdeg[j]/(1-indeg[i]*outdeg[j])) 
    }
    s_out[i] <- sum(temp_in)
  }
  return(s_out)
}

exp_in <- function(outdeg,indeg){
  N <- length(outdeg)
  s_in <- c()
  for(i in 1:N){
    temp_out <- c()
    for(j in as.vector(1:N)[-i]){
      temp_out <- c(temp_out,outdeg[i]*indeg[j]/(1-outdeg[i]*indeg[j])) 
    }
    s_in[i] <- sum(temp_out)
  }
  return(s_in)
}


