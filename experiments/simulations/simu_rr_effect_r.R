library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(rrpack)
library(corpcor)

source("experiments/simulations/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/alternative_methods/scca_chao.R')
source("experiments/evaluation.R")
source("src/reduced_rank_regression.R")


#Simulation for missing values in both X and Y

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
n <- as.numeric(args[3])
strength_theta <- args[4]
p <- as.numeric(args[5])
rs <- c(as.numeric(args[6]))
q_vals <- c(30)
overlaps <- c(1)
#props <- c(0, 0.1, 0.2)
props <- c(0)
noise = 1
seeds = 1:50
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(10)
result = c()
for (seed_n in seeds){
  #for (n in c(100, 300, 500, 1000, 10000)){
  set.seed(seed * 100 + seed_n)
  for (nnzeros in nnzero_values){
    #for(p in c(100,  200, 300,  500, 800, 80, 20)){
    #for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in q_vals){
      #for(nnzeros in c(5, 10, 15, 20, 50)){
      for (r in rs){
        if ( strength_theta == "high"){
          thetas <- diag(seq(0.9, 0.75, length.out = r))
        }else{
          if ( strength_theta == "medium"){
            thetas <- diag(seq(0.7, 0.55, length.out = r))
          }
          else{
            thetas <- diag(seq(0.5, 0.35, length.out = r))
          }
        }
        for (r_pca in c(5)){
          if ( (max(r_pca, r, nnzeros) < min(p)) & (nnzeros > max(r_pca, r)) ) {
            for (overlapping_amount in overlaps){
              for(prop_missing in props){
                cat("seed:")
                cat(seed, " ")
		nnzeros = r * 5
                print(c(n, r, r_pca, strength_theta))
                #gen = generate(n, p, q, s, prop_missing)
                start_time_creation <- system.time({
                gen = generate_example_sparse_U(n, p, q,
                                                r_pca = r_pca,
                                                nnzeros = nnzeros,
                                                theta = thetas,
                                                lambda_pca = 1,
                                                r = r,
                                                overlapping_amount = overlapping_amount,
                                                normalize_diagonal = normalize_diagonal,
                                                n_new = 5000) 
                })
                print(start_time_creation[[1]])
                X = gen$X
                Y = gen$Y
                #Xna = gen$Xna
                #Yna = gen$Yna
                Sigma0_sqrt = sqrtm(gen$Sigma)$B
                Sigma_hat_sqrt = sqrtm(gen$S)$B
                
                
                if (p < n){
                  
                  
                  start_time_rrr <- system.time({
                    rrr <- CCA_rrr(X, Y, Sx=NULL,
                                   Sy=NULL,
                                   lambda =0, Kx=NULL, r, highdim=FALSE,
                                   LW_Sy = LW_Sy)
                  })
                  
                  result = rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, rrr$U, rrr$V, gen$u, 
                                                             gen$v,
                                                             Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                             Sigma0_sqrt = Sigma0_sqrt), 
                                                            "noise" = noise, "method" = "RRR-not-highd", 
                                                            "prop_missing" = prop_missing,
                                                             "overlapping_amount" = overlapping_amount,
                                                    "nnzeros" = nnzeros,
                                                    "theta_strength" = strength_theta,
                                                    "n" = n,
                                                    "r_pca" = r_pca,
                                                    "exp" = seed * 100 + seed_n,
                                                    "normalize_diagonal" = normalize_diagonal,
                                                    "lambda_opt" = 0,
                                                    "time" = start_time_rrr[[1]]))
                }
                

                #### Oracle
                print("beginning oracle")
                set_u =  which(apply(gen$u,1, norm)>0)
                set_v =  which(apply(gen$v,1, norm)>0)
                t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                Uhat = matrix(0, p, r)
                Vhat = matrix(0, q, r)
                Uhat[set_u, ] <-  t$xcoef[, 1:r]
                Vhat[set_v, ] <-  t$ycoef[, 1:r]
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, Uhat, 
                                                            Vhat, 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  method = "Oracle",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "exp" = seed * 100 + seed_n,
                                                   "lambda_opt" = 0,
                                                   "time" = 0
                                                   
                )
                )
                
                print(paste0("Starting ", "Alt opt") )
                tryCatch({
                  start_time_alt3 <- system.time({
                    res_alt = CCA_rrr.CV(X, Y,
                                         r=r, Kx = NULL, lambda_Kx = 0,
                                         param_lambda=c(10^seq(-3, 1, length.out = 30)),
                                         kfolds=5, solver="ADMM", LW_Sy = LW_Sy, 
                                         do.scale = TRUE,
                                         rho=1, niter=2 * 1e4, thresh = 1e-6)
                  })
                  res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                  res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                  Uhat <- res_alt$ufinal[, 1:r]
                  Vhat <- res_alt$vfinal[, 1:r]
                  lambda_chosen = res_alt$lambda
                  if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                    #### Choose another lambda
                    while(sum(apply(Uhat, 1, function(x){sum(x!=0)}) >0) <r){
                      lambda_chosen = lambda_chosen / 2
                      #start_time_alt <- system.time({
                      res_alt <- CCA_rrr(X, Y, Sx = NULL, Sy=NULL,
                                         lambda =lambda_chosen, Kx=NULL, 
                                         r=r, 
                                         highdim=TRUE,
                                         solver="ADMM",
                                         LW_Sy = LW_Sy, do.scale=TRUE,
                                         thresh = 1e-6)
                      res_alt$U[which(is.na(res_alt$U))] <- 0
                      res_alt$V[which(is.na(res_alt$v))] <- 0
                      Uhat <- res_alt$U[, 1:r]
                      Vhat <- res_alt$V[, 1:r]
                      
                      #})
                      
                    }
                  }
                  result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew,
                                                              Uhat,
                                                              Vhat[, 1:r],
                                                              gen$u, gen$v,
                                                              Sigma_hat_sqrt = Sigma_hat_sqrt,
                                                              Sigma0_sqrt = Sigma0_sqrt),
                                                     "noise" = noise,
                                                     method = "RRR-ADMM-opt",
                                                     "prop_missing" = prop_missing,
                                                     "nnzeros" = nnzeros,
                                                     "theta_strength" = strength_theta,
                                                     "overlapping_amount" = overlapping_amount,
                                                     "r_pca" = r_pca,
                                                     "n" = n,
                                                     "exp" = seed * 100 + seed_n,
                                                     "normalize_diagonal" = normalize_diagonal,
                                                     "lambda_opt" = lambda_chosen,
                                                     "time" = start_time_alt3[[1]]
                  )
                  )
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
                
                for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                                 "SCCA_Parkhomenko", "Fantope", "Chao", "SGCA")){
                  
                  print(paste0("Starting ", method))
                  
                  
                  tryCatch({
                    start_time_additional_method <- system.time({
                      test1<-additional_checks(gen$X,
                                               gen$Y, S=NULL, 
                                               rank=r, kfolds=5, 
                                               method.type = method,
                                               lambdax= 10^seq(-3,1, length.out = 30),
                                               lambday = c(0))
                    })
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                test1$u[, 1:r], 
                                                                test1$v[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  method = method,  
                                                       "prop_missing" = prop_missing, 
                                                       "overlapping_amount" = overlapping_amount,
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = 0,
                                                       "time" = start_time_additional_method[[1]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                }
                
                write_csv(result, paste0("experiments/simulations/results/2024_2_newest_RRR_efficient_results", name_exp, ".csv"))
              
              }
            }
          }
       # }
      }
      #}
    }
  }
#}
}
}
