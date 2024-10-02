library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(rrpack)
library(corpcor)


source("experiments/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("src/evaluation.R")
source("src/gradient_descent.r")
source("src/reduced_rank_regression.R")
source("src/graph_reduced_rank_regression.R")


args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
n <- 500
strength_theta <- args[3]
rs <- c(as.numeric(args[4]))
p_val <- as.numeric(args[5])
q_val <- as.numeric(args[6])
overlaps <- c(1)
props <- c(0)
noise = 1
seeds = 1:100
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(5)

result = c()
for(seed_n in seeds){
  set.seed(seed * 100 + seed_n)
  print("Start loop")
  for(nnzeros in nnzero_values){
    for(p in p_val){
      for (q in q_vals){
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
            if ( (max(r_pca, r, nnzeros) < p) ) {
              for (overlapping_amount in overlaps){
                for(prop_missing in props){
                  cat("seed:")
                  cat(seed, " ")
                  print(c(n, r, r_pca, strength_theta))
                  start_time_creation <- system.time({
                    gen = generate_example_graph(n,
                                                 p1=p,
                                                 type_graph="2d-grid",
                                                 p2=q, order = 3,
                                                 r_pca = r_pca,
                                                 nnzeros = nnzeros,
                                                 do_plot = FALSE,
                                                 theta = thetas,
                                                 lambda_pca = 1,
                                                 nnzeros_pca = 20,
                                                 r = r, 
                                                 overlapping_amount = overlapping_amount,
                                                 normalize_diagonal = normalize_diagonal,
                                                 gen.using.gamma = TRUE,
                                                 n_new = 5000)
                  })
                  print("Generation time")
                  print(start_time_creation[[1]])
                  print("We're done generating")
                  X = gen$X
                  Y = gen$Y
                  groups = gen$groups
                  Sigma0_sqrt = sqrtm(gen$Sigma)$B
                  Sigma_hat_sqrt = sqrtm(gen$S)$B
                  Gamma_dagger = pinv(gen$Gamma)
                  
                  #### Try out alternative approaches
                  #### Oracle
                  print("beginning oracle")
                  tryCatch({
                  set_u =  which(apply(gen$u,1, norm)>0)
                  set_v =  which(apply(gen$v,1, norm)>0)
                  t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                  Uhat = matrix(0, dim(gen$u)[1], r)
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
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
                  
                  print(paste0("Starting ", "Alt opt") )
                  
                  tryCatch({
                    start_time_alt2 <- system.time({
                      res_alt = CCA_rrr.CV(X, Y, 
                                           r=r, Kx = NULL, lambda_Kx = 0,
                                           do.scale = TRUE,
                                           param_lambda=c(10^seq(-3, 1, length.out = 30)),
                                           kfolds=3, solver="rrr", LW_Sy = LW_Sy)
                    })
                    res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                    res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                    Uhat <- res_alt$ufinal[, 1:r]
                    Vhat <- res_alt$vfinal[, 1:r]
                    if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                      #### Choose another lambda
                      lambda_chosen = max(res_alt$resultsx$lambda[which(res_alt$resultsx$rmse > 1.05 * min(res_alt$resultsx$rmse))])
                      start_time_alt <- system.time({
                        res_alt <- CCA_rrr(X, Y, Sx = NULL, Sy=NULL,
                                           lambda =lambda_chosen, Kx=NULL, r, highdim=TRUE,
                                           solver="rrr",
                                           LW_Sy = LW_Sy, do.scale=TRUE)
                        res_alt$U[which(is.na(res_alt$U))] <- 0
                        res_alt$V[which(is.na(res_alt$v))] <- 0
                        Uhat <- res_alt$U[, 1:r]
                        Vhat <- res_alt$V[, 1:r]
                        
                      })
                    }
                    print(res_alt$ufinal)
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                Uhat, 
                                                                res_alt$vfinal[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  method = "RRR-opt",  
                                                       "prop_missing" = prop_missing, 
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "overlapping_amount" = overlapping_amount,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = res_alt$lambda,
                                                       "time" = start_time_alt2[[4]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in RRR CV:", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                  
                  
                  print(paste0("Starting ", "Alt opt2") )
                  tryCatch({
                    start_time_alt3 <- system.time({
                      res_alt = CCA_rrr.CV(X, Y,
                                           r=r, Kx = NULL, lambda_Kx = 0,
                                           param_lambda=c(10^seq(-3, 1, length.out = 30)),
                                           kfolds=3, solver="ADMM", LW_Sy = LW_Sy, do.scale = TRUE,
                                           rho=1, niter=3 * 1e4, thresh = 1e-6)
                    })
                    res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                    res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                    Uhat <- res_alt$ufinal[, 1:r]
                    Vhat <- res_alt$vfinal[, 1:r]
                    
                    
                    print(res_alt$ufinal)
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew,
                                                                Uhat,
                                                                res_alt$vfinal[, 1:r],
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
                                                       "lambda_opt" = res_alt$lambda,
                                                       "time" = start_time_alt3[[4]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                  
                  
                  print(paste0("Starting ", "CV opt-graph") )
                  
                  
                 
                  tryCatch({
                    start_time_alt4<- system.time({
                      res_alt = CCA_graph_rrr.CV(X, Y, 
                                                 Gamma = gen$Gamma,
                                                 r=r, Kx = NULL, lambda_Kx = 0,
                                                 param_lambda=c(10^seq(-3, 1, length.out = 30)),
                                                 kfolds=5, 
                                                 LW_Sy = LW_Sy, do.scale=TRUE,
                                                 rho=1, niter=2 * 1e4,
                                                 thresh=1e-4,
                                                 Gamma_dagger = Gamma_dagger)
                    })
                    res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                    res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                    Uhat <- res_alt$ufinal[, 1:r]
                    Vhat <- res_alt$vfinal[, 1:r]
                    
                    print(res_alt$ufinal)
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                Uhat, 
                                                                res_alt$vfinal[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  
                                                       method = "CVX-opt-graph",  
                                                       "prop_missing" = prop_missing, 
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "overlapping_amount" = overlapping_amount,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = res_alt$lambda,
                                                       "time" = start_time_alt4[[4]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in group CV", ":", conditionMessage(e), "\n")
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
                                                 rank=r, kfolds=3, 
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
                                                         "time" = start_time_additional_method[[4]]
                      )
                      )
                    }, error = function(e) {
                      # Print the error message
                      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                      # Skip to the next iteration
                    })
                  }
                  write_csv(result, paste0("experiments/simulations/results/graph/2024-newest_graph_RRR_efficient_results", name_exp, ".csv"))
                  print("Done loop")
      
                }
                print("Done inner most loop")
              }
            }
          }
        }
      }
    }
    
    

  }
}



