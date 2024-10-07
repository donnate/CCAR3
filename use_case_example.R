library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(corpcor)
library(rrpack)

source("experiments/simulations/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/alternative_methods/scca_chao.R')
source("experiments/evaluation.R")
source("src/reduced_rank_regression.R")



seed = 123
print(seed)
set.seed(seed)
n <- 100
strength_theta <- "high"
p <- 20
q <- 5
rs <- 2
#props <- c(0, 0.1, 0.2)
props <- c(0)
prop_missing = 0
noise = 1
seeds = 1:50
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(10)
result = c()


cat("seed:")
cat(seed, " ")
r = 2
r_pca = 2
q =10
overlapping_amount= 1
normalize_diagonal = TRUE
print(c(n, r, r_pca, strength_theta))
nnzeros= 10
thetas <- diag(seq(0.9, 0.75, length.out = r))
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
print(start_time_creation[[3]])
print(system.time({
  Sys.sleep(3)
  }))
#### check 
print(t(gen$u) %*% gen$Sigmax  %*% gen$u)
print(t(gen$v) %*% gen$Sigmay  %*% gen$v)

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
                                    "time" = start_time_rrr[[3]]))
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


prop_missing = 0
seed_n = seed
result <- rbind(result, data.frame(evaluate(X=gen$Xnew, Y=gen$Ynew, 
                                            U=Uhat, 
                                            V=Vhat, 
                                            U0 = gen$u, V0 = gen$v,
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
                                     "time" = start_time_alt3[[3]]
  )
  )
}, error = function(e) {
  # Print the error message
  cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
  # Skip to the next iteration
})

for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                 "SCCA_Parkhomenko", "Fantope", "Chao", "SGGCA")){
  
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
                                       "time" = start_time_additional_method[[3]]
    )
    )
  }, error = function(e) {
    # Print the error message
    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
    # Skip to the next iteration
  })
}
