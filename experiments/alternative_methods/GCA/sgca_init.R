# Function for Intialization via generalized fantope projection
# Inputs:
# =======
# A, B:       Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# nu:         Parameter of ADMM, default set to 1
# K:          nuclear norm constraint, equal to r
# rho:     penalty parameter on the l_1 norm of the solution, scaled by
#             sqrt(log(max(p1,p2))/n)
# epsilon:    tolerance level for convergence in ADMM
# maxiter:    maximum number of iterations in ADMM
# trace:      if set to True will print all iterations 

# Outputs:
# ========
# $Pi:     optimum of the convex program

sgca_init <-
  function(A,B,rho,K,nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE){
    p <- nrow(B)
    eigenB <- eigen(B)
    sqB <- eigenB$vectors%*%sqrt(diag(pmax(eigenB$values,0)))%*%t(eigenB$vectors)	
    tau <- 4*nu*eigenB$values[1]^2	
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    H <- Pi <- oldPi <-  diag(1,p,p)
    Gamma <- matrix(0,p,p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
      for (j in 1:20){
        Pi <- updatePi(B,sqB,A,H,Gamma,nu,rho,Pi,tau)
      }
      #Pi <- updatePi(B,sqB,A,H,Gamma,nu,lambda,Pi,tau)
      
      H <- updateH(sqB,Gamma,nu,Pi,K)
      Gamma <- Gamma + (sqB%*%Pi%*%sqB-H) * nu	
      criteria <- sqrt(sum((Pi-oldPi)^2))
      oldPi <- Pi
      i <- i+1
      if(trace==TRUE)
      {
        print(i)
        print(criteria)
      }
    }
    return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
    
  }



soft_threshold <- function(A, lambda){
  return( A * ((abs(A) >lambda) * (1-lambda* sign(A))))
}

updateH <-
  function(sqB,Gamma,nu,Pi,K){
    
    temp <- 1/nu * Gamma + sqB%*%Pi%*%sqB
    temp <- (temp+t(temp))/2
    svdtemp <- eigen(temp)
    d <- svdtemp$values
    p <- length(d)
    if(sum(pmin(1,pmax(d,0)))<=K){
      dfinal <- pmin(1,pmax(d,0))
      return(svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors))
    }
    fr <- function(x){
      sum(pmin(1,pmax(d-x,0)))
    }
    # Vincent Vu Fantope Projection
    knots <- unique(c((d-1),d))
    knots <- sort(knots,decreasing=TRUE)
    temp <- which(sapply(knots,fr)<=K)
    lentemp <- tail(temp,1)
    a=knots[lentemp]
    b=knots[lentemp+1]
    fa <- sum(pmin(pmax(d-a,0),1))
    fb <- sum(pmin(pmax(d-b,0),1))
    theta <- a+ (b-a)*(K-fa)/(fb-fa)
    dfinal <- pmin(1,pmax(d-theta,0))
    res <- svdtemp$vectors%*%diag(dfinal)%*%t(svdtemp$vectors)
    return(res)
  }


