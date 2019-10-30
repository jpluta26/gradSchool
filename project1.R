# John Pluta
# Project 1

# --- preprocessor --------------- #
# load necessary libraries, prepare data
rm(list=ls())

# install necessary packages
if("lars" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("lars") }

if("doParallel" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("doParallel") }

if("foreach" %in% rownames(installed.packages()) == FALSE)
{ install.packages("foreach") }

if("lassoshooting" %in% rownames(installed.packages()) == FALSE)
{ install.packages("lassoshooting") }

library(parallel)
library(doParallel)
library(lassoshooting)

# lars contains the diabetes dataset
library(lars)


# ----
# user-defined constants:

# lambda is chosen by some other method, usually cross-validation. this document does not
# contain any method for choosing lambda as it is beyond the scope of the project.  
# note that Hebiri & Lederer, 2013, showed that the size of the tuning parameter is 
# inversely proportional to the degree of correlation. since the simulated data is quite
# highly correlated, choose lambda < 1 (> 1 causes problems with convergence)
lambda <- .2

# epsilon is some arbitrarily small number used to compare differences in numbers. the "==" operator
# may be unreliable due to floating point errors. so instead of comparing x == y, we compare x - y < eps.
# eps >= machine precision
eps <- 1e-6

# maximum iterations for all of the algorithms
max.iter <- 10000


# --- end preprocessor -------------------- #



# ========== functions ======================= #
# ----------------------------------------------------------------------------- #
# function to create data for lasso experiments. creates an nxp matrix of data,
# with correlation induced between each of the columns.
#
# input- dimensions n and p, where n is the number of subjects and p is the number of
# covariates. two columns are created for each p, so the input should be half the desired number
# of columns.
#
# output- a matrix of simulated data, from a N(0,1) distribution with correlation induced
# 
# this function is only for simulating data and is NOT optimized. only the algorithms have been
# optimized.
#
# citation: some of this code is taken from http://comisef.wikidot.com/tutorial:correlation
create.lasso.data <- function(n,p)
{
  # variable to store the data
  Y <- NULL
  
 
  p <- as.integer(p)
  
  
  
  for(i in 1:p)
  { 
    
    # drawing exactly 1 or 0 causes a crash
    rho <- sample(seq(from=0.01, to=.875, by=0.1), 1)
    
    # create a pair of vectors, each randomly drawn from N(0,1)
    # the two vectors have correlation equal to rho
    set.seed(i)
    X <- rnorm(n * 2, mean=0, sd=1) * 0.05
    dim(X) <- c(n, 2)
    
    M <- array(rho, dim=c(2, 2))
    diag(M) <- 1
    
    cF <- chol(M)
    
    Y <- cbind(Y, X %*% cF)
  }
  
  # also center the data for LASSO
  
  return(Y)
}
# ----------------------------------------------------------------------------- #



# ----------------------------------------------------------------------------- #
# function to implement soft thresholding as described in Fu, 1998, and Friedman et al., 2007
# input:
# beta- a vector of coefficient estimates
# lambda- user-defined lambda value
#
# output:
# beta- a vector of coefficient estimates with soft thresholding applied
soft.thresh <- function(beta, lambda, x)
{
  
  if( beta > 0 & lambda < abs(beta))
  {
    beta  <- (beta - lambda) / ( t(x) %*% x )
  } else if( beta < 0 & lambda < abs(beta))
  {    
    beta <- (beta + lambda) / (  t(x) %*% x)
  } else if( lambda > abs(beta) )
  {
    beta <- 0
  }
  
  return(beta)
}
# ----------------------------------------------------------------------------- #




# ----------------------------------------------------------------------------- #
# this function implements the coordinate descent update. its identical to the update
# in the coordinate descent function, but stripped down for use in parallel processing
#
# input: 
# j, the index indicating which variable to update
# X, the design matrix
# y, the response data
# Xty, X'y precomputed for efficiency
# lambda, the lambda value
# beta, vector previous beta value
#
# output: 
# beta, updated vector of coefficient values
update <- function(j,X,Xty,lambda,beta)
{
 
    beta[j] <- Xty[j] - t(X[,j]) %*% X[,-j] %*% beta[-j]
    beta[j] <- soft.thresh(beta[j], lambda, X[,j])   
    return(beta)
}
# ----------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------- #
# function to implement the shotgun parallel stochastic coordinate descent algorithm,
# as described in Bradley et al., 2011, which in turn is a parallel version of
# stochastic coordinate descent, as described in Shwartz & Tewari 2011. 
# similar to coordinate descent but the coefficient to be updated is chosen at random
#
# input: j, the index of the coefficient to update
# X, matrix of data
# y, vector of outcome data
# Xty, precomputed X'y
# lambda, user-defined lambda value
# beta, current coefficient estimates
#
# output: beta, updated coefficient estimates
#
# note that in a regular stochastic coordinate descent, j would be randomly chosen
# inside this function. but via the shotgun method, j is chosen beforehand, and then
# each instance of scd is sent to a separate core

#
# note: Bradley et al. state that "although our analysis uses duplicate features, they
# are not needed for implementation." This refers to the indexing of 2d features, rather
# than just d, e.g. a negative and positive value of each coefficeint. this function 
# does NOT impliment duplicate features.
#
# convergence in this method is difficult to test! see code for details. this is the only
# convergence scheme i could come up with that achieves identical results to coordinate
# descent
shotgun <- function(X,y,lambda,beta)
{ 
  change.count <- 0
  n.cores <- detectCores()
  
  # the maxmium number of cores used for parallel stochastic coordinate descent depends
  # on the problem. Bradley et al. state that for rho = 0 (no correlation between features)
  # we can use p cores, whereas for rho = 1, we can only use 1. this code has only been 
  # tested on 2 cores, so 2 is enforced as a cap. change max.cores to experiment with
  # more, but it is not guaranteed to work. 4 cores is probably a sensible estimate for a 
  # maximum.
  
  # only tested on 2 cores, but should work for 4. anything greater probably violates
  # the assumptions of the problem and may cause slow down or divergence.
  max.cores <- 2
  
  if( n.cores >= 2)
  {
    n.cores <- max.cores
  }
  
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)

  # X'y is used in each calculation but is a constant. so precompute the value and just call it from memory
  Xty <- t(X) %*% y
  p <- dim(X)[2]
  
  change.count <- 0
  for(i in 1:max.iter)
  {
    # store prior results of beta to compare for convergence
    beta.prev <- beta
   
    j.list <- sample(1:p,n.cores,replace=F)
    
    
    # run the scd algorithm for each covariate, one per core
    res <- foreach(j=j.list, .export=c('update','soft.thresh', 'eps', 'max.iter')) %dopar% update(j,X,Xty,lambda,beta)
  
  
    # update beta with the results
    beta[j.list[1]] <- res[[1]][j.list[1]]
    beta[j.list[2]] <- res[[2]][j.list[2]]
  
    # for highly correlated data, LASSO will only pick a few variables and set the rest to 0.
    # as a result, its quite likely that the first two parameters randomly chosen will be set to 0
    # and trigger the convergence flag prematurely. to prevent this, don't start looking at convergence
    # criterion until at least one beta is non-zero.
    if(!all(beta == 0))
    {
      
      # has the number of betas changed?
      if( sum(beta.prev) == sum(beta))
      {
        # have the beta values changed?
        if( max(abs(beta.prev - beta)) < eps)
        {
          change.count <- change.count + 1
          
        } else( change.count <- 0)
      }
    }
    
    
    # but just because one beta weight hasn't changed doesn't mean the full model has been reached.
    # so we need to check for stability of the beta estimates as well; essentially if the number of 
    # non-zero betas hasn't changed in a certain number of iterations, AND the values of those betas hasnt
    # changed, then we have reached convergence.
    # the threshold for change.count is somewhat arbitrary, but 2*p allows us to (probably) cycle through
    # all parameters twice
    if( change.count > p )
    {
      cat("SHOTGUN: Convergence reached at iteration ", i, "\n")
      break
     }
    
  
    if( i == max.iter)
    {
      # could put stop here instead of break, but it is useful for debugging to return the beta
      # values and see whats going on
      
      cat("SHOTGUN: Maximum number of iterations reached without convergence\n")
      break
    }
  }

  stopCluster(cl)
  
  return(beta)
}
# ----------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------- #
# function to compute solution to the LASSO problem using gradient descent algorithm
# input: 
# X- a matrix of the independent variables, where each column is the data of one variable
# y- a vector of the dependent variable
# lambda- tuning parameter; user defined
# beta- initial parameters for beta estimates; usually all 0s
#
# output:
# beta- a px1 vector of beta weights. these are the coefficient estimates that minimize the
# LASSO problem
#
# note: equations refer to the implementation in Friedman et al., 2010
coord.desc <- function(X,y,lambda,beta)
{
  
  p <- dim(X)[2]
  
  # X'y is constant, so we only need to calculate it once and then call the value from memory
  Xty <- t(X) %*% y
  
  for(i in 1:max.iter)
  {
    
    beta.prev <- beta
    for(j in 1:p)
    {
      # eq 5
      beta[j] <- Xty[j] - t(X[,j]) %*% X[,-j] %*% beta[-j]
      
      # eq 6
      beta[j] <- soft.thresh(beta[j], lambda,X[,j]) 
    }
    
    # check for convergence
    if( max(abs(beta - beta.prev)) < eps)
    {
      cat("COORD: Convergence reached at iteration ", i, "\n")
      break
    } 
    
    
    
    
    # if convergence hasn't been reached after maximum iterations, something went wrong
    if( i == max.iter)
    {
      # could put a stop here but its useful for debugging to return the beta values and see whats going on
      cat("COORD: Maximum number of iterations reached without convergence")
      break
    }
    
    # if coordinate descent diverges, it will keep growing ; cut this off at some unrealistic value for beta
    if( any(beta > 100000))
    {
      cat("COORD: Beta estimates have diverged")
      break
    }
  }
  
  return(beta)
}
# ----------------------------------------------------------------------------- #





# ----------------------------------------------------------------------------- #
# function to implement the LARS algorithm with LASSO modification, as described in
# Efron et al., 2004. all equations refer to that paper.
#
# input:
# X, the matrix of independent variables
# y, the vector of outcomes
# beta, the vector of initial coefficient estimates
#
# output:
#  a list of matrices for the beta, the active set, C, and c.
get.lars <- function(X,y,beta)
{

  p <- dim(X)[2]
  
  C <- rep(0,p)
 
  
  
  # store the record of how lars proceeded to the solution. the record matrices will be created 
  # as large, max.lars.iter x p matrices and then trimmed down to only the necessary output
  active.record <- c.record <- beta.record <- matrix(0, nrow=max.iter, ncol=p)
 
  # the matrix of proportional correlations
  c <- t(X) %*% y  # eq 2.8 (with beta=0, X'B = mu = 0)

  # flag to see if we are on the last step (active set is full)
  LAST = FALSE

  # number of iterations run (to trim the record variables)
  n.iter <- 0
  
  
  
  for(k in 1:max.iter)
  {
    
    n.iter <- n.iter + 1
    C[k] <- max(abs(c)) # eq 2.9
  
  
    # using the '==' operator to compare abs(-x) == x can sometimes return an incorrect
    # answer, due to round-off error (or some other floating point error).
    # instead, an error tolerance term is used. 
  
    active <- C[k] - abs(c) < eps  # eq 2.9 (modified to add error tolerance)
    
   
    active.record[k,] <- active

    # if every element is true, then all of the covariates are in the active set
    # in which case 2.13 is undefined. need to account for this
    if( sum(active) == p ) { LAST <- TRUE }
  
    # the sign of each element of c
    s <- sign(c[active])  # eq 2.10
  
    # if there is only one variable in the active set, X[,active] will be enconded as a vector
    # and sweep will crash- so enforce that it is defined as a matrix
    # apply each element of s to each corresponding column of X 
    X.A <-  sweep(as.matrix(X[,active]), MARGIN=2, s, `*`)  # 2.4
  
   if(is.null(X.A)) { cat("X.A is null!")}
    G.A.inv <- solve(t(X.A) %*% X.A)  # eq 2.5
 
    # 1' * G.A.inv * 1 is a summation of the elements of G.A.inv
    A.A <- 1 / sqrt( sum(G.A.inv) )  # eq 2.5
  
    # G.A.inv * 1 is the sum of the rows of G.A.inv
    G.A.rsum <- as.matrix( apply(G.A.inv, 1, sum) )
    G.A.scaled <- A.A * G.A.rsum
 
    u.a <-  X.A %*% G.A.scaled # eq 2.6
  
    a <- t(X) %*% u.a # eq 2.11
  
    d <- s * G.A.scaled # eq 3.3
  
   # eq 2.13
   if( LAST == FALSE) {
    gamma <- c(
       (C[k] + c[!active])/(A.A + a[!active]), 
       (C[k] - c[!active])/(A.A - a[!active])      ) 
   
     gamma.hat <- min( gamma[gamma > 0]) } else { gamma.hat <- C[k]/A.A }
  
   
   
   
   REMOVE <- FALSE
   gamma.set <- -s * beta[active] / G.A.scaled # 3.4
   
   # if there is a sign change in beta, pick the smallest positive gamma that causes the change
   # and assign that to gamma.tilda. by convention, gamma.tilda is assigned infinity when there is no gamma > 0
   if( length( gamma.set[gamma.set > 0] ) > 0 )
   { gamma.tilda <- min(gamma.set[gamma.set > 0], gamma.hat)  } else { gamma.tilda <- Inf } # eq 3.5
   
  
    # eq 3.6
    if( gamma.tilda < gamma.hat )
    {
      REMOVE <- TRUE
      LAST <- FALSE
    
      # also need to update values w/ gamma.tilda not gamma.hat
      gamma.hat <- gamma.tilda
    }
   
    # update coefficients
    beta[active] <- beta[active] + gamma.hat * d  # eq 3.3
  
    # update proportional correlation
    c <- c - gamma.hat * a  # eq 2.15
  
    # update absolute correlation while the active set is not full
    if( LAST == FALSE )
    { C[k+1] <- C[k] - gamma.hat } # 2.16
  
    c.record[k,] <- c
    beta.record[k,] <- beta
    if(REMOVE == TRUE)
    {
      
      index <- which(gamma.tilda - gamma.set < eps)
      active[index] <- FALSE   # 3.6
      
    }
  
    # if the entire set is active, we have reached convergence. end the function.
    if( LAST == TRUE )
    {
      cat("LARS: Convergence reached at iteration", k, "\n")
      return(list(beta.record[1:n.iter,], c.record[1:n.iter,], active.record[1:n.iter,], C))
    } else
      if(LAST == FALSE & k==max.iter) 
      { 
        cat("LARS: Convergence not reached\n") 
        return(list(beta.record[1:n.iter,], c.record[1:n.iter,], active.record[1:n.iter,], C)) 
      }
  }

  
  
}
# ----------------------------------------------------------------------------- #
# ========== end functions ======================= #





# ========== main ============================== #
# first, run analysis on the diabetes dataset; extensively studied in the literature,
# so we can compare coefficient estimates to real results and make sure the implemented 
# algorithms give sensible output
data(diabetes)



# data must be centered for LASSO
y <- as.matrix(diabetes$y - mean(diabetes$y))
X <- as.matrix(diabetes$x - colMeans(diabetes$x))

# initialize betas to 0
# Fu, 1998 suggest giving beta a warm start, e.g. initializing beta values as the
# OLS estimates. However, the gains in this approach are computationally trivial
# compared to beta = 0, and warm start will fail if p > n
# if a warm start is desired, it would be implemented as:
# fit <- lm(y ~ X)
# beta <- as.matrix(fit$coef[2:(p+1)])
# out.warm.start <- grad.desc(X,y,lambda, beta)
p <- dim(X)[2]
beta <- rep(0,p)

# test results on the diabetes dataset so we have know the algorithms are working properly


cat("Estimating solution for diabetes data...\n\n")

# all three of these implementations should reach the same solution
out <- lassoshooting(X,y,lambda)

coord.time.d <- proc.time()
out.coord <- coord.desc(X,y,lambda, beta)
proc.time() - coord.time.d

shot.time.d <- proc.time()
out.shotgun <- shotgun(X,y,lambda, beta)
proc.time() - shot.time.d

# out.lars and real.out.lars should reach the same solution
real.time.d <- proc.time()
real.out.lars <- lars(X,y,type="lasso")
proc.time() - real.time.d

lars.time.d <- proc.time()
out.lars <- get.lars(X,y,beta)
proc.time() - lars.time.d

cat("Estimation complete.\n\n")

# now run the experiment
# this part of the code will take the longest; reduce max.size to 50 for quick testing
# specifies matrix sizes for simulated data
max.size <- 500
size <- seq(from=50, to=max.size, by=50)

lars.time <- coord.time <- shot.time <- real.lars.time <- rep(0, length(size))
i <- 0

cat("Estimating solution for simulated data...\n\n")
for(n in size)
{

  i <- i + 1
  
  # the second argument for this function will create 2n columns of data
  # so specify half the amount you want. this creates a square matrix.
  
  X <- create.lasso.data(n, n/2)

  p <- dim(X)[2]
  
  beta <- rep(0,p)
  y <- matrix(rnorm(n, mean=0, sd=1), nrow=n, ncol=1)
  
  # save the run time of each algorithm
  lars.time[i] <- system.time( get.lars(X,y,beta))[3]
  real.lars.time[i] <- system.time( lars(X,y,type="lasso"))[3]
  coord.time[i] <- system.time( coord.desc(X,y,lambda,beta))[3]
  shot.time[i] <- system.time( shotgun(X,y,lambda,beta))[3]
}

cat("Estimation complete.\n\n")
# create figure plotting run time of each method as a function of data size

plot(size,lars.time, type="n", ylim=c(0, max(lars.time, coord.time, shot.time) ), xlab="Size", ylab="Time")
points(size, lars.time, type="l", col="blue")
points(size, coord.time, type="l", col="black")
points(size, shot.time, type="l", col="red")
legend("topleft", inset=c(.01,.01), c("LARS", "Shooting", "Shotgun"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),
       col=c("blue", "black", "red"))

matlab.time <- c(0.0245,0.0787, 0.3008,0.4341,0.6719, 3.0569, 5.3280, 6.5276,7.3080, 9.0001)
plot(size, coord.time, type="n", ylim=c(0, max(matlab.time, coord.time)), xlab="Size", ylab="Time")
points(size, coord.time, type="l", col="blue")
points(size, matlab.time, type="l", col="red")
legend("topleft", inset=c(.01, .01), c("R", "Matlab"), lty=c(1,1), lwd=c(2.5, 2.5), col=c("blue", "red"))

# create figure of my lars function versus real lars function
# ran out of time to keep optimizing...
#plot(size,lars.time, type="n", ylim=c(0, max(lars.time, real.lars.time)), xlab="Size", ylab="Time")
#points(size, lars.time, type="l", col="red")
#points(size, real.lars.time, type="l", col="blue")
#legend("topleft", inset=c(.01,.01), c("Custom LARS", "Real LARS"), lty=c(1,1,1), lwd=c(2.5,2.5,2.5),
       col=c("red", "blue") )


# create the plot of coefficient paths from LARS output for diabetes data
# beta values plotted over increasing lambda
beta <- out.lars[[1]]
n <- dim(beta)[1]
x <- y <- rep(0,n)
for(i in 1:n)
{
  x[i] <- sum(abs(beta[i,]))
}

png("coefpath.png")
plot(x, beta[,1], type="n", ylim=c(min(beta), max(beta)), xlab="|Bj|", ylab="Coefficients")
for(i in 1:10)
{
  points(x, beta[,i], type="l", col=i)
}
dev.off()
# ========== end main ============================== #
