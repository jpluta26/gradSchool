# John Pluta
# BSTA670 - Stat Computing
# Project 2

# 12/4/14

# ------------------------------- preprocessor ------------------------------ #
rm(list=ls())

# code taken from
# http://stackoverflow.com/questions/8475102/set-default-cran-mirror-permanent-in-r
# set a CRAN mirror
local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.stat.ucla.edu"
       options(repos=r)})


# install necessary packages
if("lars" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("lars") }

library(lars)
library(MASS)
# ---------------------------------------------------------------------------- #


# ============================= begin functions ============================== #
# ---------------------------------------------------------------------------- #
# function get.n.max
# input: a design matrix X
# output: the maximum norm, as defined in Section 2, Def. 3 and Section 3.2
get.n.max <- function(X)
{
  X1 <- X[,1:q]
  X2 <- X[,(q+1):p]
  
  C11 <- (1/n) * t(X1) %*% X1
  C21 <- (1/n) * t(X2) %*% X1
  
  # take maximum norm
  n.max <- 1 - max(abs(C21 %*% solve(C11) %*% sign(B1)))
  return(n.max)
}
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# function run.sim
# for given design matrix X, this function runs a simulation n.iter times
# the simulation is:
# 1) given X and B, generate eps from N(0,sigma^2), and generate Y from Y = XB + e
# 2) run the lars algorithm to get the full coefficient path for each coefficient
# 3) examine the entire coefficient path to see if there exists a lambda for which
# lars selects the correct variable (e.g., all variables in q are non-zero, all others
# are 0), and if the selected variables are sign consistent with their true values.
# this indicates a "correct" model.
# 4) repeat step 3 n.iter times, and record the proportion of correct models found.
# return 4
run.sim <- function(X, n.iter)
{
  # keep a record of whether LARs was able to create a valid model for the generated data
  match.record <- error.record <- rep(0, n.iter)
  
  for(iter in 1:n.iter)
  {
    # randomly generate errors from N(0,sigma^2)
    set.seed(iter+100)
    eps <- rnorm(n, mean=0, sd=sigma)
    Y <- X %*% B + eps
    
    # compute coefficient path
    # data is already normalized and we don't want an intercept (mean is 0)
    out <- lars(X,Y,type="lasso", normalize=FALSE, intercept=FALSE)
    
    
    
    match.record[iter] <- check.consist2( out$beta, B )
    error.record[iter] <- get.error( out$beta)
    
  }
  
  # number of successfully matched models over the total number of models
  per.match <- sum(match.record) / length(match.record)
  avg.error <- mean(error.record)
  return(list(per.match, avg.error))
}

# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# function get.error
# input: solution path from LARS output
# output: the minimum error from the solution path
get.error <- function( beta.mat )
{
  
  end.step <- dim(beta.mat)[1]
  err <- rep(0, end.step)
  
  for( i in 1:end.step)
  {
    # at each step, find the number of mismatched signs
    err[i] <- sum(sign(beta.mat[i,]) != sign(B))
  
  }
  
  # return the minimum error; e.g. the best fitting solution
  return(min(err))
}
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# function check.consist - checks output against true beta for sign consistency
# input: beta-weights output from the LARS solution to the model; this is the entire
# coefficient path of all the beta weights
#
# output: matched, a boolean variable that is TRUE if there exists a solution
# anywhere in the path that is sign consistent with the true beta parameters
check.consist2 <- function( beta, B )
{
  # number of steps in the LARS solution path
  end.step <- dim(beta)[1]
  matched <- FALSE
  
  # iterate through each step
  for(i in 1:end.step)
  {
    if( all(sign(beta[i,]) == sign(B)))
    { 
      matched <- TRUE 
      break
    } 
  }
  
  return(matched)
}
# ---------------------------------------------------------------------------- #
# ============================= end functions ================================ #







# ----------------------------------- main ----------------------------------------- #
# set global constants- all defined in the paper
# set n.iter to 10 for testing, 1000 to replicate results
n.iter <- 1000  # number of times the simulation is run
n <- 100   # observations per variable
p <- 32    # total number of variables
q <- 5     # number of non-zero covariates
B1 <- c(7, 4, 2, 1, 1)   # beta weights
B2 <- rep(0, p - q)      # 
B <- c(B1, B2)
sigma <- sqrt(.1)        # standard deviation

# create the covariance matrix for the wishart distribution
set.seed(12)
mat <- matrix(rnorm(p*p, mean=0, sd=sigma), nrow=p, ncol=p)
mat <- (1/p)*t(mat) %*% mat

# only need one Wishart matrix
S <- rWishart(1, p, mat )
S <- S[,,1]

# variables to store output
# create 100 nxp matrices, randomly sampled from N(0,S)
X <- array(0, dim=c(n,p,100))
n.max <- rep(0,100)
per.match <- avg.error <- rep(0,100)

looptime <- proc.time()

# this part takes a while
for(i in 1:100)
{
  print(paste0("begin iteration: ", i))
  # randomly generate 100 design matrices
  set.seed(i+100)
  X[,,i] <- mvrnorm(n=n, mu=rep(0,p), Sigma = S)
  X.c <- X[,,i]
  
  # standardize the variances; typical step in LASSO
  for(j in 1:p)
  {
    X.c[,j] <- X.c[,j] / (sqrt((1/n) * t(X.c[,j]) %*% X.c[,j]))
  }
  
  
  # get the magnitude of irrepresentable condition for each matrix
  n.max[i] <- get.n.max(X.c)
  
  # get the percentage of simulations where lasso models successfully matched signs
  # and the average error
  simout <- run.sim(X.c, n.iter)
  per.match[i] <- simout[[1]]
  avg.error[i] <- simout[[2]]
  
  
}
print("done!")
proc.time() - looptime

# approximate time for a single loop; total elapsed time divided by number of loops
time <- looptime[3][[1]]/100


plot(n.max, per.match, type="p", ylim=c(0,1), cex=.5, pch=20)
abline(v=0, lty=2)

plot(n.max, avg.error, type="p", ylim=c(0,2), cex=.5, pch=20)
abline(v=0, lty=2)

# ---------------------------------------------------------------------------- #