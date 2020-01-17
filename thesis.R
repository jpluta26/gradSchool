# pluta 9/15/15
# version 1.0


# ------------------- preprocessor ----------------------- #
# clear everything
rm(list=ls())

options(error=dump.frames)

local({r <- getOption("repos")
       r["CRAN"] <- "http://cran.stat.ucla.edu"
       options(repos=r)})


# install necessary packages
if("covTest" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("covTest") }

if("glmnet" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("glmnet") }

if("ROCR" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("ROCR") }

if("pROC" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("pROC") }

if("lars" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("lars") }


if("penalized" %in% rownames(installed.packages()) == FALSE) 
{ install.packages("penalized") }

# load packages
library(covTest)
library(glmnet)
library(ROCR)
library(pROC)
library(lars)
library(penalized)






# ------------------------------------------------------------------------ #
# ----------------------------- functions -------------------------------- #
# ------------------------------------------------------------------------ #

################################################
# attach variable  names to the output of the covariance test
# so its easier to read
covTest.vals <- function( var.list , test.stat)
{
  # create the full list of variables, e.g. CA1.l, CA2.r, etc
  var.list.full <- c()
  for(side in c("l", "r"))
  {
    for(sf in var.list)
    {
      name <- paste(sf, side, sep=".")
      var.list.full <- c(var.list.full, name)
    }
  }
  
  # put the varnames in the order of the output of the covTest,
  # and attach
  order.sf <- c()
  for(i in 1:length(var.list.full))
  {
    index <- test.stat$results[i,1]
    order.sf <- c(order.sf, var.list.full[index])
  }
  p.vals <- as.matrix(test.stat$results[,3])
  p.vals <- t(p.vals)
  colnames(p.vals) <- order.sf
  
  return(p.vals)
}
################################################


################################################
plotCoefPath <- function( out )
  # input: out is the output of either lars or glmnet, it contains the full coefficient path for the LASSO
  # solution
  # output: a plot of the path for each coefficient
{
  beta <- out$beta
  n <- dim(beta)[1]
  p <- dim(beta)[2]
  x <- rep(0,n)
  for(i in 1:n)
  {
    x[i] <- sum(abs(beta[i,]))
  }
  
  plot(x, beta[,1], type="n", ylim=c(min(beta), max(beta)), xlab="|Bj|", ylab="Coefficients", main="LARS Coefficient Path")
  for(i in 1:p)
  {
    points(x, beta[,i], type="l", col=i)
  }
}
################################################


################################################
plot.qq <- function( dat, type )
  # function to create a qq-plot of the input and test for a normal distribution
  # using the shapiro-wilk test.
{
  
  # make qq plot of the data
  qqnorm(dat,  main=sprintf(paste("QQ-Plot of %s"), type))
  qqline(dat)
  
  # get p-value for shapiro-wilke test and paste to the figure
  p <- shapiro.test(dat)$p.value
  mtext(  sprintf(  paste("Shaprio-Wilk p = %5.3f"), p), cex=1.2)
  
}
################################################


################################################
# make a nice plot of age distribution
agePlot <- function( ctl.age , mci.age , type)
{
  
  age.eq <- t.test(mci.age, ctl.age)
  plot(c(1,2),c(0,0), xlim=c(0,3), ylim=c(min(mci.age,ctl.age) - 5, max(mci.age,ctl.age) + 5)
       , type="n", ylab="Age", xlab="Cohort", xaxt="n", main=paste("Age Distribution: ", type))
  points(rep(1, length(ctl.age)), ctl.age, cex=.8, pch=4)
  points(rep(2, length(mci.age)), mci.age, cex=.8, pch=4)
  points(1, mean(ctl.age), col="blue", cex=1.3, pch="o")
  points(2, mean(mci.age), col="red", cex=1.3, pch="o")
  axis(1, at=c(1,2), labels=c("CTL", "MCI"))
  mtext( paste("Ind. T-test, p= ", round(age.eq$p.value, 3)), side=3)
}
################################################


################################################
ROCstats <- function(dat, dem, y, type)
# input: X, the design matrix. assumes that the last three vectors are demographics
# y, the response (cohort)
# type: univariate or LASSO
# output: a list containing ROC performance, AUC, and coefficient estimates
{
  X <- cbind(dat, dem)
  n <- dim(X)[1]
  p <- dim(X)[2]
  pprob.vec <- rep(0,n)
  
  # type = lasso is deprecated; not used in final analysis
  # we dont look at AUC of lasso models
  if( type == "lasso")
  {
    
    # demographic data is unpenalized. get optimal lambda value through cross validation
    # this part isn't random
    lambda <- optL1(y, penalized=dat, unpenalized=dem, lambda2=0, model="logistic", 
                      standardize=FALSE, maxlambda=20, minlambda=0.005)
     
    for(i in 1:n)
    {
      fit <- penalized(y[-i], penalized=dat[-i,], unpenalized=dem[-i,], standardize=FALSE, model="logistic",
                       lambda1=lambda$lambda, lambda2=0, steps="Park")
      index <- length(fit)
      B <- c(fit[[index]]@penalized, fit[[index]]@unpenalized)
      X <- cbind(dat, demog)
      pprob.vec[i] <- X[i,] %*% B
      
     
    }
    
  } else if( type == "univar" )
  {
    for(i in 1:n)
    {
      out <- glm(y[-i] ~ X[-i,], family=binomial(link="logit"))
      B <- out$coefficients
      pprob.vec[i] <- c(1,X[i,]) %*% B
    }
    
  }
  
  # get ROC statistics
  pred <- prediction(pprob.vec, y)
  perf <- performance(pred, 'tpr', 'fpr')
  AUC <- performance(pred, 'auc')@y.values
  roc.out <- roc(response=y, predictor=pprob.vec)
  return( list( perf, AUC, roc.out ) )
}
################################################

################################################
featureSelection <- function(X, y)
{
  #n <- dim(X)[1]
  n <- 4
  p <- dim(X)[2]
  p.thresh <- 0.05
  count <- matrix(0, nrow=4, ncol=p)
  
  
  for(i in 1:n)
  {
    out <- lars.glm(X[-i,], y[-i], family="binomial")
    out.p <- covTest(out, X[-i,], y[-i])
    
    index <- which(out.p$results[,3] < p.thresh)
    count[index] <- count[index] + 1
  }
  
  return(count)
}
################################################
# ------------------------------------------------------------------------ #
# --------------------------- end functions ------------------------------ #
# ------------------------------------------------------------------------ #






# ------------------------------------------------------------------------ #
# -------------------------------- main ---------------------------------- #
# ------------------------------------------------------------------------ #
setwd("/Users/pluta/Desktop/subfield_exp")
pdf('summary_data.pdf', onefile=TRUE, useDingbats=FALSE, width=8.5, height=14)


# Set up graphics
par(mfrow=c(3,2));
par(mar=c(5,4,8,2));


# *****************************************
# ---- scrub WOLK data --------
# read in volume, demographic data, and freesurfer data
data <- read.csv("ns_fullset_rev.csv", header=TRUE)
demog <- read.csv("demog_wolk.csv")
#hbt.data <- read.csv("hvols.csv", header=TRUE)
data.fs <- read.csv("fs_aseg.csv")
var.list <- c("CA1", "CA2", "CA3", "DG", "SUB", "ERC", "BA35", "BA36")

# remove duplicate entries but retain both sides
data1 <- data[data$side == "left",]
data1 <- subset(data1, !duplicated(data1$ID))
data2 <- data[data$side == "right",]
data2 <- subset(data2, !duplicated(data2$ID))
data  <- rbind(data1, data2)


# remove subjects of the 3## class
# idk what these are
demog <- demog[demog$Subject < 300,]


# setup the cohort variable
# remove subjects that arent control or mci, recode id to 'control' or 'mci'
cohort <- c()
id.vec <- c()
for( i in 1:length(data$ID))
{
  id <- as.integer(strsplit(as.character(data$ID[i]), "DW")[[1]][[2]])
  if( id > 100 & id < 200) { type <- "control" } else if( id > 200 & id < 300)
  { type <- "mci" } else if (id > 300) { type <- "other" }
  cohort <- c(cohort, type)
  id.vec[i] <- id
}

# turn cohort into a factor
data$cohort <- cohort
data$id.num <- id.vec

# order the data by id so i can merge with cohort data
data <- data[data$cohort == "control" | data$cohort == "mci",]
data$cohort <- as.factor(data$cohort)
data <- data[order(data$id.num),]

# remove anyone whos not control or mci
demog <- demog[demog$Subject < 300,]

# remove subjects not common to both datasets; they're not in the hbt dataset because they were
# removed due to qc
for(i in c("DW116", "DW141", "DW148", "DW156", "DW201", "DW224", "DW233", "DW225"))
{
  index <- which(data.fs$ID == i)
  data.fs <- data.fs[-index,]
}

# remove 235 until i get demographic data
index <- which(data.fs$ID == "DW235")
data.fs <- data.fs[-index,]

# missing DW235 - get this!
# get rid of case DW300 and DW235 (temporarily- will add this back)
index <- which(data$id.num == 300)
data <- data[-index,]

index <- which(data$id.num == 235)
data <- data[-index,]

# remove the subjects that arent common to both lists
for( i in c(116, 141, 148, 156, 201, 224, 233, 225))
{
  index <- which(demog$Subject == i)
  demog <- demog[-index,]
}

# remove the subjects that arent common to both lists
for( i in c("DW116", "DW141", "DW148", "DW156", "DW201", "DW224", "DW233", "DW225"))
{
  index <- which(data$ID == i)
  data <- data[-index,]
}

# add the demographic data to the volume data
data1 <- data[data$side == "left",]
data1$Age <- demog$Age
data1$Education <- demog$Education
data1$Gender <- demog$Gender

data2 <- data[data$side == "right",]
data2$Age <- demog$Age
data2$Education <- demog$Education
data2$Gender <- demog$Gender



data <- rbind(data1, data2)

# voxel dimensions
xdim <- 0.4
ydim <- 0.4
zdim <- data$Slice_Thickness
vox.vol <- xdim * ydim * zdim

# cohort needs to be encoded as 0s and 1s
y <- as.numeric(data$cohort[data$side == "left"])
y <- y - rep(1, length(y))

data.fs$cohort <- y

# slice thickness correspondes to the coil used. 2.6 = 8 channel, 2 = 32 channel
# theres some weird perturbations in the slice thickness data, so round to 1 decimal. 
# should be a two-level categorical variable
data$Coil <- as.factor(round(data$Slice_Thickness, 1))
#
# end scrub
# ********************************


# ********************************
# setup vars 
var.names <- c("CA1.l", "CA2.l", "CA3.l", "DG.l", "SUB.l", "ERC.l", "BA35.l", "BA36.l",
               "CA1.r", "CA2.r", "CA3.r", "DG.r", "SUB.r", "ERC.r", "BA35.r", "BA36.l")

# create the design matrix for wolk data, and assess normality
for(sf in var.list)
{
  # nah the order is screwed up here; gotta split the variables by side first, then transform to N(0,1)
  sf.name <- paste(sf, "_vol", sep="")
  
  # which column contains the data for this subfield
  index <- which(colnames(data) == sf.name)
  
  # convert voxels to mm^3
  temp <- data[,index] * vox.vol
  
  # create new variables for each subfield by side
  # e.g. CA1.l, CA2.r, and so on
  for(side in c("l", "r"))
  {  
    name <- paste(sf, side, sep=".")
    if( side == "l") { side.full <- "left" }
    if( side == "r") { side.full <- "right" }
    assign(name, data[,index][data$side == side.full])
  }
}

# normalize cortical regions by number of slices
BA35.l <- BA35.l / data$BA35_ns[data$side == "left"]
BA36.l <- BA36.l / data$BA36_ns[data$side == "left"]
BA35.r <- BA35.r / data$BA35_ns[data$side == "right"]
BA36.r <- BA36.r / data$BA36_ns[data$side == "right"]
ERC.l  <- ERC.l / data$ERC_ns[data$side == "left"]
ERC.r <- ERC.r / data$ERC_ns[data$side == "right"]
# the design matrix X
sf.vols <- as.matrix(cbind(CA1.l, CA2.l, CA3.l, DG.l, SUB.l, ERC.l, BA35.l, BA36.l, 
                           CA1.r, CA2.r, CA3.r, DG.r, SUB.r, ERC.r, BA35.r, BA36.r))
sf.res <- sf.vols
# normalize predictors
p <- dim(sf.vols)[2]
for(i in 1:p)
{
  if(p %% 2 == 0) { SIDE <- "right" } else { SIDE <- "left"}
  
  fit <- lm(sf.res[,i] ~ Age + Gender + ICV + Coil, data=data[data$side==SIDE,])
  sf.res[,i] <- fit$res
  sf.res[,i] <- (sf.res[,i] - mean(sf.res[,i])) / sd(sf.res[,i])
  sf.vols[,i] <- (sf.vols[,i] - mean(sf.vols[,i])) / sd(sf.vols[,i])
  
}

# error check
eps <- 0.001
if( any(diag(var(sf.vols)) - 1 > eps)) { print("ERROR: variances are not 1! You need to standardize volumes!")}
# ------ WOLK data is good to go! ---------
# *****************************************


# *****************************************
# analysis

# -- marginal analysis --- #
# in this case marginal analysis refers to the series of independent linear models
# for each subfield
n.sf <- length(var.names)/2
t.l <- matrix(rep(0, n.sf * 2), nrow=n.sf, ncol=2)
colnames(t.l) <- c("z", "p")
t.r <- t.l

rownames(t.l) <- var.names[1:n.sf]
rownames(t.r) <- var.names[(n.sf + 1):(n.sf*2)]

count <- 0
for(i in var.names)
{
  # to mimic the results in the first draft, drop Coil from the model
  count <- count + 1
  fit <- glm(cohort ~ get(i) + Age + Gender + ICV + Coil, data=data[data$side=="left",], family='binomial')
  # extra z and p values from results
  
  if( count <= n.sf)
  {
    t.l[count,1] <- coef(summary(fit))[2,3]  # z values
    t.l[count,2] <- coef(summary(fit))[2,4]  # p values
  } else 
  {
    t.r[(count - n.sf),1] <- coef(summary(fit))[2,3]  # z values
    t.r[(count - n.sf),2] <- coef(summary(fit))[2,4]  # p values 
  }

}

# create tables of z and pvalues for the paper
t.l <- round(t.l,3)
t.r <- round(t.r,3)
#---------------#

# this is the correct model; now figure out how to loop through and put it
# in a table. dont forget to add coil.
fit <- glm(cohort ~ CA1.l + Age + Gender + ICV + Coil, data=data[data$side=="left",], family='binomial')
p <- coef(summary(fit))[2,4]
z <- coef(summary(fit))[2,3]

# analysis dataset as a whole:
# does it matter if we adjust for demog?
# this is the output after adjusting for demographic data (Age, Gender, ICV, Coil)
out <- lars.glm(sf.res, y, family="binomial", standardize=FALSE)
out.p <- covTestVerbose(out, sf.res, y)

# this is the output with raw volume
out2 <- lars.glm(sf.vols, y, family="binomial", standardize=FALSE)
out2.p <- covTest(out2, sf.vols, y)
# ---

# model by side
out.l<- lars.glm(sf.res[,1:8], y, family="binomial", standardize=FALSE)
out.l.p <- covTest(out.l, sf.res[,1:8], y)

out.r <- lars.glm(sf.res[,9:16], y, family="binomial", standardize=FALSE)
out.r.p <- covTest(out.r, sf.res[,9:16], y) 


# var types
# sf.vols - standard normal predictors
# sf.res - standard normal predictors, adjusted for demographic data
# CA1.l etc - raw volumes, not normalized, not adjusted
n <- dim(sf.vols)[1]
demog <- cbind( data$Age[1:n], data$ICV[1:n], data$Gender[1:n], data$Coil[1:n] )
colnames(demog) <- c("Age", "ICV", "Gender", "Coil")

# CA1.l, BA35.l are the significant predictors from the covariance test

# get AUC statistics for each method uni/multi/whole hippocampus
uniStat <- ROCstats(BA35.l, demog, y, "univar")
mvStat <- ROCstats(cbind(CA1.l, BA35.l), demog, y, "univar")
mvStat2 <- ROCstats(cbind(CA1.l, BA35.l, BA36.l, BA36.r), demog, y, "univar")
# don't need coil for whole hippocampus data, doesnt matter
fsStat  <- ROCstats(data.fs$FSH_vol.left, demog[,1:3], y, "univar")

# p-values for the two methods
uni.p <- roc.test(uniStat[[3]], fsStat[[3]], paired=TRUE, method="v")
mv.p <- roc.test(fsStat[[3]], mvStat[[3]], paired=TRUE, method="v")

# text for hte legend
text.fs <-  paste("Whole Hippocampus, AUC = ", round(fsStat[[2]][[1]][1],2), sep=" ")
text.uni <- paste("Univariate, AUC = ", round(uniStat[[2]][[1]][1],2), sep=" ")
text.mv <-  paste("Multivariate, AUC = ", round(mvStat[[2]][[1]][1],2), sep=" ")

# plot ROC curves
plot(fsStat[[1]], col="blue", main="ROC Curves")
plot(uniStat[[1]], col="green", add=TRUE)
plot(mvStat[[1]], col="red", add=TRUE)
mtext(paste("Univariate v. WH, p = ", round(uni.p$p.value,2), " MV v. WH, p = ", round(mv.p$p.value,2)), side=3 )
legend("bottomright", c(text.fs, text.uni,text.mv), col=c("blue", "green", "red"), lwd=3) 







dev.off()
# ****************************************
# ------------------------------------------------------------------------ #
# --------------------------- end --------- ------------------------------ #
# ------------------------------------------------------------------------ #
