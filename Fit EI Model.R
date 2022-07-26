rm(list = ls())
library(PoissonBinomial)
library(parallel)
setwd('~/Desktop/Repository')
source("Helper Functions.R")

# control variables
eps <- 1e-5
lr <- 1e-7

# load the data
load("trainData.rData")

#################################################
##                training loop                ##
#################################################

# prepare for Newton's method training
beta <- rep(0, ncol(data.train[[1]]$mm.precinct)) 
g <- rep(Inf, length(beta))
iter <- 0
trueLikVals <- c()
betaVals <- list()

# newton's method training
while((m <- sqrt(t(g) %*% g)) > eps) {
    
  # update
  iter <- iter + 1
  trueLikVals <- c(trueLikVals, (newLik <- lik(beta, data.train)))
  betaVals[[iter]] <- beta
    
  print(paste(iter, newLik, m))
  print(beta)
    
  # compute approximate gradient
  g <- grad.approx(beta, data = data.train)
  h <- hess.approx(beta, data = data.train)

  beta <- as.vector(beta - solve(h) %*% g)
  names(beta) <- names(g)
}
