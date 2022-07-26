
################################################
##              helper functions              ##
################################################

# sigmoid function
sigmoid <- function(x) {exp(x)/(1 + exp(x))}

# true data likelihood
lik <- function(beta, data) {
  
  precinctLogProbs <- sapply(data, FUN = function(data.p) {
    dpbinom(x = data.p$D, probs = sigmoid(data.p$mm.precinct %*% beta), log = TRUE, method = "Convolve")
  })
  
  precinctLogProbs[is.infinite(precinctLogProbs)] <- -1000#min(precinctLogProbs[!is.infinite(precinctLogProbs)])
  sum(precinctLogProbs)
  
}


# approx data likelihood
lik.approx <- function(beta, data) {
  
  precinctLogProbs <- sapply(data, FUN = function(data.p) {
    probs <- sigmoid(data.p$mm.precinct %*% beta)
    dnorm(x = data.p$D, mean = sum(probs), sd = sqrt(sum(probs*(1-probs))), log = TRUE)
  })
  
  sum(precinctLogProbs)
  
}

# true data gradient
grad <- function(beta, data) {
  
  full.mm <- do.call(rbind, lapply(data, FUN = function(x) {x$mm.precinct}))
  
  # iterate through the precincts
  conditionalProbs <- unlist(lapply(1:length(data), FUN = function(i) {
    
    cat(i)
    
    # get the data
    probs <- sigmoid(data[[i]]$mm.precinct %*% beta)
    trueProb <- dpbinom(x = data[[i]]$D, probs = probs)
    
    if(trueProb == 0)
      trueProb <- -40
    
    # get the conditional probabilities
    sapply(1:length(probs), FUN = function(j) {
      dpbinom(probs = probs[-j], x =  data[[i]]$D-1, method = 'Convolve')*probs[j]/trueProb
    })
  }))
  
  colSums((conditionalProbs - as.vector(1/(1 + exp(-full.mm %*% beta)))) * full.mm)
}


# logit swing approx
grad.ls.approx <- function(beta, data, lwr = 0, upr = 100, eps = 1e-12, subsample = NULL) {
  
  # iterate through the precincts
  if(is.null(subsample))
    indices <- 1:length(data)
  else
    indices <- sample(1:length(data.train), subsample)
  
  full.mm <- do.call(rbind, lapply(data, FUN = function(x) {x$mm.precinct}))
  
  conditionalProbs <- unlist(sapply(indices, FUN = function(i) {
    
    # get the data
    probs <- sigmoid(data[[i]]$mm.precinct %*% beta)
    D <- data[[i]]$D
    u <- upr
    l <- lwr 
    
    while(u - l > eps) {
      mid <- (l + u)/2
    
      if(sum(1/(1 + (1-probs)/probs*mid)) < D)
        u <- mid
      else
        l <- mid
    }
    
    return(1/(1 + (1-probs)/probs*mid))
  }))
  
  
  colSums((conditionalProbs - as.vector(1/(1 + exp(-full.mm %*% beta)))) * full.mm)
}

# function to compute normally approximated gradient
grad.approx <- function(beta, data) {
  
  precinctGrads <- sapply(1:length(data), FUN = function(i) {
    
    # compute the mean and variance
    dat <- data[[i]]
    probs <- sigmoid(dat$mm.precinct %*% beta)
    mean <- sum(probs)
    sigmaSq <- sum(probs*(1-probs))
    
    # compute the gradient
    term1 <- 1/sigmaSq*(dat$D-mean)*
      colSums(dat$mm.precinct*rep(probs*(1-probs), length = nrow(dat$mm.precinct)))
    term2 <- -1/2*((dat$D - mean)^2/sigmaSq^2 - 1/sigmaSq)*
      colSums(dat$mm.precinct*rep((2*probs - 1)*probs*(1-probs), length = nrow(dat$mm.precinct)))
    grad <- term1 + term2 
    
  })
  
  return(rowSums(precinctGrads))
  
}

# faster likelihood approximation
lik.approx.2 <- function(beta, data) {
  
  precinctParams <- t(sapply(data, FUN = function(data.p) {
    probs <- sigmoid(data.p$mm.precinct %*% beta)
    c(D = data.p$D, mean = sum(probs), sd = sqrt(sum(probs*(1-probs))))
  }))
  
  sum(dnorm(x = precinctParams[,1], mean = precinctParams[,2], sd = precinctParams[,3], log = TRUE))
  
}

# hessian for one precinct
computeHess.onePrecinct <- function(beta, precinct) {
  
  # get the data
  D <- precinct$D
  dm <- precinct$mm.precinct
  
  # compute the mean and variance
  probs <- sigmoid(dm %*% beta)
  mean <- sum(probs)
  sigmaSq <- sum(probs*(1-probs))
  
  # compute each term
  term1 <-  -1/sigmaSq^2*(D-mean)*colSums(dm*rep(probs*(probs-1)*(2*probs-1), length = nrow(dm))) %*%
    t(colSums(dm*rep(probs*(1-probs), length = nrow(dm)))) - 
    
    1/sigmaSq*colSums(dm*rep(probs*(1-probs), length = nrow(dm))) %*%
            t(colSums(dm*rep(probs*(1-probs), length = nrow(dm)))) + 
    
    1/sigmaSq*(D-mean)*t(dm) %*% (dm*rep(probs*(probs - 1)*(2*probs - 1), length = nrow(dm)))
  
  term2 <- 1/2/sigmaSq^(2)*colSums(dm*rep(probs*(probs-1)*(2*probs-1), length = nrow(dm))) %*% 
    t(colSums(dm*rep(probs*(probs-1)*(2*probs-1), length = nrow(dm)))) - 
    
    1/2/sigmaSq*t(dm) %*% (dm*rep(probs*(1-probs)*(1 + 6*(probs-1)*probs), length = nrow(dm)))
  
  term3 <- 1/sigmaSq^3*(D-mean)^2*
    colSums(dm*rep(probs*(probs-1)*(2*probs-1), length = nrow(dm))) %*% 
  t(colSums(dm*rep(probs*(1 - probs)*(2*probs-1), length = nrow(dm)))) + 
    
    1/sigmaSq^2*(D-mean)*colSums(dm*rep(probs*(1-probs), length = nrow(dm))) %*% 
                       t(colSums(dm*rep(probs*(1 - probs)*(2*probs-1), length = nrow(dm)))) + 
    
    1/2/sigmaSq^2*(D-mean)^2*t(dm) %*% (dm*rep(probs*(1-probs)*(1 + 6*(probs-1)*probs), length = nrow(dm)))
  
  return(list(term1, term2, term3))
}

# hessian for one precinct
computeHess.onePrecinct_v2 <- function(beta, precinct) {
  
  # get the data
  D <- precinct$D
  dm <- precinct$mm.precinct
  
  # compute the mean and variance
  probs <- sigmoid(dm %*% beta)
  mu <- sum(probs)
  sigmaSq <- sum(probs*(1-probs))
  
  # get the basic building blocks
  x <- colSums(dm*rep(probs*(1-probs), length = nrow(dm)))
  y <- colSums(dm*rep(probs*(1 - probs)*(1 - 2*probs), length = nrow(dm)))
  z <- x + (D - mu)/sigmaSq * y

  # compute each term
  hess <- -((1/2/sigmaSq - (D - mu)^2/2/sigmaSq^2)*t(dm) %*% (dm*rep(probs*(1-probs)*(1 + 6*(probs-1)*probs), length = nrow(dm))) -
              (D - mu)/sigmaSq * t(dm) %*% (dm*rep(probs*(1 - probs)*(1 - 2*probs), length = nrow(dm))) -
              1/2/sigmaSq^2 * y %*% t(y) +
              1/sigmaSq* z %*% t(z))
  
  return(hess)
}


# function to compute hessian
hess.approx <- function(beta, data) {
  
  res <- rowSums(sapply(1:length(data), FUN = function(i) {
    s <- computeHess.onePrecinct(beta, data[[i]])
    c(s[[1]], s[[2]], s[[3]])
  }))
  
  a <- matrix(res[1:length(beta)^2], nrow = length(beta))
  b <- matrix(res[(length(beta)^2 + 1):(2*length(beta)^2)], nrow = length(beta))
  c <- matrix(res[(2*length(beta)^2 +1):(3*length(beta)^2)], nrow = length(beta))
  
  return(a + b + c)
}

# better function to compute hessian
hess.approx_v2 <- function(beta, data) {
  
  res <- lapply(1:length(data), FUN = function(i) {
    computeHess.onePrecinct_v2(beta, data[[i]])
  })
  
  Reduce('+', res)
}


# hessian for one precinct
computeHess.onePrecinct_components <- function(beta, precinct) {
  
  # get the data
  D <- precinct$D
  dm <- precinct$mm.precinct
  
  # compute the mean and variance
  probs <- sigmoid(dm %*% beta)
  mu <- sum(probs)
  sigmaSq <- sum(probs*(1-probs))
  
  # get the basic building blocks
  x <- colSums(dm*rep(probs*(1-probs), length = nrow(dm)))
  y <- colSums(dm*rep(probs*(1 - probs)*(1 - 2*probs), length = nrow(dm)))
  z <- x + (D - mu)/sigmaSq * y
  
  # compute each term
  hess1 <- -((1/2/sigmaSq - (D - mu)^2/2/sigmaSq^2)*t(dm) %*% (dm*rep(probs*(1-probs)*(1 + 6*(probs-1)*probs), length = nrow(dm))) -
              (D - mu)/sigmaSq * t(dm) %*% (dm*rep(probs*(1 - probs)*(1 - 2*probs), length = nrow(dm))) -
              1/2/sigmaSq^2 * y %*% t(y))
  hess2 <- -(1/sigmaSq* z %*% t(z))
  
  return(list(hess1 = hess1, hess2 = hess2))
}

# better function to compute hessian
hess.approx_components <- function(beta, data) {
  
  res <- lapply(1:length(data), FUN = function(i) {
    computeHess.onePrecinct_components(beta, data[[i]])
  })
  
  return(res)
}


