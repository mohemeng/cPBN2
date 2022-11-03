#'Density for the continuous binomial distribution with parameters size and probability (prob)
#'
#' @param x vector of quantiles
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return density value(s) of the vector x
#' @export
#'
#' @examples dcbinom(x = 2, size = 5, prob = 0.5)

dcbinom <- function(x, size, prob){
  if(!is.numeric(x) | !is.numeric(size) | !is.numeric(prob)) stop("x, size, and prob must be numeric values")
  if(size <= 0) stop("size must be a positive integer")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")

  suppressMessages(require(cubature))
  bef <- function(y){
    g <- function(x){
         (x[1]*x[2])^(y - 1)*((1 - x[1])*(1 - x[2]))^(size - y)*(log(x[1]*(1 - x[2])) - log(x[2]*(1 - x[1])))
       }
       Int <- adaptIntegrate(g, lowerLimit = c(prob, 0), upperLimit = c(1, prob))
       return(Int$integral)
     }
    suppressWarnings(pdf.val <- sapply(x, bef))
    suppressWarnings(pdf.val <- pdf.val/((beta(x, size + 1 - x))^2))
    pdf.val[which(is.na(pdf.val) == TRUE)] = 0
    pdf.val[which(pdf.val == Inf)] = 0
    return(pdf.val)
}

#' Distribution function for the continuous binomial distribution
#'
#' @param x vector of quantiles
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return the CDF/SF values
#' @export
#'
#' @examples pcbinom(x = 4, size = 8, prob = 0.3)

pcbinom <- function(x, size, prob, lower.tail = TRUE){
  if(!is.numeric(x) | !is.numeric(size) | !is.numeric(prob)) stop("x, size, and prob must be numeric values")
  if(size <= 0) stop("size must be a positive integer")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")

  suppressWarnings(int <- 1 - pbeta(prob, x, size - x + 1, lower.tail = lower.tail))

  int[which(is.na(int) == TRUE)] = 0

  int[which(x >= size + 1)] = 1
  return(int)
}

#' Random number generation for the continuous binomial distribution
#'
#' @param n number of observations
#' @param size number of trials
#' @param prob probability of success
#'
#' @return vector of random values from the continuous Binomial distribution
#' @export
#'
#' @examples rcbinom(n = 1, size = 5, prob = 0.3)

rcbinom <- function(n, size, prob){
  if(!is.numeric(size) | !is.numeric(prob) | !is.numeric(n)) stop("inputs n, size and prob must be numeric")
  if(n <= 0 | n%%1 != 0 ) stop("input n must be a positive integer greater than 1")
  if(size < 1) stop("size >= 1")

  suppressMessages(require(rootSolve))
  bb <- 1:n
  ans.roots <- array(NA, dim = length(bb))
  gen <- runif(n)  # n uniform values

  for(i in bb){
    func <- function(x){
      f <- 1 - pbeta(prob, x, size - x + 1) - gen[i]
      return(f)
    }
    suppressWarnings(ans.roots[i] <- uniroot(func, lower =  0.001, upper = size, extendInt = "yes")$root)
  }
  return(ans.roots)
}

#' Quantile function for the continuous Binomial distribution
#'
#' @param p vector of probabilities
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return vector of quantiles
#' @export
#'
#' @examples qcbinom(p = 0.3, size = 10, prob = 0.5)

qcbinom <- function(p, size, prob){
  if(!is.numeric(size) | !is.numeric(prob) | !is.numeric(p)) stop("inputs n, size and prob must be numeric")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")
  if(p <= 0 | p >= 1) stop("prob must be in the interval (0, 1)")

  suppressMessages(require(rootSolve))

  gen <- p
  bb <- 1:length(p)
  ans.roots <- array(NA, dim = length(bb))

  for(i in bb){
    func <- function(x){
      f <- 1 - pbeta(prob, x, size - x + 1) - gen[i]
      return(f)
    }
    suppressWarnings(ans.roots[i] <- uniroot(func, lower =  1, upper = size, extendInt = "yes")$root)
  }
  return(ans.roots)
}

#' MME estimates for r and beta of the Continuous Binomial distribution
#'
#' @param data required for parameter estimation
#'
#' @return r and beta estimates by MME
#' @export
#'
#' @examples cbinomMME(data = rcbinomMME(n = 1000, r = 2, beta = 5))

cbinomMME <- function(data){
  suppressMessages(require(nleqslv)); suppressMessages(require(cubature))
  do <- function(t){
    y <- numeric(2)
    y[1] <- adaptIntegrate(function(x) (gamma(t[1] + 1)/(gamma(x[1])*gamma(t[1] - x + 1)))*(x[2]^(x[1] - 1))*(1-x[2])^(t[1]-x[1]), lowerLimit = c(0, 0), upperLimit = c(t[1] + 1, t[2]))$integral - mean(data)
    y[2] <- adaptIntegrate(function(x) 2*x[1]*(gamma(t[1] + 1)/(gamma(x[1])*gamma(t[1] - x + 1)))*(x[2]^(x[1] - 1))*(1-x[2])^(t[1]-x[1]), lowerLimit = c(0, 0), upperLimit = c(t[1] + 1, t[2]))$integral - mean(data^2)
    y
  }
  t <- c(1, 1)
  return(nleqslv(x = t, do)$x)
}


















