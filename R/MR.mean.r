# Functions for mean estimation
MREst.mean <- function(response, reg.model, mis.model, data)
{
  resp <- data[ , response]
  if (!is.numeric(resp)) stop("response variable is not numeric")
  if (all(!is.na(resp))){
    warning("no missing data; sample mean is returned")
	return(mean(resp))
  } else {
  
  J <- length(mis.model)
  K <- length(reg.model)
  g.length <- J + K

  if (g.length == 0){
    warning("no model is specified; sample mean of the complete cases is returned")
    return(mean(resp[!is.na(resp)]))
  } else {

    # names of the auxiliary variables
    aux.names <- ext.names(reg.model = reg.model, mis.model = mis.model)
    if (any(is.na(data[ , aux.names]))) stop("auxiliary variables being used need to be fully observed")
	data.sub <- data[ , c(response, aux.names)]
    data.sub$R <- 1 * !is.na(resp) # missingness indicator
    n <- NROW(data.sub)
    m <- sum(data.sub$R) # number of observed subjects
    g.hat <- matrix(0, n, g.length)

    # missingness models
    if (J > 0){
	  if (any(lapply(mis.model, function(r) all.vars(r)[1L]) != "R"))
	    stop("the response variable of all models in 'mis.model' needs to be specified as 'R'")
      for (j in 1:J){
        mis.modelj <- mis.model[[j]]
        mis.modelj$data <- data.sub
        g.hat[ , j] <- eval(mis.modelj)$fitted.values
      }
    }

    # outcome regression models
    if (K > 0){
	  if (any(lapply(reg.model, function(r) all.vars(r)[1L]) != response))
	    stop(paste("the response variable of all models in 'reg.model' needs to be \'", response, "\'", sep = ""))
      for (k in 1:K){
        reg.modelk <- reg.model[[k]]
        reg.modelk$data <- data.sub
        aug <- eval(reg.modelk)
        g.hat[ , J + k] <- predict(aug, newdata = data.sub, type = "response")
      }
    }

    g.hat <- scale(g.hat, center = TRUE, scale = FALSE)[data.sub$R == 1, ]
    g.hat <- matrix(g.hat, m, )

    # define the function to be minimized
    Fn <- function(rho, ghat){ -sum(log(1 + ghat %*% rho)) }
    Grd <- function(rho, ghat){ -colSums(ghat / c(1 + ghat %*% rho)) }
    # calculate the weights
    rho.hat <- constrOptim(theta = rep(0, g.length), f = Fn, grad = Grd, ui = g.hat, ci = rep(1 / m - 1, m), ghat = g.hat)$par
    wts <- c(1 / m / (1 + g.hat %*% rho.hat))
	wts <- wts / sum(wts)
    estimate <- sum(resp[data.sub$R == 1] * wts)
    return(estimate)
  }
  }
}

#' Multiply Robust Estimation of the Marginal Mean
#'
#' \code{MR.mean()} is used to estimate the marginal mean of a variable which is subject to missingness. Multiple missingness probability models and outcome regression models can be accommodated.
#' @param response The response variable of interest whose marginal mean is to be estimated. 
#' @param reg.model A list of outcome regression models defined by \code{\link{def.glm}}.
#' @param mis.model A list of missingness probability models defined by \code{\link{def.glm}}. The dependent variable is always specified as \code{R}.
#' @param data A data frame with missing data encoded as \code{NA}.
#' @param bootstrap Logical. Should a bootstrap method be applied to calculate the standard error of the estimator and construct a Wald confidence interval for the marginal mean.
#' @param bootstrap.size A numeric value. Number of bootstrap resamples generated if \code{bootstrap = TRUE}.
#' @param alpha Significance level used to construct the 100(1 - \code{alpha})\% Wald confidence interval.
#' @import stats
#' @return
#' \item{\code{mu}}{The estimated value of the marginal mean.}
#' \item{\code{SE}}{The bootstrap standard error of \code{mu} when \code{bootstrap = TRUE}.}
#' \item{\code{CI}}{A Wald-type confidence interval based on \code{mu} and \code{SE} when \code{bootstrap = TRUE}.}
#' @references Han, P. and Wang, L. (2013). Estimation with missing data: beyond double robustness. \emph{Biometrika}, \strong{100}(2), 417--430.
#' @references Han, P. (2014). A further study of the multiply robust estimator in missing data analysis. \emph{Journal of Statistical Planning and Inference}, \strong{148}, 101--110.
#' @examples
#' # Simulated data set
#' set.seed(123)
#' n <- 400
#' gamma0 <- c(1, 2, 3)
#' alpha0 <- c(-0.8, -0.5, 0.3)
#' X <- runif(n, min = -2.5, max = 2.5)
#' p.mis <- 1 / (1 + exp(alpha0[1] + alpha0[2] * X + alpha0[3] * X ^ 2))
#' R <- rbinom(n, size = 1, prob = 1 - p.mis)
#' a.x <- gamma0[1] + gamma0[2] * X + gamma0[3] * exp(X)
#' Y <- rnorm(n, a.x, sd = sqrt(4 * X ^ 2 + 2))
#' dat <- data.frame(X, Y)
#' dat[R == 0, 2] <- NA
#'
#' # Define the outcome regression models and missingness probability models
#' reg1 <- def.glm(formula = Y ~ X + exp(X), family = gaussian)
#' reg2 <- def.glm(formula = Y ~ X + X ^ 2, family = gaussian)
#' mis1 <- def.glm(formula = R ~ X + X ^ 2, family = binomial(link = logit))
#' mis2 <- def.glm(formula = R ~ X + exp(X), family = binomial(link = cloglog))
#' est <- MR.mean(response = Y, reg.model = list(reg1, reg2), 
#'                mis.model = list(mis1, mis2), data = dat)
#' est$mu
#' @export

# Estimating the marginal mean
MR.mean <- function(response, reg.model = NULL, mis.model = NULL, data, bootstrap = FALSE, bootstrap.size = 500, alpha = 0.05)
{
  response <- as.character(substitute(response))
  est <- MREst.mean(response = response, reg.model = reg.model, mis.model = mis.model, data = data)

  # Bootstrap method for variance estimation
  if (bootstrap == TRUE){
    set.seed(bootstrap.size)
    bs.est <- rep(0, bootstrap.size)
	n <- NROW(data)
    for (b in 1:bootstrap.size){
      bs.sample <- data[sample(1:n, n, replace = TRUE), ]
	  while (any(colSums(is.na(bs.sample)) == n)) { bs.sample <- data[sample(1:n, n, replace = TRUE), ] }
      bs.est[b] <- MREst.mean(response = response, reg.model = reg.model, mis.model = mis.model, data = bs.sample)
    }

    se <- sd(bs.est) # bootstrap standard error
    cilb <- est - qnorm(1 - alpha / 2) * se
    ciub <- est + qnorm(1 - alpha / 2) * se
	list(mu = est, SE = se, CI = c(cilb, ciub))
  } else {
    list(mu = est)
  }
}
