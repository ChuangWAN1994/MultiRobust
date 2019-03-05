# Multiple imputation
MI.quantile <- function(imp.model, L, data)
{
  n <- NROW(data)
  K <- length(imp.model) # No. of imputation models
  newdata <- vector(mode = "list", length = K)
  imp.glm <- eval.glm.K(imp.model = imp.model, data = data)

  for (k in 1:K){
    newdat <- NULL
    imp.modelk <- imp.glm[[k]]

	for (l in 1:L){
      impvar <- paste(imp.modelk$formula[[2L]])
      fam <- imp.modelk$family$family
      # link <- imp.modelk$family$link
      # coeff <- coef(imp.modelk)
      wts.tmp <- imp.modelk$prior.weights
	  if (is.null(wts.tmp)) wts.tmp <- rep(1, n)
      m <- predict(imp.modelk, newdata = data, type = "response")

      if (fam == "gaussian"){
        vars <- deviance(imp.modelk) / df.residual(imp.modelk) / wts.tmp
        imp <- rnorm(n, mean = m, sd = sqrt(vars))
      }

	  else if (fam == "binomial"){
	    if (any(m < 0) | any(m > 1))
          stop(paste("imputation model ", k, " is not appropriate, estimated probability should be within [0,1] for 'binomial' family; specify a different imputation model"))
        if (any(wts.tmp %% 1 != 0))
          stop("cannot simulate from non-integer prior.weights")

        if (!is.null(md <- imp.modelk$model)){
          y <- model.response(md)
          if(is.factor(y)){
            imp <- factor(1 + rbinom(n, size = 1, prob = m), labels = levels(y))
          } else
          imp <- rbinom(n, size = wts.tmp, prob = m)/wts.tmp
        } else imp <- rbinom(n, size = wts.tmp, prob = m)/wts.tmp

      }

      else if (fam == "poisson"){
        if (any(m <= 0))
          stop(paste("imputation model ", k, " is not appropriate, estimated mean should be non-negative for 'Poisson' family; specify a different imputation model"))
        imp <- rpois(n, lambda = m)
      }

      else if (fam == "Gamma"){
        # if(!requireNamespace("MASS", quietly = TRUE))
        #   stop("need CRAN package 'MASS' for simulation from the 'Gamma' family")
        # shape <- MASS::gamma.shape(imp.modelk)$alpha * wts.tmp
		shape <- 1 / summary(imp.modelk)$dispersion * wts.tmp
        if (any(shape <= 0) | any(shape / m <= 0))
		  stop(paste("imputation model ", k, " is not appropriate, estimated shape and rate should be positive for 'Gamma' family; specify a different imputation model"))
        imp <- rgamma(n, shape = shape, rate = shape / m)
      }

      else if (fam == "inverse.gaussian"){
        if(!requireNamespace("SuppDists", quietly = TRUE))
          stop("need CRAN package 'SuppDists' for simulation from the 'inverse.gaussian' family")
		disp <- summary(imp.modelk)$dispersion
        if (any(m <= 0) | any(disp <= 0))
		  stop(paste("imputation model ", k, " is not appropriate, estimated mean should be positive for 'inverse.gaussian' family; specify a different imputation model"))
        imp <- SuppDists::rinvGauss(n, nu = m, lambda = wts.tmp / disp)
      }

      newdatl <- data.frame(imp)
      colnames(newdatl) <- impvar
      newdat <- rbind(newdat, cbind(1:n, rep(l,n), newdatl, data[ , -match(impvar, names(data))]))
    } # end l-th imputation

    names(newdat)[1] <- "obs"
    names(newdat)[2] <- "L"
    newdata[[k]] <- data.frame(newdat)
  }
  return(newdata)
}

MREst.quantile <- function(tau, imp.model = NULL, mis.model = NULL, L = 10, data)
{
  mis.status <- (colSums(is.na(data)) != 0) # TRUE: missing; FALSE: no missing
  if (sum(mis.status) > 1)
    stop("MR.quantile is not yet implemented for multivariate missing data")
  n <- NROW(data)
  if (all(!mis.status)){ stop("no missing data; please use R function 'quantile()'") }

  response.name <- names(data)[mis.status]
  response <- data[ , response.name]
  data$R <- 1 * !is.na(response) # missingness indicator
  cc <- response[data$R == 1] # complete cases
  nobs <- sum(data$R) # number of observed subjects

  if (!is.numeric(response)) stop("response variable is not numeric")

  # define missingness and outcome regression models
  J <- length(mis.model)
  K <- length(imp.model)
  g.length <- J + K

  if (g.length == 0){
    warning("no model is specified; sample quantile of the complete cases is returned")
    return(quantile(response[!is.na(response)], probs = tau))
  } else {
    g.hat <- matrix(0, n, g.length)

    # missingness models
    if (J > 0){
      for (j in 1:J){
        mis.modelj <- mis.model[[j]]
        mis.modelj$data <- data
        g.hat[ , j] <- eval(mis.modelj)$fitted.values
      }
    }

    if (K > 0){
      MI.data <- MI.quantile(imp.model = imp.model, L = L, data = data)
      for (k in 1:K){
	    datak <- MI.data[[k]]
	    y <- datak[ , response.name]
	    quantk <- quantile(y, probs = tau)
        nidx <- rep(1:n, L)
	    eek <- cbind(nidx, (tau - (y - quantk < 0)))
	    eek <- data.frame(eek)
	    g.hat[ , J + k] <- as.matrix(aggregate(.~nidx, eek, mean))[ , -1]
      }
    }

  g.hat <- scale(g.hat, center = TRUE, scale = FALSE)[data$R == 1, ]
  g.hat <- matrix(data = g.hat, nrow = nobs, )

  # define the function to be minimized
  Fn <- function(rho, ghat){ -sum(log(1 + ghat %*% rho)) }
  Grd <- function(rho, ghat){ -colSums(ghat / c(1 + ghat %*% rho)) }
  # calculate the weights
  rho.hat <- constrOptim(theta = rep(0, NCOL(g.hat)), f = Fn, grad = Grd,
    ui = g.hat, ci = rep(1 / nobs - 1, nobs), ghat = g.hat)$par
  wts <- c(1 / nobs / (1 + g.hat %*% rho.hat))
  wts <- wts / sum(wts)
  checkf <- function(quan, tau, resp, wghts){ sum(wghts * (resp - quan) * (tau - (resp - quan <= 0))) }
  estimate <- optimize(f = checkf, interval = c(min(cc),max(cc)),
    tau = tau, resp = cc, wghts = wts)$minimum
  return(estimate)
  } # end else
}

#' Multiply Robust Estimation of the Marginal Quantile
#'
#' \code{MR.quantile()} is used to estimate the marginal quantile of a variable which is subject to missingness. Multiple missingness probability models and imputation models are allowed.
#' @param tau A numeric value in [0,1]. The quantile to be estimated. Default is 0.5.
#' @param imp.model A list of imputation models defined by \code{\link{def.glm}}.
#' @param mis.model A list of missingness probability models defined by \code{\link{def.glm}}. The dependent variable is always specified as \code{R}.
#' @param L Number of random draws from the estimated imputation model.
#' @param data A data frame with missing data encoded as \code{NA}.
#' @param bootstrap Logical. Should a bootstrap method be applied to calculate the standard error of the estimator and construct a Wald confidence interval for the estimated marginal quantile. Default is \code{FALSE}.
#' @param bootstrap.size A numeric value. Number of bootstrap resamples generated if \code{bootstrap} is \code{TRUE}. Default is 500.
#' @param alpha Significance level used to construct the 100(1 - alpha)\% Wald confidence interval.
#' @import stats
#' @return
#' \item{\code{Estimate}}{The estimated value of the marginal quantile. If \code{bootstrap} is \code{TRUE}, bootstrap standard error \code{SE} of the estimate and a Wald confidence interval are provided.}
#' @references
#' Han, P., Kong, L., Zhao, J. and Zhou, X. (2018). A general framework for quantile estimation with incomplete data. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, accepted.
#' @examples
#' # Simulated data set
#' set.seed(123)
#' n <- 400
#' gamma0 <- c(1, 2, 3)
#' alpha0 <- c(-0.8, -0.5, 0.3)
#' X <- runif(n, min = -2.5, max = 2.5)
#' X.2 <- X ^ 2
#' exp.X <- exp(X)
#' p.miss <- 1 / (1 + exp(alpha0[1] + alpha0[2] * X + alpha0[3] * X.2))
#' R <- rbinom(n, size = 1, prob = 1 - p.miss)
#' a.x <- gamma0[1] + gamma0[2] * X + gamma0[3] * exp.X
#' Y <- rnorm(n, a.x, sd = sqrt(4 * X.2 + 2))
#' dat <- data.frame(X, X.2, exp.X, Y)
#' dat[R == 0, 4] <- NA
#'
#' # Define the outcome regression models and missingness probability models
#' imp1 <- def.glm(formula = Y ~ X + exp.X, family = gaussian)
#' imp2 <- def.glm(formula = Y ~ X + X.2, family = gaussian)
#' mis1 <- def.glm(formula = R ~ X + X.2, family = binomial(link = logit))
#' mis2 <- def.glm(formula = R ~ X + exp.X, family = binomial(link = cloglog))
#' est <- MR.quantile(tau = 0.25, imp.model = list(imp1, imp2),
#'                    mis.model = list(mis1, mis2), L = 5, data = dat)
#' @export

MR.quantile <- function(tau = 0.5, imp.model = NULL, mis.model = NULL, L = 10, data, bootstrap = FALSE, bootstrap.size = 500, alpha = 0.05)
{
  est <- MREst.quantile(tau = tau, imp.model = imp.model, mis.model = mis.model, L = L, data = data)

  # Bootstrap method for variance estimation
  if (bootstrap == TRUE){
    set.seed(bootstrap.size)
    bs.est <- rep(0, bootstrap.size)
    n <- NROW(data)
    for (b in 1:bootstrap.size){
      bs.sample <- data[sample(1:n, n, replace = TRUE), ]
      while (any(colSums(is.na(bs.sample)) == n)) { bs.sample <- data[sample(1:n, n, replace = TRUE), ] }
      bs.est[b] <- MREst.quantile(tau = tau, imp.model = imp.model, mis.model = mis.model, L = L, data = bs.sample)
    }

    se <- sd(bs.est) # bootstrap standard error
    cilb <- est - qnorm(1 - alpha / 2) * se
    ciub <- est + qnorm(1 - alpha / 2) * se
    results <- c(est, se, cilb, ciub)
    names(results) <- c("Estimate", "SE", "CI Lower", "CI Upper")
    return(results)
  } else {
    names(est) <- "Estimate"
    return(est)
  }
}
