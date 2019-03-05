# Estimating the quantile regression coefficients based on the imputed datasets
beta.imp.quant <- function(quantreg.model, imp.model, L, data)
{
  if(!requireNamespace("quantreg", quietly = TRUE))
    stop("need CRAN package 'quantreg' for quantile regression")
  
  K <- length(imp.model) # No. of SETS of imputation models
  betaimp <- vector(mode = "list", length = K)
  imp.data <- MI.glm(imp.model = imp.model, L = L, data = data)

  for (k in 1:K){
    quantreg.model$data <- imp.data[[k]]
    betaimp[[k]] <- eval(quantreg.model)
  }
  return(betaimp)
}

MREst.quantreg <- function(quantreg.model, imp.model = NULL, mis.model = NULL, L = 10, data)
{
  n <- NROW(data)
  mis.status <- (colSums(is.na(data)) != 0) # TRUE: missing; FALSE: no missing
  if (length(unique(colSums(is.na(matrix(data[ , mis.status], n, ))))) > 1)
    stop("missing variables have to be simultaneously missing")

  if (all(!mis.status)){
    message("no missing data")
    quantreg.model$data <- data
    estimate <- eval(quantreg.model)
    return(estimate)
  }

  J <- length(mis.model) # No. of missingness models
  K <- length(imp.model) # No. of SETS of imputation models
  mis.names <- names(data)[mis.status] # names of the missing variables
  data$R <- 1 * complete.cases(data) # missingness indicator
  nobs <- sum(data$R) # No. of observed subjects

  if (J + K == 0){
    warning("No model is specified; a complete-case analysis is returned")
    quantreg.model$data <- data[data$R == 1, ]
    estimate <- eval(quantreg.model)
    return(estimate)
  } else {
  g.hatJ <- NULL
  if (J > 0){
    g.hatJ <- matrix(0, n, J)
    for (j in 1:J){
      mis.modelj <- mis.model[[j]]
      mis.modelj$data <- data
      g.hatJ[ , j] <- eval(mis.modelj)$fitted.values
    }
  }

  g.hatK <- NULL
  if (K > 0){
    tmp.names <- vector(mode = "list", length = K)

    for (k in 1:K){
      tmp.names[[k]] <- lapply(imp.model[[k]], function(r) r[[2]][[2]])
      if (any(duplicated(tmp.names[[k]]))) stop("only one model for each missing variable is allowed in each set of imputation models; separate multiple working models into different sets of imputation models")
      if (length(tmp.names[[k]]) != length(mis.names)){ stop("one and only one imputation model needs to be specified for each missing variable")}
    }

    betaimp <- beta.imp.quant(quantreg.model = quantreg.model, imp.model = imp.model, L = L, data = data)
    for (k in 1:K){
	  tau <- betaimp[[k]]$tau
	  xmat <- betaimp[[k]]$x
	  y <- betaimp[[k]]$y
      nidx <- rep(1:n, L)
	  eek <- cbind(nidx, xmat * (tau - (y - c(xmat %*% betaimp[[k]]$coefficients) < 0)))
	  eek <- data.frame(eek)
	  eek <- as.matrix(aggregate(.~nidx, eek, mean))[ , -1]
      g.hatK <- cbind(g.hatK, eek)
    }
  }

  g.hat <- scale(cbind(g.hatJ, g.hatK), center = TRUE, scale = FALSE)[data$R == 1, ]
  g.hat <- matrix(data = g.hat, nrow = nobs, )
  
  # define the function to be minimized
  Fn <- function(rho, ghat){ -sum(log(1 + ghat %*% rho)) }
  Grd <- function(rho, ghat){ -colSums(ghat / c(1 + ghat %*% rho)) }
  # calculate the weights
  rho.hat <- constrOptim(theta = rep(0, NCOL(g.hat)), f = Fn, grad = Grd,
    ui = g.hat, ci = rep(1 / nobs - 1, nobs), ghat = g.hat)$par
  wts <- c(1 / nobs / (1 + g.hat %*% rho.hat))
  wts <- wts / sum(wts)
  quantreg.model$weights <- wts
  quantreg.model$data <- data[data$R == 1, ]
  estimate <- eval(quantreg.model)
  return(estimate)
  } # end else
}

#' Multiply Robust Estimation of the Quantile Regression Coefficients
#'
#' \code{MR.quantreg()} is used to estimate the quantile regression coefficients. Both missing response and/or missing covariates are allowed. Multiple missingness probability models and imputation models are allowed.
#' @param quantreg.model The quantile regression model of interest, defined by \code{\link{def.quantreg}}.
#' @param imp.model A list of imputation models defined by \code{\link{def.glm}}. Within the list, each element is a list of imputation models for the missing variables. One and only one model is allowed to be specified for each missing variable. Separate multiple working models for the same variable into different lists.
#' @param mis.model A list of missingness probability models defined by \code{\link{def.glm}}. The dependent variable is always specified as \code{R}.
#' @param L Number of random draws from the estimated imputation model.
#' @param data A data frame with missing data encoded as \code{NA}.
#' @param bootstrap Logical. Should a bootstrap method be applied to calculate the standard error of the estimator and construct a Wald confidence interval for the quantile regression coefficients. Default is \code{FALSE}.
#' @param bootstrap.size A numeric value. Number of bootstrap resamples generated if \code{bootstrap} is \code{TRUE}. Default is 500.
#' @param alpha Significance level used to construct the 100(1 - alpha)\% Wald confidence interval.
#' @import stats
#' @return
#' \item{\code{Estimate}}{A table containing the estimated quantile regression coefficients. If \code{bootstrap} is \code{TRUE}, bootstrap standard errors \code{SE} of the estimates and Wald confidence intervals are also included.}
#' \item{\code{fit}}{A fitted object of class inheriting from \code{\link[quantreg]{rq}} on \code{quantreg.model}.}
#' @references
#' Han, P., Kong, L., Zhao, J. and Zhou, X. (2018). A general framework for quantile estimation with incomplete data. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, in press.
#' @examples
#' # Missing covariate and missing response simultaneously
#' # Simulated data set
#' set.seed(123)
#' n <- 400
#' gamma0 <- c(1, 2, 3)
#' alpha0 <- c(-0.8, -0.5, 0.3)
#' S <- runif(n, min = -2.5, max = 2.5) # auxiliary information
#' S.2 <- S ^ 2
#' X1 <- rbinom(n, size = 1, prob = 0.5) # covariate X1
#' X2 <- rexp(n) # covariate X2
#' p.obs <- 1 / (1 + exp(alpha0[1] + alpha0[2] * S + alpha0[3] * S.2)) # missingness probability
#' R <- rbinom(n, size = 1, prob = p.obs)
#' a.x <- gamma0[1] + gamma0[2] * X1 + gamma0[3] * X2
#' Y <- rnorm(n, a.x)
#' dat <- data.frame(S, X1, X2, Y)
#' dat[R == 0, c(2, 4)] <- NA # X1 and Y are missing
#'
#' # Estimating quantile regression coefficients
#' # quantile regression model of interest
#' reg <- def.quantreg(formula = Y ~ X1 + X2)
#' # imputation models for X1
#' impX1.1 <- def.glm(formula = X1 ~ S, family = binomial(link = logit))
#' impX1.2 <- def.glm(formula = X1 ~ S + X2, family = binomial(link = cloglog))
#' # imputation models for Y
#' impY.1 <- def.glm(formula = Y ~ S, family = gaussian)
#' impY.2 <- def.glm(formula = Y ~ S + X2, family = gaussian)
#' # missingness probability models
#' mis1 <- def.glm(formula = R ~ S + S.2, family = binomial(link = logit))
#' mis2 <- def.glm(formula = R ~ S.2, family = binomial(link = cloglog))
#'
#' imp1 <- list(impX1.1, impY.1) # 1st set of imputation models for X1 and Y
#' imp2 <- list(impX1.2, impY.2) # 2nd set of imputation models for X1 and Y
#' imp3 <- list(impX1.1, impY.2) # 3rd set of imputation models for X1 and Y
#' imp4 <- list(impX1.2, impY.1) # 4th set of imputation models for X1 and Y
#'
#' results <- MR.quantreg(quantreg.model = reg, imp.model = list(imp1, imp2, imp3, imp4),
#'                        mis.model = list(mis1, mis2), L = 5, data = dat)
#' @export

MR.quantreg <- function(quantreg.model, imp.model = NULL, mis.model = NULL, L = 10, data, bootstrap = FALSE, bootstrap.size = 500, alpha = 0.05)
{
  est <- MREst.quantreg(quantreg.model = quantreg.model, imp.model = imp.model, mis.model = mis.model, L = L, data = data)

  # Bootstrap method for variance estimation
  if (bootstrap == TRUE){
    set.seed(bootstrap.size)
    bs.est <- NULL
    n <- NROW(data)
    for (b in 1:bootstrap.size){
      bs.sample <- data[sample(1:n, n, replace = TRUE), ]
      while (any(colSums(is.na(bs.sample)) == n)) { bs.sample <- data[sample(1:n, n, replace = TRUE), ] }
      bs.est <- rbind(bs.est, MREst.quantreg(quantreg.model = quantreg.model, imp.model = imp.model,
	    mis.model = mis.model, L = L, data = bs.sample)$coefficients)
    }
    se <- apply(bs.est, 2, sd) # bootstrap standard error
    estimate <- as.matrix(cbind(est$coefficients, se))
    cilb <- estimate[ , 1] - qnorm(1 - alpha / 2) * estimate[ , 2]
    ciub <- estimate[ , 1] + qnorm(1 - alpha / 2) * estimate[ , 2]
    estimate <- as.matrix(cbind(estimate, cilb, ciub))
    colnames(estimate) <- c("Estimate", "SE", "CI Lower", "CI Upper")
    list(estimate = estimate, fit = est)
  }

  else { list(estimate = est$coefficients, fit = est) }
}
