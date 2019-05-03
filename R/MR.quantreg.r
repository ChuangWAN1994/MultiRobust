# Estimating the quantile regression coefficients based on the imputed datasets
beta.imp.quant <- function(model, imp.model, L, data)
{
  if(!requireNamespace("quantreg", quietly = TRUE))
    stop("need CRAN package 'quantreg' for quantile regression")
  
  K <- length(imp.model) # No. of SETS of imputation models
  betaimp <- vector(mode = "list", length = K)
  imp.data <- MI.glm(imp.model = imp.model, L = L, data = data)

  for (k in 1:K){
    model$data <- imp.data[[k]]
    betaimp[[k]] <- eval(model)
  }
  return(betaimp)
}

MREst.quantreg <- function(model, imp.model, mis.model, moment, order, L, data)
{
  model.names <- all.vars(model$formula)
  if (all(!is.na(data[ , model.names]))){
    warning("no missing data")
    model$data <- data
    estimate <- eval(model)
    return(estimate)
  } else {

  # names of missing variables in 'model'
  mis.var <- model.names[colSums(is.na(data[ , model.names])) != 0]
  J <- length(mis.model) # No. of missingness models
  K <- length(imp.model) # No. of LISTS of imputation models
  M <- length(moment)

  if (J + K + M == 0){
    warning("no imputation or missingness model or moment is specified; a complete-case analysis is returned")
    model$data <- data
    estimate <- eval(model)
    return(estimate)
  } else {

  var.names <- ext.names.reg(imp.model = imp.model, mis.model = mis.model)
  data.sub <- data[ , unique(c(model.names, var.names$name.cov, var.names$name.res, moment))]
  n <- NROW(data.sub)
  if (length(unique(colSums(is.na(as.matrix(data.sub[ , mis.var], n, ))))) > 1)
    stop("the current package requires missingness to be simultaneous for different variables")
  if (any(is.na(data.sub[ , unique(c(var.names$name.cov, moment))]))) stop("the covariates for imputation and missingness models and moments need to be fully observed")

  data.sub$R <- 1 * complete.cases(data.sub[ , mis.var]) # missingness indicator
  nobs <- sum(data.sub$R) # No. of observed subjects

  g.hatJ <- NULL
  if (J > 0){
	if (any(lapply(mis.model, function(r) all.vars(r)[1L]) != "R"))
      stop("the response variable of all models in 'mis.model' needs to be specified as 'R'")
    g.hatJ <- matrix(0, n, J)
    for (j in 1:J){
      mis.modelj <- mis.model[[j]]
      mis.modelj$data <- data.sub
      g.hatJ[ , j] <- eval(mis.modelj)$fitted.values
    }
  }

  g.hatK <- NULL
  if (K > 0){
    for (k in 1:K){
      tmp.names <- lapply(imp.model[[k]], function(r) all.vars(r)[1L])
	  if (any(!(mis.var %in% tmp.names))) stop(paste("imputation model", k, "needs to include a marginal imputation model for each missing variable"))
      if (length(tmp.names) != length(mis.var)) stop(paste("for imputation model", k, ", number of marginal imputation models exceeds number of missing variables"))
    }

    betaimp <- beta.imp.quant(model = model, imp.model = imp.model, L = L, data = data.sub)
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

  g.hatM <- NULL
  if (M > 0){
    for (mm in 1:M){
      g.hatM <- cbind(g.hatM, sapply(1:order, function(ord, dat){dat ^ ord}, dat = data.sub[ , moment[mm]]))
    }
  }

  g.hat <- scale(cbind(g.hatJ, g.hatK, g.hatM), center = TRUE, scale = FALSE)[data.sub$R == 1, ]
  g.hat <- matrix(data = g.hat, nrow = nobs, )
  
  # define the function to be minimized
  Fn <- function(rho, ghat){ -sum(log(1 + ghat %*% rho)) }
  Grd <- function(rho, ghat){ -colSums(ghat / c(1 + ghat %*% rho)) }
  # calculate the weights
  rho.hat <- constrOptim(theta = rep(0, NCOL(g.hat)), f = Fn, grad = Grd,
    ui = g.hat, ci = rep(1 / nobs - 1, nobs), ghat = g.hat)$par
  wts <- c(1 / nobs / (1 + g.hat %*% rho.hat))
  wts <- wts / sum(wts)
  model$weights <- wts
  model$data <- data.sub[data.sub$R == 1, ]
  estimate <- eval(model)
  return(estimate)
  } # end else
  }
}

#' Multiply Robust Estimation for Quantile Regression
#'
#' \code{MR.quantreg()} is used for quantile regression with missing responses and/or missing covariates. Multiple missingness probability models and imputation models are allowed.
#' @param model The quantile regression model of interest, defined by \code{\link{def.quantreg}}.
#' @param imp.model A list of possibly multiple lists of the form \code{list(list.1, list.2, ..., list.K)}, where \code{K} is the total number of different imputation models. For the \emph{k}-th imputation model, \code{list.k} is a list of possibly multiple models, each of which is defined by \code{\link{def.glm}} and imputes one single missing variable marginally. See details.
#' @param mis.model A list of missingness probability models defined by \code{\link{def.glm}}. The dependent variable is always specified as \code{R}.
#' @param moment A vector of auxiliary variables whose moments are to be calibrated.
#' @param order A numeric value. The order of moments to be calibrated.
#' @param L Number of imputations.
#' @param data A data frame with missing data encoded as \code{NA}.
#' @param bootstrap Logical. Should a bootstrap method be applied to calculate the standard error of the estimator and construct a Wald confidence interval for the quantile regression coefficients.
#' @param bootstrap.size A numeric value. Number of bootstrap resamples generated if \code{bootstrap = TRUE}.
#' @param alpha Significance level used to construct the 100(1 - \code{alpha})\% Wald confidence interval.
#' @import stats
#' @details The function \code{MR.quantreg()} currently deals with data with one missingness pattern. When multiple variables are subject to missingness, their values are missing simultaneously. The method in Han et al. (2019) specifies an imputation model by modeling the joint distribution of the missing variables conditional on the fully observed variables. In contrast, the function \code{MR.quantreg()} specifies an imputation model by separately modeling the marginal distribution of each missing variable conditional on the fully observed variables. These marginal distribution models for different missing variables constitute one joint imputation model. Different imputation models do not need to model the marginal distribution of each missing variable differently.
#' @return
#' \item{\code{coefficients}}{The estimated quantile regression coefficients.}
#' \item{\code{SE}}{The bootstrap standard error of \code{coefficients} when \code{bootstrap = TRUE}.}
#' \item{\code{CI}}{A Wald-type confidence interval based on \code{coefficients} and \code{SE} when \code{bootstrap = TRUE}.}
#' \item{\code{fit}}{A fitted object inheriting from class "\code{\link[quantreg]{rq}}" on \code{model}.}
#' @references
#' Han, P., Kong, L., Zhao, J. and Zhou, X. (2019). A general framework for quantile estimation with incomplete data. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}. \strong{81}(2), 305--333.
#' @examples
#' # Simulated data set
#' set.seed(123)
#' n <- 400
#' gamma0 <- c(1, 2, 3)
#' alpha0 <- c(-0.8, -0.5, 0.3)
#' S <- runif(n, min = -2.5, max = 2.5) # auxiliary variables
#' X1 <- rbinom(n, size = 1, prob = 0.5) # covariate X1
#' X2 <- rexp(n) # covariate X2
#' p.obs <- 1 / (1 + exp(alpha0[1] + alpha0[2] * S + alpha0[3] * S ^ 2)) # non-missingness probability
#' R <- rbinom(n, size = 1, prob = p.obs)
#' a.x <- gamma0[1] + gamma0[2] * X1 + gamma0[3] * X2
#' Y <- rnorm(n, a.x)
#' dat <- data.frame(S, X1, X2, Y)
#' dat[R == 0, c(2, 4)] <- NA # X1 and Y may be missing
#' 
#' # quantile regression model of interest
#' reg <- def.quantreg(formula = Y ~ X1 + X2, tau = 0.75)
#' # marginal imputation models for X1
#' impX1.1 <- def.glm(formula = X1 ~ S, family = binomial(link = logit))
#' impX1.2 <- def.glm(formula = X1 ~ S + X2, family = binomial(link = cloglog))
#' # marginal imputation models for Y
#' impY.1 <- def.glm(formula = Y ~ S, family = gaussian)
#' impY.2 <- def.glm(formula = Y ~ S + X2, family = gaussian)
#' # missingness probability models
#' mis1 <- def.glm(formula = R ~ S + I(S ^ 2), family = binomial(link = logit))
#' mis2 <- def.glm(formula = R ~ I(S ^ 2), family = binomial(link = cloglog))
#' # this example considers the following K = 3 imputation models for imputing the missing (X1, Y)
#' imp1 <- list(impX1.1, impY.1)
#' imp2 <- list(impX1.1, impY.2)
#' imp3 <- list(impX1.2, impY.1)
#' 
#' results <- MR.quantreg(model = reg, imp.model = list(imp1, imp2, imp3),
#'                        mis.model = list(mis1, mis2), L = 10, data = dat)
#' results$coefficients
#' summary(results$fit)
#' MR.quantreg(model = reg, moment = c(S, X2), order = 2, data = dat)$coefficients
#' @export

MR.quantreg <- function(model, imp.model = NULL, mis.model = NULL, moment = NULL, order = 1, 
                        L = 30, data, bootstrap = FALSE, bootstrap.size = 500, alpha = 0.05)
{
  if (!is.null(moment)) moment <- as.character(substitute(moment))[-1]
  est <- MREst.quantreg(model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = data)

  # Bootstrap method for variance estimation
  if (bootstrap == TRUE){
    set.seed(bootstrap.size)
    bbb <- function(x, model, imp.model, mis.model, moment, order, L, data){
      n <- NROW(data)
      bs.sample <- data[sample(1:n, n, replace = TRUE), ]
      while (any(colSums(is.na(bs.sample)) == n)) { bs.sample <- data[sample(1:n, n, replace = TRUE), ] }
      b.est <- MREst.quantreg(model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = bs.sample)$coefficients
      return(b.est)
    }

	bs.est <- sapply(1:bootstrap.size, bbb, model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = data)
    se <- apply(bs.est, 1, sd) # bootstrap standard error
    estimate <- est$coefficients
    cilb <- estimate - qnorm(1 - alpha / 2) * se
    ciub <- estimate + qnorm(1 - alpha / 2) * se
    list(coefficients = estimate, SE = se, CI = cbind(cilb, ciub), fit = est)
  }

  else { list(coefficients = est$coefficients, fit = est) }
}
