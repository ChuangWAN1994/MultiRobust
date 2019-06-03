def.glm <- function(formula, family, weights = NULL, ...)
{
  mf <- match.call()
  mf[[1L]] <- quote(glm)
  mf$formula <- as.formula(formula)
  mf$family <- family
  mf$weights <- weights
  return(mf)
}

# Extract score functions
scorefun <- function(x, ...)
{
  xmat <- model.matrix(x)
  xmat <- naresid(x$na.action, xmat)
  if(any(alias <- is.na(coef(x)))) xmat <- xmat[ , !alias, drop = FALSE]
  wres <- as.vector(residuals(x, "working")) * weights(x, "working")
  rval <- wres * xmat
  return(rval)
}

# Estimating the regression coefficients based on the imputed datasets
beta.imp <- function(model, imp.model, L, data)
{
  K <- length(imp.model) # No. of SETS of imputation models
  betaimp <- vector(mode = "list", length = K)
  imp.data <- MI.glm(imp.model = imp.model, L = L, data = data)

  for (k in 1:K){
    model$data <- imp.data[[k]]
    reg.impk <- eval(model)
    # record the values of the score functions
    reg.impk$sc <- data.frame(cbind(reg.impk$data$obs, scorefun(reg.impk)))
    betaimp[[k]] <- reg.impk
  }
  return(betaimp)
}

MREst.reg <- function(model, imp.model, mis.model, moment, order, L, data)
{
  model.names <- all.vars(model$formula)
  if (any(model.names == "R")) stop("variable name 'R' is reserved for the missingness indicator; please change the variable 'R' in the data set to a different name")
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
  all.names <- unique(c(model.names, var.names$name.cov, var.names$name.res, moment))
  if (any(all.names == "R")) stop("variable name 'R' is reserved for the missingness indicator; please change the variable 'R' in the data set to a different name")
  data.sub <- data[ , all.names]
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

    betaimp <- beta.imp(model = model, imp.model = imp.model, L = L, data = data.sub)
    for (k in 1:K){
      sck <- betaimp[[k]]$sc
      sck <- as.matrix(aggregate(.~V1, sck, mean))[ , -1]
      g.hatK <- cbind(g.hatK, sck)
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
  wfun <- function(w) if(any(grepl("binomial glm", w))) invokeRestart("muffleWarning")
  estimate <- withCallingHandlers(eval(model), warning = wfun)
  return(estimate)
  } # end else
  }
}

#' Multiply Robust Estimation for (Mean) Regression
#'
#' \code{MR.reg()} is used for (mean) regression under generalized linear models with missing responses and/or missing covariates. Multiple missingness probability models and imputation models are allowed.
#' @param formula The \code{\link{formula}} of the regression model of interest.
#' @param family A description of the error distribution and link function to be used for the GLM of interest.
#' @param imp.model A list of possibly multiple lists of the form \code{list(list.1, list.2, ..., list.K)}, where \code{K} is the total number of different imputation models. For the \emph{k}-th imputation model, \code{list.k} is a list of possibly multiple models, each of which is defined by \code{\link{glm.work}} and imputes one single missing variable marginally. See details.
#' @param mis.model A list of missingness probability models defined by \code{\link{glm.work}}. The dependent variable is always specified as \code{R}.
#' @param moment A vector of auxiliary variables whose moments are to be calibrated.
#' @param order A numeric value. The order of moments up to which to be calibrated.
#' @param L Number of imputations.
#' @param data A data frame with missing data encoded as \code{NA}.
#' @param bootstrap Logical. If \code{bootstrap = TRUE}, the bootstrap will be applied to calculate the standard error and construct a Wald confidence interval.
#' @param bootstrap.size A numeric value. Number of bootstrap resamples generated if \code{bootstrap = TRUE}.
#' @param alpha Significance level used to construct the 100(1 - \code{alpha})\% Wald confidence interval.
#' @param ... Addition arguments for the function \code{\link{glm}}.
#' @seealso \code{\link{glm}}.
#' @import stats
#' @details The function \code{MR.reg()} currently deals with data with one missingness pattern. When multiple variables are subject to missingness, their values are missing simultaneously. The methods in Han (2016) and Zhang and Han (2019) specify an imputation model by modeling the joint distribution of the missing variables conditional on the fully observed variables. In contrast, the function \code{MR.reg()} specifies an imputation model by separately modeling the marginal distribution of each missing variable conditional on the fully observed variables. These marginal distribution models for different missing variables constitute one joint imputation model. Different imputation models do not need to model the marginal distribution of each missing variable differently.
#' @return
#' \item{\code{coefficients}}{The estimated regression coefficients.}
#' \item{\code{SE}}{The bootstrap standard error of \code{coefficients} when \code{bootstrap = TRUE}.}
#' \item{\code{CI}}{A Wald-type confidence interval based on \code{coefficients} and \code{SE} when \code{bootstrap = TRUE}.}
#' @references
#' Han, P. (2014). Multiply robust estimation in regression analysis with missing data. \emph{Journal of the American Statistical Association}, \strong{109}(507), 1159--1173.
#'
#' Han, P. (2016). Combining inverse probability weighting and multiple imputation to improve robustness of estimation. \emph{Scandinavian Journal of Statistics}, \strong{43}, 246--260.
#'
#' Zhang, S. and Han, P. (2019). A simple implementation of multiply robust estimation for GLMs with missing data. Unpublished manuscript.
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
#' # marginal imputation models for X1
#' impX1.1 <- glm.work(formula = X1 ~ S, family = binomial(link = logit))
#' impX1.2 <- glm.work(formula = X1 ~ S + X2, family = binomial(link = cloglog))
#' # marginal imputation models for Y
#' impY.1 <- glm.work(formula = Y ~ S, family = gaussian)
#' impY.2 <- glm.work(formula = Y ~ S + X2, family = gaussian)
#' # missingness probability models
#' mis1 <- glm.work(formula = R ~ S + I(S ^ 2), family = binomial(link = logit))
#' mis2 <- glm.work(formula = R ~ I(S ^ 2), family = binomial(link = cloglog))
#' # this example considers the following K = 3 imputation models for imputing the missing (X1, Y)
#' imp1 <- list(impX1.1, impY.1)
#' imp2 <- list(impX1.1, impY.2)
#' imp3 <- list(impX1.2, impY.1)
#'
#' results <- MR.reg(formula = Y ~ X1 + X2, family = gaussian, imp.model = list(imp1, imp2, imp3),
#'                   mis.model = list(mis1, mis2), L = 10, data = dat)
#' results$coefficients
#' MR.reg(formula = Y ~ X1 + X2, family = gaussian, 
#'        moment = c(S, X2), order = 2, data = dat)$coefficients
#' @export

MR.reg <- function(formula, family = gaussian, imp.model = NULL, mis.model = NULL, moment = NULL, order = 1, 
                   L = 30, data, bootstrap = FALSE, bootstrap.size = 300, alpha = 0.05, ...)
{
  model <- def.glm(formula = formula, family = family, ...)
  if (!is.null(moment)) moment <- as.character(substitute(moment))[-1]
  est <- MREst.reg(model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = data)

  # Bootstrap method for variance estimation
  if (bootstrap == TRUE){
    set.seed(bootstrap.size)
    bbb <- function(x, model, imp.model, mis.model, moment, order, L, data){
      n <- NROW(data)
      bs.sample <- data[sample(1:n, n, replace = TRUE), ]
      while (any(colSums(is.na(bs.sample)) == n)) { bs.sample <- data[sample(1:n, n, replace = TRUE), ] }
      b.est <- MREst.reg(model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = bs.sample)$coefficients
      return(b.est)
	}

	bs.est <- sapply(1:bootstrap.size, bbb, model = model, imp.model = imp.model, mis.model = mis.model, moment = moment, order = order, L = L, data = data)
	se <- apply(bs.est, 1, sd) # bootstrap standard error
    estimate <- est$coefficients
    cilb <- estimate - qnorm(1 - alpha / 2) * se
    ciub <- estimate + qnorm(1 - alpha / 2) * se
	CI <- cbind(cilb, ciub)
	colnames(CI) <- c("lower bound", "upper bound")
    list(coefficients = estimate, SE = se, CI = CI)
  }

  else { list(coefficients = est$coefficients) }
}
