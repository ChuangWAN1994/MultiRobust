#' Define a Linear Quantile Regression Model
#'
#' Define a quantile regression model. All the arguments in \code{\link[quantreg]{rq}} of the \code{quantreg} package are allowed except \code{data}.
#' @param formula The formula of the linear quantile regression model to be fitted.
#' @param tau The quantile to be estimated, which is generally a number between 0 and 1.
#' @param weights The prior weights to be used in the model.
#' @param ... Addition arguments for the function \code{\link[quantreg]{rq}}.
#' @seealso \code{\link[quantreg]{rq}}.
#' @examples
#' # A quantile regression model with response Y and covariates X1 and X2 at the 75th percentile
#' reg <- def.quantreg(formula = Y ~ X1 + X2, tau = 0.75)
#' @export
def.quantreg <- function(formula, tau = 0.5, weights = NULL, ...)
{
  mf <- match.call()
  mf[[1L]] <- quote(quantreg::rq)  # paste "rq" in the formula
  mf$formula <- substitute(formula) # avoid using the quotation marks in the formula
  mf$weights <- weights
  return(mf)
}
