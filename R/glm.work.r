#' Define a Generalized Linear Model
#'
#' Define a generalized linear model. All the arguments in \code{\link{glm}} are allowed except for \code{data}. Supported types of \code{family} include \code{gaussian}, \code{binomial}, \code{poisson}, \code{Gamma} and \code{inverse.gaussian}.
#' @param formula The formula of the model to be fitted.
#' @param family The distribution of the response variable and the link function to be used in the model.
#' @param weights The prior weights to be used in the model.
#' @param ... Addition arguments for the function \code{\link{glm}}.
#' @seealso \code{\link{glm}}.
#' @examples
#' # A logistic regression with response R and covariates X1 and X2
#' mis1 <- glm.work(formula = R ~ X1 + X2, family = binomial(link = logit))
#' @export
glm.work <- function(formula, family = gaussian, weights = NULL, ...)
{
  mf <- match.call()
  mf[[1L]] <- quote(glm)
  mf$formula <- as.formula(formula)
  mf$family <- family
  mf$weights <- weights
  return(mf)
}
