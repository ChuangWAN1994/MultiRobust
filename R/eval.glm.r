# Fit the defined imputation models
eval.glm.K <- function(imp.model, data)
{
  K <- length(imp.model) # No. of imputation models
  imp.glm <- vector(mode = "list", length = K)
  for (k in 1:K){
    imp.modelk <- imp.model[[k]]
    imp.modelk$data <- data
	wts.tmp <- imp.modelk$weights
    imp.glm[[k]] <- eval(imp.modelk)
	imp.glm[[k]]$prior.weights <- wts.tmp
  }
  return(imp.glm)
}

# Fit the defined imputation models
eval.glm.KV <- function(imp.model, data)
{
  K <- length(imp.model) # No. of SETS of imputation models
  imp.glm <- vector(mode = "list", length = K)
  for (k in 1:K){
    imp.modelk <- imp.model[[k]]
    V <- length(imp.modelk)
    imp.glm[[k]] <- vector(mode = "list", length = V)
    for (v in 1:V){
      imp.modelkv <- imp.modelk[[v]]
      imp.modelkv$data <- data
	  wts.tmp <- imp.modelkv$weights
      imp.glm[[k]][[v]] <- eval(imp.modelkv)
	  imp.glm[[k]][[v]]$prior.weights <- wts.tmp
    }
  }
  return(imp.glm)
}
