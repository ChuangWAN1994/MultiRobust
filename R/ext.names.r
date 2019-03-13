ext.names <- function(reg.model, mis.model)
{
  J <- length(mis.model) 
  K <- length(reg.model) 
  if (J > 0) nameJ <- sapply(mis.model, function(r) all.vars(r$formula)[-1L])
  if (K > 0) nameK <- sapply(reg.model, function(r) all.vars(r$formula)[-1L])
  name <- unique(c(nameJ, nameK))
  return(name)
}

ext.names.reg <- function(imp.model, mis.model)
{
  J <- length(mis.model) 
  K <- length(imp.model) 
  
  if (J > 0) nameJ.cov <- sapply(mis.model, function(r) all.vars(r$formula)[-1L])
  if (K > 0){
    nameK.cov <- rapply(imp.model, function(r) all.vars(r$formula)[-1L], how = "unlist")
    nameK.res <- rapply(imp.model, function(r) all.vars(r$formula)[1L], how = "unlist")
  }
  
  name.cov <- unique(c(nameJ.cov, nameK.cov)) # names of covariates
  name.res <- unique(nameK.res) # names of responses
  list(name.cov = name.cov, name.res = name.res)
}
