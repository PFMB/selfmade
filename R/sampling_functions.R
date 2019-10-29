#' Generate samples based on the decomposition of Y
#'
#' @param orthdir nx1 vector; component orthogonal to the direction of interest
#' @param dir nx1 vector; component in the direction of interest
#' @param this_sd scalar; (estimated / surrogate) value for the root of error variance
#' @param nrSample integer; number of samples to be used
#' @param sampFun function; function to generate samples from
#' @param checkFun function; function to congruency with initial selection
#' @param trace logical; if \code{TRUE}, a progress bar will be printed to the console
#'
gen_samples <- function(
  orthdir,
  dir,
  this_sd, # DR: why do we need this?
  nrSample = 1000,
  sampFun = function(n) rnorm(n, mean = 0, sd = this_sd),
  checkFun,
  trace = 0,
  init_draw = FALSE,
  set_seed = FALSE,
  y_idx = NULL,
  path = NULL,
  n_cores = detectCores(logical = FALSE),
  app = "")
{
  if (set_seed == TRUE) set.seed(1)
  fac <- sampFun(nrSample)
  yb <- lapply(fac, function(tau) as.numeric(orthdir + tau*dir))
  draw <- list("yb" = yb, "fac" = fac)
  attr(draw,"seed") <- .Random.seed
  attr(draw,"time") <- Sys.time()
  attr(draw,"os_info") <- sessionInfo()
  if (is.null(path)) stop("Specify path!")
  save(draw, file = paste0(path,"PoSI/",app,"/y_draw_",c(app,as.character(Sys.time())),".RData"))

  cat("Parallel execution of checkFun active! \n")
  cluster_cl <- makeCluster(n_cores, outfile = "")
  clusterEvalQ(cluster_cl, {library(mgcv)
    library(cAIC4)})
  clusterExport(cluster_cl,"selection_function")

  if (!is.null(y_idx)) {
    load(file = paste0(path,"PoSI/",app,"/y_draw_",app,".RData"))
    fac <- draw$fac
    yb <- draw$yb
    logvals <- parSapply(cluster_cl, yb[y_idx], checkFun)
    # logvals <- sapply(yb[y_idx], checkFun)
    stopCluster(cluster_cl)
    cat("Parallel execution of checkFun and y_idx deactive! \n")
    return(list(logvals = logvals, fac = fac))
  }

  fac <- sampFun(nrSample)
  yb <- lapply(fac, function(tau) as.numeric(orthdir + tau*dir))
  logvals <- parSapply(cluster_cl, yb, checkFun)
  stopCluster(cluster_cl)
  cat("Parallel execution of checkFun deactive! \n")
  return(list(logvals = logvals, fac = fac))

}
