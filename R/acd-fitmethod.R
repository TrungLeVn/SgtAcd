#########################################
# Main Estimation: ACDfit
#' @title Model Estimation
#' @description Estimation ACD model with fitting and solver control arguments
#' @usage acdfit = function(spec, data, solver = "ucminf", out.sample = 0, solver.control = list(),
#'       fit.control = list(stationarity = 0, fixed.se = 0, scale = 0,n.sim = 2000),
#'       skew0 = NULL, shape0 = NULL, cluster = NULL, ...)
#' @param spec: ACDspec object with model specification
#' @param data: Data
#' @param solver: The chosen solver, which could be "ucminf", "optim", ...and "ms-" version that allows for multi-starts optimization strategies.
#'        refers to the documentation of Alexios's racd for more detials
#' @param out.sample: The length of data that will be used as out of sample model evaluation
#' @param solver.control: List of control values for solvers. For examples, we could specify number of restarts that will be used in multi-strarts
#'        optimization strategies
#' @param fit.control: List of fitting controls arguments, such as number of simulation to generate starting values of optimization routines
#' @param cluster: An makePSHOCKcluster object to using parallel in model estimation
#' @return An ACDfit object which contains information about model estimation
#' @export acdfit
#-------------------------
acdfit = function(spec, data, solver = "msucminf", out.sample = 0, solver.control = list(restarts = 3,trace = TRUE),
                  fit.control = list(stationarity = 0, fixed.se = 0, scale = 0,n.sim = 5000,rseed = NULL),
                  skew0 = NULL, shape10 = NULL,shape20 = NULL, cluster = NULL, ...) {
  UseMethod("acdfit")
}
.arfimaxfilteracd = function(modelinc, pars, idx, data, N, arglist) {
  # if(model[1] == 0) pars[1,1] = 0
  m = as.integer(N[1])
  T = as.integer(N[2])
  data = as.double(data)
  tmph = arglist$tmph
  garchenv = arglist$garchenv
  if(modelinc[4]>0){
    hm = as.double(tmph[,"archm"])
  } else{
    hm = double(length = T)
  }
  if(modelinc[5]>0){
    skm = as.double(tmph[,"skm"])
  } else{
    skm = double(length = T)
  }
  if(modelinc[6]>0){
    pskm = as.double(tmph[,"pskm"])
  } else{
    pskm = double(length = T)
  }
  res = double(length = T)
  # this routine is used for the mean residuals to initiate the recursion so we ignore arfima before
  zrf = double(length = T)
  constm = double(length = T)
  condm = double(length = T)
  ans = list()
  #condstm is the mu value, if we use arma, adjusted for external regressor if neccessary
  # condm is conditoinal mean, w.r.t mu, ar1, ma1 value
  # res is ret - condm
  if (sum(modelinc[2:6] > 0)) {
    ans = try(.C("armafilterC", model = as.integer(modelinc), pars = as.double(pars), idx = as.integer(idx - 1),
                 hm = hm, skm = skm, pskm = pskm,
                 x = data, res = res,
                 zrf = zrf, constm = constm, condm = condm,  m = m, T = T, PACKAGE = "SgtAcd"), silent = TRUE)
    if (inherits(ans, "try-error")) {
      assign(".csol", 1, envir = garchenv)
      assign(".filtermessage", ans, envir = garchenv)
      res = data - pars[idx[1, 1]]
      ans$res = res
      #if (model[4] > 0) {
      #  ans$zrf = rugarch:::.fracdiff(c(1, rep(0, length(data) - 1)), darfima = pars[idx[4]])
      #  ans$res = rugarch:::.fracdiff(ans$res, darfima = pars[idx[4]])
      #}
      if (any(is.na(res)))
        res[which(is.na(res))] = 0
      return(ans)
    } else {
      assign(".csol", 0, envir = garchenv)
      if (any(is.na(ans$res)))
        res[which(is.na(ans$res))] = 0
      return(ans)
    }
  } else {
    ans = list()
    ans$res = data - pars[idx[1, 1]]
    ans$zrf = zrf
    if (any(is.na(ans$res)))
      res[which(is.na(ans$res))] = 0
    return(ans)
  }
}
.extractdata = function(data){
    xdata = xts::as.xts(data)
    obj = list()
    obj$data = as.numeric(zoo::coredata(xdata))
    obj$index = zoo::index(xdata)
    obj$period = median(diff(zoo::index(xdata)))
    return(obj)
  }
#------------------------------------
# Model estimation routines
#' @importFrom zoo coredata index
#' @importFrom xts xts as.xts is.xts
#' @importFrom parallel makeForkCluster makePSOCKcluster stopCluster clusterEvalQ clusterMap parLapply parSapply clusterExport
#' @importFrom numDeriv jacobian hessian grad
#' @importFrom ucminf ucminf
#' @importFrom Rsolnp solnp gosolnp startpars
#' @importFrom Matrix Matrix crossprod t as.matrix
#' @importFrom stats arima na.omit nlminb optim runif sd var median
#' @importFrom utils head tail
#' @import rugarch
#' @import methods
#' @export acdfit
#' @useDynLib SgtAcd
#------------------------------------
.acdfit = function(spec, data, solver = "msucminf", out.sample = 0, solver.control = list(restarts =3,trace = TRUE),
                         fit.control = list(stationarity = 0, fixed.se = 0,scale = 0, n.sim = 5000,rseed = NULL),
                         skew0 = NULL, shape10 = NULL,shape20 = NULL, cluster = NULL, ...) {
  tic = Sys.time()
  vmodel = spec@model$vmodel$model
  #------------
  # Set up several default fit.control arguments
  #------------
  if(is.null(fit.control$rseed))
    rseed = 2706 else rseed = fit.control$rseed
  if (is.null(solver.control$trace))
    trace = FALSE else trace = solver.control$trace
  if (is.null(fit.control$stationarity))
    fit.control$stationarity = FALSE
  if (is.null(fit.control$fixed.se))
    fit.control$fixed.se = FALSE
  if (is.null(fit.control$scale))
    fit.control$scale = FALSE
  if (is.null(fit.control$n.sim))
    fit.control$n.sim = 5000
  #-------
  # Get the fit.control
  #-------
  mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim","rseed"))
  if (any(is.na(mm))) {
    idx = which(is.na(mm))
    enx = NULL
    for (i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in fit.control:\n", enx), sep = "", collapse = " "), call. = FALSE, domain = NULL)
  }
  #-----------
  # Set up dataset
  #-----------
  xdata = .extractdata(data)
  if (!is.numeric(out.sample))
    stop("\nacdfit-->error: out.sample must be numeric\n")
  if (as.numeric(out.sample) < 0)
    stop("\nacdfit-->error: out.sample must be positive\n")
  n.start = round(out.sample, 0)
  n = length(xdata$data)
  if ((n - n.start) < 100)
    stop("\nacdfit-->error: function requires at least 100 data\n points to run\n")
  data = xdata$data[1:(n - n.start)]
  index = xdata$index[1:(n - n.start)]
  origdata = xdata$data
  origindex = xdata$index
  period = xdata$period
  #------------
  # create a temporary environment to store values (deleted at end of function)
  #------------
  garchenv = new.env(hash = TRUE)
  arglist = list()
  arglist$garchenv <- garchenv
  arglist$sbounds = spec@model$sbounds
  arglist$pmode = 0
  model = spec@model
  modelinc = model$modelinc
  pidx = model$pidx
  arglist$index = index
  arglist$trace = trace
  m = model$maxOrder
  model$modeldata$T = T = length(as.numeric(data))
  if (fit.control$scale)
    dscale = sd(data) else dscale = 1
  zdata = data/dscale
  arglist$transform = FALSE
  arglist$fnscale = 1
  arglist$data = zdata
  arglist$dscale = dscale
  arglist$model = model
  arglist$shape10 = shape10
  arglist$shape20 = shape20
  arglist$skew0 = skew0
  ipars = model$pars
  # The temporary list arglist contains all specfications and fit arguments - to be used in estimations
  #-------------------
  # Optimization Starting Parameters Vector & Bounds
  #--------------------
  tmp = acdstart(ipars, arglist,cluster)
  arglist = tmp$arglist
  ipars = arglist$ipars = tmp$pars
  arglist$model = model
  # we now split out any fixed parameters
  # estidx is the index of parameters that will be estimated
  estidx = as.logical(ipars[, 4])
  arglist$estidx = estidx
  arglist$fit.control = fit.control
  npars = sum(estidx)
  if(as.logical(trace)){
    print("Fitting stage 1: Getting starting values and parameter bounds")
    print(ipars[estidx,])
  }
  #if (any(ipars[, 2] == 1)) {
    #if (npars == 0) {
      #if (fit.control$fixed.se == 0) {
        # if all parameters are fixed an no standard erros are to be calculated then we return a filter object
       # cat("\nacdfit-->warning: all parameters fixed...returning ACDfilter object instead\n")
       # return(acdfilter(data = data, spec = spec, out.sample = out.sample))
     # } else {
        # if all parameters are fixed but we require standard errors, we skip the solver
       # use.solver = 0
      #  ipars[ipars[, 2] == 1, 4] = 1
      #  ipars[ipars[, 2] == 1, 2] = 0
      #  arglist$ipars = ipars
      #  estidx = as.logical(ipars[, 4])
      #  arglist$estidx = estidx
     # }
  #  } else {
      # with some parameters fixed we extract them (to be rejoined at end) so that they do not enter the solver
     # use.solver = 1
    #}
  #} else {
    use.solver = 1
 # }
  #-------------------
  # start counter
  #-------------------
  assign("racd_llh", 1, envir = garchenv)
  arglist$fit.control = fit.control

  fun = switch(vmodel, sGARCH = .sacdLLH, gjrGARCH = .gjracdLLH)
  if(as.logical(trace)){
    print("Fitting stage 2: Maximum Likelihood solving")
  }
  if (use.solver) {
    parscale = rep(1, length = npars)
    names(parscale) = rownames(ipars[estidx, ])
    if (modelinc[1] > 0)
      parscale["mu"] = abs(mean(zdata))
    #if (modelinc[7] > 0)
     # parscale["omega"] = var(zdata)
    arglist$returnType = "llh"
    #arglist$returnType = "all"
    solution = .acdsolver(solver, pars = ipars[estidx, 1], fun = fun, Ifn = NULL, ILB = NULL, IUB = NULL, gr = NULL, hessian = NULL, parscale = parscale,
                          control = solver.control, LB = ipars[estidx, 5], UB = ipars[estidx, 6], cluster = cluster, arglist = arglist,rseed = rseed)
    #-----------------------------------------------------------------------
    sol = solution$sol
    hess = solution$hess
    # hess = solution$sol$hessian
    timer = Sys.time() - tic
    pars = solution$sol$pars
    if (!is.null(sol$par)) {
      ipars[estidx, 1] = sol$par
      if (modelinc[7] == 0) {
        # call it once more to get omega
        #fun here is the to get the solution pars, set it back the the filtering and calculate the omega
        tmpx = fun(sol$par, arglist)
        ipars[pidx["omega", 1], 1] = get("omega", garchenv)
      }
      if (sum(ipars[, 2]) == 0) {
        if (modelinc[1] > 0)
          ipars[pidx["mu", 1]:pidx["mu", 2], 1] = ipars[pidx["mu", 1]:pidx["mu", 2], 1] * dscale
        ipars[pidx["omega", 1], 1] = ipars[pidx["omega", 1], 1] * dscale^2
      }
    } else {
      ipars[estidx, 1] = NA
    }
    arglist$ipars = ipars
    convergence = sol$convergence
    if (convergence != 0)
      warning("\nacdfit-->warning: solver failed to converge.")
  } else {
    solution = NULL
    hess = NULL
    timer = Sys.time() - tic
    convergence = 0
    sol = list()
    sol$message = "all parameters fixed"
  }
  fit = list()
  # check convergence else write message/return create a copy of ipars in case we need to change it below to calculate standard errors which we
  # will need to reset later (because for example, infocriteria uses estimated parameters, not fixed.
  ipars2 = ipars
  if (convergence == 0) { #If the solver is converged, then go to the make fit function to calculate the tests and the standard errors and p-values
    arglist$dscale = 1
    arglist$data = data
    if (sum(ipars[, 2]) > 0 && fit.control$fixed.se == 1) {
      ipars[ipars[, 2] == 1, 4] = 1
      ipars[ipars[, 2] == 1, 2] = 0
      arglist$ipars = ipars
      estidx = as.logical(ipars[, 4])
      arglist$estidx = estidx
    }
    if(as.logical(trace)){
      print("Fitting stage 3: Make fit model")
    }
    fit = .acdmakefitmodel(f = fun, T = T, m = m, timer = timer, convergence = convergence, message = sol$message, hess, arglist = arglist)
    model$modelinc[7] = modelinc[7]
    model$modeldata$data = origdata
    model$modeldata$index = origindex
    model$modeldata$period = period
    model$pars[, 1] = fit$ipars[, 1]
    model$pars[, 5:6] = ipars2[, 5:6]
    fit$ipars[, 4] = ipars2[, 4]
    fit$ipars[, 2] = ipars2[, 2]
    fit$ipars[, 5:6] = ipars2[, 5:6]
    # make sure omega is now included (for working with object post-estimation)
    fit$ipars["omega", 3] = 1
    model$pars["omega", 3] = 1
  } else {
    fit$message = sol$message
    fit$convergence = 1
    fit$skhEst = arglist$skhEst
    model$modeldata$data = origdata
    model$modeldata$index = origindex
    model$modeldata$period = period
  }

  # make model list to return some usefule information which will be called by other functions (show, plot, sim etc)
  model = model
  model$garchLL = get("garchLL", garchenv)
  model$n.start = n.start
  fit$skhEst = arglist$skhEst
  ans = new("ACDfit", fit = fit, model = model)
  rm(garchenv)
  return(ans)
}
setMethod("acdfit", signature(spec = "ACDspec"), .acdfit)
