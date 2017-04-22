#################################################################################
##
##   R package racd by Alexios Ghalanos Copyright (C) 2012-2014
##   This file is part of the R package racd.
##
##   The R package racd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package racd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# implements nlminb, lbgs and solnp
# only solnp implements true constraint (stationarity) optimization
.acdsolver = function(solver, pars, fun, Ifn, ILB, IUB, gr, hessian, parscale,
                      control, LB, UB, cluster, arglist,rseed)
{
  if(arglist$fit.control$n.sim>0){
    ifun = function(pars, arglist){
      alpha = pars[arglist$model$pos.matrix["alpha",1]:arglist$model$pos.matrix["alpha",2]]
      beta  = pars[arglist$model$pos.matrix["beta",1]:arglist$model$pos.matrix["beta",2]]
      sum(alpha+beta)
    }
    trace = arglist$trace
    nm = names(pars)
    N = control$restarts
    if(is.null(N)) N = 1
    control$restarts = NULL
    arglist$transform = FALSE
    dopt = lapply(pars, function(x) list(mean=unname(x), sd = 2*abs(unname(x))))
    pars = Rsolnp::startpars(pars = pars, fun = fun, distr = rep(2, length(pars)), distr.opt = dopt,
                     ineqfun = ifun, ineqLB = 1e-12, ineqUB = 0.99,
                     LB = LB, UB = UB, bestN = N, n.sim = arglist$fit.control$n.sim,
                     cluster = cluster, arglist = arglist,rseed = rseed)
    pars = pars[,-NCOL(pars), drop=FALSE]
    colnames(pars) = nm
  } else{
    pars = matrix(pars, nrow = 1)
  }
  retval = switch(solver,
                  nlminb = .nlminbsolver(pars, fun, gr, hessian, parscale, control, LB, UB, arglist),
                  solnp = .solnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist),
                  cmaes = .cmaessolver(pars, fun, control, LB, UB, arglist),
                  ucminf = .ucminfsolver(pars, fun, gr, hessian, parscale, control, LB, UB, arglist),
                  optim = .optimsolver(pars, fun, gr, control, LB, UB, arglist),
                  msoptim = .msoptimsolver(pars, fun, gr, control, LB, UB, arglist, cluster),
                  msucminf = .msucminfsolver(pars, fun, gr, hessian, parscale, control, LB, UB, arglist, cluster),
                  msnlminb = .msnlminbsolver(pars, fun, gr, hessian, parscale, control, LB, UB, arglist, cluster),
                  mssolnp = .mssolnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist, cluster)
  )
  return(retval)
}

.solnpsolver = function(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist)
{
  control = .solnp.ctrl(control)
  ans = try(solnp(pars, fun = fun, eqfun = NULL, eqB = NULL,
                  ineqfun = Ifn, ineqLB = ILB, ineqUB = IUB, LB = LB, UB = UB,
                  control = control, arglist), silent = TRUE)
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$message = ans
    sol$pars = rep(NA, length(pars))
    names(sol$par) = names(pars)
  }
  else {
    sol = ans
  }
  hess = NULL
  if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")
  return(list(sol = sol, hess = hess))
}

.cmaessolver = function(pars, fun, control, LB, UB, arglist){
  control = .cmaes.ctrl(control)
  ans = try(cmaes(pars, fun, lower = LB, upper = UB, insigma = 1, ctrl = control, arglist = arglist), silent=TRUE)
  if(inherits(ans, "try-error")){
    sol = list()
    sol$convergence = 1
    sol$message = ans
    sol$pars = rep(NA, length(pars))
    names(sol$pars) = names(pars)
  } else{
    ansf = try(nlminb(start = ans$par, fun, lower = LB, upper = UB, arglist = arglist,
                      control = list(trace= ifelse(control$options$DispFinal, 1, 0),
                                     eval.max = 2000, iter.max = 1000, step.min = 0.1)), silent = TRUE)
    if (inherits(ansf, "try-error") || ansf$convergence>0){
      # revert to cmaes solution
      sol = ans
      sol$pars = ans$par
      sol$par = NULL
      sol$convergence = 0
      sol$message = ans$stopflag
    }
    else {
      sol = ansf
    }
  }
  if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")
  return(list(sol = sol, hess = NULL))
}

# for use with cmaes vectorized + parallel application
pfun = function(pars, arglist, fun, cluster){
  ans = parallel::parRapply(cluster, pars, function(x) fun(x, arglist))
  return(ans)
}

.nlminbsolver = function (pars, fun, gr, hessian, parscale, control, LB, UB, arglist)
{
  control = .nlminb.ctrl(control)
  ans = try(nlminb(start = pars, objective = fun, gradient = gr,
                   hessian = hessian, arglist = arglist, scale = 1/parscale,
                   control = control, lower = LB, upper = UB), silent = TRUE)
  pscale = rep(1, length(pars))
  smin = 0.1
  maxtries = 1
  while(ans$convergence!=0 && maxtries<15) {
    x_pars = get("x_pars", envir = arglist$garchenv)
    control$step.min = smin*0.1
    smin = smin*0.1
    pscale = 0.25*pscale
    ans = try(nlminb(start = x_pars, objective = fun, gradient = gr,
                     hessian = hessian, arglist = arglist, scale = pscale,
                     control = control, lower = LB, upper = UB), silent = TRUE)
    maxtries = maxtries+1
  }
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$message = ans
    sol$pars = rep(NA, length(pars))
    names(sol$par) = names(pars)
  }
  else{
    sol = ans
    sol$pars = ans$par
    names(sol$pars) = names(pars)
    sol$par = NULL
  }
  hess = NULL
  if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")
  return(list(sol = sol, hess = hess))
}

#.mlsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB, arglist){
#	arglist$LB = LB
#	arglist$UB = UB
#	arglist$transform = TRUE
#	arglist$fnscale = -1
#	xfun = function(pars){
#		return(fun(pars, arglist))
#	}
#	ans = maxLik(xfun, grad = NULL, hess = NULL, start = logtransform(pars, LB, UB, inverse = TRUE),
#			method = "Newton-Raphson", constraints=NULL, print.level = 3,
#			tol = 1e-08, reltol=sqrt(.Machine$double.eps), gradtol = 1e-06,
#			steptol = 1e-10, lambdatol = 1e-06, qrtol = 1e-10, iterlim = 150,
#			finalHessian = FALSE, bhhhHessian=FALSE)
#	if (inherits(ans, "try-error")) {
#		sol = list()
#		sol$convergence = 1
#		sol$message = ans
#		sol$pars = rep(NA, length(pars))
#		names(sol$par) = names(pars)
#		hess = NULL
#	}
#	else{
#		sol = ans
#		sol$pars = logtransform(ans$estimate, LB, UB)
#		names(sol$pars) = names(pars)
#		sol$par = NULL
#		hess = NULL
#	}
#	return(list(sol = sol, hess = hess))
#}

.ucminfsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB, arglist){
  control = .ucminf.ctrl(control)
  arglist$LB = LB
  arglist$UB = UB
  arglist$transform = TRUE
  ans = ucminf(fn = fun, gr = NULL, par = .invlogtransform(pars, LB, UB), arglist = arglist,
               control = control, hessian=0)
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$message = ans
    sol$pars = rep(NA, length(pars))
    names(sol$par) = names(pars)
    hess = NULL
  }
  else{
    sol = ans
    sol$pars = logtransform(ans$par, LB, UB)
    names(sol$pars) = names(pars)
    sol$par = NULL
    if(ans$convergence>0) sol$convergence = 0 else sol$convergence = 1
    hess = NULL
  }
  if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")

  return(list(sol = sol, hess = hess))
}

.optimsolver = function(pars, fun, gr, control, LB, UB, arglist){
  # optim will not check for default control parameters
  arglist$LB = LB
  arglist$UB = UB
  arglist$transform = TRUE
  if(is.null(control$method)) method = "BFGS" else method = control$method
  control$method = NULL
  ans = optim(fn = fun, gr = NULL, par = .invlogtransform(pars, LB, UB), arglist = arglist,
              control = control, method = method)
  if (inherits(ans, "try-error")) {
    sol = list()
    sol$convergence = 1
    sol$message = ans$message
    sol$pars = rep(NA, length(pars))
    names(sol$par) = names(pars)
    sol$lik = NA
    hess = NULL
  }
  else{
    sol = ans
    sol$pars = logtransform(ans$par, LB, UB)
    names(sol$pars) = names(pars)
    sol$par = NULL
    sol$convergence = ans$convergence
    hess = NULL
  }
  if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")
  return(list(sol = sol, hess = hess))
}

# Multi-Start Optim
.msoptimsolver = function(pars, fun, gr, control, LB, UB, arglist, cluster){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, require(racd))
    parallel::clusterExport(cluster, c("fun", "LB", "UB", "control", "pars", "arglist"),
                            envir = environment())
    xsol = parallel::parLapplyLB(cluster, 1:N, function(i){
      return( .optimsolver(pars[i,], fun, gr = NULL, control, LB, UB, arglist))
    })
  } else{
    for(i in 1:N){
      xsol[[i]] = .optimsolver(pars[i,], fun, gr, control, LB, UB, arglist)
    }
  }
  best = sapply(xsol, function(x) x$sol$value)
  best = which(best == min(best, na.rm=TRUE))[1]
  return(xsol[[best]])
}

.msucminfsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB, arglist, cluster){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, require(SgtAcd))
    parallel::clusterExport(cluster, c("fun", "gr", "hessian","LB", "UB", "control", "pars", "arglist", "parscale"),
                            envir = environment())
    xsol = parallel::parLapplyLB(cluster, 1:N, function(i){
      return( .ucminfsolver(pars[i,], fun, gr = NULL, hessian = NULL,
                            parscale = parscale, control, LB, UB, arglist) )
    })
  } else{
    for(i in 1:N){
      xsol[[i]] = .ucminfsolver(pars[i,], fun, gr = NULL, hessian = NULL,
                                parscale = parscale, control, LB, UB, arglist)
    }
  }
  best = sapply(xsol, function(x) x$sol$value)
  best = which(best == min(best, na.rm=TRUE))
  return(xsol[[best]])
}


.msnlminbsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB, arglist, cluster){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, require(SgtAcd))
    parallel::clusterExport(cluster, c("fun", "LB", "UB", "control", "pars", "arglist", "parscale"),
                            envir = environment())
    xsol = parallel::parLapplyLB(cluster, 1:N, function(i){
      return( .nlminbsolver(pars, fun, gr = NULL, hessian = NULL,
                            parscale, control, LB, UB, arglist) )
    })
  } else{
    for(i in 1:N){
      xsol[[i]] = .nlminbsolver(pars, fun, gr = NULL, hessian = NULL,
                                parscale, control, LB, UB, arglist)
    }
  }
  best = sapply(xsol, function(x) x$sol$objective)
  best = which(best == min(best, na.rm=TRUE))
  return(xsol[[best]])
}


.mssolnpsolver = function(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist, cluster){
  N = NROW(pars)
  xsol = vector(mode="list", length = N)
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, require(SgtAcd))
    parallel::clusterExport(cluster, c("fun", "LB", "UB", "control", "pars", "arglist"),
                            envir = environment())
    xsol = parallel::parLapplyLB(cluster, 1:N, function(i){
      return( .solnpsolver(pars[i,], fun, Ifn = NULL, ILB = NULL, IUB = NULL, control, LB, UB, arglist) )
    })
  } else{
    for(i in 1:N){
      xsol[[i]] = .solnpsolver(pars[i,], fun, Ifn = NULL, ILB = NULL, IUB = NULL, control, LB, UB, arglist)
    }
  }
  best = sapply(xsol, function(x) tail(x$sol$value,1))
  best = which(best == min(best, na.rm=TRUE))
  return(xsol[[best]])
}


#.newuoasolver = function(pars, fun, control, LB, UB, arglist){
#	control = .minqa.ctrl(control, pars)
#	arglist$LB = LB
#	arglist$UB = UB
#	arglist$transform = TRUE
#	ans = newuoa(fn = fun, par = .invlogtransform(pars, LB, UB), arglist = arglist,
#			control = control)
#	if (inherits(ans, "try-error") | ans$ierr!=0){
#		sol = list()
#		sol$convergence = 1
#		sol$message = ans$msg
#		sol$pars = rep(NA, length(pars))
#		names(sol$par) = names(pars)
#		hess = NULL
#	}
#	else{
#		sol = ans
#		sol$pars = logtransform(ans$par, LB, UB)
#		names(sol$pars) = names(pars)
#		sol$par = NULL
#		sol$convergence = ans$ierr
#		hess = NULL
#	}
#	if(sol$convergence!=0) warning("\nracd-->warning: no convergence...\n")
#	return(list(sol = sol, hess = hess))
#}

################################################################################
# Solver control parameters
.ucminf.ctrl = function(control){
  if(is.null(control$trace)) control$trace = 0
  if(is.null(control$grtol)) control$grtol = 1e-8
  if(is.null(control$xtol)) control$xtol = 1e-12
  if(is.null(control$stepmax)) control$stepmax = 1
  if(is.null(control$maxeval)) control$maxeval = 1500
  if(is.null(control$grad)) control$grad = "forward"
  if(is.null(control$gradstep)) control$gradstep = c(1e-6, 1e-8)
  mm = match(names(control), c("trace", "grtol", "xtol", "stepmax", "maxeval", "grad", "gradstep"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}

.nlminb.ctrl = function(control){
  if(is.null(control$trace)) control$trace = 0
  if(is.null(control$eval.max)) control$eval.max = 1500
  if(is.null(control$iter.max)) control$iter.max = 500
  if(is.null(control$abs.tol)) control$abs.tol = 0
  if(is.null(control$rel.tol)) control$rel.tol = 1e-10
  if(is.null(control$x.tol)) control$x.tol = 2.2e-8
  if(is.null(control$xf.tol)) control$xf.tol = 2.2e-14
  if(is.null(control$step.min)) control$step.min = 0.1
  if(is.null(control$step.max)) control$step.max = 1
  if(is.null(control$sing.tol)) control$sing.tol = control$rel.tol
  mm = match(names(control), c("trace", "eval.max", "iter.max", "abs.tol", "rel.tol", "x.tol", "xf.tol",
                               "step.min", "step.max", "sing.sing"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}
.solnp.ctrl = function(control){
  if(is.null(control$trace)) control$trace = 0
  if(is.null(control$rho)) control$rho = 1
  if(is.null(control$outer.iter)) control$outer.iter = 50
  if(is.null(control$inner.iter)) control$inner.iter = 1800
  if(is.null(control$delta)) control$delta = 1e-9
  if(is.null(control$tol)) control$tol = 1e-8
  mm = match(names(control), c("trace", "rho", "outer.iter", "inner.iter", "delta", "tol"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}
.cmaes.ctrl = function(control){
  cmaes.control(options = control$options, CMA = control$CMA)
}

.minqa.ctrl = function(control, pars){
  n = length(pars)
  if(is.null(control$npt)) control$npt = min(n*2, n+2)
  if(is.null(control$iprint)) control$iprint = 0
  if(is.null(control$maxfun)) control$maxfun = 10000
  mm = match(names(control), c("npt", "rhobeg", "rhoend", "iprint", "maxfun"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    wrong_opts = NULL
    for(i in 1:length(idx)) wrong_opts = c(wrong_opts, names(control)[idx[i]])
    warning(paste(c("\nunidentified option(s) in solver.control:\n", wrong_opts), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  return(control)
}
################################################################################
