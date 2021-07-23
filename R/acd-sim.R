#------------------
# ACDsim method
#' @title ACD simulation
#' @description Produce simulation for ACD model.
#' @usage acdsim = function(fit, n.sim = 1000, n.start = 0, endpoint = NA, m.sim = 1,
#'            presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA,
#'           preshape = NA, rseed = NA,...)
#' @param fit ACDfit object
#' @param n.sim Simulation horizon
#' @param n.start Burn-in period in the simulation
#' @param endpoint This is my additional argument, specify the endpoint in the fitted series that will be used as
#'                 pre-values to start the simulation
#' @param m.sim Number of simulation serires
#' @param pre.. The user-specified pre-values of sigma, returns, residuals, skew, shape that will be used to start the simulation
#' @param rseed The seed number that must be suplied at the same legth as m.sim arguments to ensure the reproducibility
#' @return An ACDsim object with simulated values for returns, sigma, skew, shape
#' @export acdsim

#-------------------
# simulation
#-----------------

acdsim = function(fit, n.sim = 1000, n.start = 0, endpoint = NA, m.sim = 1,
                  presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA,
                  preshape1 = NA,preshape2 = NA, rseed = NA, cluster = NULL,...)
{
  UseMethod("acdsim")
}

.acdsim = function(fit, n.sim = 1000, n.start = 0,endpoint = NA, m.sim = 1,
                   presigma = NA, prereturns = NA, preresiduals = NA, preskew = NA,
                   preshape1 = NA, preshape2 = NA, rseed = NA, cluster = NULL,...){
  ans = switch(fit@model$vmodel$model,
               sGARCH = .sacdsim(fit = fit, n.sim = n.sim,endpoint = endpoint, n.start = n.start,
                                 m.sim = m.sim, presigma = presigma, prereturns = prereturns,
                                 preresiduals = preresiduals, preskew = preskew,
                                 preshape1 = preshape1, preshape2 = preshape2, rseed = rseed, cluster = cluster,...),
               gjrGARCH = .gjracdsim(fit = fit, n.sim = n.sim, endpoint = endpoint,n.start = n.start,
                                     m.sim = m.sim, presigma = presigma, prereturns = prereturns,
                                     preresiduals = preresiduals, preskew = preskew,
                                     preshape1 = preshape1, preshape2 = preshape2, rseed = rseed,cluster = cluster, ...))
  return(ans)
}

setMethod(f = "acdsim", signature(fit = "ACDfit"), .acdsim)
#--------------------------------------------
# acdsim for gjrGARCH process
#--------------------------------------------
.gjracdsim = function(fit, n.sim = 1000, n.start = 0,endpoint = NA, m.sim = 1, presigma = NA,
                      prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA,preshape2 = NA, rseed = NA,
                      cluster = NULL,  ...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  n = n.sim + n.start
  model = fit@model
  if(is.na(endpoint)){
    endpoint = model$modeldata$T
  }
  modelinc = model$modelinc
  idx = model$pidx
  ipars = model$pars
  sbounds = model$sbounds
  m = model$maxOrder
  distribution = model$dmodel$model
  # check if necessary the external regressor forecasts provided first

  if(!is.na(presigma[1])){
    presigma = as.vector(presigma)
    if(length(presigma)<m) stop(paste("\nacdsim-->error: presigma must be of length ", m, sep=""))
  } else{
    presigma = as.vector(tail(as.numeric(sigmaAcd(fit)[1:endpoint]), m))
  }
  if(!is.na(prereturns[1])){
    prereturns = as.vector(prereturns)
    if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
  } else{
    prereturns = as.vector(tail(model$modeldata$data[1:endpoint], m))
  }
  if(!is.na(preresiduals[1])){
    preresiduals = as.vector(preresiduals)
    if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
    preres = matrix(preresiduals, nrow = m, ncol = m.sim)
  } else{
    preresiduals = as.vector(tail(residualsAcd(fit)[1:endpoint], m))
    preres = matrix(preresiduals, nrow = m, ncol = m.sim)
  }
  # Random Samples from the Distribution are calculated at every recursion in the
  # c-code as they depend on the actual time-varying skew & shape
  z = matrix(0, ncol = m.sim, nrow = n.sim + n.start)
  z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
  # z = matrix(0, ncol = m.sim, nrow = n.sim+n.start)
  pretskew = pretempskew = rep(0, m)

  if(model$modelinc[14]>0)
  {
    if( is.na(preskew[1]) ){
      # The tempskew[1] is the transformed skew parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempskew = as.vector(tail(fit@fit$tempskew[1:endpoint], m))
      pretskew = as.vector(tail(fit@fit$tskew[1:endpoint], m))
    } else{
      # preskew is provided un-transformed
      pretempskew = logtransform(tail(as.vector(preskew), m), sbounds[1], sbounds[2], inverse = TRUE)
      pretskew = tail(as.vector(preskew), m)
    }
  }
  if(model$modelinc[11]>0){
    tskew = rep(ipars["skew", 1], n+m)
    pretskew = ipars["skew",1]
  } else{
    tskew = c(pretskew, rep(0, n))
  }

  pretshape1 = pretempshape1 = pretshape2 = pretempshape2 = rep(0, m)
  if(model$modelinc[19]>0)
  {
    if( is.na(preshape1[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape1 = as.vector(tail(fit@fit$tempshape1[1:endpoint], m))
      pretshape1 = as.vector(tail(fit@fit$tshape1[1:endpoint], m))
    } else{
      pretempshape1 = exptransform1(tail(as.vector(preshape1), m), sbounds[3], sbounds[7], inverse = TRUE)
      pretshape1 = tail(as.vector(preshape1), m)
    }
  }
  if(model$modelinc[12]>0){
    tshape1 = rep(ipars["shape1", 1], n+m)
    pretshape1 = ipars["shape1",1]
  } else if(sum(model$modelinc[11:13])==2){
    tshape1 = rep(0,n+m)
    pretshape1 = 0
  } else{
    tshape1 = c(as.vector(pretshape1), rep(0, n))
  }
  if(model$modelinc[24]>0)
  {
    if( is.na(preshape2[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape2 = as.vector(tail(fit@fit$tempshape2[1:endpoint], m))
      pretshape2 = as.vector(tail(fit@fit$tshape2[1:endpoint], m))
    } else{
      pretempshape2 = logtransform(tail(as.vector(preshape2), m), sbounds[5], sbounds[6], inverse = TRUE)
      pretshape2 = tail(as.vector(preshape2), m)
    }
  }
  if(model$modelinc[13]>0){
    tshape2 = rep(ipars["shape2", 1], n+m)
    pretshape2 = ipars["shape2",1]
  } else if(sum(model$modelinc[11:12]) == 2){
    tshape2 = rep(0,n+m)
    pretshape2 = 0
  } else{
    tshape2 = c(as.vector(pretshape2), rep(0, n))
  }
  # input vectors/matrices
  prePskew = rep(0,m)
  prePskew = Pskew(pretskew,pretshape1,pretshape2,distribution)
  h = c(presigma^2, rep(0, n))
  x = c(prereturns, rep(0, n))
  tmpskew = c(pretempskew, rep(0, n))
  tmpshape1 = c(pretempshape1, rep(0, n))
  tmpshape2 = c(pretempshape2, rep(0, n))
  pskew = c(prePskew,rep(0,n))
  constmean = ipars[idx["mu",1]:idx["mu",2], 1]
  # Taking into account the initial m periods with prePskew and Presigma. From m+1, the varying part in conditional mean
  # depends on varying Pskew and Sigma
  if(modelinc[36]>0){
    constm = c(rep(constmean+prePskew*presigma,m),rep(constmean,n))
  } else{
    constm = c(rep(constmean,n+m))
  }
  constm = matrix(constm, ncol = m.sim, nrow = n + m)
  # MATRIX
  zz = preres[1:m]/presigma[1:m]
  for(j in 1:m.sim){
    z[1:m, j] = zz
  }
  pre_nres = matrix(NA,ncol = m.sim,nrow = m)
  for(i in 1:m){
    for(ii in 1:m.sim){
      if(preres[i,ii] < 0){
        pre_nres[i,ii] = preres[i,ii]*preres[i,ii]
      }else{
        pre_nres[i,ii] =0
      }
    }
  }
  res = rbind( preres, matrix(0, ncol = m.sim, nrow = n) )
  nres = rbind(pre_nres,matrix(0,ncol = m.sim,nrow = n))
  # outpus matrices
  sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
  seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
  residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
  skewSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
  shape1Sim 		= matrix(0, ncol = m.sim, nrow = n.sim)
  shape2Sim     = matrix(0, ncol = m.sim, nrow = n.sim)
  zSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
  pskewSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, library(SgtAcd))
    parallel::clusterExport(cluster, c("modelinc", "ipars", "idx","x", "constm", "h","z", "res","nres",
                                       "tmpskew", "tmpshape1","tmpshape2", "tskew", "tshape1","tshape2", "sbounds","pskew", "sseed",
                                       "n", "m","n.sim"), envir = environment())
    S = parallel::parLapply(cluster, 1:m.sim, function(i){
      set.seed(sseed[i])
      tmp = try(.C("gjracdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]),
                   idx = as.integer(idx[,1]-1),x = as.double(x),constm = as.double(constm[,i]), h = as.double(h), z = as.double(z[,i]),
                   res = as.double(res[,i]), e = as.double(res[,i]*res[,i]),nres = as.double(nres[,i]),
                   tempskew = as.double(tmpskew), tempshape1 = as.double(tmpshape1),tempshape2 = as.double(tmpshape2),
                   tskew = as.double(tskew), tshape1 = as.double(tshape1),tshape2 = as.double(tshape2),
                   sbounds = as.double(sbounds),pskew = as.double(pskew), T = as.integer(n+m), m = as.integer(m),
                   PACKAGE = "SgtAcd"), silent = TRUE)
      #tmpPskew = Pskew(lambda = tmp$tskew,kappa = tmp$tshape1,nu = tmp$tshape2,distribution = distribution)
      #ans2 = rugarch:::.armaxsim(modelinc[1:3], ipars = ipars, idx = idx, constm = constm[,i],
      #                           x = x, res = tmp$res, T = n + m, m)
      #seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
      ret = cbind(tail(tmp$x,n.sim),tail(sqrt(tmp$h), n.sim), tail(tmp$res, n.sim), tail(tmp$tskew, n.sim),
                  tail(tmp$tshape1, n.sim),tail(tmp$tshape2,n.sim), tail(tmp$z, n.sim),tail(tmp$pskew,n.sim))
      return(ret)
    })
    seriesSim = sapply(S, function(x) x[,1])
    sigmaSim = sapply(S, function(x) x[,2])
    residSim = sapply(S, function(x) x[,3])
    skewSim = sapply(S, function(x) x[,4])
    shape1Sim = sapply(S, function(x) x[,5])
    shape2Sim = sapply(S, function(x) x[,6])
    zSim = sapply(S, function(x) x[,7])
    pskewSim = sapply(S, function(x) x[,8])
  } else{
    for(i in 1:m.sim){
      set.seed(sseed[i])
      tmp = try(.C("gjracdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]),
                   idx = as.integer(idx[,1]-1),x = as.double(x),constm = as.double(constm[,i]), h = as.double(h), z = as.double(z[,i]),
                   res = as.double(res[,i]), e = as.double(res[,i]*res[,i]),nres = as.double(nres[,i]),
                   tempskew = as.double(tmpskew), tempshape1 = as.double(tmpshape1),tempshape2 = as.double(tmpshape2),
                   tskew = as.double(tskew), tshape1 = as.double(tshape1),tshape2 = as.double(tshape2),
                   sbounds = as.double(sbounds),pskew = as.double(pskew), T = as.integer(n+m), m = as.integer(m),
                   PACKAGE = "SgtAcd"), silent = TRUE)
      #ans2 = rugarch:::.armaxsim(modelinc[1:3], ipars = ipars, idx = idx, constm = constm[,i],
      #                           x = x, res = tmp$res, T = n + m, m)
      seriesSim[,i] = tail(tmp$x, n.sim)
      sigmaSim[,i] 	 = tail(sqrt(tmp$h), n.sim)
      residSim[,i] 	 = tail(tmp$res, n.sim)
      skewSim[,i] 	 = tail(tmp$tskew, n.sim)
      shape1Sim[,i] 	 = tail(tmp$tshape1, n.sim)
      shape2Sim[,i]    = tail(tmp$tshape2, n.sim)
      zSim[,i] 		 = tail(tmp$z, n.sim)
      pskewSim[,i]   = tail(tmp$pskew, n.sim)
    }
  }
  sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, skewSim = skewSim,
             shape1Sim = shape1Sim,shape2Sim = shape2Sim, zSim = zSim,pskewSim = pskewSim)
  sim$n.sim  = n.sim
  sim$m.sim  = m.sim

  sol = new("ACDsim",
            simulation = sim,
            model = model,
            seed = as.integer(sseed))
  return(sol)
}


#--------------------------------------------
# acdsim for sGARCH process
#--------------------------------------------

.sacdsim = function(fit, n.sim = 1000, n.start = 0,endpoint = NA, m.sim = 1, presigma = NA,
                    prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA,preshape2 = NA, rseed = NA,
                    cluster = NULL,...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  n = n.sim + n.start
  model = fit@model
  if(is.na(endpoint)){
    endpoint = model$modeldata$T
  }
  modelinc = model$modelinc
  idx = model$pidx
  ipars = model$pars
  sbounds = model$sbounds
  m = model$maxOrder
  distribution = model$dmodel$model
  # check if necessary the external regressor forecasts provided first

  if(!is.na(presigma[1])){
    presigma = as.vector(presigma)
    if(length(presigma)<m) stop(paste("\nacdpath-->error: presigma must be of length ", m, sep=""))
  } else{
    presigma = as.vector(tail(as.numeric(sigmaAcd(fit)[1:endpoint]), m)) #take the last sigma as pre-sigma
  }
  if(!is.na(prereturns[1])){
    prereturns = as.vector(prereturns)
    if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
  } else{
    prereturns = as.vector(tail(model$modeldata$data[1:endpoint], m)) #take the last return (r_{t=1}) as pre-sigma
  }
  if(!is.na(preresiduals[1])){
    preresiduals = as.vector(preresiduals)
    if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
    preres = matrix(preresiduals, nrow = m, ncol = m.sim)
  } else{
    preresiduals = as.vector(tail(residualsAcd(fit)[1:endpoint], m))
    preres = matrix(preresiduals, nrow = m, ncol = m.sim)
  }

  # Random Samples from the Distribution are calculated at every recursion in the
  # c-code as they depend on the actual time-varying skew & shape
  z = matrix(0, ncol = m.sim, nrow = n.sim + n.start)
  z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
  # z = matrix(0, ncol = m.sim, nrow = n.sim+n.start)
  pretskew = pretempskew = rep(0, m)

  if(model$modelinc[14]>0)
  {
    if( is.na(preskew[1]) ){
      # The tempskew[1] is the transformed skew parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempskew = as.vector(tail(fit@fit$tempskew[1:endpoint], m))
      pretskew = as.vector(tail(fit@fit$tskew[1:endpoint], m))
    } else{
      # preskew is provided un-transformed
      pretempskew = logtransform(tail(as.vector(preskew), m), sbounds[1], sbounds[2], inverse = TRUE)
      pretskew = tail(as.vector(preskew), m)
    }
  }
  if(model$modelinc[11]>0){
    tskew = rep(ipars["skew", 1], n+m)
  } else{
    tskew = c(pretskew, rep(0, n))
  }

  pretshape1 = pretempshape1 = pretshape2 = pretempshape2 = rep(0, m)
  if(model$modelinc[19]>0)
  {
    if( is.na(preshape1[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape1 = as.vector(tail(fit@fit$tempshape1[1:endpoint], m))
      pretshape1 = as.vector(tail(fit@fit$tshape1[1:endpoint], m))
    } else{
      pretempshape1 = exptransform1(tail(as.vector(preshape1), m), sbounds[3], sbounds[7], inverse = TRUE)
      pretshape1 = tail(as.vector(preshape1), m)
    }
  }
  if(model$modelinc[12]>0){
    tshape1 = rep(ipars["shape1", 1], n+m)
  } else{
    tshape1 = c(as.vector(pretshape1), rep(0, n))
  }
  if(model$modelinc[24]>0)
  {
    if( is.na(preshape2[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape2 = as.vector(tail(fit@fit$tempshape2[1:endpoint], m))
      pretshape2 = as.vector(tail(fit@fit$tshape2[1:endpoint], m))
    } else{
      pretempshape2 = logtransform(tail(as.vector(preshape2), m), sbounds[5], sbounds[6], inverse = TRUE)
      pretshape2 = tail(as.vector(preshape2), m)
    }
  }
  if(model$modelinc[13]>0){
    tshape2 = rep(ipars["shape2", 1], n+m)
  } else{
    tshape2 = c(as.vector(pretshape2), rep(0, n))
  }
  # input vectors/matrices
  prePskew = rep(0,m)
  if(pretskew != 0){
    prePskew = Pskew(pretskew,pretshape1,pretshape2,distribution)
  }
  h = c(presigma^2, rep(0, n))
  x = c(prereturns, rep(0, n))
  tmpskew = c(pretempskew, rep(0, n))
  tmpshape1 = c(pretempshape1, rep(0, n))
  tmpshape2 = c(pretempshape2, rep(0, n))
  pskew = c(prePskew,rep(0,n))
  constmean = ipars[idx["mu",1]:idx["mu",2], 1]
  # Taking into account the initial m periods with prePskew and Presigma. From m+1, the varying part in conditional mean
  # depends on varying Pskew and Sigma
  if(modelinc[36]>0){
    constm = c(rep(constmean+prePskew*presigma,m),rep(constmean,n))
  } else{
    constm = c(rep(constmean,n+m))
  }
  constm = matrix(constm, ncol = m.sim, nrow = n + m)
  # MATRIX
  zz = preres[1:m]/presigma[1:m] #z_{t-1}
  for(j in 1:m.sim){
    z[1:m, j] = zz
  }
  res = rbind( preres, matrix(0, ncol = m.sim, nrow = n) )
  # outpus matrices
  sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
  seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
  residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
  skewSim 		= matrix(0, ncol = m.sim, nrow = n.sim)
  shape1Sim 		= matrix(0, ncol = m.sim, nrow = n.sim)
  shape2Sim     = matrix(0, ncol = m.sim, nrow = n.sim)
  zSim 			= matrix(0, ncol = m.sim, nrow = n.sim)
  pskewSim 			= matrix(0, ncol = m.sim, nrow = n.sim)

  if(!is.null(cluster)){
    parallel::clusterEvalQ(cluster, library(SgtAcd))
    parallel::clusterExport(cluster, c("modelinc", "ipars", "idx","x", "constm", "h","z", "res",
                                       "tmpskew", "tmpshape1","tmpshape2", "tskew", "tshape1","tshape2", "sbounds","pskew", "sseed",
                                       "n", "m","n.sim"), envir = environment())
    S = parallel::parLapply(cluster, 1:m.sim, function(i){
      set.seed(sseed[i])
      tmp = try(.C("sacdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]),
                   idx = as.integer(idx[,1]-1),x = as.double(x),constm = as.double(constm[,i]), h = as.double(h), z = as.double(z[,i]),
                   res = as.double(res[,i]), e = as.double(res[,i]*res[,i]),
                   tempskew = as.double(tmpskew), tempshape1 = as.double(tmpshape1),tempshape2 = as.double(tmpshape2),
                   tskew = as.double(tskew), tshape1 = as.double(tshape1),tshape2 = as.double(tshape2),
                   sbounds = as.double(sbounds),pskew = as.double(pskew), T = as.integer(n+m), m = as.integer(m),
                   PACKAGE = "SgtAcd"), silent = TRUE)
      #tmpPskew = Pskew(lambda = tmp$tskew,kappa = tmp$tshape1,nu = tmp$tshape2,distribution = distribution)
      #ans2 = rugarch:::.armaxsim(modelinc[1:3], ipars = ipars, idx = idx, constm = constm[,i],
      #                           x = x, res = tmp$res, T = n + m, m)
      #seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
      ret = cbind(tail(tmp$x,n.sim),tail(sqrt(tmp$h), n.sim), tail(tmp$res, n.sim), tail(tmp$tskew, n.sim),
                  tail(tmp$tshape1, n.sim),tail(tmp$tshape2,n.sim), tail(tmp$z, n.sim), tail(tmp$pskewSim,n.sim))
      return(ret)
    })
    seriesSim = sapply(S, function(x) x[,1])
    sigmaSim = sapply(S, function(x) x[,2])
    residSim = sapply(S, function(x) x[,3])
    skewSim = sapply(S, function(x) x[,4])
    shape1Sim = sapply(S, function(x) x[,5])
    shape2Sim = sapply(S, function(x) x[,6])
    zSim = sapply(S, function(x) x[,7])
    pskewSim = sapply(S, function(x) x[,8])
  } else{
    for(i in 1:m.sim){
      set.seed(sseed[i])
      tmp = try(.C("sacdsimC", model = as.integer(modelinc), pars = as.double(ipars[,1]),
                   idx = as.integer(idx[,1]-1),x = as.double(x),constm = as.double(constm[,i]), h = as.double(h), z = as.double(z[,i]),
                   res = as.double(res[,i]), e = as.double(res[,i]*res[,i]),
                   tempskew = as.double(tmpskew), tempshape1 = as.double(tmpshape1),tempshape2 = as.double(tmpshape2),
                   tskew = as.double(tskew), tshape1 = as.double(tshape1),tshape2 = as.double(tshape2),
                   sbounds = as.double(sbounds),pskew = as.double(pskew), T = as.integer(n+m), m = as.integer(m),
                   PACKAGE = "SgtAcd"), silent = TRUE)
      #ans2 = rugarch:::.armaxsim(modelinc[1:3], ipars = ipars, idx = idx, constm = constm[,i],
      #                           x = x, res = tmp$res, T = n + m, m)
      seriesSim[,i] = tail(tmp$x, n.sim)
      sigmaSim[,i] 	 = tail(sqrt(tmp$h), n.sim)
      residSim[,i] 	 = tail(tmp$res, n.sim)
      skewSim[,i] 	 = tail(tmp$tskew, n.sim)
      shape1Sim[,i] 	 = tail(tmp$tshape1, n.sim)
      shape2Sim[,i]    = tail(tmp$tshape2, n.sim)
      zSim[,i] 		 = tail(tmp$z, n.sim)
      pskewSim[,i] = tail(tmp$pskew,n.sim)
    }
  }
  sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, skewSim = skewSim,
             shape1Sim = shape1Sim,shape2Sim = shape2Sim, zSim = zSim,pskewSim = pskewSim)
  sim$n.sim  = n.sim
  sim$m.sim  = m.sim

  sol = new("ACDsim",
            simulation = sim,
            model = model,
            seed = as.integer(sseed))
  return(sol)
}
