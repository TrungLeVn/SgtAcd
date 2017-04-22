#----------------------------------------------------------------------------------
# path
#' @importFrom sgt rsgt
#' @export acdpath
#' @title ACD simulated path
#' @usage acdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA,
#'        prereturns = NA, preresiduals = NA, preskew = NA, preshape = NA, rseed = NA,
#'        cluster = NULL, ...)
#' @param spec A spec object with fixed parameters
#' @param pre.. Must be supplied. The values of returns, residuals, skew, shape that will be used to start the simulation process
#' @return A ACDpath object that will be a simulated process for an ACD process.
#---------------------------------------------------------------------
acdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA,
                   prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA, preshape2 = NA, rseed = NA,
                   cluster = NULL, ...)
{
  UseMethod("acdpath")
}

.acdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA,
                    prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA,preshape2 = NA, rseed = NA,
                    cluster = NULL, ...){
  ans = switch(spec@model$vmodel$model,
               sGARCH = .sacdpath(spec = spec, n.sim = n.sim, n.start = n.start,
                                  m.sim = m.sim, presigma = presigma, prereturns = prereturns,
                                  preresiduals = preresiduals, preskew = preskew,
                                  preshape1 = preshape1,preshape2 = preshape2, rseed = rseed, cluster = cluster, ...),
               gjrGARCH = .gjracdpath(spec = spec, n.sim = n.sim, n.start = n.start,
                                      m.sim = m.sim, presigma = presigma, prereturns = prereturns,
                                      preresiduals = preresiduals, preskew = preskew,
                                      preshape1 = preshape1, preshape2 = preshape2, rseed = rseed, cluster = cluster, ...))
  return(ans)
}
setMethod("acdpath", signature(spec = "ACDspec"), .acdpath)
#----------------------------------------------------------------------------------
# Get Pskewness for the conditional mean equation
#----------------------------------------------------------------------------------
#' @export  Pskew
Pskew = function(lambda,kappa,nu,distribution)
{
    if(distribution == "sgt"){
      beta1 = beta(1/kappa,nu/kappa)
      beta2 = beta(2/kappa,(nu-1)/kappa)
      beta3 = beta(3/kappa,(nu-2)/kappa)
      A = beta2 * 1/sqrt(beta1) * 1/sqrt(beta3)
      S = sqrt(1 + 3*lambda^2 - 4 * A^2 * lambda^2)
      ans = 2 * lambda * A * 1/S
    }
    if(distribution == "sged"){
      gamma1 = gamma(1.0/kappa)
      gamma2 = gamma(2.0/kappa)
      gamma3 = gamma(3.0/kappa)
      A = gamma2/sqrt(gamma1 * gamma3);
      S = sqrt(1 + 3*lambda^2 - 4 * A^2 * lambda^2)
      ans = 2 * lambda * A/S;
    }
    if(distribution == "sst"){
      kappa = 2.0;
      beta1 = beta(1/kappa,nu/kappa)
      beta2 = beta(2/kappa,(nu-1)/kappa)
      beta3 = beta(3/kappa,(nu-2)/kappa)
      A = beta2 * 1/sqrt(beta1) * 1/sqrt(beta3)
      S = sqrt(1 + 3*lambda^2 - 4 * A^2 * lambda^2)
      ans = 2 * lambda * A * 1/S
    }
    return(ans);
}

#--------------------------------------
# ACD Path for GJR-GARCH process
#-------------------------------------
.gjracdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA,
                       prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA,preshape2 = NA, rseed = NA,
                       cluster = NULL, ...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  n = n.sim + n.start
  model = spec@model
  modelinc = model$modelinc
  idx = model$pidx
  ipars = model$pars
  sbounds = model$sbounds
  m = model$maxOrder
  distribution = model$dmodel$model

  if(!is.na(presigma[1])){
    presigma = as.vector(presigma)
    if(length(presigma)<m) stop(paste("\nacdpath-->error: presigma must be of length ", m, sep=""))
  } else{
    stop("\nacdpath-->error: presigma cannot be NA.")
  }
  if(!is.na(prereturns[1])){
    prereturns = as.vector(prereturns)
    if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
  } else{
    prereturns = as.numeric(tail(data, m))
  }
  if(!is.na(preresiduals[1])){
    preresiduals = as.vector(preresiduals)
    if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
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
      pretempskew = rep(ipars[idx["skcons",1],1], m)
      pretskew = logtransform(pretempskew, sbounds[1], sbounds[2])
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
# First shape1 parameter
  pretshape1 = pretempshape1 = rep(0, m)
  if(model$modelinc[19]>0)
  {
    if( is.na(preshape1[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape1 = rep(ipars[idx["sh1cons",1],1], m)
      pretshape1 = exptransform1(pretempshape1, sbounds[3], sbounds[7])
    } else{
      pretempshape1 = exptransform1(tail(as.vector(preshape1), m), sbounds[3], sbounds[7], inverse = TRUE)
      pretshape1 = tail(as.vector(preshape1), m)
    }
  }
  if(model$modelinc[12]>0){
    tshape1 = rep(ipars["shape1", 1], n+m)
  } else{
    tshape1 = c(pretshape1, rep(0, n))
  }
 # Set Pre-shape2 parameters
  pretshape2 = pretempshape2 = rep(0, m)
  if(model$modelinc[24]>0)
  {
    if( is.na(preshape2[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape2 = rep(ipars[idx["sh2cons",1],1], m)
      pretshape2 = logtransform(pretempshape2, sbounds[5], sbounds[6])
    } else{
      pretempshape2 = logtransform(tail(as.vector(preshape2), m), sbounds[5], sbounds[6], inverse = TRUE)
      pretshape2 = tail(as.vector(preshape2), m)
    }
  }
  if(model$modelinc[13]>0){
    tshape2 = rep(ipars["shape2", 1], n+m)
  } else{
    tshape2 = c(pretshape2, rep(0, n))
  }
  prePskew = rep(0,m)
  prePskew = Pskew(pretskew,pretshape1,pretshape2,distribution)
  # input vectors/matrices
  h = c(presigma^2, rep(0, n))
  x = c(prereturns, rep(0, n))
  tmpskew = c(pretempskew, rep(0, n))
  tmpshape1 = c(pretempshape1, rep(0, n))
  tmpshape2 = c(pretempshape2,rep(0,n))
  pskew = c(prePskew,rep(0,n))
  constmean = ipars[idx["mu",1]:idx["mu",2], 1]
  # Taking into account the initial m periods with prePskew and Presigma. From m+1, the varying part in conditional mean
  # depends on varying Pskew and Sigma
  constm = c(rep(constmean+prePskew*presigma,m),rep(constmean,n))
  constm = matrix(constm, ncol = m.sim, nrow = n + m)
  # TRUNG: From the simulated value of z_t, based on conditional mean, conditoinal sigma, conditional skew, conditional shape, we will start the path
  # MATRIX

  # If we do not provide preresiduals then the uncertainty will be start from n.ahead = 1
  # If we do provide the preresiduals, then the uncertainty will only start fomr n.ahead = 2
  if( !is.na(preresiduals) && !is.na(presigma) ){
    zz = preres[1:m]/presigma[1:m]
    for(j in 1:m.sim){
      z[1:m, j] = zz
    }
  } else{
    # ? Do we want the same for all m.sim? If yes, there is no uncertainty for
    # the n.sim = 1 for sigma, and higher moment (equal to their forecast value).
    # If no, then uncertainty is introduced in the n.sim=1 values.
    set.seed(sseed[1])
    for(k in 1:m){
      if(distribution == "sgt"){
        z[k,] = sgt::rsgt(n = m.sim,lambda = tskew[k],p = tshape1[k],q = tshape2[k]/tshape1[k])
      } else if(distribution == "sged"){
        z[k,] = sgt::rsgt(n = m.sim,lambda = tskew[k],p = tshape1[k],q = Inf)
      } else {
        z[k,] = sgt::rsgt(n = m.sim,lambda = tskew[k],p = 2, q = tshape2[k]/tshape1[k])
      }
    }
  }
  if(is.na(preresiduals[1])){
    preres = z[1:m, , drop = FALSE]*matrix(presigma, ncol = m.sim, nrow = m)
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
#Doing the simulation. Each period, an simulated inovation z will be generated
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

  sol = new("ACDpath",
            path = sim,
            model = model,
            seed = as.integer(sseed))
  return(sol)
}
#-------------------------------------
# ACDpath for sGARCH process
#-------------------------------------
.sacdpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA,
                     prereturns = NA, preresiduals = NA, preskew = NA, preshape1 = NA,preshape2 = NA, rseed = NA,
                     cluster = NULL, ...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  n = n.sim + n.start
  model = spec@model
  modelinc = model$modelinc
  idx = model$pidx
  ipars = model$pars
  sbounds = model$sbounds
  m = model$maxOrder
  distribution = model$dmodel$model

  if(!is.na(presigma[1])){
    presigma = as.vector(presigma)
    if(length(presigma)<m) stop(paste("\nacdpath-->error: presigma must be of length ", m, sep=""))
  } else{
    stop("\nacdpath-->error: presigma cannot be NA.")
  }
  if(!is.na(prereturns[1])){
    prereturns = as.vector(prereturns)
    if(length(prereturns)<m) stop(paste("\nuacdsim-->error: prereturns must be of length ", m, sep=""))
  } else{
    prereturns = as.numeric(tail(data, m))
  }
  if(!is.na(preresiduals[1])){
    preresiduals = as.vector(preresiduals)
    if(length(preresiduals)<m) stop(paste("\nuacdsim-->error: preresiduals must be of length ", m, sep=""))
    preres = matrix(preresiduals, nrow = m, ncol = m.sim)
  } else {
    stop("\acdpath --> error: Preresiduals cannot be NA")
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
      pretempskew = rep(ipars[idx["skcons",1],1], m)
      pretskew = logtransform(pretempskew, sbounds[1], sbounds[2])
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
  # First shape1 parameter
  pretshape1 = pretempshape1 = rep(0, m)
  if(model$modelinc[19]>0)
  {
    if( is.na(preshape1[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape1 = rep(ipars[idx["sh1cons",1],1], m)
      pretshape1 = exptransform1(pretempshape1, sbounds[3], sbounds[7])
    } else{
      pretempshape1 = exptransform1(tail(as.vector(preshape1), m), sbounds[3], sbounds[7], inverse = TRUE)
      pretshape1 = tail(as.vector(preshape1), m)
    }
  }
  if(model$modelinc[12]>0){
    tshape1 = rep(ipars["shape1", 1], n+m)
  } else{
    tshape1 = c(pretshape1, rep(0, n))
  }
  # Set Pre-shape2 parameters
  pretshape2 = pretempshape2 = rep(0, m)
  if(model$modelinc[19]>0)
  {
    if( is.na(preshape2[1]) ){
      # The tempshape[1] is the transformed shape parameter of the
      # non-time varying model from which we initiated the original fit.
      pretempshape2 = rep(ipars[idx["sh2cons",1],1], m)
      pretshape2 = logtransform(pretempshape2, sbounds[5], sbounds[6])
    } else{
      pretempshape2 = logtransform(tail(as.vector(preshape2), m), sbounds[5], sbounds[6], inverse = TRUE)
      pretshape2 = tail(as.vector(preshape2), m)
    }
  }
  if(model$modelinc[12]>0){
    tshape2 = rep(ipars["shape2", 1], n+m)
  } else{
    tshape2 = c(pretshape2, rep(0, n))
  }
  prePskew = rep(0,m)
  prePskew = Pskew(pretskew,pretshape1,pretshape2,distribution)
  # input vectors/matrices
  h = c(presigma^2, rep(0, n))
  x = c(prereturns, rep(0, n))
  tmpskew = c(pretempskew, rep(0, n))
  tmpshape1 = c(pretempshape1, rep(0, n))
  tmpshape2 = c(pretempshape2,rep(0,n))
  constmean = ipars[idx["mu",1]:idx["mu",2], 1]
  # Taking into account the initial m periods with prePskew and Presigma. From m+1, the varying part in conditional mean
  # depends on varying Pskew and Sigma
  constm = c(rep(constmean+prePskew*presigma,m),rep(constmean,n))
  constm = matrix(constm, ncol = m.sim, nrow = n + m)
  # TRUNG: From the simulated value of z_t, based on conditional mean, conditoinal sigma, conditional skew, conditional shape, we will start the path
  # MATRIX
  zz = preres[1:m]/presigma[1:m]
  for(j in 1:m.sim){
    z[1:m, j] = zz
  }
  if(is.na(preresiduals[1])){
    preres = z[1:m, , drop = FALSE]*matrix(presigma, ncol = m.sim, nrow = m)
    #preres = z[1:m, , drop = FALSE] * presigma
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
  pskewSim = matrix(0, ncol = m.sim, nrow = n.sim)

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

  sol = new("ACDpath",
            path = sim,
            model = model,
            seed = as.integer(sseed))
  return(sol)
}
