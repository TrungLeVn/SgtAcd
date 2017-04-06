##--------------------
# Generate starting values for parameteres.
# The starting parameters for mean equation, if not specified, is the fitted value of ARMA process
# The strating parameters for variance equation, if not specified, is the fitted value of GARCh process
#------
TinY = 1e-08
#-------
# Starting parameters for mean equation using ARMA fit
#-------
.meqstartpars = function(pars, arglist) {
  data = zoo::zoo(arglist$data,order.by = arglist$index);
  N = length(as.numeric(unlist(data)))
  model = arglist$model
  modelinc = model$modelinc
  modeldesc = model$modeldesc
  idx = model$pidx
  tmph = 0
  # Get the estimation for ARMA specfication in mean equation.  To be used to find the starting parameters in variance equation
  if (sum(modelinc[4:6]) > 0){
    tempspec = acdspec(mean.model = list(armaOrder = FALSE, skm = FALSE, shm = FALSE),
                       distribution.model = model$dmodel, variance.model = model$vmodel)
    tempfit = acdfit(tempspec,data = data)
    if(tempfit@fit$convergence!=0){
      tempfit = acdfit(tempspec, data = data,solver = "mssolnp")
      if(tempfit@fit$convergence!=0){
        stop("\nacdfit-->error: could not find appropriate starting values for recursion\n")
      }
    }
    tmph = cbind(archm = tempfit$sigma,skm = tempfit$skewness,kum = tempfit$kurtosis)
    mexdata = NULL
    if(modelinc[4] > 0){
      mexdata = cbind(mexdata,tmph[,"archm"])
    }
    if(modelin[5]>0){
      mexdata = cbind(mexdata,tmph[,"skm"])
    }
    if(modelinc[6]>0){
      mexdata = cbind(mexdata,tmph[,"kum"])
    }
    y = coredata(data);
    fit.mean = lm(y ~ mexdata)
    pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean$coef["(Intercept)"]
    if(modelinc[4]>0){
      pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean$coef["archm"]
    }
    if(modelinc[5]>0){
      pars[idx["skm", 1]:idx["skm", 2], 1] = fit.mean$coef["skm"]
    }
    if(modelinc[6]>0){
      pars[idx["kum", 1]:idx["kum", 2], 1] = fit.mean$coef["kum"]
    }
  } else if (modelinc[2] > 0 | modelinc[3] > 0) {
    ttemp = arima(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], method = "CSS")
    fit.mean = ttemp$coef
    if (modelinc[1] > 0) {
      pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean["intercept"]
    }
    if (modelinc[2] > 0){
      pars[idx["ar", 1]:idx["ar", 2], 1] = fit.mean[c(paste("ar", 1:modelinc[2], sep = ""))]
    }
    if (modelinc[3] > 0){
      pars[idx["ma", 1]:idx["ma", 2], 1] = fit.mean[c(paste("ma", 1:modelinc[3], sep = ""))]
    }
   } else{
     pars[idx["mu", 1]:idx["mu", 2], 1] = 0
   }
  arglist$tmph = tmph
  return(list(pars = pars, arglist = arglist))
}
#------
# Dealing with fixed parameters in specification, starting values and define the bounds for mean equation parameters
#------
.meqstart = function(pars, arglist) {
  dscale = arglist$dscale
  model = arglist$model
  start.pars = model$start.pars
  start.names = names(start.pars)
  fixed.pars = model$fixed.pars
  fixed.names = names(fixed.pars)
  idx = model$pidx
  modelinc = model$modelinc
  data = arglist$data
  # mu
  if (modelinc[1] > 0) {
    pars[idx["mu", 1]:idx["mu", 2], 5] = -10 * abs(mean(data))
    pars[idx["mu", 1]:idx["mu", 2], 6] = 10 * abs(mean(data))
    if (!is.null(start.pars$mu))
      pars[idx["mu", 1]:idx["mu", 2], 1] = start.pars$mu[1]/dscale
    if (any(substr(fixed.names, 1, 2) == "mu")) {
      pars[idx["mu", 1]:idx["mu", 2], 1] = as.numeric(fixed.pars$mu)
      pars[idx["mu", 1]:idx["mu", 2], 5] = fixed.pars$mu
      pars[idx["mu", 1]:idx["mu", 2], 6] = fixed.pars$mu
    }
  }

  if (modelinc[2] > 0) {
    arnames = paste("ar", 1:modelinc[2], sep = "")
    pars[idx["ar", 1]:idx["ar", 2], 5] = -1 + TinY
    pars[idx["ar", 1]:idx["ar", 2], 6] = 1 - TinY
    if (any(substr(start.names, 1, 2) == "ar")) {
      j = which(substr(start.names, 1, 2) == "ar")
      armatch = charmatch(start.names[j], arnames)
      pars[arnames[armatch], 1] = as.numeric(start.pars[j])
    }
    if (any(substr(fixed.names, 1, 2) == "ar")) {
      j = which(substr(fixed.names, 1, 2) == "ar")
      armatch = charmatch(fixed.names[j], arnames)
      pars[arnames[armatch], 1] = as.numeric(fixed.pars[j])
      pars[arnames[armatch], 5] = as.numeric(fixed.pars[j])
      pars[arnames[armatch], 6] = as.numeric(fixed.pars[j])
    }
  }
  # ma
  if (modelinc[3] > 0) {
    manames = paste("ma", 1:modelinc[3], sep = "")
    pars[idx["ma", 1]:idx["ma", 2], 5] = -1 + TinY
    pars[idx["ma", 1]:idx["ma", 2], 6] = 1 - TinY
    if (any(substr(start.names, 1, 2) == "ma")) {
      j = which(substr(start.names, 1, 2) == "ma")
      mamatch = charmatch(start.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(start.pars[j])
    }
    if (any(substr(fixed.names, 1, 2) == "ma")) {
      j = which(substr(fixed.names, 1, 2) == "ma")
      mamatch = charmatch(fixed.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 5] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 6] = as.numeric(fixed.pars[j])
    }
  }
  # archm
  if (modelinc[4] > 0) {
    manames = paste("archm", 1:modelinc[4], sep = "")
    pars[idx["archm", 1]:idx["archm", 2], 5] = -1 + TinY
    pars[idx["archm", 1]:idx["archm", 2], 6] = 1 - TinY
    if (any(substr(start.names, 1, 2) == "archm")) {
      j = which(substr(start.names, 1, 2) == "archm")
      mamatch = charmatch(start.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(start.pars[j])
    }
    if (any(substr(fixed.names, 1, 2) == "archm")) {
      j = which(substr(fixed.names, 1, 2) == "archm")
      mamatch = charmatch(fixed.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 5] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 6] = as.numeric(fixed.pars[j])
    }
  }
  # skm
  if (modelinc[5] > 0) {
    manames = paste("skm", 1:modelinc[5], sep = "")
    pars[idx["skm", 1]:idx["skm", 2], 5] = -5 + TinY
    pars[idx["skm", 1]:idx["skm", 2], 6] = 5 - TinY
    if (any(substr(start.names, 1, 2) == "skm")) {
      j = which(substr(start.names, 1, 2) == "skm")
      mamatch = charmatch(start.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(start.pars[j])
    }
    if (any(substr(fixed.names, 1, 2) == "skm")) {
      j = which(substr(fixed.names, 1, 2) == "skm")
      mamatch = charmatch(fixed.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 5] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 6] = as.numeric(fixed.pars[j])
    }
  }
  # kum
  if (modelinc[6] > 0) {
    manames = paste("kum", 1:modelinc[6], sep = "")
    pars[idx["kum", 1]:idx["kum", 2], 5] = -5 + TinY
    pars[idx["kum", 1]:idx["kum", 2], 6] = 5 - TinY
    if (any(substr(start.names, 1, 2) == "kum")) {
      j = which(substr(start.names, 1, 2) == "kum")
      mamatch = charmatch(start.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(start.pars[j])
    }
    if (any(substr(fixed.names, 1, 2) == "kum")) {
      j = which(substr(fixed.names, 1, 2) == "kum")
      mamatch = charmatch(fixed.names[j], manames)
      pars[manames[mamatch], 1] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 5] = as.numeric(fixed.pars[j])
      pars[manames[mamatch], 6] = as.numeric(fixed.pars[j])
    }
  }
  return(pars)
}
#-----
# starting parameters w.r.t model specification, getting from the GARCH fit
#-----
acdstart = function(pars, arglist) {
  tmp = .meqstartpars(pars, arglist)
  pars = tmp$pars
  arglist = tmp$arglist
  model = arglist$model
  dscale = arglist$dscale
  modelinc = model$modelinc
  start.pars = model$start.pars
  fixed.pars = model$fixed.pars
  idx = model$pidx
  fixed.names = names(fixed.pars)
  start.names = names(start.pars)
  pars = .meqstart(pars, arglist)
  shape10 = arglist$shape10
  shape20 = arglist$shape20
  skew0 = arglist$skew0
  data = arglist$data
  garchenv = arglist$garchenv
  assign("garchLL", NA, envir = garchenv)
  dbounds = arglist$sbounds
  bounds = .DistributionBounds(model$dmodel$model)
  numdispar = sum(modelinc[c(14,19,24)])
  arglist$skhEst = rep(0,1,numdispar)
  # to be used for starting the recursion in case of problems
  if (modelinc[7] > 0) {
    pars[idx["omega", 1]:idx["omega", 2], 5] = var(data)/1e+05
    pars[idx["omega", 1]:idx["omega", 2], 6] = var(data) * 1e+05
    if (is.null(start.pars$omega))
      pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
    if (any(substr(fixed.names, 1, 5) == "omega")) {
      pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
      pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
      pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
    }
  }
  if (modelinc[8] > 0) {
    pars[idx["alpha", 1]:idx["alpha", 2], 5] = 0 + TinY
    pars[idx["alpha", 1]:idx["alpha", 2], 6] = 0.6 -TinY
    if (is.null(start.pars$alpha))
      pars[idx["alpha", 1]:idx["alpha", 2], 1] = 0.05 else pars[idx["alpha", 1]:idx["alpha", 2], 1] = start.pars$alpha[1]
    if (any(substr(fixed.names, 1, 5) == "alpha")) {
      pars[idx["alpha", 1]:idx["alpha", 2], 1] = as.numeric(fixed.pars$alpha)
      pars[idx["alpha", 1]:idx["alpha", 2], 5] = fixed.pars$alpha
      pars[idx["alpha", 1]:idx["alpha", 2], 6] = fixed.pars$alpha
    }
  }
  if (modelinc[9] > 0) {
    pars[idx["beta", 1]:idx["beta", 2], 5] = 0.2 + TinY
    pars[idx["beta", 1]:idx["beta", 2], 6] = 1-TinY
    if (is.null(start.pars$beta))
      pars[idx["beta", 1]:idx["beta", 2], 1] = 0.7 else pars[idx["beta", 1]:idx["beta", 2], 1] = start.pars$beta[1]
    if (any(substr(fixed.names, 1, 5) == "beta")) {
      pars[idx["beta", 1]:idx["beta", 2], 1] = as.numeric(fixed.pars$beta)
      pars[idx["beta", 1]:idx["beta", 2], 5] = fixed.pars$beta
      pars[idx["beta", 1]:idx["beta", 2], 6] = fixed.pars$beta
    }
  }
  if (modelinc[10] > 0) {
    pars[idx["gamma", 1]:idx["gamma", 2], 5] = -.5+ TinY
    pars[idx["gamma", 1]:idx["gamma", 2], 6] = 1-TinY
    if (is.null(start.pars$gamma))
      pars[idx["gamma", 1]:idx["gamma", 2], 1] = 0.5 else pars[idx["gamma", 1]:idx["gamma", 2], 1] = start.pars$gamma[1]
    if (any(substr(fixed.names, 1, 5) == "gamma")) {
      pars[idx["gamma", 1]:idx["gamma", 2], 1] = as.numeric(fixed.pars$gamma)
      pars[idx["gamma", 1]:idx["gamma", 2], 5] = fixed.pars$gamma
      pars[idx["gamma", 1]:idx["gamma", 2], 6] = fixed.pars$gamma
    }
  }
  if (modelinc[11] > 0) {
    if (is.na(pars[idx["skew", 1]:idx["skew", 2], 5]))
      pars[idx["skew", 1]:idx["skew", 2], 5] = dbounds[1]
    if (is.na(pars[idx["skew", 1]:idx["skew", 2], 6]))
      pars[idx["skew", 1]:idx["skew", 2], 6] = dbounds[2]
    if (is.null(start.pars$skew))
      pars[idx["skew", 1]:idx["skew", 2], 1] = bounds$skew else pars[idx["skew", 1]:idx["skew", 2], 1] = start.pars$skew[1]
      if (any(substr(fixed.names, 1, 4) == "skew")) {
        pars[idx["skew", 1]:idx["skew", 2], 1] = as.numeric(fixed.pars$skew)
        pars[idx["skew", 1]:idx["skew", 2], 5] = as.numeric(fixed.pars$skew)
        pars[idx["skew", 1]:idx["skew", 2], 6] = as.numeric(fixed.pars$skew)
      }
  }
  if (modelinc[12] > 0) {
    if (is.na(pars[idx["shape1", 1]:idx["shape1", 2], 5]))
      pars[idx["shape1", 1]:idx["shape1", 2], 5] = dbounds[3]
    if (is.na(pars[idx["shape1", 1]:idx["shape1", 2], 6]))
      pars[idx["shape1", 1]:idx["shape1", 2], 6] = dbounds[4]
    if (is.null(start.pars$shape1))
      pars[idx["shape1", 1]:idx["shape1", 2], 1] = bounds$shape1 else pars[idx["shape1", 1]:idx["shape1", 2], 1] = start.pars$shape1[1]
      if (any(substr(fixed.names, 1, 6) == "shape1")) {
        pars[idx["shape1", 1]:idx["shape1", 2], 1] = as.numeric(fixed.pars$shape1)
        pars[idx["shape1", 1]:idx["shape1", 2], 5] = NA
        pars[idx["shape1", 1]:idx["shape1", 2], 6] = NA
      }
  }
  if (modelinc[13] > 0) {
    if (is.na(pars[idx["shape2", 1]:idx["shape2", 2], 5]))
      pars[idx["shape2", 1]:idx["shape2", 2], 5] = dbounds[5]
    if (is.na(pars[idx["shape2", 1]:idx["shape2", 2], 6]))
      pars[idx["shape2", 1]:idx["shape2", 2], 6] = dbounds[6]
    if (is.null(start.pars$shape2))
      pars[idx["shape2", 1]:idx["shape2", 2], 1] = bounds$shape2 else pars[idx["shape2", 1]:idx["shape2", 2], 1] = start.pars$shape2[1]
      if (any(substr(fixed.names, 1, 6) == "shape2")) {
        pars[idx["shape2", 1]:idx["shape2", 2], 1] = as.numeric(fixed.pars$shape2)
        pars[idx["shape2", 1]:idx["shape2", 2], 5] = NA
        pars[idx["shape2", 1]:idx["shape2", 2], 6] = NA
      }
  }
  if (modelinc[14] > 0) {
    uncskew = bounds$skew
    xskew = .acdskewbounds(modelinc[15:17], uncskew, model$dmodel$model, dbounds[1:2])
    arglist$skhEst[1] = xskew$sk0
    if (is.na(pars[idx["skcons", 1]:idx["skcons", 2], 5]))
      pars[idx["skcons", 1]:idx["skcons", 2], 5] = xskew$skewpar.LB[1]
    if (is.na(pars[idx["skcons", 1]:idx["skcons", 2], 6]))
      pars[idx["skcons", 1]:idx["skcons", 2], 6] = xskew$skewpar.UB[1]
    if (is.null(start.pars$skcons))
      pars[idx["skcons", 1]:idx["skcons", 2], 1] = xskew$skewpars[1] else pars[idx["skcons", 1]:idx["skcons", 2], 1] = start.pars$skcons[1]
    if (any(substr(fixed.names, 1, 6) == "skcons")) {
      pars[idx["skcons", 1]:idx["skcons", 2], 1] = as.numeric(fixed.pars$skcons)
      pars[idx["skcons", 1]:idx["skcons", 2], 5] = as.numeric(fixed.pars$skcons)
      pars[idx["skcons", 1]:idx["skcons", 2], 6] = as.numeric(fixed.pars$skcons)
    }
    if (modelinc[15] > 0) {
      if (is.na(pars[idx["skalpha", 1]:idx["skalpha", 2], 5]))
        pars[idx["skalpha", 1]:idx["skalpha", 2], 5] = xskew$skewpar.LB[2]
      if (is.na(pars[idx["skalpha", 1]:idx["skalpha", 2], 6]))
        pars[idx["skalpha", 1]:idx["skalpha", 2], 6] = xskew$skewpar.UB[2]
      if (is.null(start.pars$skalpha))
        pars[idx["skalpha", 1]:idx["skalpha", 2], 1] = xskew$skewpars[2] else pars[idx["skalpha", 1]:idx["skalpha", 2], 1] = start.pars$skalpha[1]
        if (any(substr(fixed.names, 1, 6) == "skalpha")) {
          pars[idx["skalpha", 1]:idx["skalpha", 2], 1] = as.numeric(fixed.pars$skalpha)
          pars[idx["skalpha", 1]:idx["skalpha", 2], 5] = as.numeric(fixed.pars$skalpha)
          pars[idx["skalpha", 1]:idx["skalpha", 2], 6] = as.numeric(fixed.pars$skalpha)
        }
    }
    if (modelinc[16] > 0) {
      if (is.na(pars[idx["skgamma", 1]:idx["skgamma", 2], 5]))
        pars[idx["skgamma", 1]:idx["skgamma", 2], 5] = xskew$skewpar.LB[3]
      if (is.na(pars[idx["skgamma", 1]:idx["skgamma", 2], 6]))
        pars[idx["skgamma", 1]:idx["skgamma", 2], 6] = xskew$skewpar.UB[3]
      if (is.null(start.pars$skgamma))
        pars[idx["skgamma", 1]:idx["skgamma", 2], 1] = xskew$skewpars[3] else pars[idx["skgamma", 1]:idx["skgamma", 2], 1] = start.pars$skgamma[1]
        if (any(substr(fixed.names, 1, 6) == "skgamma")) {
          pars[idx["skgamma", 1]:idx["skgamma", 2], 1] = as.numeric(fixed.pars$skgamma)
          pars[idx["skgamma", 1]:idx["skgamma", 2], 5] = as.numeric(fixed.pars$skgamma)
          pars[idx["skgamma", 1]:idx["skgamma", 2], 6] = as.numeric(fixed.pars$skgamma)
        }
    }
    if (modelinc[17] > 0) {
      if (is.na(pars[idx["skbeta", 1]:idx["skbeta", 2], 5]))
        pars[idx["skbeta", 1]:idx["skbeta", 2], 5] = xskew$skewpar.LB[4]
      if (is.na(pars[idx["skbeta", 1]:idx["skbeta", 2], 6]))
        pars[idx["skbeta", 1]:idx["skbeta", 2], 6] = xskew$skewpar.UB[4]
      if (is.null(start.pars$skbeta))
        pars[idx["skbeta", 1]:idx["skbeta", 2], 1] = xskew$skewpars[4] else pars[idx["skbeta", 1]:idx["skbeta", 2], 1] = start.pars$skbeta[1]
        if (any(substr(fixed.names, 1, 6) == "skbeta")) {
          pars[idx["skbeta", 1]:idx["skbeta", 2], 1] = as.numeric(fixed.pars$skbeta)
          pars[idx["skbeta", 1]:idx["skbeta", 2], 5] = as.numeric(fixed.pars$skbeta)
          pars[idx["skbeta", 1]:idx["skbeta", 2], 6] = as.numeric(fixed.pars$skbeta)
        }
    }
    if (modelinc[18] > 0) {
    if (is.na(pars[idx["volsk", 1]:idx["volsk", 2], 5]))
      pars[idx["volsk", 1]:idx["volsk", 2], 5] = -2
    if (is.na(pars[idx["volsk", 1]:idx["volsk", 2], 6]))
      pars[idx["volsk", 1]:idx["volsk", 2], 6] = 2
    if (is.null(start.pars$volsk))
      pars[idx["volsk", 1]:idx["volsk", 2], 1] = 0 else pars[idx["volsk", 1]:idx["volsk", 2], 1] = start.pars$volsk[1]
      if (any(substr(fixed.names, 1, 6) == "volsk")) {
        pars[idx["volsk", 1]:idx["volsk", 2], 1] = as.numeric(fixed.pars$volsk)
        pars[idx["volsk", 1]:idx["volsk", 2], 5] = as.numeric(fixed.pars$volsk)
        pars[idx["volsk", 1]:idx["volsk", 2], 6] = as.numeric(fixed.pars$volsk)
      }
  }
}
  if (modelinc[19] > 0) {
    uncshape1 = bounds$shape1
    xshape1 = .acdshape1bounds(modelinc[20:22], uncshape1, model$dmodel$model, dbounds[c(3,4,7)])
    arglist$skhEst[2] = xshape1$sh0
    pxd = which(is.na(pars[idx["sh1cons", 1]:idx["sh1cons", 2], 5]))
    if (length(pxd) > 0)
      pars[(idx["sh1cons", 1]:idx["sh1cons", 2])[pxd], 5] = xshape1$shapepar.LB[1]
    pxd = which(is.na(pars[idx["sh1cons", 1]:idx["sh1cons", 2], 6]))
    if (length(pxd) > 0)
      pars[(idx["sh1cons", 1]:idx["sh1cons", 2])[pxd], 6] = xshape1$shapepar.UB[1]
    if (is.null(start.pars$sh1cons))
      pars[idx["sh1cons", 1]:idx["sh1cons", 2], 1] = xshape1$shapepars[1] else pars[idx["sh1cons", 1]:idx["sh1cons", 2], 1] = start.pars$sh1cons[1]
    if (any(substr(fixed.names, 1, 6) == "sh1cons")) {
      pars[idx["sh1cons", 1]:idx["sh1cons", 2], 1] = as.numeric(fixed.pars$sh1cons)
      pars[idx["sh1cons", 1]:idx["sh1cons", 2], 5] = as.numeric(fixed.pars$sh1cons)
      pars[idx["sh1cons", 1]:idx["sh1cons", 2], 6] = as.numeric(fixed.pars$sh1cons)
    }

    if (modelinc[20] > 0) {
      if (is.na(pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 5]))
        pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 5] = xshape1$shapepar.LB[2]
      if (is.na(pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 6]))
        pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 6] = xshape1$shapepar.UB[2]
      if (is.null(start.pars$sh1alpha))
        pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 1] = xshape1$shapepars[2] else pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 1] = start.pars$sh1alpha[1]
        if (any(substr(fixed.names, 1, 6) == "sh1alpha")) {
          pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 1] = as.numeric(fixed.pars$sh1alpha)
          pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 5] = as.numeric(fixed.pars$sh1alpha)
          pars[idx["sh1alpha", 1]:idx["sh1alpha", 2], 6] = as.numeric(fixed.pars$sh1alpha)
        }
    }
    if (modelinc[21] > 0) {
      if (is.na(pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 5]))
        pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 5] = xshape1$shapepar.LB[3]
      if (is.na(pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 6]))
        pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 6] = xshape1$shapepar.UB[3]
      if (is.null(start.pars$sh1gamma))
        pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 1] = xshape1$shapepars[3] else pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 1] = start.pars$sh1gamma[1]
        if (any(substr(fixed.names, 1, 6) == "sh1gamma")) {
          pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 1] = as.numeric(fixed.pars$sh1gamma)
          pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 5] = as.numeric(fixed.pars$sh1gamma)
          pars[idx["sh1gamma", 1]:idx["sh1gamma", 2], 6] = as.numeric(fixed.pars$sh1gamma)
        }
    }
    if (modelinc[22] > 0) {
      if (is.na(pars[idx["sh1beta", 1]:idx["sh1beta", 2], 5]))
        pars[idx["sh1beta", 1]:idx["sh1beta", 2], 5] = xshape1$shapepar.LB[4]
      if (is.na(pars[idx["sh1beta", 1]:idx["sh1beta", 2], 6]))
        pars[idx["sh1beta", 1]:idx["sh1beta", 2], 6] = xshape1$shapepar.UB[4]
      if (is.null(start.pars$sh1beta))
        pars[idx["sh1beta", 1]:idx["sh1beta", 2], 1] = xshape1$shapepars[4] else pars[idx["sh1beta", 1]:idx["sh1beta", 2], 1] = start.pars$sh1beta[1]
        if (any(substr(fixed.names, 1, 6) == "sh1beta")) {
          pars[idx["sh1beta", 1]:idx["sh1beta", 2], 1] = as.numeric(fixed.pars$sh1beta)
          pars[idx["sh1beta", 1]:idx["sh1beta", 2], 5] = as.numeric(fixed.pars$sh1beta)
          pars[idx["sh1beta", 1]:idx["sh1beta", 2], 6] = as.numeric(fixed.pars$sh1beta)
        }
    }
    if (modelinc[23] > 0) {
  if (is.na(pars[idx["volsh1", 1]:idx["volsh1", 2], 5]))
    pars[idx["volsh1", 1]:idx["volsh1", 2], 5] = -2
  if (is.na(pars[idx["volsh1", 1]:idx["volsh1", 2], 6]))
    pars[idx["volsh1", 1]:idx["volsh1", 2], 6] = 2
  if (is.null(start.pars$volsh1))
    pars[idx["volsh1", 1]:idx["volsh1", 2], 1] = 0 else pars[idx["volsh1", 1]:idx["volsh1", 2], 1] = start.pars$volsh1[1]
    if (any(substr(fixed.names, 1, 6) == "volsh1")) {
      pars[idx["volsh1", 1]:idx["volsh1", 2], 1] = as.numeric(fixed.pars$volsh1)
      pars[idx["volsh1", 1]:idx["volsh1", 2], 5] = as.numeric(fixed.pars$volsh1)
      pars[idx["volsh1", 1]:idx["volsh1", 2], 6] = as.numeric(fixed.pars$volsh1)
    }
}
}
  if (modelinc[24] > 0) {
      uncshape2 = bounds$shape2
      xshape2 = .acdshape2bounds(modelinc[25:27], uncshape2, model$dmodel$model, dbounds[5:7])
      arglist$skhEst[3] = xshape2$sh0
      pxd = which(is.na(pars[idx["sh2cons", 1]:idx["sh2cons", 2], 5]))
      if (length(pxd) > 0)
        pars[(idx["sh2cons", 1]:idx["sh2cons", 2])[pxd], 5] = xshape2$shapepar.LB[1]
      pxd = which(is.na(pars[idx["sh2cons", 1]:idx["sh2cons", 2], 6]))
      if (length(pxd) > 0)
        pars[(idx["sh2cons", 1]:idx["sh2cons", 2])[pxd], 6] = xshape2$shapepar.UB[1]
      #pars[(idx["sh2cons", 1]:idx["sh2cons", 2])[pxd], 6] = -3

      if (is.null(start.pars$sh2cons))
        pars[idx["sh2cons", 1]:idx["sh2cons", 2], 1] = xshape2$shapepars[1] else pars[idx["sh2cons", 1]:idx["sh2cons", 2], 1] = start.pars$sh2cons[1]
      if (any(substr(fixed.names, 1, 6) == "sh2cons")) {
        pars[idx["sh2cons", 1]:idx["sh2cons", 2], 1] = as.numeric(fixed.pars$sh2cons)
        pars[idx["sh2cons", 1]:idx["sh2cons", 2], 5] = as.numeric(fixed.pars$sh2cons)
        pars[idx["sh2cons", 1]:idx["sh2cons", 2], 6] = as.numeric(fixed.pars$sh2cons)
      }

      if (modelinc[25] > 0) {
        if (is.na(pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 5]))
          pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 5] = xshape2$shapepar.LB[2]
        if (is.na(pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 6]))
          pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 6] = xshape2$shapepar.UB[2]
        if (is.null(start.pars$sh2alpha))
          pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 1] = xshape2$shapepars[2] else pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 1] = start.pars$sh2alpha[1]
          if (any(substr(fixed.names, 1, 6) == "sh2alpha")) {
            pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 1] = as.numeric(fixed.pars$sh2alpha)
            pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 5] = as.numeric(fixed.pars$sh2alpha)
            pars[idx["sh2alpha", 1]:idx["sh2alpha", 2], 6] = as.numeric(fixed.pars$sh2alpha)
          }
      }
      if (modelinc[26] > 0) {
        if (is.na(pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 5]))
          pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 5] = xshape2$shapepar.LB[3]
        if (is.na(pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 6]))
          pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 6] = xshape2$shapepar.UB[3]
        if (is.null(start.pars$sh2gamma))
          pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 1] = xshape2$shapepars[3] else pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 1] = start.pars$sh2gamma[1]
          if (any(substr(fixed.names, 1, 6) == "sh2gamma")) {
            pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 1] = as.numeric(fixed.pars$sh2gamma)
            pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 5] = as.numeric(fixed.pars$sh2gamma)
            pars[idx["sh2gamma", 1]:idx["sh2gamma", 2], 6] = as.numeric(fixed.pars$sh2gamma)
          }
      }
      if (modelinc[27] > 0) {
        if (is.na(pars[idx["sh2beta", 1]:idx["sh2beta", 2], 5]))
          pars[idx["sh2beta", 1]:idx["sh2beta", 2], 5] = xshape2$shapepar.LB[4]
        if (is.na(pars[idx["sh2beta", 1]:idx["sh2beta", 2], 6]))
          pars[idx["sh2beta", 1]:idx["sh2beta", 2], 6] = xshape2$shapepar.UB[4]
        if (is.null(start.pars$sh2beta))
          pars[idx["sh2beta", 1]:idx["sh2beta", 2], 1] = xshape2$shapepars[4] else pars[idx["sh2beta", 1]:idx["sh2beta", 2], 1] = start.pars$sh2beta[1]
          if (any(substr(fixed.names, 1, 6) == "sh2beta")) {
            pars[idx["sh2beta", 1]:idx["sh2beta", 2], 1] = as.numeric(fixed.pars$sh2beta)
            pars[idx["sh2beta", 1]:idx["sh2beta", 2], 5] = as.numeric(fixed.pars$sh2beta)
            pars[idx["sh2beta", 1]:idx["sh2beta", 2], 6] = as.numeric(fixed.pars$sh2beta)
          }
      }
      if (modelinc[28] > 0) {
        if (is.na(pars[idx["volsh2", 1]:idx["volsh2", 2], 5]))
          pars[idx["volsh2", 1]:idx["volsh2", 2], 5] = -2
        if (is.na(pars[idx["volsh2", 1]:idx["volsh2", 2], 6]))
          pars[idx["volsh2", 1]:idx["volsh2", 2], 6] = 2
        if (is.null(start.pars$volsh2))
          pars[idx["volsh2", 1]:idx["volsh2", 2], 1] = 0 else pars[idx["volsh2", 1]:idx["volsh2", 2], 1] = start.pars$volsh2[1]
          if (any(substr(fixed.names, 1, 6) == "volsh2")) {
            pars[idx["volsh2", 1]:idx["volsh2", 2], 1] = as.numeric(fixed.pars$volsh2)
            pars[idx["volsh2", 1]:idx["volsh2", 2], 5] = as.numeric(fixed.pars$volsh2)
            pars[idx["volsh2", 1]:idx["volsh2", 2], 6] = as.numeric(fixed.pars$volsh2)
          }
      }
  }
  ans = list(pars = pars, arglist = arglist)
  return(ans)
}
