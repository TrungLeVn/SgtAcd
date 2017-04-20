#--------------------
# ACDforecast method
#' @export acdforecast
#' @title Autoregressive Conditional Density model forecast
#' @description Produce forecast for ACD models, basing on the simulation
#' @param fitORspec Either ACDfit or ACDspec object is allowed
#--------------------
# forecast
acdforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
                       m.sim = 1000, cluster = NULL,
                       skew0 = NULL, shape10 = NULL,skew20 = NULL, ...)
{
  UseMethod("acdforecast")
}

.acdforecastFit = function(fitORspec, n.ahead = 10, n.roll = 0,
                               m.sim = 1000, cluster = NULL, ...)
{
  fit = fitORspec
  ans = .acdforecast(fit = fit, n.ahead = n.ahead,
                     n.roll = n.roll,
                     m.sim = m.sim,
                     cluster = cluster, ...)
 return(ans)
}
# switch if we supply model specification to the acdforecast function
.acdforecastSpec = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
                               m.sim = 1000,
                               cluster = NULL, skew0 = NULL, shape0 = NULL, ...)
{
  spec = fitORspec
  ans = .acdforecast2(spec = spec, data = data, n.ahead = n.ahead,
                      n.roll = n.roll, out.sample = out.sample,
                      m.sim = m.sim,
                      cluster = cluster, skew0 = skew0, shape10 = shape10,shape20 = shape20, ...)
return(ans)
}
setMethod("acdforecast", signature(fitORspec = "ACDfit"), .acdforecastFit)
setMethod("acdforecast", signature(fitORspec = "ACDspec"), .acdforecastSpec)
#------------------
#---------------------------
.acdforecast = function(fit, n.ahead = 10, n.roll = 0,
                        m.sim = 1000, cluster = NULL, rseed = NA, ...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdsim-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  data = fit@model$modeldata$data
  Nor = length(as.numeric(data))
  index = fit@model$modeldata$index
  period = fit@model$modeldata$period
  ns = fit@model$n.start ##ns is the number of periods before the end of dataset that is used to fit the model
  if(n.roll>0 && ns ==0) stop("\uacdsim-->error: Rolling forecast is only available when we have out-of-sample data!\n")
  N = Nor - ns #N is the length of data that is used to fit model or data[N] is the last point of fitted data
  model = fit@model
  ipars = fit@fit$ipars
  modelinc = model$modelinc
  idx = model$pidx
  if( n.roll > ns ) stop("\nacdforecast-->error: n.roll must not be greater than out.sample!")
  pars = fit@fit$coef
  ipars = fit@fit$ipars
  # filter data (check external regressor data - must equal length of origData)
  # Generate the 1 extra ahead forecast using filter process with the last residuals and z as based values
  if(n.roll > 0){ # if we have out-of-sample in the fit object and we wish to rolling forecast
  fcreq = ifelse(ns >= (n.ahead+n.roll+1), n.ahead+n.roll + 1, ns)
  tmp =  xts::xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
  fspec = acdspec(variance.model = list(model = model$vmodel$model,
                                        garchOrder = model$vmodel$garchOrder,
                                        variance.targeting = FALSE),
                  mean.model = list(armaOrder = model$mmodel$armaOrder,
                                    include.mean = model$mmodel$include.mean,
                                    archm = model$mmodel$archm, skm = model$mmodel$skm, pskm = model$mmodel$pskm),
                  distribution.model = list(model = model$dmodel$model,
                                            skewOrder = model$dmodel$skewOrder, skewshock = model$dmodel$skewshock,
                                            skewmodel = model$dmodel$skewmodel, volsk = model$dmodel$volsk,
                                            shape1Order = model$dmodel$shape1Order, shape1shock = model$dmodel$shape1shock,
                                            shape1model = model$dmodel$shape1model, volsh1 = model$dmodel$volsh1,
                                            shape2Order = model$dmodel$shape2Order, shape2shock = model$dmodel$shape2shock,
                                            shape2model = model$dmodel$shape2model, volsh2 = model$dmodel$volsh2,
                                            exp.rate=model$sbounds[7]))
  fspec = setfixedacd(fspec, as.list(coefacd(fit)))
  fspec = setboundsacd(fspec,list(shape = fit@model$sbounds[3:4], skew = fit@model$sbounds[1:2]))
  flt = acdfilter(spec = fspec, data = tmp, n.old = N, skew0 = fit@fit$tskew[1], shape0 = fit@fit$tshape[1])
  sigmafilter 	= flt@filter$sigma
  resfilter 		= flt@filter$residuals
  zfilter 		= flt@filter$z
  tskewfilter 	= flt@filter$tskew
  tshape1filter 	= flt@filter$tshape1
  tshape2filter   =flt@filter$tshape2
  tempskewfilter 	= flt@filter$tempskew
  tempshape1filter = flt@filter$tempshape1
  tempshape2filter = flt@filter$tempshape2
  pskewnessfilter = flt@filter$Pskewness
  seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
  tskewfor = tshape1for = tshape2for = pskewnessfor= matrix(NA, ncol = n.roll+1, nrow = n.ahead)
  seriesfor[1,] = fittedAcd(flt)[(N+1):(N+n.roll+1)]
  sigmafor[1,]  =  sigmaAcd(flt)[(N+1):(N+n.roll+1)]
  tskewfor[1,]  =  skew(flt)[(N+1):(N+n.roll+1)]
  tshape1for[1,] = shape1(flt)[(N+1):(N+n.roll+1)]
  tshape2for[1,] = shape2(flt)[(N+1):(N+n.roll+1)]
  pskewnessfor[1,] = Pskewness(flt)[(N+1):(N+n.roll+1)]
  # n.roll x n.ahead (n.ahead=1 generted by model)
  # n.ahead>1 by simulation
  colnames(seriesfor) = colnames(sigmafor) = as.character(index[N:(N+n.roll)])
  colnames(tskewfor) = colnames(tshape1for) = colnames(tshape2for) = colnames(pskewnessfor) = as.character(index[N:(N+n.roll)])
  rownames(seriesfor) = rownames(sigmafor) = paste("T+", 1:n.ahead, sep="")
  rownames(tskewfor) = rownames(tshape2for) = rownames(tshape2for) = paste("T+", 1:n.ahead, sep="")
  mx = model$maxOrder #At the moment, the package only allows for maximum lag of 1 in both mean, sigma and higher moment equations.
  # Hence mx is always 1
  #constmean = unname(coefacd(fit)["mu"])
  for(i in 1:(n.roll+1)){
    if(n.ahead>1){
      spec = fspec
      presig     = tail(sigmafilter[1:(N+i-1)],  mx)
      preskew    = tail(tskewfilter[1:(N+i-1)],  mx)
      preshape1   = tail(tshape1filter[1:(N+i-1)], mx)
      preshape2  = tail(tshape2filter[1:(N+i-1)], mx)
      prereturns = tail(data[1:(N+i-1)],         mx)
      prePskewness = tail(pskewnessfilter[1:(N+i-1)], mx)
      for(j in 2:n.ahead){
        sim = acdpath(spec, n.sim = 1, n.start = 0, m.sim = m.sim,
                      presigma = presig, preskew = preskew,
                      preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                      preresiduals = NA, rseed = rseed,
                      cluster = cluster)
        seriesfor[j, i] = apply(sim@path$seriesSim, 1, "mean")
        sigmafor[j,  i] = sqrt(apply(sim@path$sigmaSim^2, 1, "mean"))
        tskewfor[j,  i]	= apply(sim@path$skewSim, 1, "mean")
        tshape1for[j, i]	= apply(sim@path$shape1Sim,  1, "mean")
        tshape2for[j, i] = apply(sim@path$shape2Sim, 1, "mean")
        pskewnessfor[j,i] = apply(sim@path$pskewSim,1,"mean")
        # update the previous values:
        prereturns[mx] = seriesfor[j, i]
        presig[mx]   = sigmafor[j, i]
        preskew[mx]  = tskewfor[j, i]
        preshape1[mx] = tshape1for[j, i]
        preshape2[mx] = tshape2for[j, i]
        prePskewness[mx] = pskewnessfor[j,i]
      }
    }
  }
  } else{
    seriesfor = sigmafor = matrix(NA, ncol = 1, nrow = n.ahead)
    tskewfor = tshape1for = tshape2for = tempskewfor = pskewnessfor= matrix(NA, ncol = 1, nrow = n.ahead)
    colnames(seriesfor) = colnames(sigmafor) = as.character(index[N])
    colnames(tskewfor) = colnames(tshape1for) = colnames(tshape2for) = colnames(pskewnessfor) = as.character(index[N])
    rownames(seriesfor) = rownames(sigmafor) = paste("T+", 1:n.ahead, sep="")
    rownames(tskewfor) = rownames(tshape2for) = rownames(tshape2for) = paste("T+", 1:n.ahead, sep="")
    mx = model$maxOrder
    presig = as.numeric(tail(sigmaAcd(fit),mx))
    preskew = as.numeric(tail(skew(fit),mx))
    preshape1 = as.numeric(tail(shape1(fit),mx))
    preshape2 = as.numeric(tail(shape2(fit),mx))
    preresiduals = as.numeric(tail(residualsAcd(fit),mx))
    prereturns = as.numeric(tail(fit@model$modeldata$data[N],mx))
    prePskewness = as.numeric(tail(Pskewness(fit),mx))
    sim = acdsim(fit, n.sim = n.ahead, m.sim = m.sim,
                    presigma = presig, preskew = preskew,
                    preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                    preresiduals = preresiduals, rseed = rseed,
                    cluster = cluster)
    for(j in 1 : n.ahead){
      seriesfor[j,1] = mean(sim@simulation$seriesSim[j,])
      sigmafor[j,1] = sqrt(mean(sim@simulation$sigmaSim[j,]^2))
      tskewfor[j,1]	= mean(sim@simulation$skewSim[j,])
      tshape1for[j,1]	= mean(sim@simulation$shape1Sim[j,])
      tshape2for[j,1] = mean(sim@simulation$shape2Sim[j,])
      pskewnessfor[j,1] = mean(sim@simulation$pskewSim[j,])
    }
    }
  fcst = list()
  fcst$n.ahead = n.ahead
  fcst$n.roll = n.roll
  fcst$N = N+ns
  fcst$n.start = ns
  fcst$seriesFor = seriesfor
  fcst$sigmaFor  = sigmafor
  fcst$tskewFor  = tskewfor
  fcst$tshape1For = tshape1for
  fcst$tshape2For = tshape2for
  fcst$pskewnessfor = pskewnessfor
  if(n.roll>0){
    model$modeldata$sigma = flt@filter$sigma
    model$modeldata$residuals = flt@filter$residuals
    model$modeldata$tskew  = flt@filter$tskew
    model$modeldata$tshape1 = flt@filter$tshape1
    model$modeldata$tshape2 = flt@filter$tshape2
    model$modeldata$Pskewness = flt@filter$Pskewness
  } else{
    model = fit@model
  }
  ans = new("ACDforecast",
            forecast = fcst,
            model = model)
  return(ans)
}

# CONTINUE HERE:
.acdforecast2 = function(spec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
                         m.sim = 1000, cluster = NULL, skew0 = NULL,
                         shape0 = NULL, ...)
{
  vmodel = spec@model$vmodel$model
  xdata = rugarch:::.extractdata(data)
  data = xdata$data
  index = xdata$index
  period = xdata$period
  Nor = length(as.numeric(data))
  ns = out.sample
  N = Nor - ns
  model = spec@model
  modelinc = model$modelinc
  idx = model$pidx
  ipars = model$pars
  pars = unlist(model$fixed.pars)
  parnames = names(pars)
  modelnames = rugarch:::.checkallfixed(spec)
  if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
    cat("\nacdfilter-->error: parameters names do not match specification\n")
    cat("Expected Parameters are: ")
    cat(paste(modelnames))
    cat("\n")
    stop("Exiting", call. = FALSE)
  }
  fspec = spec
  setfixed(fspec)<-as.list(pars)
  if( n.roll > ns ) stop("\nacdforecast-->error: n.roll must not be greater than out.sample!")
  # filter data (check external regressor data - must equal length of origData)
  fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
  # Generate the 1 extra ahead forecast
  if((n.ahead+n.roll)<=ns){
    tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
  } else{
    tmp =  xts(c(data[1:(N + fcreq)],0), c(index[1:(N + fcreq)], index[(N + fcreq)]+1))
  }
  flt = acdfilter(spec = fspec, data = tmp, n.old = N, skew0 = skew0, shape0 = shape0)
  sigmafilter 	= flt@filter$sigma
  resfilter 		= flt@filter$residuals
  zfilter 		= flt@filter$z
  tskewfilter 	= flt@filter$tskew
  tshapefilter 	= flt@filter$tshape
  tempskewfilter 	= flt@filter$tempskew
  tempshapefilter = flt@filter$tempshape

  seriesfor = sigmafor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
  tskewfor = tshapefor = tempshafor = tempskewfor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
  seriesfor[1,] = fitted(flt)[(N+1):(N+n.roll+1)]
  sigmafor[1,]  =  sigma(flt)[(N+1):(N+n.roll+1)]
  tskewfor[1,]  =  skew(flt)[(N+1):(N+n.roll+1)]
  tshapefor[1,] = shape(flt)[(N+1):(N+n.roll+1)]
  # n.roll x n.ahead (n.ahead=1 generted by model)
  # n.ahead>1 by simulation
  colnames(seriesfor) = colnames(sigmafor) = as.character(index[N:(N+n.roll)])
  colnames(tskewfor) = colnames(tshapefor) = as.character(index[N:(N+n.roll)])
  rownames(seriesfor) = rownames(sigmafor) = paste("T+", 1:n.ahead, sep="")
  rownames(tskewfor) = rownames(tshapefor) = paste("T+", 1:n.ahead, sep="")
  mx = model$maxOrder
  for(i in 1:(n.roll+1)){
    if(n.ahead>1){
      spec = fspec
      presig     = tail(sigmafilter[1:(N+1+i-1)],  mx)
      preskew    = tail(tskewfilter[1:(N+1+i-1)],  mx)
      preshape   = tail(tshapefilter[1:(N+1+i-1)], mx)
      prereturns = tail(data[1:(N+1+i-1)],         mx)
      for(j in 2:n.ahead){
        sim = acdpath(fspec, n.sim = 1, n.start = 0, m.sim = m.sim,
                      presigma = presig, preskew = preskew,
                      preshape = preshape, prereturns = prereturns,
                      preresiduals = NA, rseed = NA,
                      mexsimdata = mxsim, vexsimdata = vxsim,
                      skxsimdata = skxsim, shxsimdata = shxsim,
                      cluster = cluster)
        seriesfor[j, i] = apply(sim@path$seriesSim, 1, "mean")
        sigmafor[j,  i] = sqrt(apply(sim@path$sigmaSim^2, 1, "mean"))
        tskewfor[j,  i]	= apply(sim@path$skewSim, 1, "mean")
        tshapefor[j, i]	= apply(sim@path$shapeSim,  1, "mean")
        # update the previous values:
        prereturns[mx] = seriesfor[j, i]
        presig[mx]   = sigmafor[j, i]
        preskew[mx]  = tskewfor[j, i]
        preshape[mx] = tshapefor[j, i]
      }
    }
  }

  fcst = list()
  fcst$n.ahead = n.ahead
  fcst$n.roll = n.roll
  fcst$N = N+ns
  fcst$n.start = ns
  fcst$seriesFor = seriesfor
  fcst$sigmaFor  = sigmafor
  fcst$tskewFor  = tskewfor
  fcst$tshapeFor = tshapefor
  model$modeldata$sigma = flt@filter$sigma
  model$modeldata$residuals = flt@filter$residuals
  model$modeldata$tskew  = flt@filter$tskew
  model$modeldata$tshape = flt@filter$tshape
  ans = new("ACDforecast",
            forecast = fcst,
            model = model)
  return(ans)
}
