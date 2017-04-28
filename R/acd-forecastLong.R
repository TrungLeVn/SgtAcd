#--------------------
# ACDforecast method
#' @export acdforecastLong
#' @title Autoregressive Conditional Density model forecast
#' @description Produce forecast for ACD models, basing on the simulation
#' @param fitORspec Either ACDfit or ACDspec object is allowed
#--------------------
# forecast
acdforecastLong = function(fitORspec, data = NULL, horizon = 10, n.roll = 0, out.sample = 0,
                       m.sim = 1000, cluster = NULL,
                       skew0 = NULL, shape10 = NULL,skew20 = NULL, ...)
{
  UseMethod("acdforecastLong")
}

.acdforecastFitL = function(fitORspec, horizon = 10, n.roll = 0,
                           m.sim = 1000, cluster = NULL, ...)
{
  fit = fitORspec
  ans = .acdforecastL(fit = fit, horizon = horizon,
                     n.roll = n.roll,
                     m.sim = m.sim,
                     cluster = cluster, ...)
  return(ans)
}
# switch if we supply model specification to the acdforecast function
.acdforecastSpecL = function(fitORspec, data = NULL, horizon = 10, n.roll = 0, out.sample = 0,
                            m.sim = 1000,
                            cluster = NULL, skew0 = NULL, shape0 = NULL, ...)
{
  spec = fitORspec
  ans = .acdforecast2L(spec = spec, data = data, horizon = horizon,
                      n.roll = n.roll, out.sample = out.sample,
                      m.sim = m.sim,
                      cluster = cluster, skew0 = skew0, shape10 = shape10,shape20 = shape20, ...)
  return(ans)
}
setMethod("acdforecastLong", signature(fitORspec = "ACDfit"), .acdforecastFitL)
setMethod("acdforecastLong", signature(fitORspec = "ACDspec"), .acdforecastSpecL)
#------------------
#---------------------------
.acdforecastL = function(fit, horizon = 10, n.roll = 0,
                        m.sim = 1000, cluster = NULL, rseed = NA, ...)
{
  if(is.na(rseed[1])){
    sseed = 2706:(2706+m.sim-1)
  } else{
    if(length(rseed) != m.sim) stop("\uacdforecast-->error: rseed must be of length m.sim!\n")
    sseed = rseed[1:m.sim]
  }
  data = fit@model$modeldata$data
  Nor = length(fit@model$modeldata$data)
  index = fit@model$modeldata$index
  period = fit@model$modeldata$period
  n.ahead = 1 #Only allow for one-steap ahead as the long-horizon conditional density is not additive
  ns = fit@model$n.start ##ns is the number of periods before the end of dataset that is used to fit the model
  if(n.roll>0 && ns ==0) stop("\uacdforecast-->error: Rolling forecast is only available when we have out-of-sample data!\n")
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
                                      archm = model$mmodel$archm, skm = model$mmodel$skm, pskm = model$mmodel$pskm, adjm = model$mmodel$adjm),
                    distribution.model = list(model = model$dmodel$model,
                                              skewOrder = model$dmodel$skewOrder, skewshock = model$dmodel$skewshock,
                                              skewmodel = model$dmodel$skewmodel, volsk = model$dmodel$volsk,
                                              shape1Order = model$dmodel$shape1Order, shape1shock = model$dmodel$shape1shock,
                                              shape1model = model$dmodel$shape1model, volsh1 = model$dmodel$volsh1,
                                              shape2Order = model$dmodel$shape2Order, shape2shock = model$dmodel$shape2shock,
                                              shape2model = model$dmodel$shape2model, volsh2 = model$dmodel$volsh2,
                                              exp.rate=model$sbounds[7]))
    fspec = setfixedacd(fspec, as.list(coefacd(fit)))
    fspec = setboundsacd(fspec,list(shape1 = fit@model$sbounds[3:4],shape2 = fit@model$sbounds[5:6], skew = fit@model$sbounds[1:2]))
    flt = acdfilter(spec = fspec, data = tmp, n.old = N, skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape[1],shape20 = fit@fit$tshape2[1])
    sigmafilter 	= flt@filter$sigma
    resfilter 		= flt@filter$residuals
    zfilter 		= flt@filter$z
    tskewfilter 	= flt@filter$tskew
    tshape1filter 	= flt@filter$tshape1
    tshape2filter   =flt@filter$tshape2
    tempskewfilter 	= flt@filter$tempskew
    tempshape1filter = flt@filter$tempshape1
    tempshape2filter = flt@filter$tempshape2
    seriesfor = sigmafor = matrix(NA, nrow = n.roll+1, ncol = n.ahead)
    skewnessFor = kurtosisFor = matrix(NA,nrow = n.roll + 1, ncol = n.ahead)
    # n.roll x n.ahead (n.ahead=1 generted by model)
    # n.ahead>1 by simulation
    rownames(seriesfor) = rownames(sigmafor) = rownames(skewnessFor) = rownames(kurtosisFor)= as.character(index[N:(N+n.roll)])
    colnames(seriesfor) = colnames(sigmafor) = colnames(skewnessFor) = colnames(kurtosisFor) = paste("T+", 1:n.ahead, sep="")
    mx = model$maxOrder #At the moment, the package only allows for maximum lag of 1 in both mean, sigma and higher moment equations.
    # Hence mx is always 1
    #constmean = unname(coefacd(fit)["mu"])
    for(i in 1:(n.roll+1)){
        spec = fspec
        presig     = tail(sigmafilter[1:(N+i-1)],  mx)
        preskew    = tail(tskewfilter[1:(N+i-1)],  mx)
        preshape1   = tail(tshape1filter[1:(N+i-1)], mx)
        preshape2  = tail(tshape2filter[1:(N+i-1)], mx)
        prereturns = tail(data[1:(N+i-1)],         mx)
          sim = acdpath(spec, n.sim = horizon, n.start = 0, m.sim = m.sim,
                        presigma = presig, preskew = preskew,
                        preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                        preresiduals = NA, rseed = rseed,
                        cluster = cluster)
          seriesfor[i, ] = mean(colSums(sim@path$seriesSim))
          sigmafor[i,  ] = sqrt(mean(colSums(sim@path$sigmaSim^2)))
          skewnessFor[i, ] = PerformanceAnalytics::skewness(colSums(sim@path$seriesSim))
          kurtosisFor[i, ] = PerformanceAnalytics::kurtosis(colSums(sim@path$seriesSim),method = "excess")
    }
  } else{
    seriesfor = sigmafor = matrix(NA, nrow = 1, ncol = n.ahead)
    skewnessFor = kurtosisFor = matrix(NA,nrow = 1, ncol = n.ahead)
    # n.roll x n.ahead (n.ahead=1 generted by model)
    # n.ahead>1 by simulation
    rownames(seriesfor) = rownames(sigmafor) = rownames(skewnessFor) = rownames(kurtosisFor)= as.character(index[N])
    colnames(seriesfor) = colnames(sigmafor) = colnames(skewnessFor) = colnames(kurtosisFor) = paste("T+", 1:n.ahead, sep="")
    mx = model$maxOrder
    presig = as.numeric(tail(sigmaAcd(fit),mx))
    preskew = as.numeric(tail(skew(fit),mx))
    preshape1 = as.numeric(tail(shape1(fit),mx))
    preshape2 = as.numeric(tail(shape2(fit),mx))
    preresiduals = as.numeric(tail(residualsAcd(fit),mx))
    prereturns = as.numeric(tail(fit@model$modeldata$data[N],mx))
    sim = acdsim(fit, n.sim = horizon, m.sim = m.sim,
                 presigma = presig, preskew = preskew,
                 preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                 preresiduals = preresiduals, rseed = rseed,
                 cluster = cluster)
    seriesfor[,1] = mean(colSums(sim@simulation$seriesSim))
    sigmafor[,1] = sqrt(mean(colSums(sim@simulation$sigmaSim^2)))
    skewnessFor[,1] = PerformanceAnalytics::skewness(colSums(sim@simulation$seriesSim))
    }
  fcst = list()
  fcst$n.ahead = n.ahead
  fcst$n.roll = n.roll
  fcst$N = N+ns
  fcst$n.start = ns
  fcst$seriesFor = seriesfor
  fcst$sigmaFor  = sigmafor
  fcst$skewnessFor  = skewnessFor
  fcst$kurtosisFor = kurtosisFor
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
.acdforecast2L = function(spec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
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

  seriesfor = sigmafor = matrix(NA, nrow = n.roll+1, ncol = n.ahead)
  skewnessFor = kurtosisFor = matrix(NA,nrow = n.roll + 1, ncol = n.ahead)
  # n.roll x n.ahead (n.ahead=1 generted by model)
  # n.ahead>1 by simulation
  rownames(seriesfor) = rownames(sigmafor) = rownames(skewnessFor) = rownames(kurtosisFor)= as.character(index[N:(N+n.roll)])
  colnames(seriesfor) = colnames(sigmafor) = colnames(skewnessFor) = colnames(kurtosisFor) = paste("T+", 1:n.ahead, sep="")
  mx = model$maxOrder #At the moment, the package only allows for maximum lag of 1 in both mean, sigma and higher moment equations.
  # Hence mx is always 1
  #constmean = unname(coefacd(fit)["mu"])
  for(i in 1:(n.roll+1)){
    spec = fspec
    presig     = tail(sigmafilter[1:(N+i-1)],  mx)
    preskew    = tail(tskewfilter[1:(N+i-1)],  mx)
    preshape1   = tail(tshape1filter[1:(N+i-1)], mx)
    preshape2  = tail(tshape2filter[1:(N+i-1)], mx)
    prereturns = tail(data[1:(N+i-1)],         mx)
    sim = acdpath(spec, n.sim = horizon, n.start = 0, m.sim = m.sim,
                  presigma = presig, preskew = preskew,
                  preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                  preresiduals = NA, rseed = rseed,
                  cluster = cluster)
    seriesfor[i, ] = mean(colSums(sim@path$seriesSim))
    sigmafor[i,  ] = sqrt(mean(colSums(sim@path$sigmaSim^2)))
    skewnessFor[i, ] = PerformanceAnalytics::skewness(colSums(sim@path$seriesSim))
    kurtosisFor[i, ] = PerformanceAnalytics::kurtosis(colSums(sim@path$seriesSim),method = "excess")
  }
  fcst = list()
  fcst$n.ahead = n.ahead
  fcst$n.roll = n.roll
  fcst$N = N+ns
  fcst$n.start = ns
  fcst$seriesFor = seriesfor
  fcst$sigmaFor  = sigmafor
  fcst$skewnessFor  = skewnessfor
  fcst$kurtosisFor = kurtosisFor
  model$modeldata$sigma = flt@filter$sigma
  model$modeldata$residuals = flt@filter$residuals
  model$modeldata$tskew  = flt@filter$tskew
  model$modeldata$tshape = flt@filter$tshape
  ans = new("ACDforecast",
            forecast = fcst,
            model = model)
  return(ans)
}
