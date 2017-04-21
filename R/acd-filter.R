#---------------
# ACDfilter method
#' @title Filtering new data
#' @description Filtering new data with regards to the model specification with fixed parameters
#' @usage acdfilter = function(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, shape0 = NULL, ...)
#' @param spec An ACDspec object with fixed parameters
#' @param data dataset
#' @param n.old Number of old data that will be used to fit the model specification. The last value will be
#'          used to start the filtering
#' @return An ACDfiltering object
#' @export acdfilter
#---------------
acdfilter = function(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, shape10 = NULL,shape20 = NULL, ...)
{
  UseMethod("acdfilter")
}

.acdfilterswitch = function(spec, data, out.sample = 0,  n.old = NULL, skew0 = NULL, shape10 = NULL,shape20 = NULL, ...)
{
  switch(spec@model$vmodel$model,
         sGARCH = .acdfilter(spec = spec, data = data, out.sample = out.sample,
                             n.old = n.old, skew0 = skew0, shape0 = shape0, ...),
         gjrGARCH = .acdfilter(spec = spec, data = data, out.sample = out.sample,
                               n.old = n.old, skew0 = skew0, shape0 = shape0, ...))
}
setMethod("acdfilter", signature(spec = "ACDspec"), .acdfilterswitch)
#--------------------
.acdfilter = function(spec, data, out.sample = 0, n.old = NULL, skew0 = NULL, shape10 = NULL,shape20 = NULL, ...)
{
  # n.old is optional and indicates the length of the original dataseries (in
  # cases when this represents a dataseries augmented by newer data). The reason
  # for using this is so that the old and new datasets agree since the original
  # recursion uses the sum of the residuals to start the recursion and therefore
  # is influenced by new data. For a small augmentation the values converge after
  # x periods, but it is sometimes preferable to have this option so that there is
  # no forward looking information contaminating the study.
  # TRUNG: The spec object in the function is setfixed with the fitted object, using old data.
  tic = Sys.time()
  vmodel = spec@model$vmodel$model
  xdata = .extractdata(data)
  data = xdata$data
  index = xdata$index
  period = xdata$period
  origdata = data
  origindex = xdata$index
  T = length(origdata)  - out.sample
  # T is the length of oldata or  in other words, the length of original data that is used to compute the estimated parameters.

  if(!is.null(n.old) && n.old>T) stop("\nn.old cannot be greater than length data - out.sample!")

  data = origdata[1:T]
  index = origindex[1:T]
  model = spec@model
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
  # once more into the spec
  # NB Any changes made to the spec are not preserved once we apply set fixed
  #setfixed(spec)<-as.list(pars)
  spec <- setfixedacd(spec,as.list(pars))
  garchenv = new.env(hash = TRUE)
  #racd_llh = 1
  #assign("racd_llh", 1, envir = garchenv)
  arglist = list()
  arglist$n.old = n.old
  arglist$transform = FALSE
  arglist$garchenv <- garchenv
  arglist$pmode = 0
  modelinc = model$modelinc
  pidx = model$pidx
  arglist$index = index
  arglist$trace = 0
  # store length of data for easy retrieval
  arglist$data = data
  arglist$ipars  = ipars
  # starting parameter for skew+shape1+shape2
  arglist$skhEst = rep(0,3)
  if(!is.null(shape10)){
    if(modelinc[19]==1) arglist$skhEst[2] = exptransform1(shape10, model$sbounds[3], rate = model$sbounds[7], inverse=TRUE)
  } else{
    if(modelinc[19]==1) arglist$skhEst[2] = pars["sh1cons"] #modelinc[27] is shcons
  }
  if(!is.null(shape20)){
    if(modelinc[24]==1) arglist$skhEst[3] = logtransform(shape20, model$sbounds[5], model$sbounds[6], inverse=TRUE)
  } else{
    if(modelinc[24]==1) arglist$skhEst[3] = pars["sh2cons"] #modelinc[27] is shcons
  }

  if(!is.null(skew0)){
    if(modelinc[14]==1) arglist$skhEst[1] = logtransform(skew0, model$sbounds[1], model$sbounds[2], inverse=TRUE)
  } else{
    if(modelinc[14]==1) arglist$skhEst[1] = pars["skcons"] #modelinc[21] is skcons
  }

  arglist$tmph  = 0
  arglist$model = model
  # we now split out any fixed parameters
  estidx = as.logical( ipars[,3] )
  arglist$estidx = estidx
  arglist$returnType = "all"
  arglist$fit.control=list(stationarity = 0)
  ans  = switch(vmodel,
                sGARCH   = .sacdLLH(pars, arglist),
                gjrGARCH  = .gjracdLLH(pars, arglist))

  filter = list()
  filter$z = ans$z
  filter$sigma = sqrt(ans$h)
  filter$residuals = ans$res
  filter$fitted.values = data - ans$res
  filter$tskew = ans$tskew
  filter$tshape1 = ans$tshape1
  filter$tshape2 = ans$tshape2
  filter$tempskew = ans$tempskew
  filter$tempshape1 = ans$tempshape1
  filter$tempshape2 = ans$tempshape2
  filter$Pskewness = ans$Pskewness
  filter$skewness = ans$skewness
  filter$kurtosis = ans$kurtosis
  filter$LLH = -ans$llh
  filter$log.likelihoods = ans$LHT
  filter$ipars = ipars
  filter$skhEst = arglist$skhEst
  model$modeldata$data = origdata
  model$modeldata$index = origindex
  model$modeldata$period = period
  model$modeldata$T = T
  model$n.start = out.sample
  filter$timer = Sys.time() - tic
  sol = new("ACDfilter",
            filter = filter,
            model = model)
  return(sol)
}
