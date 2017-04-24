#' @title Get model specification
#' @exportMethod getspecacd
getspecacd = function(object){
  UseMethod("getspecacd")
}
.getspecacd = function(object)
{
  vmodel = object@model$vmodel
  mmodel = object@model$mmodel
  dmodel = object@model$dmodel
  spec = acdspec( variance.model = vmodel, mean.model = mmodel,
                  distribution.model = dmodel, start.pars = object@model$start.pars,
                  fixed.pars = object@model$fixed.pars)
  return(spec)
}

setMethod(f = "getspecacd", signature(object = "ACDfit"), definition = .getspecacd)
#----------------------------------------------------------------------------------
# fixed parameters
#' @title Set fixed parameters
#' @description Set fixed parameters to a model specification
#' @param object Model sepecification
#' @param value A named list of parameters to be fixed
#' @exportMethod setfixedacd
setfixedacd = function(object, value){
  UseMethod("setfixedacd")
}
# get parameter values

.setfixedacd <- function(object,value){
  oldfixed = unlist(object@model$fixed.pars)
  model = object@model
  ipars = model$pars
  pars = unlist(value)
  names(pars) = parnames = tolower(names(pars))
  # included parameters in model
  modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
  inc = NULL
  for(i in seq_along(parnames)){
    if(is.na(match(parnames[i], modelnames))){
      warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
    } else{
      inc = c(inc, i)
    }
  }
  fixed.pars = pars[inc]
  names(fixed.pars) = tolower(names(pars[inc]))
  # need to check for duplicates
  xidx = match(names(fixed.pars), names(oldfixed))
  if(any(!is.na(xidx))){
    oldfixed=oldfixed[-na.omit(xidx)]
  }
  if(!is.null(oldfixed) && length(oldfixed)>0) fixed.pars = c(fixed.pars, oldfixed)
  # set parameter values
  vmodel = object@model$vmodel
  mmodel = object@model$mmodel
  dmodel = object@model$dmodel
  tmp = acdspec( variance.model = vmodel, mean.model = mmodel,
                 distribution.model = dmodel, start.pars  = model$start.pars,
                 fixed.pars = as.list(fixed.pars))
  tmp@model$pars[tmp@model$pars[,2]==0,5:6] = object@model$pars[tmp@model$pars[,2]==0,5:6]
  tmp <- setboundsacd(tmp,list(shape=object@model$sbounds[3:4], skew = object@model$sbounds[1:2]))
  return(tmp)
}
setMethod(f="setfixedacd", signature= c(object = "ACDspec", value = "vector"), definition = .setfixedacd)
#----------------------------------------------------------------------------------
# starting parameters
#' @title Set starting values for a model specification
#' @description Set the starting values for estimation routine of a model specification objects
#' @param object Model specification class
#' @param Value A named list of starting values for parameters, which are then put to the estimation procedures
#' @exportMethod setstartacd
setstartacd <- function(object,value){
  UseMethod("setstartacd")
}
.setstartacd = function(object, value){
  # get parameter values
  model = object@model
  ipars = model$pars
  pars = unlist(value)
  names(pars) = parnames = tolower(names(pars))
  # included parameters in model
  modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
  inc = NULL
  for(i in seq_along(parnames)){
    if(is.na(match(parnames[i], modelnames))){
      warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
    } else{
      inc = c(inc, i)
    }
  }
  start.pars = pars[inc]
  names(start.pars) = tolower(names(pars[inc]))
  # set parameter values
  # set parameter values
  vmodel = object@model$vmodel
  mmodel = object@model$mmodel
  dmodel = object@model$dmodel

  tmp = acdspec( variance.model = vmodel, mean.model = mmodel,
                 distribution.model = dmodel, fixed.pars  = model$fixed.pars,
                 start.pars = as.list(start.pars))
  tmp@model$pars[tmp@model$pars[,2]==0,5:6] = object@model$pars[tmp@model$pars[,2]==0,5:6]
  tmp = setboundsacd(tmp,object = list(skew=object@model$sbounds[1:2], shape1 = object@model$sbounds[3:4],shape2 = object@model$sbounds[5,6]))
  return(tmp)
}

setMethod(f="setstartacd", signature= c(object = "ACDspec", value = "vector"), definition = .setstartacd)
#----------------------------------------------------------------------------------
.checkallfixed = function( spec ){
  # check that a given spec with fixed parameters
  model = spec@model
  pars = model$pars
  pnames = rownames(pars)
  estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
  return( estpars )
}
#----------------------------------------------------------------------------------
# set parameters bounds
# Set the lower and upper bounds
# value is a list with names parameters taking 2 values (lower and upper)
# e.g. value  = list(alpha1 = c(0, 0.1), beta1 = c(0.9, 0.99))
# special attention given to skew/shape in case of dynamics
#' @title Set bounds for conditional density distribution
#' @description Set the bounds to the shape parameters of conditional distribution to which the density is exist
#' @param object The model specification object
#' @exportMethod setboundsacd
setboundsacd <- function(object,value){
  UseMethod("setboundsacd")
}
.acdsetbounds = function(object, value){
  model = object@model
  ipars = model$pars
  parnames = tolower(names(value))
  if(model$modelinc[14]>0 && any(parnames=="skew")){
    idx = which(parnames=="skew")
    if(!is.na(value[[idx]][1])) object@model$sbounds[1] = value[[idx]][1]
    if(!is.na(value[[idx]][2])) object@model$sbounds[2] = value[[idx]][2]
  }
  if(model$modelinc[19]>0 && any(parnames=="shape1")){
    idx = which(parnames=="shape1")
    if(!is.na(value[[idx]][1])) object@model$sbounds[3] = value[[idx]][1]
    if(!is.na(value[[idx]][2])) object@model$sbounds[4] = value[[idx]][2]
  }
  if(model$modelinc[24]>0 && any(parnames=="shape2")){
    idx = which(parnames=="shape2")
    if(!is.na(value[[idx]][1])) object@model$sbounds[5] = value[[idx]][1]
    if(!is.na(value[[idx]][2])) object@model$sbounds[6] = value[[idx]][2]
  }
  # included parameters in model
  modelnames = rownames(ipars[which(ipars[,4] == 1), ])
  sp = na.omit(match(parnames, modelnames))
  if(length(sp)>0){
    for(i in 1:length(sp)){
      #if(length(value[[modelnames[sp[i]]]])!=2)
      ipars[modelnames[sp[i]], 5] = as.numeric(value[[modelnames[sp[i]]]][1])
      ipars[modelnames[sp[i]], 6] = as.numeric(value[[modelnames[sp[i]]]][2])
    }
  }
  object@model$pars = ipars
  return(object)
}
setMethod(f="setboundsacd", signature= c(object = "ACDspec", value = "vector"), definition = .acdsetbounds)
#----------------------------------------------------------------
# Convergence check method
#-----------------------------------------------------------------
# convergence method
#' @exportMethod acdconvergence
acdconvergence = function(object){
  UseMethod("acdconvergence")
}
.acdconvergence = function(object)
{
  return( object@fit$convergence )
}
setMethod("acdconvergence", signature(object = "ACDfit"),  .acdconvergence)

.acdconvergenceroll = function(object){
  nonc = object@model$noncidx
  if(is.null(nonc)){
    ans = 0
  } else{
    ans = 1
    attr(ans, 'nonconverged')<-nonc
  }
  return(ans)
}

setMethod("acdconvergence", signature(object = "ACDroll"),  definition = .acdconvergenceroll)
#-----------------------------------------------------------------
# conditional mean (fitted method)
#-----------------------------------------------------------------
fittedAcd <- function(object)
{
  UseMethod("fittedAcd")
}
.acdfitted = function(object)
{
  if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
    D = object@model$modeldata$index[1:object@model$modeldata$T]
  }
  ans = switch(class(object)[1],
               ACDfit = xts(object@fit$fitted.values, D),
               ACDfilter = xts(object@filter$fitted.values, D),
               ACDforecast = object@forecast$seriesFor,
               ACDsim = {
                 ans = object@simulation$seriesSim
                 rownames(ans) = paste("T+",1:NROW(object@simulation$seriesSim), sep="")
                 return(ans)
               },
               ACDpath ={
                 ans = object@path$seriesSim
                 rownames(ans) = paste("T+",1:NROW(object@path$seriesSim), sep="")
                 return(ans)
               },
               ACDroll = as.xts(object@forecast$density[,"Mu",drop=FALSE]))
  return(ans)
}
setMethod("fittedAcd", signature(object = "ACDfit"), .acdfitted)
setMethod("fittedAcd", signature(object = "ACDfilter"), .acdfitted)
setMethod("fittedAcd", signature(object = "ACDforecast"), .acdfitted)
setMethod("fittedAcd", signature(object = "ACDpath"), .acdfitted)
setMethod("fittedAcd", signature(object = "ACDsim"), .acdfitted)
setMethod("fittedAcd", signature(object = "ACDroll"), .acdfitted)
#' @exportMethod fittedAcd

#------
#-------------------------------------------------------------------------------
# residuals method
residualsAcd <- function(object, standardize = FALSE)
{
  UseMethod("residualsAcd")
}
.acdfitresids = function(object, standardize = FALSE)
{
  if(is(object, "ACDfit")){
    if(standardize) ans = object@fit$residuals/object@fit$sigma else  ans = object@fit$residuals
  } else{
    if(standardize) ans = object@filter$residuals/object@filter$sigma else  ans = object@filter$residuals
  }
  return(ans)
}


setMethod("residualsAcd", signature(object = "ACDfit"), .acdfitresids)
setMethod("residualsAcd", signature(object = "ACDfilter"), .acdfitresids)
#' @exportMethod residualsAcd

#-------------------------------------------------------------------------------
# conditional sigma
sigmaAcd <- function(object)
{
  UseMethod("sigmaAcd")
}
.acdsigma = function(object)
{
  if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
    D = object@model$modeldata$index[1:object@model$modeldata$T]
  }
  ans = switch(class(object)[1],
               ACDfit = xts(object@fit$sigma, D),
               ACDfilter = xts(object@filter$sigma, D),
               ACDforecast = object@forecast$sigmaFor,
               ACDsim = {
                 ans = object@simulation$sigmaSim
                 rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
                 return(ans)
               },
               ACDpath ={
                 ans = object@path$sigmaSim
                 rownames(ans) = paste("T+",1:NROW(object@path$sigmaSim), sep="")
                 return(ans)
               },
               ACDroll = as.xts(object@forecast$density[,"Sigma",drop=FALSE]))
  return(ans)
}

setMethod("sigmaAcd", signature(object = "ACDfit"), .acdsigma)
setMethod("sigmaAcd", signature(object = "ACDfilter"), .acdsigma)
setMethod("sigmaAcd", signature(object = "ACDforecast"), .acdsigma)
setMethod("sigmaAcd", signature(object = "ACDpath"), .acdsigma)
setMethod("sigmaAcd", signature(object = "ACDsim"), .acdsigma)
setMethod("sigmaAcd", signature(object = "ACDroll"), .acdsigma)
#' @exportMethod sigmaAcd

#-------------------------------------------------------------------------------
# conditional skew
skew = function(object, transformed = TRUE, ...)
{
  UseMethod("skew")
}

.acdskew = function(object, transformed = TRUE)
{
  if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
    D = object@model$modeldata$index[1:object@model$modeldata$T]
  }
  if(is(object, "ACDfit")){
    if(transformed) ans = xts(object@fit$tskew, D) else ans = xts(object@fit$tempskew, D)
  } else if(is(object, "ACDfilter")){
    if(transformed) ans = xts(object@filter$tskew, D) else ans = xts(object@filter$tempskew, D)
  } else if(is(object, "ACDforecast")){
    if(transformed){
      ans = object@forecast$tskewFor
    } else{
      ans = apply(object@forecast$tskewFor, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
    }
  } else if(is(object, "ACDpath")){
    if(transformed){
      ans = object@path$skewSim
      rownames(ans) = paste("T+",1:NROW(object@path$skewSim), sep="")
    } else{
      ans = apply(object@path$skewSim, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@path$skewSim), sep="")

    }
  } else if(is(object, "ACDsim")){
    if(transformed){
      ans = object@simulation$skewSim
      rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
    } else{
      ans = apply(object@simulation$skewSim, 2, function(x) logtransform(x, object@model$sbounds[1], object@model$sbounds[2], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
    }
  } else if(is(object, "ACDroll")){
    if(transformed){
      ans = as.xts(object@forecast$density[,"Skew",drop=FALSE])
    } else{
      stop("\nACDroll object does return the (possibly) dynamic bounds needed for the transformation.")
    }
  } else{
    ans = NA
  }
  return(ans)
}


setMethod("skew", signature(object = "ACDfit"), .acdskew)
setMethod("skew", signature(object = "ACDfilter"), .acdskew)
setMethod("skew", signature(object = "ACDforecast"), .acdskew)
setMethod("skew", signature(object = "ACDpath"), .acdskew)
setMethod("skew", signature(object = "ACDsim"), .acdskew)
setMethod("skew", signature(object = "ACDroll"), .acdskew)
#' @exportMethod skew
#-------------------------------------------------------------------------------
# conditional shape1
shape1 = function(object, transformed = TRUE, ...)
{
  UseMethod("shape")
}

.acdshape1 = function(object, transformed = TRUE)
{
  if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
    D = object@model$modeldata$index[1:object@model$modeldata$T]
  }
  if(is(object, "ACDfit")){
    if(transformed) ans = xts(object@fit$tshape1, D) else ans = xts(object@fit$tempshape1, D)
  } else if(is(object, "ACDfilter")){
    if(transformed) ans = xts(object@filter$tshape1, D) else ans = xts(object@filter$tempshape1, D)
  } else if(is(object, "ACDforecast")){
    if(transformed){
      ans = object@forecast$tshape1For
    } else{
      ans = apply(object@forecast$tshape1For, 2, function(x) exptransform1(x, object@model$sbounds[3], object@model$sbounds[7], inverse = TRUE))
    }
  } else if(is(object, "ACDpath")){
    if(transformed){
      ans = object@path$shape1Sim
      rownames(ans) = paste("T+",1:NROW(object@path$shape1Sim), sep="")
    } else{
      ans = apply(object@path$shape1Sim, 2, function(x) exptransform1(x, object@model$sbounds[3], object@model$sbounds[7], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@path$shape1Sim), sep="")
    }
  }
  else if(is(object, "ACDsim")){
    if(transformed){
      ans = object@simulation$shape1Sim
      rownames(ans) = paste("T+",1:NROW(object@simulation$shape1Sim), sep="")
    } else{
      ans = apply(object@simulation$shape1Sim, 2, function(x) exptransform1(x, object@model$sbounds[3], object@model$sbounds[7], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@simulation$shape1Sim), sep="")
    }
  } else if(is(object, "ACDroll")){
    if(transformed){
      ans = as.xts(object@forecast$density[,"Shape1",drop=FALSE])
    } else{
      stop("\nACDroll object does return the (possibly) dynamic bounds needed for the transformation.")
    }
  } else{
    ans = NA
  }
  return(ans)
}

setMethod("shape1", signature(object = "ACDfit"), .acdshape1)
setMethod("shape1", signature(object = "ACDfilter"), .acdshape1)
setMethod("shape1", signature(object = "ACDforecast"), .acdshape1)
setMethod("shape1", signature(object = "ACDpath"), .acdshape1)
setMethod("shape1", signature(object = "ACDsim"), .acdshape1)
setMethod("shape1", signature(object = "ACDroll"), .acdshape1)

#' @exportMethod shape1
#------------------------------------------------------------------------------
shape2 = function(object, transformed = TRUE, ...)
{
  UseMethod("shape2")
}

.acdshape2 = function(object, transformed = TRUE)
{
  if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
    D = object@model$modeldata$index[1:object@model$modeldata$T]
  }
  if(is(object, "ACDfit")){
    if(transformed) ans = xts(object@fit$tshape2, D) else ans = xts(object@fit$tempshape2, D)
  } else if(is(object, "ACDfilter")){
    if(transformed) ans = xts(object@filter$tshape2, D) else ans = xts(object@filter$tempshape2, D)
  } else if(is(object, "ACDforecast")){
    if(transformed){
      ans = object@forecast$tshape2For
    } else{
      ans = apply(object@forecast$tshape2For, 2, function(x) logtransform(x, object@model$sbounds[5], object@model$sbounds[6], inverse = TRUE))
    }
  } else if(is(object, "ACDpath")){
    if(transformed){
      ans = object@path$shape2Sim
      rownames(ans) = paste("T+",1:NROW(object@path$shape2Sim), sep="")
    } else{
      ans = apply(object@path$shape2Sim, 2, function(x) logtransform(x, object@model$sbounds[5], object@model$sbounds[6], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@path$shape2Sim), sep="")
    }
  }
  else if(is(object, "ACDsim")){
    if(transformed){
      ans = object@simulation$shape2Sim
      rownames(ans) = paste("T+",1:NROW(object@simulation$shape2Sim), sep="")
    } else{
      ans = apply(object@simulation$shape2Sim, 2, function(x) logtransform(x, object@model$sbounds[5], object@model$sbounds[6], inverse = TRUE))
      rownames(ans) = paste("T+",1:NROW(object@simulation$shape2Sim), sep="")
    }
  } else if(is(object, "ACDroll")){
    if(transformed){
      ans = as.xts(object@forecast$density[,"shape2",drop=FALSE])
    } else{
      stop("\nACDroll object does return the (possibly) dynamic bounds needed for the transformation.")
    }
  } else{
    ans = NA
  }
  return(ans)
}

setMethod("shape2", signature(object = "ACDfit"), .acdshape2)
setMethod("shape2", signature(object = "ACDfilter"), .acdshape2)
setMethod("shape2", signature(object = "ACDforecast"), .acdshape2)
setMethod("shape2", signature(object = "ACDpath"), .acdshape2)
setMethod("shape2", signature(object = "ACDsim"), .acdshape2)
setMethod("shape2", signature(object = "ACDroll"), .acdshape2)

#' @exportMethod shape2
#-------------------------------------------------------------------------------
# conditional skewness
skewness = function(object, ...)
{
  UseMethod("skewness")
}

.acdskewness = function(object, ...)
{
  if(is(object, "ACDfit") | is(object, "ACDfilter") | is(object, "ACDroll")){
    S = object@fit$skewness
    if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
      D = object@model$modeldata$index[1:object@model$modeldata$T]
      S = xts(S, D)
    } else{
      S = xts(S, as.POSIXct(rownames(object@forecast$density)))
    }
  } #else {
    #skew = skew(object)
    #shape = shape(object)
    #lambda = object@model$pars["ghlambda",1]
    #m = NCOL(skew)
    #n = NROW(skew)
    #S = matrix(NA, ncol = m, nrow = n)
    #for(i in 1:m){ S[,i] = dskewness(distribution = dist, skew = skew[,i], shape = shape[,i], lambda = lambda) }
  #}
  return(S)
}

setMethod("skewness", signature(object = "ACDfit"), .acdskewness)
setMethod("skewness", signature(object = "ACDfilter"), .acdskewness)
#setMethod("skewness", signature(object = "ACDforecast"), .acdskewness)
#setMethod("skewness", signature(object = "ACDpath"), .acdskewness)
#setMethod("skewness", signature(object = "ACDsim"), .acdskewness)
setMethod("skewness", signature(object = "ACDroll"), .acdskewness)
#' @exportMethod skewness
#-------------------------------------------------------------------------------
# conditional excess kurtosis
kurtosis = function(object, ...)
{
  UseMethod("kurtosis")
}

.acdkurtosis = function(object, ...)
{
  if(is(object, "ACDfit") | is(object, "ACDfilter") | is(object, "ACDroll")){
    S = object@fit$kurtosis
    if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
      D = object@model$modeldata$index[1:object@model$modeldata$T]
      S = xts(S, D)
    } else{
      S = xts(S, as.POSIXct(rownames(object@forecast$density)))
    }
  } #else {
  #skew = skew(object)
  #shape = shape(object)
  #lambda = object@model$pars["ghlambda",1]
  #m = NCOL(skew)
  #n = NROW(skew)
  #S = matrix(NA, ncol = m, nrow = n)
  #for(i in 1:m){ S[,i] = dkurtosis(distribution = dist, skew = skew[,i], shape = shape[,i], lambda = lambda) }
  #}
  return(S)
}

setMethod("kurtosis", signature(object = "ACDfit"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDfilter"), .acdkurtosis)
#setMethod("kurtosis", signature(object = "ACDforecast"), .acdkurtosis)
#setMethod("kurtosis", signature(object = "ACDpath"), .acdkurtosis)
#setMethod("kurtosis", signature(object = "ACDsim"), .acdkurtosis)
setMethod("kurtosis", signature(object = "ACDroll"), .acdkurtosis)
#' @exportMethod kurtosis
#-------------------------------------------------------------------------------
# conditional Pskew
Pskewness = function(object, ...)
{
  UseMethod("Pskewness")
}

.acdPskew = function(object, ...)
{
  if(is(object, "ACDfit")){
    S = object@fit$Pskewness
  }else if(is(object, "ACDfilter")){
      S = object@filter$Pskewness
    }
    if(class(object)[1] == "ACDfit" | class(object)[1] == "ACDfilter"){
      D = object@model$modeldata$index[1:object@model$modeldata$T]
      S = xts(S, D)
    } else{
      S = xts(S, as.POSIXct(rownames(object@forecast$density)))
    }
  #} else {
  #skew = skew(object)
  #shape = shape(object)
  #lambda = object@model$pars["ghlambda",1]
  #m = NCOL(skew)
  #n = NROW(skew)
  #S = matrix(NA, ncol = m, nrow = n)
  #for(i in 1:m){ S[,i] = dPskew(distribution = dist, skew = skew[,i], shape = shape[,i], lambda = lambda) }
  #}
  return(S)
}

setMethod("Pskewness", signature(object = "ACDfit"), .acdPskew)
setMethod("Pskewness", signature(object = "ACDfilter"), .acdPskew)
#setMethod("Pskew", signature(object = "ACDforecast"), .acdPskew)
#setMethod("Pskew", signature(object = "ACDpath"), .acdPskew)
#setMethod("Pskew", signature(object = "ACDsim"), .acdPskew)
#setMethod("Pskewness", signature(object = "ACDroll"), .acdPskew)
#' @exportMethod Pskewness
#-------------------------------------------------------------------------------
# likelihood method
#' @exportMethod acdlikelihood
acdlikelihood = function(object){
  UseMethod("acdlikelihood")
}
.acdfitLikelihood = function(object)
{
  if(is(object, "ACDfit")){
    return(c("ACD"=object@fit$LLH, "GARCH" = object@model$garchLL))
  } else if(is(object, "ACDfilter")){
    return(c("ACD"=object@filter$LLH, "GARCH" = object@model$garchLL))
  } else{
    if(!is.null(object@model$noncidx)) stop("\nroll object contains non-converged windows...use resume method first.\n")
    tmp = sapply(object@model$LL, function(x) x$log.lik)
    colnames(tmp) = sapply(object@model$LL, function(x) as.character(x$date))
    rownames(tmp) = c("ACD", "GARCH")
    return(tmp)
  }
}

setMethod("acdlikelihood", signature(object = "ACDfit"), .acdfitLikelihood)
setMethod("acdlikelihood", signature(object = "ACDfilter"), .acdfitLikelihood)
#------
# infocriteria method
.acdinfocriteria = function(object)
{
  if(is(object, "ACDfilter")){
    # np = sum(object@filter$ipars[,2])
    # all parameters fixed
    acdnp = 0
  } else{
    acdnp = sum(object@fit$ipars[,4])
    garchnp = sum(object@fit$ipars[1:13,4])
  }
  acditest = rugarch:::.information.test(acdlikelihood(object)[1], nObs = length(fittedAcd(object)),
                                         nPars = acdnp)
  garchitest = rugarch:::.information.test(acdlikelihood(object)[2], nObs = length(fittedAcd(object)),
                                           nPars = garchnp)
  itestm = matrix(0, ncol = 2, nrow = 4)
  itestm[1,1] = acditest$AIC
  itestm[2,1] = acditest$BIC
  itestm[3,1] = acditest$SIC
  itestm[4,1] = acditest$HQIC

  itestm[1,2] = garchitest$AIC
  itestm[2,2] = garchitest$BIC
  itestm[3,2] = garchitest$SIC
  itestm[4,2] = garchitest$HQIC
  colnames(itestm) = c("ACD", "GARCH")
  rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
  return(itestm)
}

setMethod("infocriteria", signature(object = "ACDfit"), .acdinfocriteria)
#--------------------------
# show methods

setMethod("show",
          signature(object = "ACDspec"),
          function(object){
            mmodel = object@model$mmodel
            vmodel = object@model$vmodel$model
            model = object@model
            modelinc = object@model$modelinc
            cat(paste("\n*---------------------------------*", sep = ""))
            cat(paste("\n*          ACD Model Spec         *", sep = ""))
            cat(paste("\n*---------------------------------*", sep = ""))
            cat("\n\nConditional mean Dynamics \t")
            cat(paste("\n-----------------------------------", sep = ""))
            cat(paste("\nARMA Order\t:","(",modelinc[2],",",modelinc[3],")",sep=""))
            cat(paste("\nGARCH-in-mean effect\t:",as.logical(modelinc[4]),sep = ""))
            cat(paste("\nSkew-in-mean effect\t:",as.logical(modelinc[6]),sep = ""))
            cat(paste("\nAdjust conditional mean\t:",as.logical(modelinc[36]),sep= ""))
            cat("\n\nConditional Variance Dynamics \t")
            cat(paste("\n-----------------------------------", sep = ""))
            cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
            if(vmodel == "fGARCH"){
              cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
            }
            cat("Mean Model\t: ARFIMA(", modelinc[2],",",",",modelinc[3],")\n", sep = "")
            cat("Distribution\t:", model$dmodel$model,"\n")
            if(sum(object@model$dmodel$skewOrder)>0){
              cat("\nConditional Skew Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[15], ",", modelinc[16], ",", modelinc[17],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
            }
            if(sum(object@model$dmodel$shape1Order)>0){
              cat("\nConditional Shape1 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape1 Model\t: ", object@model$dmodel$shape1model,"(", modelinc[20], ",", modelinc[21], ",", modelinc[22],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape1shock],"\n",sep=""))
            }
            if(sum(object@model$dmodel$shape2Order)>0){
              cat("\nConditional Shape2 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape2 Model\t: ", object@model$dmodel$shape2model,"(", modelinc[25], ",", modelinc[26], ",", modelinc[27],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape2shock],sep=""))
            }
            cat("\n")
            invisible(object)
          })

setMethod("show",
          signature(object = "ACDfit"),
          function(object){
            vmodel = object@model$vmodel$model
            model = object@model
            modelinc = object@model$modelinc
            cat(paste("\n*---------------------------------*", sep = ""))
            cat(paste("\n*          ACD Model Fit          *", sep = ""))
            cat(paste("\n*---------------------------------*", sep = ""))
            cat("\n\nConditional Variance Dynamics \t")
            cat(paste("\n-----------------------------------", sep = ""))
            cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
            if(vmodel == "fGARCH"){
              cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
            }
            cat("Mean Model\t: ARFIMA(", modelinc[2],",",",",modelinc[3],")\n", sep = "")
            cat("Adjust mean\t:", as.logical(modelinc[36]),"\n",sep = "")
            cat("Distribution\t:", model$dmodel$model,"\n")
            if(sum(object@model$dmodel$skewOrder)>0){
              cat("\nConditional Skew Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[15], ",", modelinc[16], ",", modelinc[17],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
            }
            if(sum(object@model$dmodel$shape1Order)>0){
              cat("\nConditional Shape1 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape1 Model\t: ", object@model$dmodel$shape1model,"(", modelinc[20], ",", modelinc[21], ",", modelinc[22],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape1shock],sep=""))
            }
            if(sum(object@model$dmodel$shape2Order)>0){
              cat("\nConditional Shape2 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape2 Model\t: ", object@model$dmodel$shape2model,"(", modelinc[25], ",", modelinc[26], ",", modelinc[27],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape2shock],sep=""))
            }
            if(object@fit$convergence == 0){
              cat("\n\nOptimal Parameters")
              cat(paste("\n------------------------------------\n",sep=""))
              print(round(object@fit$matcoef,6), digits = 5)
              cat("\nRobust Standard Errors:\n")
              print(round(object@fit$robust.matcoef,6), digits = 5)
              if(!is.null(object@fit$hessian.message)){
                cat(paste("\n", object@fit$hessian.message))
              }
              cat("\nLogLikelihood :", unname(object@fit$LLH[1]), "\n")
              stdresid = object@fit$residuals/object@fit$sigma
              itestm = infocriteria(object)
              cat("\nInformation Criteria")
              cat(paste("\n------------------------------------\n",sep=""))
              print(itestm,digits=5)
              cat("\nWeighted Ljung-Box Test on Standardized Residuals")
              cat(paste("\n------------------------------------\n",sep=""))
              tmp1 = rugarch:::.weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
              print(tmp1, digits = 4)
              cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
              cat("\nH0 : No serial correlation\n")
              cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
              cat(paste("\n------------------------------------\n",sep=""))
              tmp2 = rugarch:::.weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
              print(tmp2, digits = 4)
              cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
              cat("\n\nWeighted ARCH LM Tests")
              cat(paste("\n------------------------------------\n",sep=""))
              gdf = sum(modelinc[8:9])
              L2 = rugarch:::.weightedarchlmtest(residualsAcd(object), sigmaAcd(object), lags  = gdf+1, fitdf=gdf)
              L5 = rugarch:::.weightedarchlmtest(residualsAcd(object), sigmaAcd(object), lags  = gdf+3, fitdf=gdf)
              L10 = rugarch:::.weightedarchlmtest(residualsAcd(object), sigmaAcd(object), lags = gdf+5, fitdf=gdf)
              alm = matrix(0,ncol = 4,nrow = 3)
              alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
              alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
              alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
              colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
              rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
              print(alm,digits = 4)
              nyb = rugarch:::.nyblomTest(object)
              if(is.character(nyb$JointCritical)){
                colnames(nyb$IndividualStat)<-""
                cat("\nNyblom stability test")
                cat(paste("\n------------------------------------\n",sep=""))
                cat("Joint Statistic: ", "no.parameters>20 (not available)")
                cat("\nIndividual Statistics:")
                print(nyb$IndividualStat, digits = 4)
                cat("\nAsymptotic Critical Values (10% 5% 1%)")
                cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
                cat("\n\n")
              } else{
                colnames(nyb$IndividualStat)<-""
                cat("\nNyblom stability test")
                cat(paste("\n------------------------------------\n",sep=""))
                cat("Joint Statistic: ", round(nyb$JointStat,4))
                cat("\nIndividual Statistics:")
                print(nyb$IndividualStat, digits = 4)
                cat("\nAsymptotic Critical Values (10% 5% 1%)")
                cat("\nJoint Statistic:     \t", round(nyb$JointCritical, 3))
                cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
                cat("\n\n")
              }
              #cat("Sign Bias Test")
              #cat(paste("\n------------------------------------\n",sep=""))
              #sgtest = signbias(object)
              #print(sgtest, digits = 4)
              #cat("\n")
              #cat("\nAdjusted Pearson Goodness-of-Fit Test:")
              #cat(paste("\n------------------------------------\n",sep=""))
              #gofm = gof(object,c(20, 30, 40, 50))
              #print(gofm, digits = 4)
              #cat("\n")
              cat("\nElapsed time :", object@fit$timer,"\n\n")
            } else{
              cat("\nConvergence Problem:")
              cat("\nSolver Message:", object@fit$message,"\n\n")

            }
            invisible(object)
          })
setMethod("show",
          signature(object = "ACDfilter"),
          function(object){
            vmodel = object@model$vmodel$model
            model = object@model
            modelinc = object@model$modelinc
            cat(paste("\n*---------------------------------*", sep = ""))
            cat(paste("\n*          ACD Model Filter       *", sep = ""))
            cat(paste("\n*---------------------------------*", sep = ""))
            cat("\n\nConditional Variance Dynamics \t")
            cat(paste("\n-----------------------------------", sep = ""))
            cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
            if(vmodel == "fGARCH"){
              cat(paste("fGARCH Sub-Model\t: ", model$vmodel$vsubmodel, "\n", sep = ""))
            }
            cat("Mean Model\t: ARFIMA(", modelinc[2],",",",",modelinc[3],")\n", sep = "")
            cat("Distribution\t:", model$dmodel$model,"\n")
            if(sum(object@model$dmodel$skewOrder)>0){
              cat("\nConditional Skew Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Skew Model\t: ", object@model$dmodel$skewmodel,"(", modelinc[15], ",", modelinc[16], ",", modelinc[17],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$skewshock],"\n",sep=""))
            }
            if(sum(object@model$dmodel$shape1Order)>0){
              cat("\nConditional Shape1 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape1 Model\t: ", object@model$dmodel$shape1model,"(", modelinc[20], ",", modelinc[21], ",", modelinc[22],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape1shock],sep=""))
            }
            if(sum(object@model$dmodel$shape2Order)>0){
              cat("\nConditional Shape2 Dynamics \t")
              cat(paste("\n-----------------------------------", sep = ""))
              cat(paste("\nACD Shape2 Model\t: ", object@model$dmodel$shape2model,"(", modelinc[25], ",", modelinc[26], ",", modelinc[27],")", sep=""))
              cat(paste("\nShock type:\t", c("Std.Residuals", "Residuals")[object@model$dmodel$shape2shock],sep=""))
            }
            cat("\n\nFilter Parameters")
            cat(paste("\n------------------------------------\n",sep=""))
            print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
            cat("\nLogLikelihood :", unname(likelihood(object)), "\n")
            stdresid = object@filter$residuals/object@filter$sigma
            cat("\nQ-Statistics on Standardized Residuals")
            cat(paste("\n---------------------------------------\n",sep=""))
            tmp1 = rugarch:::.weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
            print(tmp1, digits = 4)
            cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
            cat("\nH0 : No serial correlation\n")
            cat("\nQ-Statistics on Standardized Squared Residuals")
            cat(paste("\n---------------------------------------\n",sep=""))
            tmp2 = rugarch:::.weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
            print(tmp2, digits = 4)
            cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
            cat("\n\nARCH LM Tests")
            cat(paste("\n---------------------------------------\n",sep=""))
            gdf = sum(modelinc[8:9])
            L2 = rugarch:::.weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+1, fitdf=gdf)
            L5 = rugarch:::.weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+3, fitdf=gdf)
            L10 = rugarch:::.weightedarchlmtest(residuals(object), sigma(object), lags = gdf+5, fitdf=gdf)
            alm = matrix(0,ncol = 4,nrow = 3)
            alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
            alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
            alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
            colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
            rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
            print(alm,digits = 4)
            cat("\nElapsed time :", object@filter$timer,"\n\n")
            invisible(object)
          })

setMethod("show",
          signature(object = "ACDforecast"),
          function(object){
            vmodel = object@model$vmodel$model
            model = object@model
            cat(paste("\n*------------------------------------*", sep = ""))
            cat(paste("\n*       ACD Model Forecast           *", sep = ""))
            cat(paste("\n*------------------------------------*", sep = ""))
            cat(paste("\nGARCH Model  : ", vmodel, sep = ""))
            n.ahead = object@forecast$n.ahead
            cat(paste("\nHorizon      : ", n.ahead, sep = ""))
            cat(paste("\nRoll Steps   : ", object@forecast$n.roll, sep = ""))
            n.start = object@forecast$n.start
            if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
            cat(paste("\nOut of Sample: ", infor, "\n", sep = ""))
            cat(paste("\n0-roll forecast [T0=", as.character(object@model$modeldata$index[length(object@model$modeldata$data)]), "]:\n", sep=""))
            zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1],
                       object@forecast$tskewFor[,1], object@forecast$tshape1For[,1],object@forecast$tshape2For[,1] )
            colnames(zz) = c("series", "sigma", "skew", "shape1", "shape2")
            print(zz, digits = 4)
            cat("\n\n")
          })

setMethod("show",
          signature(object = "ACDroll"),
          function(object){
            if(!is.null(object@model$noncidx)){
              cat("\nObject contains non-converged estimation windows. Use resume method to re-estimate.\n")
              invisible(object)
            } else{
              cat(paste("\n*-------------------------------------*", sep = ""))
              cat(paste("\n*              ACD Roll               *", sep = ""))
              cat(paste("\n*-------------------------------------*", sep = ""))
              N = object@model$n.refits
              gmodel = object@model$spec@model$vmodel$model
              model = object@model$spec@model
              modelinc = object@model$spec@model$modelinc
              cat("\nNo.Refits\t\t:", N)
              cat("\nRefit Horizon\t:", object@model$refit.every)
              cat("\nNo.Forecasts\t:", NROW(object@forecast$density))
              cat(paste("\nGARCH Model\t\t: ", gmodel, "(",modelinc[8],",",modelinc[9],")\n", sep = ""))
              cat("Mean Model\t\t: ARFIMA(", modelinc[2],",",",",modelinc[3],")\n", sep = "")
              if(sum(model$dmodel$skewOrder)>0){
                cat(paste("\nACD Skew Model\t: ", model$dmodel$skewmodel,"(", modelinc[15], ",", modelinc[16], ",", modelinc[17],")", sep=""))
              }
              if(sum(model$dmodel$shape1Order)>0){
                cat(paste("\nACD Shape1 Model\t: ", model$dmodel$shape1model,"(", modelinc[20], ",", modelinc[21], ",", modelinc[22],")", sep=""))
              }
              if(sum(model$dmodel$shape2Order)>0){
                cat(paste("\nACD Shape2 Model\t: ", model$dmodel$shape2model,"(", modelinc[25], ",", modelinc[26], ",", modelinc[27],")", sep=""))
              }
              cat("\nDistribution\t:", model$dmodel$model,"\n")
              cat("\nForecast Density:\n")
              print(round(head(object@forecast$density),4))
              cat("\n..........................\n")
              print(round(tail(object@forecast$density),4))
              cat("\nElapsed:", format(object@model$elapsed))
              cat("\n")
              invisible(object)
            }
          })
#' @exportMethod show
#----
# coef method
#' @exportMethod coefacd

coefacd <- function(object){
  UseMethod("coefacd")
}
.acdfitcoef = function(object)
{
  if(is(object, "ACDfit")){
    return(object@fit$coef)
  } else{
    return(object@model$pars[object@model$pars[,3]==1, 1])
  }
}

setMethod("coefacd", signature(object = "ACDfit"), .acdfitcoef)
#setMethod("coefacd", signature(object = "ACDfilter"), .acdfitcoef)

