#----------------------------------------------------------------------------------
# quantile method
#' @exportMethod acdquantile
acdquantile <- function(x, probs = c(0.01,0.05)){
  UseMethod("acdquantile")
}
.acdfitquantile = function(x, probs = c(0.01, 0.05))
{
  if(class(x)=="ACDroll"){
    d = x@model$spec@model$dmodel$model
    skew = x@forecast$density[,"Skew"]
    shape1 = x@forecast$density[,"Shape1"]
    shape2 = x@forecast$density[,"Shape2"]
    s = x@forecast$density[,"Sigma"]
    m = x@forecast$density[,"Mu"]
    Q = matrix(NA, ncol = length(probs), nrow = length(s))
    for(i in 1:length(probs)) Q[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape,
                                                            lambda = lambda) * as.numeric(s)
    colnames(Q) = paste("q[", probs,"]", sep="")
    Q = xts(Q, as.POSIXct(rownames(x@forecast$density)))
  } else{
    d = x@model$dmodel$model
    mu = fittedAcd(x)
    sig = sigmaAcd(x)
    tskew = skew(x)
    tshape1 = shape1(x)
    tshape2 = shape2(x)
    Q = matrix(NA, nrow = length(sig), ncol = length(probs))
    for(i in seq_along(probs)){
      if(d == "sgt") quant = qsgt(probs[i],lambda = tskew,p = tshape1,q = tshape2/tshape1)
      if(d == "sged") quant = qsgt(probs[i],lambda = tskew,p = tshape1)
      if(d == "sst") quant = qsgt(probs[i],lambda = tskew,p = 2,q = tshape2/2)
      Q[,i] = mu + sig*quant
    }
    colnames(Q) = paste("q[", probs,"]", sep="")
    Q = xts(Q, index(sig))
  }
  return(Q)
}

setMethod("acdquantile", signature(x = "ACDfit"),  .acdfitquantile)
setMethod("acdquantile", signature(x = "ACDfilter"),  .acdfitquantile)
setMethod("acdquantile", signature(x = "ACDroll"),  .acdfitquantile)

