#---
#' @description Returns the bounds to which a distribution density exist
#' @title Distribution bounds
#' @usage distbounds(distribution)
#' @param distribution: a string input with three supported conditional distribution: sgt, sged and sst
#' @return a vector with the Lowerbounds(LB), Upperbounds(UB) for skew and shape parameters to which the denstiy exist
#-----
## Set distribution bounds
#-----
TinY = 1e-08

distbounds = function(distribution) {
  distribution = tolower(distribution[1])  #Translate all distribution names to lower case
  distribution = match.arg(distribution, c("sgt", "sged","sst"))
  ans = .DistributionBounds(distribution)
  ans = c(ans$skew.LB, ans$skew.UB, ans$shape1.LB, ans$shape1.UB, ans$shape2.LB, ans$shape2.UB)
  names(ans) = c("skew(LB)", "skew(UB)", "shape1(LB)", "shape1(UB)", "shape2(LB)", "shape2(UB)")
  return(ans)
}
#----
#' @description logtransform or Inverse-logtransform w.r.t specify LB and UB bounds
#' @title Log Transformation
#' @param x value to be transfromed
#' @param lower: Lower bound
#' @param upper: Uper bound
#' @param inverse: Logical value. The defaul is FALSE
#' @usage logtransform(x,lower,upper,inverse = FALSE)
#' @export logtransform
#-----
logtransform = function(x, lower, upper, inverse = FALSE) {
  if (!inverse) {
    ans = lower + (upper - lower)/(1 + exp(-1 * x))
  } else {
    ans = -1 * log(-(upper - x)/(-x + lower))
  }
  return(ans)
}
#----
#' @description Exponential-transform or Inverse Exponential transform w.r.t specify LB and UB bounds
#' @title Exponential Transformation
#' @param x value to be transfromed
#' @param lower: Lower bound
#' @param rate: exponential rate. The default is 1
#' @param inverse: Logical value. The default is FALSE
#' @usage exptransform(x,lower,rate,inverse = FALSE)
#' @export exptransform1
exptransform1 = function(x, lower, rate = 1, inverse = FALSE) {
  if (!inverse) {
    ans = lower + exp(rate * x)
  }
  else {
    ans = (1/rate)*log((x - lower))
  }
  return(ans)
}
#-----
#' @description Exponential-transform or Inverse Exponential transform w.r.t specify LB and UB bounds
#' @title Exponential Transformation
#' @param x value to be transfromed
#' @param lower: Lower bound
#' @param upper: Upper bound
#' @param rate: exponential rate. The default is 1
#' @param inverse: Logical value. The default is FALSE
#' @usage exptransform2(x,lower,upper, rate,inverse = FALSE)
#' @export exptransform2
exptransform2 = function(x, lower, upper, rate = 1, inverse = FALSE) {
  if (!inverse) {
    ans = lower + upper * exp(-rate * x)
  }
    else {
    ans = -(1/rate)*log((x - lower)/upper)
  }
  return(ans)
}


#--------------------------
# Set distribution bounds
#-------------------------
.DistributionBounds = function(distribution) {
   if (distribution == "sgt") {
    skew = 0
    skew.LB = -0.9999
    skew.UB = 0.9999
    shape1 = 2
    shape1.LB = 0.5
    #shape1.LB = 1
    shape1.UB = 10
    shape2 = 8
    shape2.LB = 2
    shape2.UB = 50
    include.shape1 = TRUE
    include.shape2 = TRUE
  }
  if (distribution == "sged") {
    skew = 0
    skew.LB = -0.999
    skew.UB = 0.999
    shape1 = 2
    shape1.LB = 0.5
    #shape1.LB = 1
    shape1.UB = 10
    shape2 = 999
    shape2.LB = 2
    shape2.UB = 50
    include.shape1 = TRUE
    include.shape2 = FALSE
  }
  if (distribution == "sst") {
    skew = 0
    skew.LB = -0.999
    skew.UB = 0.999
    shape1 = 2
    shape1.LB = 0.5
    #shape1.LB = 1
    shape1.UB = 10
    shape2 = 8
    shape2.LB = 2
    shape2.UB = 50
    include.shape1 = FALSE
    include.shape2 = TRUE
  }
  skew0 = 0
  shape10 = 0
  shape20 = 0
  include.skew = TRUE
   ans = list(shape1 = shape1, shape1.LB = shape1.LB, shape1.UB = shape1.UB,
              shape2 = shape2, shape2.LB = shape2.LB, shape2.UB = shape2.UB,
              skew = skew, skew.LB = skew.LB, skew.UB = skew.UB,
              include.skew = include.skew,include.shape1 = include.shape1,
              include.shape2 = include.shape2,
              skew0 = skew0, shape10 = shape10, shape20 = shape20)
  return(ans)
}

.acdskewbounds = function(acdOrder, unconpar, distribution, dbounds) {
  .eps = .Machine$double.eps
  skew.LB = dbounds[1]
  skew.UB = dbounds[2]
  par = par.LB = par.UB = numeric()
  # intercept
  par[1] = logtransform(unconpar, skew.LB, skew.UB, inverse = TRUE)
  par.LB[1] = logtransform(skew.LB * 0.99, skew.LB, skew.UB,inverse = TRUE)
  par.UB[1] = logtransform(skew.UB * 0.99, skew.LB, skew.UB,inverse = TRUE)
  # acdOrder
  par = c(par,0.1,0.1,0.5)
  par.LB = c(par.LB, -1, -1,-1+.eps)
  par.UB = c(par.UB,1,1,1-.eps)
  return(list(skewpars = par, skewpar.LB = par.LB, skewpar.UB = par.UB, sk0 = par[1]))
}
.acdshape1bounds = function(acdOrder, unconpar, distribution, dbounds) {
  .eps = .Machine$double.eps
  shape1.LB = dbounds[1]
  shape1.UB = dbounds[2]
  par = par.LB = par.UB = numeric()
   # intercept
  par[1] = exptransform1(unconpar, shape1.LB,inverse = TRUE)
  par.LB[1] = exptransform1(shape1.LB*1.001,shape1.LB,inverse = TRUE)
  par.UB[1] = exptransform1(shape1.UB*0.99,shape1.LB,inverse = TRUE)
  #par[1] = logtransform(unconpar,shape1.LB,shape1.UB,inverse = TRUE)
  #par.LB[1] = logtransform(shape1.LB * 1.0001, shape1.LB, shape1.UB,inverse = TRUE)
  #par.UB[1] = logtransform(shape1.UB * 0.999, shape1.LB, shape1.UB,inverse = TRUE)
  # alpha1 and alpha2 with lower/upper bounds
  par = c(par,0,0,0.5)
  par.LB = c(par.LB, -1, -1,-0.9+.eps)
  #par.LB = c(par.LB, -1, -1,0+.eps)
  par.UB = c(par.UB,1,1,0.9-.eps)
  return(list(shapepars = par, shapepar.LB = par.LB, shapepar.UB = par.UB, sh0 = par[1]))
}
.acdshape2bounds = function(acdOrder, unconpar, distribution, dbounds) {
  .eps = .Machine$double.eps
  shape2.LB = dbounds[1]
  shape2.UB = dbounds[2]
  par = par.LB = par.UB = numeric()
     # intercept
  par[1] = logtransform(unconpar,shape2.LB,shape2.UB,inverse = TRUE)
  #par[1] = -5
  par.LB[1] = logtransform(shape2.LB * 1.001, shape2.LB, shape2.UB,inverse = TRUE)
  par.UB[1] = logtransform(shape2.UB * 0.99, shape2.LB, shape2.UB,inverse = TRUE)
  #par[1] = exptransform2(unconpar,shape2.LB,shape2.UB,inverse = TRUE)
  #par.LB[1] = exptransform2(shape2.LB * 1.0001,shape2.LB,shape2.UB,inverse = TRUE)
  #par.UB[2] = exptransform2(shape2.UB * 0.999,shape2.LB,shape2.UB,inverse = TRUE)
    # alpha1 and alpha2 with lower/upper bounds
    par = c(par,0,0,0.5)
    par.LB = c(par.LB, -2, -2,-1+.eps)
    #par.LB = c(par.LB, -2, -2,0+.eps)
    par.UB = c(par.UB,2,2,1-.eps)
  return(list(shapepars = par, shapepar.LB = par.LB, shapepar.UB = par.UB, sh0 = par[1]))
}
# logistic transformation and inverse transformation
.invlogtransform = function(y, LB, UB)
{
  x = -1*log(-(UB-y)/(-y+LB))
  return(x)
}

.logtransform = function(x, LB, UB)
{
  y = LB + (UB-LB)/(1+exp(-1*x))
  return(y)
}


.exptransform2 = function(x, lower,upper, rate=1){
  lower+upper*exp(-rate*x)
}
.exptransform1 = function(x, lower, rate=1){
  lower+exp(rate*x)
}

.invexptransform2 = function(x, lower, upper, rate=1)
{
  -(1/rate)*log((x - lower)/upper)
}
.invexptransform1 = function(x, lower, rate=1)
{
  (1/rate)*log(x - lower)
}
