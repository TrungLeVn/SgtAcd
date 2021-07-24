#------------------------------------
# Model estimation routines with sGARCH specification for variance process
#' @useDynLib SgtAcd
#' @import rugarch
#------------------------------------
.sacdLLH = function(pars, arglist) {
  if (arglist$transform) {
    pars = logtransform(pars, arglist$LB, arglist$UB)
  }
  # prepare inputs
  eps = .Machine$double.eps
  data = arglist$data
  assign("x_pars", pars, envir = arglist$garchenv)
  if (!is.null(arglist$n.old)){
    Nx = arglist$n.old} else {
      Nx = length(data)}
  model = arglist$model
  estidx = arglist$estidx
  idx = model$pidx
  ipars = arglist$ipars
  ipars[estidx, 1] = pars
  T = length(data)
  fit.control = arglist$fit.control
  m = model$maxOrder
  N = c(m, T)
  modelinc = model$modelinc
  rx = .arfimaxfilteracd(modelinc, ipars[, 1], idx, data = data, N = N, arglist)
  # res is the residuals of mean process. zrf is the standardized residuals
  res = rx$res
  zrf = rx$zrf
  res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
  if( !is.null(arglist$n.old) ){
    rx = .arfimaxfilteracd(modelinc, ipars[,1], idx,data = data[1:Nx], N = c(m, Nx), arglist)
    res2 = rx$res
    res2[is.na(res2) | !is.finite(res2) | is.nan(res2)] = 0
    mvar = mean(res2*res2)
  } else{
    mvar = mean(res*res)
  }
  persist = (sum(ipars[idx["alpha",1]:idx["alpha",2],1]) + sum(ipars[idx["beta",1]:idx["beta",2],1]))
  # unconditional sigma value
  #mvar = mean(res*res)
  if(modelinc[7]>0){
    ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1])
    hEst = mvar
  } else{
    mv = 0
    persist = sum(ipars[idx["alpha",1]:idx["alpha",2],1])+sum(ipars[idx["beta",1]:idx["beta",2],1])
    ipars[idx["omega",1],1] = mvar * (1 - persist) - mv
    hEst = mvar
    assign("omega", ipars[idx["omega",1],1], arglist$garchenv)
  }
  if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data) * (1 - persist)
  assign("racd_ipars", ipars, envir = arglist$garchenv)
  if(fit.control$stationarity == 1){
    if(!is.na(persist) && persist >= 1) return(llh = get("racd_llh", arglist$garchenv) + 2*(abs(get("racd_llh", arglist$garchenv))))
  }

  sbounds = model$sbounds
  skhEst = arglist$skhEst
  tempskew  = double(length = T)
  tempshape1 = double(length = T)
  tempshape2 = double(length = T)
  tskew 	= double(length = T)
  tshape1 	= double(length = T)
  tshape2 	= double(length = T)
  h 		= double(length = T)
  z 		= double(length = T)
  constm 	= double(length = T)
  condm 	= double(length = T)
  llh 	= double(length = 1)
  LHT 	= double(length = T)
  skew = double(length = T)
  kurt = double(length = T)
  pskew = double(length = T)
  ans = try(.C("sacd",
               model = as.integer(modelinc),
               pars = as.double(ipars[,1]),
               idx = as.integer(idx[,1]-1),
               hEst = as.double(hEst),
               x = as.double(data),
               res = as.double(res),
               e = double(T),
               zrf = as.double(zrf),
               constm = double(T),
               condm = double(T),
               m = as.integer(m),
               T = as.integer(T),
               h = double(T),
               z = double(T),
               tempskew = double(T),
               tempshape1 = double(T),
               tempshape2 = double(T),
               skhEst = as.double(skhEst),
               tskew = double(T),
               tshape1 = double(T),
               tshape2 = double(T),
               sbounds = as.double(sbounds),
               llh = double(1),
               LHT = double(T),
               skew = double(T),
               kurt = double(T),
               pskew = double(T),
               PACKAGE = "SgtAcd"), silent = TRUE )
  if( inherits(ans, "try-error") ){
    cat(paste("\nacdfit-->warning: ", ans,"\n", sep=""))
    return( llh = get("racd_llh", arglist$garchenv) + 0.1*abs( get("racd_llh", arglist$garchenv) ) )
  }

  z = ans$z
  h = ans$h
  res = ans$res
  llh  = ans$llh
  tskew  = ans$tskew
  tshape1 = ans$tshape1
  tshape2 = ans$tshape2
  tempskew = ans$tempskew
  tempshape1 = ans$tempshape1
  tempshape2 = ans$tempshape2
  skew = ans$skew;
  kurt = ans$kurt;
  pskew = ans$pskew;
  if( is.finite(llh) && !is.na(llh) && !is.nan(llh) ){
    assign("racd_llh", llh, envir = arglist$garchenv)
  } else {
    llh = (get("racd_llh", arglist$garchenv) + 0.1*(abs(get("racd_llh",arglist$garchenv))))
  }
  # LHT = raw scores
  LHT = -ans$LHT
  ans = switch(arglist$returnType,
               llh = arglist$fnscale*llh,
               LHT = LHT,
               all = list(llh = llh, h = h, res = res, z = z, kappa = kappa,
                          tskew = tskew, tshape1 = tshape1, tempshape1 = tempshape1,
                          tshape2 = tshape2, tempshape2 = tempshape2,
                          tempskew = tempskew, LHT = LHT, skewness = skew, kurtosis = kurt,Pskewness = pskew))
  return( ans )
}
