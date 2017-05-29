# rolling estimation/forecast
#---------------------
#' @importFrom PerformanceAnalytics skewness kurtosis

#' @export acdrollsim
acdrollsim = function(spec, data, horizon = 22,m.sim = 10000, forecast.length = 500,n.start = NULL, burn = 0,
                      refit.every = 25, refit.window = c("recursive", "moving"),
                      window.size = NULL, solver = "msucminf", fit.control = list(), solver.control = list(trace = FALSE),
                      calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
                      keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE,compareGARCH = c("LL","none"), ...)
{
  UseMethod("acdrollsim")
}
.acdrollsim = function(spec, data, horizon = 22,m.sim = 10000, forecast.length = 500, n.start = NULL,burn = 0,
                    refit.every = 25, refit.window = c("recursive", "moving"),
                    window.size = NULL, solver = "msucminf", fit.control = list(),solver.control = list(trace = FALSE),
                    calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05), cluster = NULL,
                    keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE,compareGARCH = c("LL","none"),...)
{
  tic = Sys.time()
  if(is.null(solver.control$trace)) trace = FALSE else trace = solver.control$trace
  if(is.null(solver.control$restarts)) solver.control$restarts = 3
  if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
  if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
  if(is.null(fit.control$scale)) fit.control$scale = FALSE
  if(is.null(fit.control$n.sim)) fit.control$n.sim = 5000
  compareGARCH = compareGARCH[1]
  mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim"))
  if(any(is.na(mm))){
    idx = which(is.na(mm))
    enx = NULL
    for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
    warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
  }
  xdata = .extractdata(data)
  data = xdata$data
  index = xdata$index
  period = xdata$period
  T = NROW(data)
  model = spec@model
    #If we want to expand the acdrolling forecast, we should change from here
  if(is.null(n.start)){
    if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
    n.start = T - forecast.length
  } else{
    forecast.length = T - n.start
  }
  if(T<=n.start) stop("\nacdrollsim:--> start cannot be greater than length of data")
    # the ending points of the estimation window - Set up the data for each estimation and rolling forecast
  s = seq(n.start+refit.every, T, by = refit.every)
  m = length(s)
  # the rolling forecast length
  out.sample = rep(refit.every, m)
  # adjustment to include all the datapoints from the end
  if(s[m]<T){
    s = c(s,T)
    m = length(s)
    out.sample = c(out.sample, s[m]-s[m-1])
  }
  if(refit.window == "recursive"){
    rollind = lapply(1:m, FUN = function(i) 1:s[i])
  } else{
    if(!is.null(window.size)){
      if(window.size<100) stop("\nacdroll:--> window size must be greater than 100.")
      rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])+1):s[i])
    } else{
      rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
    }
  }
  gspec = .spec2GARCH(spec)
  if( !is.null(cluster) ){
    parallel::clusterEvalQ(cl = cluster, library(SgtAcd))
    parallel::clusterExport(cluster, c("data", "index", "s","refit.every","trace","acdlikelihood","calculate.VaR","VaR.alpha",
                                       "keep.coef",  "gspec", "fixARMA","horizon","m.sim","burn",
                                       "fixGARCH", "rollind", "spec", "out.sample", "solver","acdconvergence",
                                       "solver.control", "fit.control"), envir = environment())
    tmp = parallel::parLapplyLB(cl = cluster, 1:m, fun = function(i){
      print(paste("Now estimating window:",i,sep = " "))
      zspec = spec
      xspec = gspec
      if(as.logical(trace)) print("Start the GARCH fitting procedure")
      if(sum(spec@model$modelinc[c(14,19,24)])==0){
        fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control), silent=TRUE)
      }else{
        gfit = acdfit(xspec,  zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                      solver = "msucminf",solver.control = list(trace = FALSE))
        if(acdconvergence(gfit)==0){
          if(fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          } else if(fixARMA && !fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
          } else if(!fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
          } else{
            zspec <- setstartacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          }
          if(xspec@model$modelinc[11]>0) skew0 = coefacd(gfit)["skew"] else skew0 = NULL
          if(xspec@model$modelinc[12]>0) shape10 = coefacd(gfit)["shape1"] else shape10 = NULL
          if(xspec@model$modelinc[13]>0) shape20 = coefacd(gfit)["shape2"] else shape20 = NULL
          glik = unname(acdlikelihood(gfit)[1])
        } else{
          shape0 = NULL
          skew0 = NULL
          glik = NA
        }
        solver.control$trace = FALSE
        fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
      }
      if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
        if(sum(spec@model$modelinc[c(14,19,24)])==0) {
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, NA))}else{
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
        }
      } else{
        # compare GARCH likelihood with ACD model and reject if lik less than
        clik = acdlikelihood(fit)[1]
        if(sum(spec@model$modelinc[c(14,19,24)])==0){glik = 0}
        if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
        } else{
          fspec <- getspecacd(fit)
          fspec <- setfixedacd(fspec,as.list(coefacd(fit)))
          fspec <- setboundsacd(fspec,list(shape1 = fit@model$sbounds[3:4],shape2 = fit@model$sbounds[5:6], skew = fit@model$sbounds[1:2]))
          n.old = fit@model$modeldata$T
          data = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
          flt <- acdfilter(fspec,data,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape[1],shape20 = fit@fit$tshape2[1])
          sigmafilter 	= flt@filter$sigma
          resfilter 		= flt@filter$residuals
          zfilter 		= flt@filter$z
          tskewfilter 	= flt@filter$tskew
          tshape1filter 	= flt@filter$tshape1
          tshape2filter   =flt@filter$tshape2
          tempskewfilter 	= flt@filter$tempskew
          tempshape1filter = flt@filter$tempshape1
          tempshape2filter = flt@filter$tempshape2
          mx = fspec@model$maxOrder
          sig = matrix(NA,ncol = 1, nrow = out.sample[i])
          ret = matrix(NA,ncol = 1, nrow = out.sample[i])
          skewness = matrix(NA,ncol = 1, nrow = out.sample[i])
          kurtosis = matrix(NA,ncol = 1, nrow = out.sample[i])
          if(calculate.VaR) VaR.matrix = matrix(NA,ncol = length(VaR.alpha),nrow = out.sample[i])
          if(as.logical(trace)) print("Start the simulation procedures")
          for(ii in 0:(out.sample[i]-1)){
            presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
            preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
            preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
            preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
            prereturns = tail(data[1:(n.old+ii)],         mx)
            tempPath <- acdpath(fspec,n.sim = horizon,m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                preresiduals = NA, rseed = NA)
            return = as.numeric(colSums(tempPath@path$seriesSim))
            sig[ii+1,] = sqrt(mean(colSums(tempPath@path$sigmaSim^2)))
            ret[ii+1,] = mean(return)
            Lskewness = PerformanceAnalytics::skewness(return)
            Lkurtosis = PerformanceAnalytics::kurtosis(return,method = "excess")
            skewness[ii+1,] = Lskewness
            kurtosis[ii+1,] = Lkurtosis
            if(calculate.VaR) VaR.matrix[ii+1,] = quantile(return,VaR.alpha)
            rm(list = c("tempPath","return"))
          }
          if(calculate.VaR){
            y = as.data.frame(cbind(ret, sig, skewness, kurtosis, VaR.matrix))
          } else {
            y = as.data.frame(cbind(ret, sig, skewness, kurtosis))
          }
          rownames(y) = tail(as.character(index[rollind[[i]]]),out.sample[i])
          if(calculate.VaR) {
            colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis",as.character(VaR.alpha))
          } else{
            colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis")
          }
          if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
          ans = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
        }
      }
      return(ans)})
  } else{
    tmp = vector(mode = "list", length = m)
    for(i in 1:m){
      if(as.logical(trace)) print(paste("Now estimating window:",i,sep = " "))
      zspec = spec
      xspec = gspec
      if(as.logical(trace)) print("Start the GARCH fitting procedure")
      if(sum(spec@model$modelinc[c(14,19,24)])==0){
        fit = try(acdfit(spec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control), silent=TRUE)
      }else{
        gfit = acdfit(xspec,  zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                      solver = "msucminf",solver.control = list(trace = FALSE))
        if(acdconvergence(gfit)==0){
          if(fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          } else if(fixARMA && !fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
          } else if(!fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
          } else{
            zspec <- setstartacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          }
          if(xspec@model$modelinc[11]>0) skew0 = coefacd(gfit)["skew"] else skew0 = NULL
          if(xspec@model$modelinc[12]>0) shape10 = coefacd(gfit)["shape1"] else shape10 = NULL
          if(xspec@model$modelinc[13]>0) shape20 = coefacd(gfit)["shape2"] else shape20 = NULL
          glik = unname(acdlikelihood(gfit)[1])
        } else{
          shape0 = NULL
          skew0 = NULL
          glik = NA
        }
        solver.control$trace = FALSE
        fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
      }
      if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
        if(sum(spec@model$modelinc[c(14,19,24)])==0) {
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, NA))}else{
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
          }
      } else{
        # compare GARCH likelihood with ACD model and reject if lik less than
        clik = acdlikelihood(fit)[1]
        if(sum(spec@model$modelinc[c(14,19,24)])==0){glik = 0}
        if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
        } else{
          fspec <- getspecacd(fit)
          fspec <- setfixedacd(fspec,as.list(coefacd(fit)))
          n.old = fit@model$modeldata$T
          ffilter <- acdfilter(fspec,zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index),n.old = n.old)
          sig = matrix(NA,ncol = 1, nrow = out.sample[i])
          ret = matrix(NA,ncol = 1, nrow = out.sample[i])
          skewness = matrix(NA,ncol = 1, nrow = out.sample[i])
          kurtosis = matrix(NA,ncol = 1, nrow = out.sample[i])
          if(calculate.VaR) VaR.matrix = matrix(NA,ncol = length(VaR.alpha),nrow = out.sample[i])
          if(as.logical(trace)) print("Start the simulation procedures")
          for(ii in 0:(out.sample[i]-1)){
            tempPath <- acdpath(fspec,n.sim = horizon,m.sim = m.sim,n.start = burn, presigma = sigmaAcd(ffilter)[n.old+ii],
                                prereturns = ffilter@model$modeldata$data[n.old +ii],preresiduals = residualsAcd(ffilter)[n.old+ii],
                                preskew = skew(ffilter)[n.old+ii],preshape1 = shape1(ffilter)[n.old+ii],preshape2 = shape2(ffilter)[n.old +ii])
            sigma = as.numeric(sqrt(colSums(tempPath@path$sigmaSim^2)))
            return = as.numeric(colSums(tempPath@path$seriesSim))
            sig[ii+1,] = mean(sigma)
            ret[ii+1,] = mean(return)
            Lskewness = PerformanceAnalytics::skewness(return)
            Lkurtosis = PerformanceAnalytics::kurtosis(return,method = "excess")
            skewness[ii+1,] = Lskewness
            kurtosis[ii+1,] = Lkurtosis
            if(calculate.VaR) VaR.matrix[ii+1,] = quantile(return,VaR.alpha)
            rm(list = c("tempPath","sigma","return"))
          }
          if(calculate.VaR){
            y = as.data.frame(cbind(ret, sig, skewness, kurtosis, VaR.matrix))
          } else {
            y = as.data.frame(cbind(ret, sig, skewness, kurtosis))
          }
          rownames(y) = tail(as.character(index[rollind[[i]]]),out.sample[i])
          if(calculate.VaR) {
            colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis",as.character(VaR.alpha))
          } else{
            colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis")
          }
          if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
          tmp[[i]] = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
        }
      }
    }
  }
  conv = sapply(tmp, FUN = function(x) x$converge)
  if(any(!conv)){
    warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
    noncidx = which(!conv)
    model = list()
    model$fixARMA = fixARMA
    model$fixGARCH = fixGARCH
    model$compareGARCH = compareGARCH
    model$spec = spec
    model$data = data
    model$index = index
    model$period = period
    model$horizon = horizon
    model$m.sim = m.sim
    model$burn = burn
    model$forecast.length = forecast.length
    model$n.start = n.start
    model$n.refits = m
    model$refit.every = refit.every
    model$refit.window = refit.window
    model$window.size = window.size
    model$calculate.VaR = calculate.VaR
    model$VaR.alpha = VaR.alpha
    model$keep.coef = keep.coef
    model$noncidx = noncidx
    model$rollind = rollind
    model$out.sample = out.sample
    forecast = tmp
    toc = Sys.time()-tic
    model$elapsed = toc
    ans = new("ACDroll",
              model = model,
              forecast = forecast)
    return(ans)
  } else{
    noncidx = NULL
    forc = tmp[[1]]$y
    if(m>1){
      for(i in 2:m){
        forc = rbind(forc, tmp[[i]]$y)
      }
    }
    if(keep.coef){
      cf = vector(mode = "list", length = m)
      for(i in 1:m){
        cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
        cf[[i]]$coef = tmp[[i]]$cf
      }
    } else{
      cf = NULL
    }
    LL = vector(mode = "list", length = m)
    for(i in 1:m){
      LL[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
      LL[[i]]$log.likelihood = tmp[[i]]$lik
    }
    if(calculate.VaR){
      VaR.matrix = forc[,5:NCOL(forc)]
    } else{
      VaR.matrix = NULL
    }
    model = list()
    model$spec = spec
    model$data = data
    model$index = index
    model$period = period
    model$n.ahead = horizon
    model$m.sim = m.sim
    model$burn = burn
    model$forecast.length = forecast.length
    model$n.start = n.start
    model$refit.every = refit.every
    model$n.refits = m
    model$refit.window = refit.window
    model$window.size = window.size
    model$calculate.VaR = calculate.VaR
    model$VaR.alpha = VaR.alpha
    model$keep.coef = keep.coef
    model$noncidx = noncidx
    model$coef = cf
    model$LL = LL
    model$rollind = rollind
    model$out.sample = out.sample
    forecast = list(VaR = VaR.matrix, density = forc)
  }
  toc = Sys.time()-tic
  model$elapsed = toc
  model$fixARMA = fixARMA
  model$fixGARCH = fixGARCH
  model$compareGARCH = compareGARCH
  ans = new("ACDroll",
            model = model,
            forecast = forecast)
  return( ans )
}

setMethod("acdrollsim", signature(spec = "ACDspec"), .acdrollsim)
#---------------------------------------------------
# Resume method
#----------------------------------------------------
#' @export acdresumeSim
acdresumeSim = function(object, spec = NULL, solver = "mssolnp", fit.control = list(),
                     solver.control = list(), cluster = NULL, fixARMA = NULL, fixGARCH = NULL,
                     compareGARCH = NULL)
{
  UseMethod("acdresumeSim")
}
.acdresumeSim = function(object, spec = NULL, solver = "mssolnp", fit.control = list(),
                         solver.control = list(), cluster = NULL, fixARMA = NULL, fixGARCH = NULL,
                         compareGARCH = NULL)
{
  if(!is.null(object@model$noncidx)){
    noncidx = object@model$noncidx
    tic = Sys.time()
    if(is.null(solver.control$trace)) trace = FALSE else trace = solver.control$trace
    if(is.null(solver.control$restarts)) solver.control$restarts = 3
    if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
    if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
    if(is.null(fit.control$scale)) fit.control$scale = FALSE
    if(is.null(fit.control$n.sim)) fit.control$n.sim = 5000
    mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "n.sim"))
    if(any(is.na(mm))){
      idx = which(is.na(mm))
      enx = NULL
      for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
      warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
    }
    model = object@model
    if(is.null(fixARMA))  fixARMA = model$fixARMA
    if(is.null(fixGARCH)) fixGARCH = model$fixGARCH
    if(is.null(compareGARCH)) compareGARCH = model$compareGARCH
    keep.coef = model$keep.coef
    if(is.null(spec)) spec = model$spec
    gspec = .spec2GARCH(spec)
    data = model$data
    index = model$index
    period = model$period
    T = NROW(data)
    modelinc = spec@model$modelinc
    calculate.VaR = model$calculate.VaR
    VaR.alpha = model$VaR.alpha
    n.ahead = model$n.ahead
    n.start = model$n.start
    horizon = model$horizon
    m.sim = model$m.sim
    burn = model$burn
    forecast.length = model$forecast.length
    refit.every = model$refit.every
    refit.window = model$refit.window
    window.size = model$window.size
    #If we want to expand the acdrolling forecast, we should change from here
    if(is.null(n.start)){
      if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
      n.start = T - forecast.length
    } else{
      forecast.length = T - n.start
    }
    if(T<=n.start) stop("\nacdrollsim:--> start cannot be greater than length of data")
    # the ending points of the estimation window - Set up the data for each estimation and rolling forecast
    s = seq(n.start+refit.every, T, by = refit.every)
    m = length(s)
    # the rolling forecast length
    out.sample = rep(refit.every, m)
    # adjustment to include all the datapoints from the end
    if(s[m]<T){
      s = c(s,T)
      m = length(s)
      out.sample = c(out.sample, s[m]-s[m-1])
    }
    if(refit.window == "recursive"){
      rollind = lapply(1:m, FUN = function(i) 1:s[i])
    } else{
      if(!is.null(window.size)){
        if(window.size<100) stop("\nacdroll:--> window size must be greater than 100.")
        rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])+1):s[i])
      } else{
        rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
      }
    }
    if( !is.null(cluster)){
      parallel::clusterEvalQ(cl = cluster, library(SgtAcd))
      parallel::clusterExport(cluster, c("data", "index", "s","refit.every","trace","acdlikelihood","calculate.VaR","VaR.alpha",
                                         "keep.coef",  "gspec", "fixARMA","horizon","m.sim","burn",
                                         "fixGARCH", "rollind", "spec", "out.sample", "solver","acdconvergence",
                                         "solver.control", "fit.control"), envir = environment())
      tmp = parallel::parLapplyLB(cl = cluster, as.list(noncidx), fun = function(i){
        print(paste("Now estimating window:",i,sep = " "))
        zspec = spec
        xspec = gspec
        if(as.logical(trace)) print("Start the GARCH fitting procedure")
        if(sum(spec@model$modelinc[c(14,19,24)])==0){
          fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                           solver = solver, solver.control = solver.control,
                           fit.control = fit.control), silent=TRUE)
        }else{
          gfit = acdfit(xspec,  zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                        solver = "msucminf",solver.control = list(trace = FALSE))
          if(acdconvergence(gfit)==0){
            if(fixARMA && fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
            } else if(fixARMA && !fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
            } else if(!fixARMA && fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
            } else{
              zspec <- setstartacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
            }
            if(xspec@model$modelinc[11]>0) skew0 = coefacd(gfit)["skew"] else skew0 = NULL
            if(xspec@model$modelinc[12]>0) shape10 = coefacd(gfit)["shape1"] else shape10 = NULL
            if(xspec@model$modelinc[13]>0) shape20 = coefacd(gfit)["shape2"] else shape20 = NULL
            glik = unname(acdlikelihood(gfit)[1])
          } else{
            shape0 = NULL
            skew0 = NULL
            glik = NA
          }
          solver.control$trace = FALSE
          fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                           solver = solver, solver.control = solver.control,
                           fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
        }
        if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
          if(sum(spec@model$modelinc[c(14,19,24)])==0) {
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, NA))}else{
              ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
            }
        } else{
          # compare GARCH likelihood with ACD model and reject if lik less than
          clik = acdlikelihood(fit)[1]
          if(sum(spec@model$modelinc[c(14,19,24)])==0){glik = 0}
          if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
          } else{
            fspec <- getspecacd(fit)
            fspec <- setfixedacd(fspec,as.list(coefacd(fit)))
            fspec <- setboundsacd(fspec,list(shape1 = fit@model$sbounds[3:4],shape2 = fit@model$sbounds[5:6], skew = fit@model$sbounds[1:2]))
            n.old = fit@model$modeldata$T
            data = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
            flt <- acdfilter(fspec,data,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape[1],shape20 = fit@fit$tshape2[1])
            sigmafilter 	= flt@filter$sigma
            resfilter 		= flt@filter$residuals
            zfilter 		= flt@filter$z
            tskewfilter 	= flt@filter$tskew
            tshape1filter 	= flt@filter$tshape1
            tshape2filter   =flt@filter$tshape2
            tempskewfilter 	= flt@filter$tempskew
            tempshape1filter = flt@filter$tempshape1
            tempshape2filter = flt@filter$tempshape2
            mx = fspec@model$maxOrder
            sig = matrix(NA,ncol = 1, nrow = out.sample[i])
            ret = matrix(NA,ncol = 1, nrow = out.sample[i])
            skewness = matrix(NA,ncol = 1, nrow = out.sample[i])
            kurtosis = matrix(NA,ncol = 1, nrow = out.sample[i])
            if(calculate.VaR) VaR.matrix = matrix(NA,ncol = length(VaR.alpha),nrow = out.sample[i])
            if(as.logical(trace)) print("Start the simulation procedures")
            for(ii in 0:(out.sample[i]-1)){
              presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
              preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
              preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
              preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
              prereturns = tail(data[1:(n.old+ii)],         mx)
              tempPath <- acdpath(fspec,n.sim = horizon,m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                  preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                  preresiduals = NA, rseed = NA)
              return = as.numeric(colSums(tempPath@path$seriesSim))
              sig[ii+1,] = sqrt(mean(colSums(tempPath@path$sigmaSim^2)))
              ret[ii+1,] = mean(return)
              Lskewness = PerformanceAnalytics::skewness(return)
              Lkurtosis = PerformanceAnalytics::kurtosis(return,method = "excess")
              skewness[ii+1,] = Lskewness
              kurtosis[ii+1,] = Lkurtosis
              if(calculate.VaR) VaR.matrix[ii+1,] = quantile(return,VaR.alpha)
              rm(list = c("tempPath","return"))
            }
            if(calculate.VaR){
              y = as.data.frame(cbind(ret, sig, skewness, kurtosis, VaR.matrix))
            } else {
              y = as.data.frame(cbind(ret, sig, skewness, kurtosis))
            }
            rownames(y) = tail(as.character(index[rollind[[i]]]),out.sample[i])
            if(calculate.VaR) {
              colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis",as.character(VaR.alpha))
            } else{
              colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis")
            }
            if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
            ans = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
          }
        }
        return(ans)})
    } else{
      tmp = lapply(as.list(noncidx), FUN = function(i){
        zspec = spec
        xspec = gspec
        if(sum(spec@model$modelinc[c(14,19,24)])==0){
          fit = try(acdfit(spec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                           solver = solver, solver.control = solver.control,
                           fit.control = fit.control), silent=TRUE)
        }else{
          gfit = acdfit(xspec,  zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                        solver = "msucminf",solver.control = list(trace = FALSE))
          if(acdconvergence(gfit)==0){
            if(fixARMA && fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
            } else if(fixARMA && !fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
            } else if(!fixARMA && fixGARCH){
              zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
            } else{
              zspec <- setstartacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
            }
            if(xspec@model$modelinc[11]>0) skew0 = coefacd(gfit)["skew"] else skew0 = NULL
            if(xspec@model$modelinc[12]>0) shape10 = coefacd(gfit)["shape1"] else shape10 = NULL
            if(xspec@model$modelinc[13]>0) shape20 = coefacd(gfit)["shape2"] else shape20 = NULL
            glik = unname(acdlikelihood(gfit)[1])
          } else{
            shape0 = NULL
            skew0 = NULL
            glik = NA
          }
          solver.control$trace = FALSE
          fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                           solver = solver, solver.control = solver.control,
                           fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
        }
        if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
          if(sum(spec@model$modelinc[c(14,19,24)])==0) {
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, NA))}else{
              ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
            }
        } else{
          # compare GARCH likelihood with ACD model and reject if lik less than
          clik = acdlikelihood(fit)[1]
          if(sum(spec@model$modelinc[c(14,19,24)])==0){glik = 0}
          if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
          } else{
            fspec <- getspecacd(fit)
            fspec <- setfixedacd(fspec,as.list(coefacd(fit)))
            n.old = fit@model$modeldata$T
            ffilter <- acdfilter(fspec,zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index),n.old = n.old)
            sig = matrix(NA,ncol = 1, nrow = out.sample[i])
            ret = matrix(NA,ncol = 1, nrow = out.sample[i])
            skewness = matrix(NA,ncol = 1, nrow = out.sample[i])
            kurtosis = matrix(NA,ncol = 1, nrow = out.sample[i])
            if(calculate.VaR) VaR.matrix = matrix(NA,ncol = length(VaR.alpha),nrow = out.sample[i])
            if(as.logical(trace)) print("Start the simulation procedures")
            for(ii in 0:(out.sample[i]-1)){
              tempPath <- acdpath(fspec,n.sim = horizon,m.sim = m.sim,n.start = burn, presigma = sigmaAcd(ffilter)[n.old+ii],
                                  prereturns = ffilter@model$modeldata$data[n.old +ii],preresiduals = residualsAcd(ffilter)[n.old+ii],
                                  preskew = skew(ffilter)[n.old+ii],preshape1 = shape1(ffilter)[n.old+ii],preshape2 = shape2(ffilter)[n.old +ii])
              sigma = as.numeric(sqrt(colSums(tempPath@path$sigmaSim^2)))
              return = as.numeric(colSums(tempPath@path$seriesSim))
              sig[ii+1,] = mean(sigma)
              ret[ii+1,] = mean(return)
              Lskewness = PerformanceAnalytics::skewness(return)
              Lkurtosis = PerformanceAnalytics::kurtosis(return,method = "excess")
              skewness[ii+1,] = Lskewness
              kurtosis[ii+1,] = Lkurtosis
              if(calculate.VaR) VaR.matrix[ii+1,] = quantile(return,VaR.alpha)
              rm(list = c("tempPath","sigma","return"))
            }
            if(calculate.VaR){
              y = as.data.frame(cbind(ret, sig, skewness, kurtosis, VaR.matrix))
            } else {
              y = as.data.frame(cbind(ret, sig, skewness, kurtosis))
            }
            rownames(y) = tail(as.character(index[rollind[[i]]]),out.sample[i])
            if(calculate.VaR) {
              colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis",as.character(VaR.alpha))
            } else{
              colnames(y) = c("Mu", "Sigma", "skewness", "Kurtosis")
            }
            if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
            ans = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
          }
        }
      })
    }
    forecast = object@forecast
    conv = sapply(tmp, FUN = function(x) x$converge)
    for(i in 1:length(noncidx)){
      if(conv[i]) forecast[[noncidx[i]]] = tmp[[i]]
    }
    if(any(!conv)){
      warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
      noncidx = which(!conv)
      model = list()
      model$fixARMA = fixARMA
      model$fixGARCH = fixGARCH
      model$compareGARCH = compareGARCH
      model$spec = spec
      model$data = data
      model$index = index
      model$period = period
      model$horizon = horizon
      model$m.sim = m.sim
      model$burn = burn
      model$forecast.length = forecast.length
      model$n.start = n.start
      model$n.refits = m
      model$refit.every = refit.every
      model$refit.window = refit.window
      model$window.size = window.size
      model$calculate.VaR = calculate.VaR
      model$VaR.alpha = VaR.alpha
      model$keep.coef = keep.coef
      model$noncidx = noncidx
      model$rollind = rollind
      model$out.sample = out.sample
      forecast = tmp
      toc = Sys.time()-tic
      model$elapsed = toc
      ans = new("ACDroll",
                model = model,
                forecast = forecast)
      return(ans)
    } else{
      noncidx = NULL
      forc = forecast[[1]]$y
      if(m>1){
        for(i in 2:m){
          forc = rbind(forc, forecast[[i]]$y)
        }
      }
      if(keep.coef){
        cf = vector(mode = "list", length = m)
        for(i in 1:m){
          cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
          cf[[i]]$coef = forecast[[i]]$cf
        }
      } else{
        cf = NULL
      }
      LL = vector(mode = "list", length = m)
      for(i in 1:m){
        LL[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
        LL[[i]]$log.likelihood = forecast[[i]]$lik
      }
      if(calculate.VaR){
        VaR.matrix = forc[,5:NCOL(forc)]
      } else{
        VaR.matrix = NULL
      }
      model = list()
      model$spec = spec
      model$data = data
      model$index = index
      model$period = period
      model$n.ahead = horizon
      model$m.sim = m.sim
      model$burn = burn
      model$forecast.length = forecast.length
      model$n.start = n.start
      model$refit.every = refit.every
      model$n.refits = m
      model$refit.window = refit.window
      model$window.size = window.size
      model$calculate.VaR = calculate.VaR
      model$VaR.alpha = VaR.alpha
      model$keep.coef = keep.coef
      model$noncidx = noncidx
      model$coef = cf
      model$LL = LL
      model$rollind = rollind
      model$out.sample = out.sample
      forecast = list(VaR = VaR.matrix, density = forc)
    }
    toc = Sys.time()-tic
    model$elapsed = toc
    model$fixARMA = fixARMA
    model$fixGARCH = fixGARCH
    model$compareGARCH = compareGARCH
    ans = new("ACDroll",
              model = model,
              forecast = forecast)
  } else{
    # do nothing...all converged
    ans = object
  }
  return( ans )
}
setMethod("acdresumeSim", signature(object = "ACDroll"),  .acdresumeSim)
#-------------------------------------------------------------------------------
