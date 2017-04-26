# rolling estimation/forecast
#---------------------
#' @export acdrollfor
acdrollfor = function(spec, data, n.ahead = 1, forecast.length = 500,
                   n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"),
                   window.size = NULL, solver = "msucminf", fit.control = list(), solver.control = list(),
                   calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
                   keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE,compareGARCH = c("LL","none"), ...)
{
  UseMethod("acdrollfor")
}
.acdrollfor = function(spec, data, n.ahead = 1, forecast.length = 500, n.start = NULL,
                    refit.every = 25, refit.window = c("recursive", "moving"),
                    window.size = NULL, solver = "msucminf", fit.control = list(),
                    solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05), cluster = NULL,
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
  modelinc = spec@model$modelinc
  if(n.ahead>1) stop("\nacdroll:--> n.ahead>1 not supported...try again.")
  #If we want to expand the acdrolling forecast, we should change from here
  if(is.null(n.start)){
    if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
    n.start = T - forecast.length
  } else{
    forecast.length = T - n.start
  }
  if(T<=n.start) stop("\nacdroll:--> start cannot be greater than length of data")
  distribution = spec@model$dmodel$model
  if(distribution=="sgt"){
    InShape1 = TRUE
    InShape2 = TRUE
  }
  if(distribution=="sged"){
    InShape1 = TRUE
    InShape2 = FALSE
  }
  if(distribution=="sst"){
    InShape1 = FALSE
    InShape2 = TRUE
  }
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
      rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
    } else{
      rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
    }
  }
  gspec = .spec2GARCH(spec)
  if( !is.null(cluster) ){
    parallel::clusterEvalQ(cl = cluster, library(SgtAcd))
    parallel::clusterExport(cluster, c("data", "index", "s","refit.every","trace",
                                       "keep.coef",  "gspec", "fixARMA",
                                       "fixGARCH", "rollind", "spec", "out.sample", "solver",
                                       "solver.control", "fit.control"), envir = environment())
    tmp = parallel::parLapplyLB(cl = cluster, 1:m, fun = function(i){
      if(as.logical(trace)) print(paste("Now estimating window:",i,sep = " "))
      zspec = spec
      xspec = gspec
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
          zspec <- setstartacd(spec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:15])]))
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
      fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                       solver = solver, solver.control = solver.control,
                       fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
      if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
        ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
      } else{
        # compare GARCH likelihood with ACD model and reject if lik less than
        clik = acdlikelihood(fit)[1]
        if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
        } else{
          f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1)
          sig = as.numeric(sigmaAcd(f))
          ret = as.numeric(fittedAcd(f))
          if(InShape1) sh1p = as.numeric(shape1(f))
          if(InShape2) sh2p = as.numeric(shape2(f))
          skw = as.numeric(skew(f))
          rlz = zoo::coredata(tail(data[rollind[[i]]], out.sample[i]))
          # use xts for indexing the forecasts
          if(InShape1&&InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, sh2p, rlz))
          if(!InShape1) y = as.data.frame(cbind(ret, sig, skw, sh2p, rlz))
          if(!InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, rlz))
          rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
          if(InShape1&&InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Shape2", "Realized")
          if(!InShape1) colnames(y) = c("Mu", "Sigma", "Skew", "Shape2", "Realized")
          if(!InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Realized")
          if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
          ans = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
        }
      }
      return(ans)})
  } else{
    tmp = vector(mode = "list", length = m)
    for(i in 1:m){
      if(trace) paste("Now estimating window:",i,sep = " ")
      zspec = spec
      xspec = gspec
      gfit = acdfit(xspec,  zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]),order.by = , out.sample = out.sample[i],
                    solver = "msucminf",solver.control = list(trace = FALSE))
      if(acdconvergence(gfit)==0){
        if(fixARMA && fixGARCH){
          zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
        } else if(fixARMA && !fixGARCH){
          zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
        } else if(!fixARMA && fixGARCH){
          zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
        } else{
          zspec <- setstartacd(spec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:15])]))
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
      fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                       solver = solver, solver.control = solver.control,
                       fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
      if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
        ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
      } else{
        # compare GARCH likelihood with ACD model and reject if lik less than
        clik = acdlikelihood(fit)[1]
        if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
        } else{
          f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1)
          sig = as.numeric(sigmaAcd(f))
          ret = as.numeric(fittedAcd(f))
          if(InShape1) sh1p = as.numeric(shape1(f))
          if(InShape2) sh2p = as.numeric(shape2(f))
          skw = as.numeric(skew(f))
          rlz = zoo::coredata(tail(data[rollind[[i]]], out.sample[i]))
          # use xts for indexing the forecasts
          if(InShape1&&InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, sh2p, rlz))
          if(!InShape1) y = as.data.frame(cbind(ret, sig, skw, sh2p, rlz))
          if(!InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, rlz))
          rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
          if(InShape1&&InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Shape2", "Realized")
          if(!InShape1) colnames(y) = c("Mu", "Sigma", "Skew", "Shape2", "Realized")
          if(!InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Realized")
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
    model$n.ahead = n.ahead
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
      if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
      VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
      for(i in 1:length(VaR.alpha)){
        if(distribution == "sgt"){
          VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                 lambda = forc[,3], p = forc[,4], q = forc[,5]/forc[,4])
        } else if(distribution == "sged"){
          VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                lambda = forc[,3], p = forc[,4])
        } else{
          VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                lambda = forc[,3], p = 2, q = forc[,4]/2)
        }
      }
      if(distribution == "sgt"){
        VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
      } else{
        VaR.matrix[,length(VaR.alpha)+1] = forc[,5]
      }
      colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
      VaR.matrix = as.data.frame(VaR.matrix)
      rownames(VaR.matrix) = rownames(forc)
    } else{
      VaR.matrix = NULL
    }
    model = list()
    model$spec = spec
    model$data = data
    model$index = index
    model$period = period
    model$n.ahead = n.ahead
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

setMethod("acdrollfor", signature(spec = "ACDspec"),  .acdrollfor)
#---------------------------------------------------
# Resume method
#----------------------------------------------------
#' @export acdresumeFor
acdresumeFor = function(object, spec = NULL, solver = "mssolnp", fit.control = list(),
                  solver.control = list(), cluster = NULL, fixARMA = NULL, fixGARCH = NULL,
                  compareGARCH = NULL)
{
  UseMethod("acdresumeFor")
}
.acdresumeFor = function(object, spec = NULL, solver = "mssolnp", fit.control = list(),
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
    compareGARCH = compareGARCH[1]
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
    datanames = model$datanames
    data = model$data
    index = model$index
    period = model$period
    T = NROW(data)
    modelinc = spec@model$modelinc
    calculate.VaR = model$calculate.VaR
    VaR.alpha = model$VaR.alpha
    n.ahead = model$n.ahead
    n.start = model$n.start
    forecast.length = model$forecast.length
    refit.every = model$refit.every
    refit.window = model$refit.window
    window.size = model$window.size
    if(n.ahead>1) stop("\nacdroll:--> n.ahead>1 not supported...try again.")
    #If we want to expand the acdrolling forecast, we should change from here
    if(is.null(n.start)){
      if(is.null(forecast.length)) stop("\nacdroll:--> forecast.length amd n.start are both NULL....try again.")
      n.start = T - forecast.length
    } else{
      forecast.length = T - n.start
    }
    if(T<=n.start) stop("\nacdroll:--> start cannot be greater than length of data")
    distribution = spec@model$dmodel$model
    if(distribution=="sgt"){
      InShape1 = TRUE
      InShape2 = TRUE
    }
    if(distribution=="sged"){
      InShape1 = TRUE
      InShape2 = FALSE
    }
    if(distribution=="sst"){
      InShape1 = FALSE
      InShape2 = TRUE
    }
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
        rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
      } else{
        rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
      }
    }
    gspec = .spec2GARCH(spec)
    if( !is.null(cluster) ){
      parallel::clusterEvalQ(cl = cluster, library(SgtAcd))
      parallel::clusterExport(cluster, c("data", "index", "s","refit.every","trace",
                                         "keep.coef", "gspec", "fixARMA",
                                         "fixGARCH", "rollind", "spec", "out.sample", "solver",
                                         "solver.control", "fit.control"), envir = environment())
      tmp = parallel::parLapplyLB(cl = cluster, 1:m, fun = function(i){
        zspec = spec
        xspec = gspec
        gfit = acdfit(xspec, data[rollind[[i]]], out.sample = out.sample[i],
                      solver = "msucminf",solver.control = list(trace = FALSE))
        if(acdconvergence(gfit)==0){
          if(fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          } else if(fixARMA && !fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
          } else if(!fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
          } else{
            zspec <- setstartacd(spec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:15])]))
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
        fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
        if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
        } else{
          # compare GARCH likelihood with ACD model and reject if lik less than
          clik = acdlikelihood(fit)[1]
          if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
          } else{
            f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1)
            sig = as.numeric(sigmaAcd(f))
            ret = as.numeric(fittedAcd(f))
            if(InShape1) sh1p = as.numeric(shape1(f))
            if(InShape2) sh2p = as.numeric(shape2(f))
            skw = as.numeric(skew(f))
            rlz = zoo::coredata(tail(data[rollind[[i]]], out.sample[i]))
            # use xts for indexing the forecasts
            if(InShape1&&InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, sh2p, rlz))
            if(!InShape1) y = as.data.frame(cbind(ret, sig, skw, sh2p, rlz))
            if(!InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, rlz))
            rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
            if(InShape1&&InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Shape2", "Realized")
            if(!InShape1) colnames(y) = c("Mu", "Sigma", "Skew", "Shape2", "Realized")
            if(!InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Realized")
            if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
            ans = list(y = y, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
          }
        }
        return(ans)})
    } else{
      tmp = vector(mode = "list", length = m)
      for(i in 1:m){
        zspec = spec
        xspec = gspec
        gfit = acdfit(xspec, data[rollind[[i]]], out.sample = out.sample[i],
                      solver = "msucminf",solver.control = list(trace = FALSE))
        if(acdconvergence(gfit)==0){
          if(fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:10])]))
          } else if(fixARMA && !fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:6])]))
          } else if(!fixARMA && fixGARCH){
            zspec <- setfixedacd(zspec,as.list(coefacd(gfit)[(sum(gspec@model$modelinc[1:6])+1):sum(gspec@model$modelinc[7:10])]))
          } else{
            zspec <- setstartacd(spec,as.list(coefacd(gfit)[1:sum(gspec@model$modelinc[1:15])]))
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
        fit = try(acdfit(zspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control, shape10 = shape10, shape20 = shape20,skew0 = skew0), silent=TRUE)
        if(inherits(fit, 'try-error') || acdconvergence(fit)!=0 || is.null(fit@fit$cvar)){
          ans = list(y = NA, cf = NA, converge = FALSE, lik = c(NA, glik))
        } else{
          # compare GARCH likelihood with ACD model and reject if lik less than
          clik = acdlikelihood(fit)[1]
          if(!is.na(glik) && clik[1]<glik[1] && compareGARCH=="LL"){
            ans = list(y = NA, cf = NA, converge = FALSE, lik = c(clik, glik))
          } else{
            f = acdforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1)
            sig = as.numeric(sigmaAcd(f))
            ret = as.numeric(fittedAcd(f))
            if(InShape1) sh1p = as.numeric(shape1(f))
            if(InShape2) sh2p = as.numeric(shape2(f))
            skw = as.numeric(skew(f))
            rlz = zoo::coredata(tail(data[rollind[[i]]], out.sample[i]))
            # use xts for indexing the forecasts
            if(InShape1&&InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, sh2p, rlz))
            if(!InShape1) y = as.data.frame(cbind(ret, sig, skw, sh2p, rlz))
            if(!InShape2) y = as.data.frame(cbind(ret, sig, skw, sh1p, rlz))
            rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
            if(InShape1&&InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Shape2", "Realized")
            if(!InShape1) colnames(y) = c("Mu", "Sigma", "Skew", "Shape2", "Realized")
            if(!InShape2) colnames(y) = c("Mu", "Sigma", "Skew", "Shape1", "Realized")
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
      model$n.ahead = n.ahead
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
        if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
        VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
        for(i in 1:length(VaR.alpha)){
          if(distribution == "sgt"){
            VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                  lambda = forc[,3], p = forc[,4], q = forc[,5]/forc[,4])
          } else if(distribution == "sged"){
            VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                  lambda = forc[,3], p = forc[,4])
          } else{
            VaR.matrix[,i] = qsgt(prob = VaR.alpha[i], mu = forc[,1], sigma = forc[,2],
                                  lambda = forc[,3], p = 2, q = forc[,5]/2)
          }
        }
        VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
        colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
        VaR.matrix = as.data.frame(VaR.matrix)
        rownames(VaR.matrix) = rownames(forc)
      } else{
        VaR.matrix = NULL
      }
      model = list()
      model$spec = spec
      model$data = data
      model$index = index
      model$period = period
      model$n.ahead = n.ahead
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
setMethod("acdresumeFor", signature(object = "ACDroll"),  .acdresumeFor)
#-------------------------------------------------------------------------------
.spec2GARCH = function(spec){
  modelinc = spec@model$modelinc
  gspec = acdspec(mean.model=list(include.mean = as.logical(modelinc[1]),
                                     armaOrder = modelinc[2:3], archm = spec@model$mmodel$archm,
                                     skm = spec@model$mmodel$skm,pskm = spec@model$mmodel$pskm,adjm = spec@model$mmodel$adjm),
                     variance.model=list(model = spec@model$vmodel$model, garchOrder = modelinc[8:9],
                                         variance.targeting = spec@model$vmodel$variance.targeting),
                     distribution.model = list(model = spec@model$dmodel$model,skewOrder = NULL,shape1Order = NULL, shape2Order = NULL))
  return(gspec)
}
