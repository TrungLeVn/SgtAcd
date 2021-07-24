# rolling estimation/forecast
#---------------------
#' @importFrom PerformanceAnalytics skewness kurtosis

#' @export acdrollsim
acdrollsim = function(spec, data, horizon = c(22,66),m.sim = 10000, forecast.length = 500,n.start = NULL, burn = 0,
                      refit.every = 22, refit.window = c("recursive", "moving"),
                      window.size = NULL, solver = "msucminf", fit.control = list(), solver.control = list(trace = FALSE),
                    calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
                      keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE,compareGARCH = c("LL","none"), ...)
{
  UseMethod("acdrollsim")
}
.acdrollsim = function(spec, data, horizon = c(22,66),m.sim = 10000, forecast.length = 500, n.start = NULL,burn = 0,
                    refit.every = 25, refit.window = c("recursive", "moving"),
                    window.size = NULL, solver = "msucminf", fit.control = list(),solver.control = list(trace = FALSE),
                    calculate.VaR = TRUE, VaR.alpha = c(0.01,0.05), cluster = NULL,
                    keep.coef = TRUE, fixARMA = TRUE, fixGARCH = TRUE,compareGARCH = c("LL","none"),...)
{
  tic = Sys.time()
  if(is.null(solver.control$trace)) trace = FALSE else trace = solver.control$trace
  if(is.null(solver.control$restarts)) solver.control$restarts = 2
  if(is.null(fit.control$stationarity)) fit.control$stationarity = FALSE
  if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
  if(is.null(fit.control$scale)) fit.control$scale = FALSE
  if(is.null(fit.control$n.sim)) fit.control$n.sim = 1000
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
                                       "solver.control", "fit.control","print"), envir = environment())
    tmp = parallel::parLapply(cl = cluster, 1:m, fun = function(i){
      print(paste("doing window",i,sep = ""))
      zspec = spec
      xspec = gspec
      print("Start the GARCH fitting procedure")
      if(sum(spec@model$modelinc[c(14,19,24)])==0){
        fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control), silent=TRUE)
      }else{
        gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                      solver = solver, solver.control = solver.control,
                      fit.control = fit.control), silent=TRUE)
        if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
          gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                            solver = "mssolnp", solver.control = list(restarts = 10),
                            fit.control = list(n.sim = 10000)), silent=TRUE)
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
          }else if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
            stop("\nacdroll:--> GARCh model can not be converged.\n")
          }
        }else{
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
        }
        print(paste('Now Start the main fit at',i))
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
          dataflt = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
          flt <- acdfilter(fspec,dataflt,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape1[1],shape20 = fit@fit$tshape2[1])
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
          tm = matrix(NA,ncol = length(horizon),nrow = out.sample[i])
          sig = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          ret = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          skewness = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          kurtosis = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          if(calculate.VaR){
            VaR = list()
            for(h in 1:length(horizon)){
              VaR[[as.character(horizon[h])]] = matrix(NA,ncol = length(VaR.alpha)*2,nrow = out.sample[i]) #*2 to compute ES as well
            }
          }
          for(h in 1:length(horizon)){
            # Number of simulation should be equal to out.sample/horizon
            for(ii in seq(0,(out.sample[i]-1),by = horizon[h])){
              presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
              preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
              preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
              preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
              prereturns = tail(dataflt[1:(n.old+ii)],         mx)
              preresiduals = tail(resfilter[1:(n.old+ii)],mx)
              tempPath <- acdpath(fspec,n.sim = horizon[h],m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                  preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                  preresiduals = preresiduals, rseed = NA)
              if(horizon[h] == 1){
                return = as.numeric(tempPath@path$seriesSim[horizon[h],])
                sig[ii+1,h] = sqrt(mean(tempPath@path$sigmaSim[horizon[h],]^2))
              } else{
                return = as.numeric(colSums(tempPath@path$seriesSim[1:horizon[h],]))
                sig[ii+1,h] = sqrt(mean(colSums(tempPath@path$sigmaSim[1:horizon[h],]^2)))
              }
              tm[ii+1,h] = sum(return^3)/m.sim
              ret[ii+1,h] = mean(return)
              skewness[ii+1,h] = PerformanceAnalytics::skewness(return,method = "sample")
              kurtosis[ii+1,h] = PerformanceAnalytics::kurtosis(return,method = "sample_excess")
              if(calculate.VaR){
                VaRest <- quantile(return,VaR.alpha)
                ESest <- VaRest
                for(a in 1:length(VaRest)){
                  ESest[a] <- mean(return[return < VaRest[a]])
                }
                VaR[[as.character(horizon[h])]][ii+1,1:length(VaRest)] = VaRest
                VaR[[as.character(horizon[h])]][ii+1,(length(VaRest)+1):(2*length(VaRest))] =ESest
              }
            }
            rm(list = c("tempPath","return"))
          }
          forecast = list()
          if(calculate.VaR){
            for(h in 1:length(horizon)){
              out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h],matrix(VaR[[h]][seq(1,(out.sample[i]),by = horizon[h]),],ncol = length(VaR.alpha)*2)))
              rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
              colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis",c(paste("VaR",VaR.alpha),paste("ES",VaR.alpha)))
              forecast[[as.character(horizon[h])]] = out
            }
          } else {
            for(h in 1:length(horizon)){
              out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h]))
              rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
              colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis")
              forecast[[as.character(horizon[h])]] = out
            }
          }
          if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
          print(paste("Finish at window",i))
          ans = list(forecast = forecast, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
        }
      }
      return(ans)
      })
  } else{
    tmp = list()
    for(i in 1:m){
      print(paste("doing window",i,sep = ""))
      zspec = spec
      xspec = gspec
      if(as.logical(trace)) print("Start the GARCH fitting procedure")
      if(sum(spec@model$modelinc[c(14,19,24)])==0){
        fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                         solver = solver, solver.control = solver.control,
                         fit.control = fit.control), silent=TRUE)
      }else{
        gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                          solver = solver, solver.control = solver.control,
                          fit.control = fit.control), silent=TRUE)
        if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
          gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                            solver = "mssolnp", solver.control = list(restarts = 10),
                            fit.control = list(n.sim = 10000)), silent=TRUE)
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
          }else if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
            stop("\nacdroll:--> GARCh model can not be converged.\n")
          }
        }else{
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
        }
        print(paste('Now Start the main fit at',i))
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
          dataflt = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
          flt <- acdfilter(fspec,dataflt,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape1[1],shape20 = fit@fit$tshape2[1])
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
          tm = matrix(NA,ncol = length(horizon),nrow = out.sample[i])
          sig = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          ret = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          skewness = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          kurtosis = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
          if(calculate.VaR){
            VaR = list()
            for(h in 1:length(horizon)){
              VaR[[as.character(horizon[h])]] = matrix(NA,ncol = length(VaR.alpha)*2,nrow = out.sample[i]) #*2 to compute ES as well
            }
          }
          for(h in 1:length(horizon)){
            # Number of simulation should be equal to out.sample/horizon
            for(ii in seq(0,(out.sample[i]-1),by = horizon[h])){
              presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
              preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
              preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
              preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
              prereturns = tail(dataflt[1:(n.old+ii)],         mx)
              preresiduals = tail(resfilter[1:(n.old+ii)],mx)
              tempPath <- acdpath(fspec,n.sim = horizon[h],m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                  preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                  preresiduals = preresiduals, rseed = NA)
              if(horizon[h] == 1){
                return = as.numeric(tempPath@path$seriesSim[horizon[h],])
                sig[ii+1,h] = sqrt(mean(tempPath@path$sigmaSim[horizon[h],]^2))
              } else{
                return = as.numeric(colSums(tempPath@path$seriesSim[1:horizon[h],]))
                sig[ii+1,h] = sqrt(mean(colSums(tempPath@path$sigmaSim[1:horizon[h],]^2)))
              }
              tm[ii+1,h] = sum(return^3)/m.sim
              ret[ii+1,h] = mean(return)
              skewness[ii+1,h] = PerformanceAnalytics::skewness(return,method = "sample")
              kurtosis[ii+1,h] = PerformanceAnalytics::kurtosis(return,method = "sample_excess")
              if(calculate.VaR){
                VaRest <- quantile(return,VaR.alpha)
                ESest <- VaRest
                for(a in 1:length(VaRest)){
                  ESest[a] <- mean(return[return < VaRest[a]])
                }
                VaR[[as.character(horizon[h])]][ii+1,1:length(VaRest)] = VaRest
                VaR[[as.character(horizon[h])]][ii+1,(length(VaRest)+1):(2*length(VaRest))] =ESest
              }
            }
            rm(list = c("tempPath","return"))
          }
          forecast = list()
          if(calculate.VaR){
            for(h in 1:length(horizon)){
              out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h],matrix(VaR[[h]][seq(1,(out.sample[i]),by = horizon[h]),],ncol = length(VaR.alpha)*2)))
              rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
              colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis",c(paste("VaR",VaR.alpha),paste("ES",VaR.alpha)))
              forecast[[as.character(horizon[h])]] = out
            }
          } else {
            for(h in 1:length(horizon)){
              out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                        kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h]))
              rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
              colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis")
              forecast[[as.character(horizon[h])]] = out
            }
          }
          if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
          tmp[[i]] = list(forecast = forecast, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
        }
      }
    }
  }
  conv = vector()
  for(i in 1:m){conv[i] = tmp[[i]]$converge}
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
    forecast = list()
    forecast[[1]] = tmp[[1]]$forecast
    if(m>1){
      for(i in 2:m){
        forecast[[i]] = tmp[[i]]$forecast
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
    rollind = model$rollind
    out.sample = model$out.sample
    m = length(noncidx)
    cf = list()
    for(i in 1:length(rollind)){
      cf[[i]]= object@forecast[[i]]$cf
    }
    if( !is.null(cluster)){
      tmp = foreach::foreach(i = 1:m,.packages = "SgtAcd",.export = c("data", "index","refit.every","trace","acdlikelihood",
                                                                      "calculate.VaR","VaR.alpha",
                                                                      "keep.coef",  "gspec", "fixARMA",
                                                                      "horizon","m.sim","burn","fixGARCH", "rollind",
                                                                      "spec", "out.sample", "solver","acdconvergence",
                                                                      "solver.control", "fit.control","print","noncidx","cf"),
                             .combine = list,.multicombine = TRUE)%dopar%{
                               print(paste("doing window",i,sep = ""))
                               zspec = spec
                               xspec = gspec
                               print("Start the GARCH fitting procedure")
                               if(sum(spec@model$modelinc[c(14,19,24)])==0){
                                 fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                                                  solver = solver, solver.control = solver.control,
                                                  fit.control = fit.control), silent=TRUE)
                               }else{
                                 gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                                                   solver = solver, solver.control = solver.control,
                                                   fit.control = fit.control), silent=TRUE)
                                 if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
                                   gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                                                     solver = "mssolnp", solver.control = list(restarts = 10),
                                                     fit.control = list(n.sim = 10000)), silent=TRUE)
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
                                   }else if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
                                     stop("\nacdroll:--> GARCh model can not be converged.\n")
                                   }
                                 }else{
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
                                 }
                                 print(paste('Now Start the main fit at',i))
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
                                   dataflt = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
                                   flt <- acdfilter(fspec,dataflt,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape1[1],shape20 = fit@fit$tshape2[1])
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
                                   tm = matrix(NA,ncol = length(horizon),nrow = out.sample[i])
                                   sig = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
                                   ret = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
                                   skewness = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
                                   kurtosis = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
                                   if(calculate.VaR){
                                     VaR = list()
                                     for(h in 1:length(horizon)){
                                       VaR[[as.character(horizon[h])]] = matrix(NA,ncol = length(VaR.alpha)*2,nrow = out.sample[i]) #*2 to compute ES as well
                                     }
                                   }
                                   for(h in 1:length(horizon)){
                                     # Number of simulation should be equal to out.sample/horizon
                                     for(ii in seq(0,(out.sample[i]-1),by = horizon[h])){
                                       presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
                                       preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
                                       preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
                                       preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
                                       prereturns = tail(dataflt[1:(n.old+ii)],         mx)
                                       preresiduals = tail(resfilter[1:(n.old+ii)],mx)
                                       tempPath <- acdpath(fspec,n.sim = horizon[h],m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                                           preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                                           preresiduals = preresiduals, rseed = NA)
                                       if(horizon[h] == 1){
                                         return = as.numeric(tempPath@path$seriesSim[horizon[h],])
                                         sig[ii+1,h] = sqrt(mean(tempPath@path$sigmaSim[horizon[h],]^2))
                                       } else{
                                         return = as.numeric(colSums(tempPath@path$seriesSim[1:horizon[h],]))
                                         sig[ii+1,h] = sqrt(mean(colSums(tempPath@path$sigmaSim[1:horizon[h],]^2)))
                                       }
                                       tm[ii+1,h] = sum(return^3)/m.sim
                                       ret[ii+1,h] = mean(return)
                                       skewness[ii+1,h] = PerformanceAnalytics::skewness(return,method = "sample")
                                       kurtosis[ii+1,h] = PerformanceAnalytics::kurtosis(return,method = "sample_excess")
                                       if(calculate.VaR){
                                         VaRest <- quantile(return,VaR.alpha)
                                         ESest <- VaRest
                                         for(a in 1:length(VaRest)){
                                           ESest[a] <- mean(return[return < VaRest[a]])
                                         }
                                         VaR[[as.character(horizon[h])]][ii+1,1:length(VaRest)] = VaRest
                                         VaR[[as.character(horizon[h])]][ii+1,(length(VaRest)+1):(2*length(VaRest))] =ESest
                                       }
                                     }
                                     rm(list = c("tempPath","return"))
                                   }
                                   forecast = list()
                                   if(calculate.VaR){
                                     for(h in 1:length(horizon)){
                                       out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                                                 tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                                                 kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h],matrix(VaR[[h]][seq(1,(out.sample[i]),by = horizon[h]),],ncol = length(VaR.alpha)*2)))
                                       rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
                                       colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis",c(paste("VaR",VaR.alpha),paste("ES",VaR.alpha)))
                                       forecast[[as.character(horizon[h])]] = out
                                     }
                                   } else {
                                     for(h in 1:length(horizon)){
                                       out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                                                 tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                                                 kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h]))
                                       rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
                                       colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis")
                                       forecast[[as.character(horizon[h])]] = out
                                     }
                                   }
                                   if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
                                   print(paste("Finish at window",i))
                                   ans = list(forecast = forecast, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
                                 }
                               }
        }
    } else{
      tmp = list()
      for(i in 1:m){
        print(paste("doing window",i,sep = ""))
        zspec = spec
        xspec = gspec
        if(as.logical(trace)) print("Start the GARCH fitting procedure")
        if(sum(spec@model$modelinc[c(14,19,24)])==0){
          fit = try(acdfit(gspec, zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                           solver = solver, solver.control = solver.control,
                           fit.control = fit.control), silent=TRUE)
        }else{
          gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                            solver = solver, solver.control = solver.control,
                            fit.control = fit.control), silent=TRUE)
          if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
            gfit = try(acdfit(xspec,zoo::zoo(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i],
                              solver = "mssolnp", solver.control = list(restarts = 10),
                              fit.control = list(n.sim = 10000)), silent=TRUE)
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
            }else if(inherits(gfit, 'try-error') || acdconvergence(gfit)!=0){
              stop("\nacdroll:--> GARCh model can not be converged.\n")
            }
          }else{
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
          }
          print(paste('Now Start the main fit at',i))
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
            dataflt = zoo::zoo(fit@model$modeldata$data,order.by = fit@model$modeldata$index)
            flt <- acdfilter(fspec,dataflt,n.old = n.old,skew0 = fit@fit$tskew[1], shape10 = fit@fit$tshape1[1],shape20 = fit@fit$tshape2[1])
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
            tm = matrix(NA,ncol = length(horizon),nrow = out.sample[i])
            sig = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
            ret = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
            skewness = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
            kurtosis = matrix(NA,ncol = length(horizon), nrow = out.sample[i])
            if(calculate.VaR){
              VaR = list()
              for(h in 1:length(horizon)){
                VaR[[as.character(horizon[h])]] = matrix(NA,ncol = length(VaR.alpha)*2,nrow = out.sample[i]) #*2 to compute ES as well
              }
            }
            for(h in 1:length(horizon)){
              # Number of simulation should be equal to out.sample/horizon
              for(ii in seq(0,(out.sample[i]-1),by = horizon[h])){
                presig     = tail(sigmafilter[1:(n.old+ii)],  mx)
                preskew    = tail(tskewfilter[1:(n.old+ii)],  mx)
                preshape1   = tail(tshape1filter[1:(n.old+ii)], mx)
                preshape2  = tail(tshape2filter[1:(n.old+ii)], mx)
                prereturns = tail(dataflt[1:(n.old+ii)],         mx)
                preresiduals = tail(resfilter[1:(n.old+ii)],mx)
                tempPath <- acdpath(fspec,n.sim = horizon[h],m.sim = m.sim,n.start = burn,  presigma = presig, preskew = preskew,
                                    preshape1 = preshape1,preshape2 = preshape2, prereturns = prereturns,
                                    preresiduals = preresiduals, rseed = NA)
                if(horizon[h] == 1){
                  return = as.numeric(tempPath@path$seriesSim[horizon[h],])
                  sig[ii+1,h] = sqrt(mean(tempPath@path$sigmaSim[horizon[h],]^2))
                } else{
                  return = as.numeric(colSums(tempPath@path$seriesSim[1:horizon[h],]))
                  sig[ii+1,h] = sqrt(mean(colSums(tempPath@path$sigmaSim[1:horizon[h],]^2)))
                }
                tm[ii+1,h] = sum(return^3)/m.sim
                ret[ii+1,h] = mean(return)
                skewness[ii+1,h] = PerformanceAnalytics::skewness(return,method = "sample")
                kurtosis[ii+1,h] = PerformanceAnalytics::kurtosis(return,method = "sample_excess")
                if(calculate.VaR){
                  VaRest <- quantile(return,VaR.alpha)
                  ESest <- VaRest
                  for(a in 1:length(VaRest)){
                    ESest[a] <- mean(return[return < VaRest[a]])
                  }
                  VaR[[as.character(horizon[h])]][ii+1,1:length(VaRest)] = VaRest
                  VaR[[as.character(horizon[h])]][ii+1,(length(VaRest)+1):(2*length(VaRest))] =ESest
                }
              }
              rm(list = c("tempPath","return"))
            }
            forecast = list()
            if(calculate.VaR){
              for(h in 1:length(horizon)){
                out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                          tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                          kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h],matrix(VaR[[h]][seq(1,(out.sample[i]),by = horizon[h]),],ncol = length(VaR.alpha)*2)))
                rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
                colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis",c(paste("VaR",VaR.alpha),paste("ES",VaR.alpha)))
                forecast[[as.character(horizon[h])]] = out
              }
            } else {
              for(h in 1:length(horizon)){
                out = as.data.frame(cbind(ret[seq(1,(out.sample[i]),by = horizon[h]),h], sig[seq(1,(out.sample[i]),by = horizon[h]),h],
                                          tm[seq(1,(out.sample[i]),by = horizon[h]),h],skewness[seq(1,(out.sample[i]),by = horizon[h]),h],
                                          kurtosis[seq(1,(out.sample[i]),by = horizon[h]),h]))
                rownames(out) = tail(as.character(index[rollind[[i]]]),out.sample[i])[seq(1,(out.sample[i]),by = horizon[h])]
                colnames(out) = c("Mu", "Sigma","TM", "Skewness", "Kurtosis")
                forecast[[as.character(horizon[h])]] = out
              }
            }
            if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
            tmp[[i]] = list(forecast = forecast, cf = cf, converge = TRUE, lik = c(acdlikelihood(fit)[1], glik))
          }
        }
      }
    }
    forecast = object@forecast
    for(i in 1:length(rollind)){
      if(!is.element(i,noncidx)) forecast[[i]] = forecast[[i]]$forecast
    }
    conv = vector()
    for(i in 1:m){
      conv[i] = tmp[[i]]$converge
      if(tmp[[i]]$converge) forecast[[noncidx[i]]] = tmp[[i]]$forecast
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
      model$rollind = rollind
      model$out.sample = out.sample
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
