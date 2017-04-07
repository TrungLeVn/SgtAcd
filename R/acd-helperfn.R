#################################################################################
##
##   R package racd by Alexios Ghalanos Copyright (C) 2012-2014
##   This file is part of the R package racd.
##
##   The R package racd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package racd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
###############################################################

.acdmakefitmodel = function(f, T, m, timer, convergence, message, hess, arglist)
{
	# Turn Stationarity Off for numerical derivative calculation
	fit.control = arglist$fit.control
	fit.control$stationarity = arglist$fit.control$stationarity = 0
	data = arglist$data
	ipars = arglist$ipars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	arglist$returnType = "llh"
	fit = vector(mode = "list")
	if(is.null(hess)){
		#fit$hessian = rugarch:::.hessian2sidedcpp(f, ipars[estidx, 1], arglist = arglist)
		fit$hessian = hessian(f, x = ipars[estidx, 1], arglist = arglist, method.args=list(zero.tol=.Machine$double.eps))
		E = eigen(fit$hessian)$values
		# approx. number of decimal places lost to roundoff/numerical estimation error
		condH = log10(max(E)/min(E))
	} else{
		fit$hessian = hess
		E = eigen(fit$hessian)$values
		condH = log10(max(E)/min(E))
	}
	fit$cvar = try(solve(fit$hessian), silent = TRUE)
	# (might also fail in which case user will see an error about the inversion failure)
	if(inherits(fit$cvar, "try-error")){
		warning("\nracd-->warning: failed to invert hessian\n")
		fit$cvar = NULL
	}
	arglist$returnType = "all"
	temp = f(pars = ipars[estidx, 1],	arglist = arglist)
	fit$var = abs(temp$h)
	fit$sigma = sqrt(abs(temp$h))
	fit$condH = condH
	fit$z = temp$z
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$LHT
	fit$residuals = temp$res
	fit$tskew = temp$tskew
	fit$tshape1 = temp$tshape1
	fit$tshape2 = temp$tshape2
	fit$tempskew = temp$tempskew
	fit$tempshape1 = temp$tempshape1
	fit$tempshape2 = temp$tempshape2
	fit$skewness = temp$skewness
	fit$kurtosis = temp$kurtosis
	fit$Pskew = temp$Pskewness

	if(sum(ipars[,2])>0){
		pall = ipars[estidx | as.logical(ipars[,2]==1), 1]
		fixed = match(rownames(ipars[ipars[,2]==1, , drop = FALSE]), names(pall))
		fixedn = length(fixed)
		fNA = rep(NA, fixedn)
		nfixedn = length(pall) - fixedn
		fit$coef = pall
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA, nfixedn)
			fit$tval = rep(NA, nfixedn)
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(pall))
			fit$matcoef[-fixed,] = cbind(ipars[estidx, 1], fit$se.coef,
			fit$tval, rep(NA, nfixedn))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$robust.se.coef = rep(NA, nfixedn)
			fit$robust.tval = rep(NA, nfixedn)
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef,
					fit$robust.tval, rep(NA, nfixedn))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = "failed to invert hessian"
		} else{
			arglist$returnType = "LHT"
			tmp = rugarch:::robustvcv(fun = f, pars = ipars[estidx, 1], nlag = 0, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(zero.tol=.Machine$double.eps), arglist = arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef[-fixed]/fit$se.coef
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$se.coef,
			fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA,fNA, fNA)
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef[-fixed]/fit$robust.se.coef
			# change here
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef,
			fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = NULL
		}
		if(model$modelinc[7] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))

			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	} else{
		fit$coef = ipars[estidx, 1]
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA,length(fit$coef))
			fit$tval = rep(NA,length(fit$coef))
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
			fit$tval, rep(NA,length(fit$coef)))
			fit$robust.se.coef = rep(NA,length(fit$coef))
			fit$robust.tval = rep(NA,length(fit$coef))
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval,
			rep(NA,length(fit$coef)))
			fit$hessian.message = "failed to invert hessian"
		} else{
			nlag=min(floor(1.2*(T)^(1/3)),(T))
			arglist$returnType = "LHT"
			tmp = rugarch:::robustvcv(fun = f, pars = ipars[estidx,1], nlag = nlag, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(zero.tol=.Machine$double.eps), arglist = arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef/fit$se.coef
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
			fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef/fit$robust.se.coef
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval,2*(1-pnorm(abs(fit$robust.tval))))
			fit$hessian.message = NULL
		}
		# variance targeting case
		if(model$modelinc[7] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))
			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	}
	# return the correct indicators to the ipars (changed in the case of fixed pars and fixed.se = TRUE)
	ipars[,2:4] = model$pars[,2:4]
	dimnames(fit$matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	dimnames(fit$robust.matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	fit$fitted.values = data-fit$residuals
	fit$convergence = convergence
	fit$message = message
	fit$timer = timer
	fit$ipars = ipars
	return(fit)
}
