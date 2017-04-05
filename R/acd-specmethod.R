#----------------------------------------------------------------------------------
## acdspec method
#' @description  Specify models to estimate with mean, variance and higher moment equations
#' @title Model specification
#' @usage function(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1),variance.targeting = FALSE),
#'       mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
#'       distribution.model = list(model = "sged", skewOrder = c(1,0, 1),skewshock = 1, skewshocktype = 1, skewmodel = "quad",
#'                                 shape1Order = c(1, 0, 1), shape1shock = 1, shape1shocktype = 1,shape1model = "quad",
#'                                  shape2Order = c(1, 0, 1), shape2shock = 1, shape2shocktype = 1,shape2model = "quad",exp.rate = 1),
#'       start.pars = list(), fixed.pars = list())
#' @param variance.model Specification of variance process. Only support
#'         sGARCH and gjrGARCH. Do not support external regessor.
#' @param mean.model Specification of mean process. Only support maximum ARMA(1,1)
#' @param distribution.model Specification of skew/shape distribution's
#'            parameter dynamids. Only support "sgt", "sged" and "nig"
#'         as conditional distribution.
#' @param start.pars List of starting values for specific parameters using in optimization routines
#' @param fixed.pars List of parameters that we want to fixed, i.e, drop out of optimization routines
#' @return Return an ACDspec,i.e., specification objects.
#' @export acdspec
#-----------------------------------------------------------------------------------

acdspec <- function(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), variance.targeting = FALSE),
                   mean.model = list(armaOrder = c(0,0), include.mean = TRUE,archm = TRUE, skm = FALSE, kum = FALSE),
                   distribution.model = list(model = "sgt", skewOrder = c(1,0, 1),skewshock = 1, skewshocktype = 1, skewmodel = "quad",volsk = FALSE,
                                             shape1Order = c(1, 0, 1), shape1shock = 1, shape1shocktype = 1,shape1model = "quad",volsh1 = FALSE,
                                             shape2Order = c(1, 0, 1), shape2shock = 1, shape2shocktype = 1,shape2model = "quad",volsh2 = FALSE,exp.rate = 1),
                   start.pars = list(), fixed.pars = list()) {
  UseMethod("acdspec") #acdspec here is the name of the generic function acdspec, that will be used to create an object of ACDspec class
}
## TRUNG: shock types: 1 in z^2 2, in resids^2 3 in abs(z) 4 in abs(resid)
.expand.model <- function(model) {
  modelnames = NULL
  for (i in 1:28) {
    if (model[i] > 0) {
      if (any(c(2,3,8,9,10,15,16,17,20,21,22,25,26,27) == i)) {
        modelnames = c(modelnames, paste(names(model)[i]))
      } else {
        modelnames = c(modelnames, names(model)[i])
      }
    }
  }
  return(modelnames)
}

.acdspec <- function(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1), variance.targeting = FALSE),
                     mean.model = list(armaOrder = c(0,0), include.mean = TRUE,archm = FALSE, skm = FALSE, kum = FALSE),
                     distribution.model = list(model = "sgt", skewOrder = c(1,0, 1 ),skewshock = 1, skewshocktype = 1, skewmodel = "quad", volsk = FALSE,
                                               shape1Order = c(1, 0, 1), shape1shock = 1, shape1shocktype = 1,shape1model = "quad", volsh1 = FALSE,
                                               shape2Order = c(1, 0, 1), shape2shock = 1, shape2shocktype = 1,shape2model = "quad", volsh2 = FALSE, exp.rate = 1),
                     start.pars = list(), fixed.pars = list()) {
  # specify the parameter list
  modelinc = rep(0, 35)
  names(modelinc) = c("mu", "ar", "ma", "archm", "skm", "kum",
                      "omega", "alpha", "beta", "gamma",
                      "skew", "shape1", "shape2",
                      "skcons", "skalpha", "skgamma", "skbeta", "volsk",
                      "sh1cons", "sh1alpha", "sh1gamma", "sh1beta", "volsh1",
                      "sh2cons", "sh2alpha", "sh2gamma", "sh2beta", "volsh2",
                      "skewmodel", "shape1model","shape2model", "skshock", "sh1shock","sh2shock",
                      "aux")

  modeldesc = list()

  #-----
  # Check the input of arguments in acdspec
  #----
  mm = match(names(mean.model), c("armaOrder", "include.mean","archm","skm","kum"))
  if (any(is.na(mm))) {
    idx = which(is.na(mm))
    enx = NULL
    for (i in 1:length(idx)) enx = c(enx, names(mean.model)[idx[i]])
    warning(paste(c("unidentified option(s) in mean.model:\n", enx), sep = "", collapse = " "), call. = FALSE, domain = NULL)
  }
  vm = match(names(variance.model), c("model", "garchOrder", "variance.targeting"))
  if (any(is.na(vm))) {
    idx = which(is.na(vm))
    enx = NULL
    for (i in 1:length(idx)) enx = c(enx, names(variance.model)[idx[i]])
    warning(paste(c("unidentified option(s) in variance.model:\n", enx), sep = "", collapse = " "), call. = FALSE, domain = NULL)
  }
  dm = match(names(distribution.model), c("model",
                                          "skewOrder", "skewshock", "skewshocktype", "skewmodel", "volsk",
                                          "shape1Order", "shape1shock","shape1shocktype", "shape1model", "volsh1",
                                          "shape2Order", "shape2shock", "shape2shocktype", "shape2model", "volsh2",
                                          "exp.rate"))
  if (any(is.na(dm))) {
    idx = which(is.na(dm))
    enx = NULL
    for (i in 1:length(idx)) enx = c(enx, names(distribution.model)[idx[i]])
    warning(paste(c("unidentified option(s) in distribution.model:\n", enx), sep = "", collapse = " "), call. = FALSE, domain = NULL)
  }

  #------
  # Specify the variance equation. The undefined argument will goes with  default
  #-----
  vmodel = list(model = "gjrGARCH", garchOrder = c(1, 1), variance.targeting = FALSE)
  idx = na.omit(match(names(variance.model), names(vmodel)))

  if (length(idx) > 0)
    for (i in 1:length(idx)) vmodel[idx[i]] = variance.model[i]
  valid.model = c("sGARCH", "gjrGARCH")
  if (!any(vmodel$model == valid.model))
    stop("\nacdpec-->error: the garch model does not appear to be a valid choice.\n", call. = FALSE)
  modelinc[8] = vmodel$garchOrder[1]
  modelinc[9] = vmodel$garchOrder[2]
  if (vmodel$model == "gjrGARCH" &&modelinc[8]!=0) {
    modelinc[10] = 1
  }
  if (is.null(vmodel$variance.targeting))
    modelinc[7] = 1 else modelinc[7] = as.integer(1 - vmodel$variance.targeting)
  #---------------------------------------
  # Specify the mean equation. The undefined argument will take default specifications
  #---------------------------------------
  mmodel = list(armaOrder = c(0, 0), include.mean = TRUE, archm = FALSE, skm = FALSE, shm = FALSE)
  idx = na.omit(match(names(mean.model), names(mmodel)))
  if (length(idx) > 0)
    for (i in 1:length(idx)) mmodel[idx[i]] = mean.model[i]
  if(as.logical(mmodel$archm) || as.logical(mmodel$skm) || as.logical(mmodel$shm))
  {
    if(sum(mmodel$armaOrder) > 0){
      stop("\nacdpec-->error: ARMA specification is not allowed if we allow for either archm, skm or shm.\n", call. = FALSE)
    }
    if( as.logical(mmodel$archm) ){
      modelinc[4] = 1
    }
    if( as.logical(mmodel$skm) ){
      modelinc[5] = 1
    }
    if( as.logical(mmodel$shm) ){
      modelinc[6] = 1
    }
  } else{
    modelinc[2] = mmodel$armaOrder[1]
    modelinc[3] = mmodel$armaOrder[2]
  }
  if (is.null(mmodel$include.mean))
    modelinc[1] = 1 else modelinc[1] = as.integer(mmodel$include.mean)
  #--------
  # Specify the higher moment equations. The undefined argument will goes with default
  #-------
  dmodel = list(model = "sgt", skewOrder = c(1, 0, 1), skewshock = 1, skewshocktype = 1, skewmodel = "quad", volsk = FALSE,
                shape1Order = c(1 ,0, 1), shape1shock = 1, shape1shocktype = 1, shape1model = "quad", volsh1 = FALSE,
                shape2Order = c(1 ,0, 1), shape2shock = 1, shape2shocktype = 1, shape2model = "quad", volsh2 = FALSE, exp.rate = 1)
  idx = na.omit(match(names(distribution.model), names(dmodel)))
  if (length(idx) > 0)
    for (i in 1:length(idx)) {
      dmodel[idx[i]] = distribution.model[i]
    }
  valid.distribution = c("sgt","sged","sst")
  if (!any(dmodel$model == valid.distribution))
    stop("\nacdspec-->error: the cond.distribution does not appear to be a valid choice.")
  if(dmodel$model == "sged" && !is.null(dmodel$shape2Order))
    stop("\nacdspec-->error: for the sged distribution, shape2 parameter is Inf.")
  if(dmodel$model == "sst" && !is.null(dmodel$shape1Order))
    stop("\nacdspec-->error: for the sst distribution, shape1 parameter is always 2.")
  if(dmodel$model == "sged"){
    fixed.pars = list(shape2 = 999)
  }
  if(dmodel$model == "sst"){
    fixed.pars = list(shape1 = 2)
  }
  #Skewshocktype ==1 using the squared values, else using the absolute value
  # modelinc = 0: quad with squared values
  # modelinc = 1: quad with abs values
  # modelinc = 2: pwl with squared values
  if (dmodel$skewmodel == "pwl") {
    modelinc[29] = 2
  } else {
    # default (quad)
    if (dmodel$skewshocktype == 1) #Skewshocktype ==1 using the squared values, else using the absolute value
      modelinc[29] = 0 else modelinc[29] = 1
  }
  # Shapeshocktype ==1 using the squared values, else using the absolute value
  # modelinc = 0: quad with squared values
  # modelinc = 1: quad with abs values
  # modelinc = 2: pwl with squared values
  # modelinc = 3: pwl with abs value
  if (dmodel$shape1model == "pwl") {
    if (dmodel$shape1shocktype == 1)
      modelinc[30] = 2 else modelinc[30] = 3
  } else {
    if (dmodel$shape1shocktype == 1)
      modelinc[30] = 0 else modelinc[30] = 1
  }
  if (dmodel$shape2model == "pwl") {
    if (dmodel$shape2shocktype == 1)
      modelinc[31] = 2 else modelinc[31] = 3
  } else {
    if (dmodel$shape2shocktype == 1)
      modelinc[31] = 0 else modelinc[31] = 1
  }
  if (dmodel$skewshock == 1)
    modelinc[32] = 1
  if (dmodel$shape1shock == 1&&dmodel$model != "sst")
    modelinc[33] = 1
  if (dmodel$shape2shock == 1&&dmodel$model != "sged")
    modelinc[34] = 1
  #------
  # Set the distribution bounds and specify the parameters to be estimated (which will have value 1 in modelinc)
  #-----
  # because we exlucde the normal, we add 1 to the value (for c code)
  modeldesc$distno = which(dmodel$model == valid.distribution)
  # 1 is sgt; 2 is sged; 3 is sst
  di = .DistributionBounds(dmodel$model)
  sbounds = rep(0, 7)
  # check if time-varying first
  if (is.null(dmodel$skewOrder)) { # If we do not have time-varying skew, then only calculate the unconditional skew
    modelinc[14:18] = 0
    modelinc[11] = di$include.skew
    sbounds[1] = di$skew.LB
    sbounds[2] = di$skew.UB
  } else { # if we specify time-varying skew parameters
    modelinc[14] = 1
    modelinc[15] = dmodel$skewOrder[1]
    modelinc[16] = dmodel$skewOrder[2]
    modelinc[17] = dmodel$skewOrder[3]
    modelinc[11] = 0
    sbounds[1] = di$skew.LB
    sbounds[2] = di$skew.UB
    if(as.logical(dmodel$volsk)) modelinc[18] = 1
  }
  if(dmodel$model == "sgt"){
    if (is.null(dmodel$shape1Order)) {
      modelinc[19:23] = 0
      modelinc[12] = di$include.shape1
      sbounds[3] = di$shape1.LB
      sbounds[4] = di$shape1.UB
    } else {
      modelinc[19] = 1
      modelinc[20] = dmodel$shape1Order[1]
      modelinc[21] = dmodel$shape1Order[2]
      modelinc[22] = dmodel$shape1Order[3]
      modelinc[12] = 0
      sbounds[3] = di$shape1.LB
      sbounds[4] = di$shape1.UB
      if(as.logical(dmodel$volsh1)) modelinc[23] = 1
    }
    if (is.null(dmodel$shape2Order)) {
      modelinc[24:28] = 0
      modelinc[13] = di$include.shape2
      sbounds[5] = di$shape2.LB
      sbounds[6] = di$shape2.UB
    } else {
      modelinc[24] = 1
      modelinc[25] = dmodel$shape2Order[1]
      modelinc[26] = dmodel$shape2Order[2]
      modelinc[27] = dmodel$shape2Order[3]
      modelinc[13] = 0
      sbounds[5] = di$shape2.LB
      sbounds[6] = di$shape2.UB
      if(as.logical(dmodel$volsh2)) modelinc[28]= 1
    }
  }
  if(dmodel$model == "sged"){
    if (is.null(dmodel$shape1Order)) {
      modelinc[19:23] = 0
      modelinc[12] = di$include.shape1
      sbounds[3] = di$shape1.LB
      sbounds[4] = di$shape1.UB
    } else {
      modelinc[19] = 1
      modelinc[20] = dmodel$shape1Order[1]
      modelinc[21] = dmodel$shape1Order[2]
      modelinc[22] = dmodel$shape1Order[3]
      modelinc[12] = 0
      sbounds[3] = di$shape1.LB
      sbounds[4] = di$shape1.UB
      if(as.logical(dmodel$volsh1)) modelinc[23] = 1
    }
    modelinc[24:28] = 0
    modelinc[13] = 1
    sbounds[5] = di$shape2.LB
    sbounds[6] = di$shape2.UB
  }
  if(dmodel$model == "sst"){
    modelinc[19:22] = 0
    modelinc[12] = 1
    sbounds[3] = di$shape1.LB
    sbounds[4] = di$shape1.UB
    if (is.null(dmodel$shape2Order)) {
      modelinc[24:28] = 0
      modelinc[13] = di$include.shape2
      sbounds[5] = di$shape2.LB
      sbounds[6] = di$shape2.UB
    } else {
      modelinc[24] = 1
      modelinc[25] = dmodel$shape2Order[1]
      modelinc[26] = dmodel$shape2Order[2]
      modelinc[27] = dmodel$shape2Order[3]
      modelinc[13] = 0
      sbounds[5] = di$shape2.LB
      sbounds[6] = di$shape2.UB
      if(as.logical(dmodel$volsh2)) modelinc[28]= 1
    }
  }

  #-------
  # Create a matrix to trace the estimation
  #------
  modelinc[35] = modeldesc$distno
  # 1 = sgt  ; 2 = sged; 3 = sst;
  maxOrder = 1
  modelnames = .expand.model(modelinc)

  pos = 1
  pos.matrix = matrix(0, ncol = 3, nrow = 28)
  colnames(pos.matrix) = c("start", "stop", "include")
  rownames(pos.matrix) = c("mu", "ar", "ma", "archm", "skm", "kum",
                           "omega", "alpha", "beta", "gamma",
                           "skew", "shape1","shape2",
                           "skcons", "skalpha", "skgamma", "skbeta", "volsk",
                           "sh1cons", "sh1alpha", "sh1gamma", "sh1beta","volsh1",
                           "sh2cons", "sh2alpha", "sh2gamma", "sh2beta", "volsh2")
  for (i in 1:28) {
    if (modelinc[i] > 0) {
      pos.matrix[i, 1:3] = c(pos, pos + modelinc[i] - 1, 1)
      pos = max(pos.matrix[1:i, 2] + 1)
    }
  }
  nn = length(modelnames)
  modelmatrix = matrix(0, ncol = 3, nrow = nn)
  rownames(modelmatrix) = modelnames
  colnames(modelmatrix) = c("opt", "fixed", "start")
  # Opt column will have value 1 for the parameters that need to be optimized
  # Fixed column will have value 1 for the parameters that will be fixed through out the esimation
  # Start column will have the value 1 for the parameters that have been specified the starting values
  fixed.names = names(fixed.pars)
  fp = charmatch(fixed.names, modelnames)

  if (!is.null(fixed.names) && any(!is.na(fp))) {
    fixed = fp[!is.na(fp)]
    modelmatrix[fixed, 2] = 1
    fz = charmatch(modelnames, fixed.names)
    fz = fz[!is.na(fz)]
    fixed.pars = fixed.pars[fz]
    names(fixed.pars) = fixed.names[fz]
  } else {
    fixed.pars = NULL
  }
  modelmatrix[, 1] = 1 - modelmatrix[, 2]
  start.names = names(start.pars)
  sp = charmatch(start.names, modelnames)
  if (!is.null(start.names) && any(!is.na(sp))) {
    start = sp[!is.na(sp)]
    modelmatrix[start, 3] = 1
    sz = charmatch(modelnames, start.names)
    sz = sz[!is.na(sz)]
    start.pars = start.pars[sz]
  } else {
    start.pars = NULL
  }
  #-----
  # Parameter Matrix
  #------
  mm = sum(modelinc[c(2,3,8,9,10,15,16,17,20,21,22,25,26,27)])
  mm = mm - length(which(modelinc[c(2,3,8,9,10,15,16,17,20,21,22,25,26,27)] > 0)) #Incase we have order larger than 1
  pars = matrix(0, ncol = 6, nrow = 28 + mm)
  colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
  pidx = matrix(NA, nrow = 28, ncol = 2)
  colnames(pidx) = c("begin", "end")
  rownames(pidx) = c("mu", "ar", "ma", "archm", "skm", "kum",
                     "omega", "alpha", "beta", "gamma",
                     "skew", "shape1","shape2",
                     "skcons", "skalpha", "skgamma", "skbeta", "volsk",
                     "sh1cons", "sh1alpha", "sh1gamma", "sh1beta","volsh1",
                     "sh2cons", "sh2alpha", "sh2gamma", "sh2beta", "volsh2")
  fixed.names = names(fixed.pars)
  pnames = NULL #parameter names
  nx = 0
  # post.matrx[,3] having value 1 if we need to include that parameters in the model
  if (pos.matrix[1, 3] == 1) {
    pars[1, 3] = 1
    pars[1, 1] = 0
    if (any(substr(fixed.names, 1, 2) == "mu"))
      pars[1, 2] = 1 else pars[1, 4] = 1
  }
  pidx[1, 1] = 1
  pidx[1, 2] = 1
  pnames = c(pnames, "mu")

  nx = 1
  pn = 1
  pidx[2, 1] = 2
  if (pos.matrix[2, 3] == 1) {
    pn = length(seq(pos.matrix[2, 1], pos.matrix[2, 2], by = 1))
    for (i in 1:pn) {
      pars[(nx + i), 1] = 0
      pars[(nx + i), 3] = 1
      nnx = paste("ar", i, sep = "")
      sp = na.omit(match(fixed.names, nnx))
      if (length(sp) > 0)
        pars[(nx + i), 2] = 1 else pars[(nx + i), 4] = 1
      pnames = c(pnames, nnx)
    }
  } else {
    pnames = c(pnames, "ar")
  }
  pidx[2, 2] = 1 + pn

  nx = nx + pn
  pn = 1
  pidx[3, 1] = nx + 1
  if (pos.matrix[3, 3] == 1) {
    pn = length(seq(pos.matrix[3, 1], pos.matrix[3, 2], by = 1))
    for (i in 1:pn) {
      pars[(nx + i), 1] = 0
      pars[(nx + i), 3] = 1
      nnx = paste("ma", i, sep = "")
      sp = na.omit(match(fixed.names, nnx))
      if (length(sp) > 0)
        pars[(nx + i), 2] = 1 else pars[(nx + i), 4] = 1
      pnames = c(pnames, nnx)
    }
  } else {
    pnames = c(pnames, "ma")
  }
  pidx[3, 2] = nx + pn

  nx = nx + pn
  pn = 1
  pidx[4, 1] = nx + 1
  if (pos.matrix[4, 3] == 1) {
    pars[nx + pn, 3] = 1
    pars[nx + pn, 1] = 0
    if (any(!is.na(match(fixed.names, "archm"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames, "archm")
  pidx[4, 2] = nx + pn

  nx = nx + pn
  pn = 1
  pidx[5, 1] = nx + 1
  if (pos.matrix[5, 3] == 1) {
    pars[nx + pn, 3] = 1
    pars[nx + pn, 1] = 0
    if (any(!is.na(match(fixed.names, "skm"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames, "skm")
  pidx[5, 2] = nx + pn

  nx = nx + pn
  pn = 1
  pidx[6, 1] = nx + 1
  if (pos.matrix[6, 3] == 1) {
    pars[nx + pn, 3] = 1
    pars[nx + pn, 1] = 0
    if (any(!is.na(match(fixed.names, "kum"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames, "kum")
  pidx[6, 2] = nx + pn

  nx = nx + pn
  pn = 1
  pidx[7, 1] = nx + 1
  if (pos.matrix[7, 3] == 1) {
    pars[nx + pn, 3] = 1
    pars[nx + pn, 1] = 0
    if (any(!is.na(match(fixed.names, "omega"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames, "omega")
  pidx[7, 2] = nx + pn

  nx = nx + pn
  pn = 1
  pidx[8, 1] = nx + 1
  if (pos.matrix[8, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "alpha"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"alpha")
  pidx[8,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[9, 1] = nx + 1
  if (pos.matrix[9, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "beta"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"beta")
  pidx[9,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[10, 1] = nx + 1
  if (pos.matrix[10, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "gamma"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"gamma")
  pidx[10,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[11, 1] = nx + 1
  if (pos.matrix[11, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "skew"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"skew")
  pidx[11,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[12, 1] = nx + 1
  if (pos.matrix[12, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "shape1"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"shape1")
  pidx[12,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[13, 1] = nx + 1
  if (pos.matrix[13, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "shape2"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"shape2")
  pidx[13,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[14, 1] = nx + 1
  if (pos.matrix[14, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "skcons"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"skcons")
  pidx[14,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[15, 1] = nx + 1
  if (pos.matrix[15, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "skalpha"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"skalpha")
  pidx[15,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[16, 1] = nx + 1
  if (pos.matrix[16, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "skgamma"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"skgamma")
  pidx[16,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[17, 1] = nx + 1
  if (pos.matrix[17, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "skbeta"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"skbeta")
  pidx[17,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[18, 1] = nx + 1
  if (pos.matrix[18, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "volsk"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"volsk")
  pidx[18,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[19, 1] = nx + 1
  if (pos.matrix[19, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh1cons"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh1cons")
  pidx[19,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[20, 1] = nx + 1
  if (pos.matrix[20, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh1alpha"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh1alpha")
  pidx[20,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[21, 1] = nx + 1
  if (pos.matrix[21, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh1gamma"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh1gamma")
  pidx[21,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[22, 1] = nx + 1
  if (pos.matrix[22,3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh1beta"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh1beta")
  pidx[22,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[23, 1] = nx + 1
  if (pos.matrix[23,3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "volsh1"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"volsh1")
  pidx[23,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[24, 1] = nx + 1
  if (pos.matrix[24, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh2cons"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh2cons")
  pidx[24,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[25, 1] = nx + 1
  if (pos.matrix[25, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh2alpha"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh2alpha")
  pidx[25,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[26, 1] = nx + 1
  if (pos.matrix[26, 3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh2gamma"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh2gamma")
  pidx[26,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[27, 1] = nx + 1
  if (pos.matrix[27,3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "sh2beta"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"sh2beta")
  pidx[27,2] = nx+pn

  nx = nx + pn
  pn = 1
  pidx[28, 1] = nx + 1
  if (pos.matrix[28,3] == 1) {
    pars[nx+pn,3] =1
    pars[nx+pn,1] =0
    if (any(!is.na(match(fixed.names, "volsh2"))))
      pars[nx + pn, 2] = 1 else pars[nx + pn, 4] = 1
  }
  pnames = c(pnames,"volsh2")
  pidx[28,2] = nx+pn

  rownames(pars) = pnames
  zf = match(fixed.names, rownames(pars))
  if (length(zf) > 0)
    pars[zf, 1] = unlist(fixed.pars)
  pars[, "LB"] = NA
  pars[, "UB"] = NA
  sbounds[7] = dmodel$exp.rate
  #-----
  # Output
  # modelinc: A nammed vector which specify which parameter to be included in the model
  # pars is a matrix, which sepecify all parameters, with level and fixed columns for parameters that is fixed, includes or not, should be estimate or not and the bounds
  # sbouds is the bounds of distribution parameters w.r.t the conditoinal distribution
  # start.pars is the starting values of parameters that is specified
  # fixed.pars is the fixed values of parameters that is specified
  # maxorder is the maximum order of all equations
  # pos.matrix is a matrix that sepecify the position of parameters, starting from the mean equation
  # pidx is the index value of parameters
  # vmodel is variance equation specification
  # mmmodel is the mean equation specification
  # dmodel is the higher moment equations specification
  #-----
  model = list(modelinc = modelinc, pars = pars, sbounds = sbounds, start.pars = start.pars, fixed.pars = fixed.pars, maxOrder = maxOrder,
               pos.matrix = pos.matrix, pidx = pidx, vmodel = vmodel, mmodel = mmodel, dmodel = dmodel)
  # set the output as an ACDspec class
  ans = new("ACDspec", model = model)
  return(ans)  ##Return an object with class ACDspec, which is a list, specified in line 878
}

setMethod(f = "acdspec", definition = .acdspec)
## Here we specify that the function .acdspec the definition function of the method "acdspec".
