#' @exportClass ACDspec
#' @exportClass ACDfit
#' @exportClass ACDsim
#' @exportClass ACDpath
#----
setClass("rACD", "VIRTUAL")
#-----------
# ACDspec objects
#----------
setClass("ACDspec", representation(model = "vector"), contains = "rACD")
# If we want to create a new oject of class ACDspec, use new("ACDspec",model = objecttobeACDspecclass)
# Specify ACDspec as a S4 class
# contains specify the classes that ACDspec class inherit from.
#-------------------------------------------------------------------------------
setClass("ACDfit", representation(fit = "vector", model = "vector"), contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDfilter", representation(filter = "vector", model = "vector"), contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDforecast", representation(forecast = "vector", model = "vector"), contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDsim", representation(simulation = "vector", model = "vector", seed = "integer"), contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDpath", representation(path = "vector", model = "vector", seed = "integer"), contains = "rACD")
#-------------------------------------------------------------------------------
setClass("ACDroll", representation(model = "vector", forecast = "vector"), contains = "rACD")
