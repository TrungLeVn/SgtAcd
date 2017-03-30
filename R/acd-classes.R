#' @exportClass ACDspec
#' @exportClass ACDfit
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
