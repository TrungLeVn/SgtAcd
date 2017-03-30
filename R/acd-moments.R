momentsT <- function(x, skew,shape,opt){
  A = gamma(2/shape) * (gamma(1/shape)^(-1/2)) * (gamma(3/shape))^(-1/2)
  S = sqrt(1 + 3*skew^2 - 4 * A^2 * skew^2)
  m = 2 * skew * A * (1/S)
  theta = (gamma(1/shape))^(1/2) * (gamma(3/shape))^(-1/2) * (1/S)
  C = shape/(2*theta) * (1/gamma(1/shape))
  A3 = 4*skew*(1+skew^2)*gamma(4/shape)*(1/gamma(1/shape))*(theta^3)
  A4 = (1+10*(skew^2) + 5*(skew^4)) * gamma(5/shape) * (1/gamma(1/shape)) * (theta^4)
  pdf = C * exp(-((abs(x+m))^shape)/(((1+sign(x + m) * skew)^shape)*(theta^shape)))
  skewness = A3 - 3 * m - m^3
  kurtosis = A4 - 4*A3*m + 6*(m^2) + 3*(m^4) -3
  if(opt == 1){
    return(skewness)
  }else if(opt ==2){
    return(kurtosis)
  } else{
    return(pdf)
  }
}
