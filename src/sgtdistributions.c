# include "sgtdistributions.h"
# include <R.h>
# include <limits.h>
# include <math.h>
# include <Rmath.h>
double sgn(double x)
{
  double ans;
  if(x >= 0){
    ans = 1;
  }
  if(x < 0){
    ans = -1;
  }
  return(ans);
}
/*
 * Skewd Generalized Error Distribution
 */
double dsged(const double value, const double mean, const double sig, const double sk, const double ku)
{
  double x = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku;
  double ans;

  sigma = sigma/sqrt((PI*(1.0+3.0*pow(lambda,2.0))*gammafn(3.0/kappa)- pow(16,1.0/kappa)
                        *pow(lambda,2)*pow(gammafn(0.5+(1.0/kappa)),2)*gammafn(1.0/kappa))/(PI*gammafn(1.0/kappa)));

  x = x + (pow(2.0,2.0/kappa)*sigma*lambda*gammafn(0.5+(1.0/kappa)))/sqrt(PI);

  ans =  kappa/(2.0*sigma*gammafn(1.0/kappa)*exp(pow(fabs(x-mu)/(sigma* (1.0 + lambda*sgn(x-mu))),kappa)));
  return(ans);
}

double qsged(const double value, const double mean, const double sig, const double sk, const double ku)
{
  double prob = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku;
  sigma = sigma/sqrt((PI*(1.0+3.0*pow(lambda,2.0))*gammafn(3.0/kappa)- pow(16,1.0/kappa)*pow(lambda,2)*pow(gammafn(0.5+(1.0/kappa)),2)*gammafn(1.0/kappa))/(PI*gammafn(1.0/kappa)));

  double lam;
  double out;
  double ans;

  if(prob < (1.0-lambda)/2){
    prob = 1 - prob;
    lam = -1.0*lambda;
    out = mu+ (-1.0)*sigma*(1+lam)*pow(qgamma(2.0*prob/(1.0 + lam) + (lam-1)/(lam+1),1.0/kappa,1.0,1,0),1.0/kappa);
  }
  else{
    out = mu+ sigma*(1+lambda)*pow(qgamma(2.0*prob/(1.0 + lambda) + (lambda-1)/(lambda+1),1.0/kappa,1.0,1,0),1.0/kappa);
  }
  ans = out - (pow(2,2.0/kappa)*sigma*lambda*gammafn(0.5+(1.0/kappa)))/(sqrt(PI));
  return(ans);
}
double psged(const double value, const double mean,const double sig,const double sk,const double ku)
{
  double x = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku;
  double ans;
  sigma = sigma/sqrt((PI*(1.0+3.0*pow(lambda,2.0))*gammafn(3.0/kappa)- pow(16,1.0/kappa)*pow(lambda,2)*pow(gammafn(0.5+(1.0/kappa)),2)*gammafn(1.0/kappa))/(PI*gammafn(1.0/kappa)));
  x = x + (pow(2.0,2.0/kappa)*sigma*lambda*gammafn(0.5+(1.0/kappa)))/sqrt(PI);
  x = x - mu;
  if(x < 0){
    lambda = -1.0 * lambda;
    x = -1.0 * x;
    ans = 1.0 - ((1.0-lambda)/2.0 + ((1.0+lambda)/2)*pgamma(pow((x/(sigma*(1+lambda))),kappa),1.0/kappa,1.0,1,0));
  } else{
    ans = (1.0-lambda)/2.0 + ((1.0+lambda)/2)*pgamma(pow((x/(sigma*(1+lambda))),kappa),1.0/kappa,1.0,1,0);
  }
  return(ans);
}
double rsged(const double mean,const double sig, const double sk, const double ku)
{
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku;
  double z, ans;
  z = runif(0.0,1.0);
  ans = qsged(z,mu,sigma,lambda,kappa);
  return(ans);
}

/*
* Skewed Student-t distribution
*/

double dsgt(const double value,const double mean,const double sig,const double sk,const double ku1,const double Rku2)
{
  double x = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku1;
  double nu = Rku2/ku1;
  double ans;
/* sigma here is the v*sigma in the sgt R documents */
  sigma = sigma/(pow(nu,(1/kappa)) * sqrt((3.0 * pow(lambda,2.0) + 1.0) * (beta(3.0/kappa,
                                          nu - 2.0/kappa)/beta(1.0/kappa, nu)) - 4.0 * pow(lambda,2.0) * pow((beta(2.0/kappa,
                                                               nu - 1.0/kappa)/beta(1.0/kappa, nu)),2)));
  x = x + (2.0 * sigma * lambda * pow(nu,(1/kappa)) * beta(2.0/kappa, nu - 1.0/kappa))/beta(1.0/kappa, nu);
  ans =  kappa/(2 * sigma * pow(nu,(1.0/kappa)) * beta(1.0/kappa, nu) * pow((1.0 + pow(fabs(x - mu),kappa)/(nu * pow(sigma,kappa)* pow((1.0 + lambda * sgn(x - mu)),kappa))),(nu + 1/kappa)));
  return(ans);
}
/* double dsgt(const double value,const double mean,const double sig,const double sk,const double ku1,const double ku2)
{
  double x = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku1;
  double nu = ku2;
  double beta1, beta2, beta3, nuAdj, varScale, A1, A2, p, u,m,C, ans;
  beta1 = beta(1/kappa,nu/kappa);
  beta2 = beta(2/kappa,(nu-1)/kappa);
  beta3 = beta(3/kappa,(nu-2)/kappa);
  nuAdj = (nu+1)/kappa;
  A2 = (1 + 3*pow(lambda,2)) * 1/beta1 * pow(nuAdj,2/kappa) * beta3;
  A1 = 2 * lambda * 1/beta1 * pow(nuAdj,1/kappa) * beta2;
  p = A1/sqrt(A2 - pow(A1,2));
  varScale = sigma/sqrt(A2-pow(A1,2));
  m = mu - p*sigma;
  u = x - m;
  C = 0.5 * kappa * pow(nuAdj,-1/kappa) * 1/varScale * 1/beta1;
  ans = C * pow((1 + (pow(fabs(u),kappa))/(nuAdj * pow(1 + sign(u)*lambda,kappa) * pow(varScale,kappa))),(-nuAdj));
  return(ans);
} */
double qsgt(const double value,const double mean,const double sig,const double sk,const double ku1,const double Rku2)
{
  double prob = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku1;
  double nu = Rku2/ku1;
  double ans;
  sigma = sigma/(pow(nu,(1/kappa)) * sqrt((3.0 * pow(lambda,2.0) + 1.0) * (beta(3.0/kappa,
                                          nu - 2.0/kappa)/beta(1.0/kappa, nu)) - 4.0 * pow(lambda,2.0) * pow((beta(2.0/kappa,
                                                               nu - 1.0/kappa)/beta(1.0/kappa, nu)),2)));
  double lam;
  double out;
  if(prob > (1.0-lambda)/2){
    prob = 1 - prob;
    lam = -1.0*lambda;
    out = mu+ (-1.0)*sigma * (lam - 1.0) * pow((1.0/(nu * qbeta(1.0 - 2.0 *(prob/(1.0 - lam)),1.0/kappa,nu,1,0)) - 1.0/nu),(-1.0/kappa));
  }
  else{
    out = mu+ sigma * (lambda - 1.0) * pow((1.0/(nu * qbeta(1.0 - 2.0 *(prob/(1.0 - lambda)),1.0/kappa,nu,1,0)) - 1.0/nu),(-1.0/kappa));
  }
  ans = out - (2.0 * sigma * lambda * pow(nu,(1.0/kappa)) * beta(2.0/kappa,nu - 1.0/kappa))/beta(1.0/kappa, nu);
  return(ans);
}
double psgt(const double value,const double mean,const double sig,const double sk,const double ku1,const double Rku2)
{
  double x = value;
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku1;
  double nu = Rku2/ku1;
  double ans;
  sigma = sigma/(pow(nu,(1/kappa)) * sqrt((3.0 * pow(lambda,2.0) + 1.0) * (beta(3.0/kappa,
                                          nu - 2.0/kappa)/beta(1.0/kappa, nu)) - 4.0 * pow(lambda,2.0) * pow((beta(2.0/kappa,
                                                               nu - 1.0/kappa)/beta(1.0/kappa, nu)),2)));
  x = x + (2.0 * sigma * lambda * pow(nu,(1/kappa)) * beta(2.0/kappa, nu - 1.0/kappa))/beta(1.0/kappa, nu);
  x = x - mu;
  double out;
  if(x > 0){
    lambda = -1.0 * lambda;
    x = -1.0 * x;
    out = (1.0 - lambda)/2.0 + (lambda - 1.0)/2.0 * pbeta(1.0/(1.0 + nu * pow((sigma * (1.0 - lambda)/(-1.0*x)),kappa)), 1.0/kappa, nu,1,0);
    ans = 1.0 - out;
  } else{
    out = (1.0 - lambda)/2.0 + (lambda - 1.0)/2.0 * pbeta(1.0/(1.0 + nu * pow((sigma * (1.0 - lambda)/(-1.0*x)),kappa)), 1.0/kappa, nu,1,0);
    ans = out;
  }
  return(ans);
}
double rsgt(const double mean,const double sig,const double sk,const double ku1,const double Rku2)
{
  double mu = mean;
  double sigma = sig;
  double lambda = sk;
  double kappa = ku1;
  double nu = Rku2;
  double z,ans;
  z = runif(0.0,1.0);
  ans = qsgt(z,mu,sigma,lambda,kappa,nu);
  return(ans);
}
/*
 * SGT - Bali et al 2008

double dsgtB(const double value,const double sk,const double ku1,const double ku2)
{
  double x = value;
  double lambda = sk;
  double kappa = ku1;
  double nu = ku2;
  double g,rho,theta,C,sig,ans,beta1,beta2,beta3,C1;
  beta1 = beta(nu/kappa,1.0/kappa);
  beta2 = beta((nu-1.0)/kappa,2.0/kappa);
  beta3 = beta((nu-2.0)/kappa,3.0/kappa);
  g = (1.0 + 3.0*pow(lambda,2.0)) * pow(beta1,(-1.0)) * pow((nu+1.0)/kappa,2.0/kappa) * beta3;
  rho = 2.0 * lambda * pow(beta1,(-1.0)) * pow((nu+1.0)/kappa,1.0/kappa) * beta2;
  theta = 1.0/(sqrt(g - pow(rho,2)));
  sig = rho * theta;
  C = 0.5 * kappa * pow((nu+1.0)/kappa,(-1.0/kappa)) * pow(beta1,(-1.0)) * pow(theta,(-1.0));
  C1 = 1.0 + (pow(fabs(x + sig),kappa))/(((nu+1.0)/kappa) * pow((1 + sgn(x + sig) * lambda),kappa) * pow(theta,kappa));
  ans = C*pow(C1,-1.0 * (nu + 1.0)/kappa);
  return(ans);
}
*/
/*
 * wrapper function
 */
double ddist(const double zz, const double hh, const double sk, const double ku1, const double ku2, const int ndis)
{
  /*ndist: 1- Skewed Generalized Student-t
  2 - Skewed Generalized Erorr distribution
  */
  double pdf=0;
  if(ndis==1)
  {
    pdf=dsgt(zz,0,1,sk, ku1 ,ku2)/hh;
  }
  if(ndis==2)
  {
    pdf=dsged(zz,0,1,sk,ku1)/hh;
  }
  if(ndis==3)
  {
    pdf=dsgt(zz,0,1,sk,2.0,ku2)/hh;
  }
  return pdf;
}

double rdist(const double sk, const double ku1, const double ku2, const int ndis)
{
  /*ndist: 1- Skewed Generalized Student-t
   2 - Skewed Generalized Erorr distribution
  3 - Skew-Student-t
  */
  double ans=0.0;
  if(ndis==1)
  {
    ans=rsgt(0, 1,sk,ku1,ku2);
  }
  if(ndis==2)
  {
    ans=rsged(0,1,sk,ku1);
  }
  if(ndis==3)
  {
    ans=rsgt(0,1,sk,2,ku2);
  }
  return(ans);
}
double pdist(const double q, const double mu, const double sigma, const double sk, const double ku1, const double ku2, const int ndis)
{
  /*ndist: 1- Skewed Generalized Student-t
   2 - Skewed Generalized Erorr distribution
  3 - Skew-Student-t
  */
  double ans=0.0;
  if(ndis==1)
  {
    ans=psgt(q, mu, sigma, sk, ku1,ku2);
  }
  if(ndis==2)
  {
    ans=psged(q, mu, sigma, sk,ku1);
  }
  if(ndis==3)
  {
    ans=psgt(q, mu, sigma, sk, 2, ku2);
  }
  return(ans);
}
double qdist(const double prob,const double mu, const double sigma, const double sk, const double ku1, const double ku2, const int ndis)
{
  /*ndist: 1- Skewed Generalized Student-t
   2 - Skewed Generalized Erorr distribution
  3 - Skew-Student-t
  */
  double ans=0.0;
  if(ndis==1)
  {
    ans=qsgt(prob, mu, sigma, sk, ku1,ku2);
  }
  if(ndis==2)
  {
    ans=psged(prob, mu, sigma, sk,ku1);
  }
  if(ndis==3)
  {
    ans=psgt(prob, mu, sigma, sk, 2, ku2);
  }
  return(ans);
}
/*
 * Higher moment calculation
 */
double Mr(const int r, const double theta, const double lambda, const double kappa, const double nu)
{
  double ans = 0.0;
  ans = 0.5 * (pow((1.0 + lambda),(r+1)) + pow(-1.0,r) * pow((1.0- lambda),(r+1)))*
    pow((nu+1.0)/kappa,r/kappa)*
    1/beta(nu/kappa,1.0/kappa)*
    beta((nu-r)/kappa,(r+1)/kappa)*
    pow(theta,r);
  return(ans);
}
double skewness(const double sk, const double ku1, const double ku2, const int ndis)
{
  double lambda = sk;
  double kappa = ku1;
  double nu = ku2;
  double ans = 0.0;
  if(ndis ==1){
    double beta1, beta2, beta3, nuAdj, g, rho, theta, delta, M3;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1.0)/kappa);
    beta3 = beta(3.0/kappa,(nu-2.0)/kappa);
    nuAdj = (nu+1.0)/kappa;
    g = (1.0 + 3.0*pow(lambda,2.0)) * 1.0/beta1 * pow(nuAdj,2.0/kappa) * beta3;
    rho = 2.0 * lambda * 1.0/beta1 * pow(nuAdj,1.0/kappa) * beta2;
    theta = 1.0/sqrt(g - pow(rho,2.0));
    delta = rho * theta;
    M3 = Mr(3,theta,lambda,kappa,nu);
    ans = M3 - 3.0*delta - pow(delta,3.0);
  }
  if(ndis == 2){
  double theta, S, m, A, A3;
  A = gammafn(2.0/kappa)/sqrt(gammafn(1.0/kappa) * gammafn(3.0/kappa));
  S = sqrt(1.0 + 3.0*pow(lambda,2) - 4.0 * pow(A,2) * pow(lambda,2));
  theta = sqrt(gammafn(1.0/kappa)/gammafn(3.0/kappa)) * 1.0/S;
  m = 2.0 * lambda * (A/S);
  A3 = 4 * lambda * (1 + pow(lambda,2)) * gammafn(4.0/kappa)* (1/gammafn(1.0/kappa)) * pow(theta,3);
  ans = A3 - 3.0*m - pow(m,3);
  }
  if(ndis == 3){
    double beta1, beta2, beta3, nuAdj, g, rho, theta, delta, M3;
    kappa = 2.0;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1)/kappa);
    beta3 = beta(3.0/kappa,(nu-2)/kappa);
    nuAdj = (nu+1.0)/kappa;
    g = (1.0 + 3.0*pow(lambda,2.0)) * 1.0/beta1 * pow(nuAdj,(2.0/kappa)) * beta3;
    rho = 2.0 * lambda * 1.0/beta1 * pow(nuAdj,(1.0/kappa)) * beta2;
    theta = 1.0/sqrt(g - pow(rho,2));
    delta = rho * theta;
    M3 = Mr(3,theta,lambda,kappa,nu);
    ans = M3 - 3*delta - pow(delta,3);
  }
  return(ans);
}

double kurtosis(const double sk, const double ku1, const double ku2, const int ndis)
{
  double lambda = sk;
  double kappa = ku1;
  double nu = ku2;
  double ans = 0.0;
  if(ndis ==1){
    double beta1, beta2, beta3, nuAdj, g, rho, theta, delta, M3, M4;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1)/kappa);
    beta3 = beta(3.0/kappa,(nu-2)/kappa);
    nuAdj = (nu+1.0)/kappa;
    g = (1.0 + 3.0*pow(lambda,2.0)) * 1.0/beta1 * pow(nuAdj,(2.0/kappa)) * beta3;
    rho = 2.0 * lambda * 1.0/beta1 * pow(nuAdj,(1.0/kappa)) * beta2;
    theta = 1.0/sqrt(g - pow(rho,2));
    delta = rho * theta;
    M3 = Mr(3,theta,lambda,kappa,nu);
    M4 = Mr(4,theta,lambda,kappa,nu);
    ans = M4 - 4*M3*delta +6*pow(delta,2) + 3*pow(delta,4);
  }
  if(ndis == 2){
    double theta, S, m, A, A3, A4;
    A = gammafn(2.0/kappa)/sqrt(gammafn(1.0/kappa) * gammafn(3.0/kappa));
    S = sqrt(1.0 + 3.0*pow(lambda,2.0) - 4.0 * pow(A,2) * pow(lambda,2));
    theta = sqrt(gammafn(1.0/kappa)/gammafn(3.0/kappa)) * 1.0/S;
    m = 2.0 * lambda * (A/S);
    A3 = 4.0 * lambda * (1.0 + pow(lambda,2.0)) * gammafn(4.0/kappa)* (1.0/gammafn(1.0/kappa)) * pow(theta,3.0);
    A4 = (1.0 + 10.0 * pow(lambda,2.0) + 5.0* pow(lambda,4.0)) * gammafn(5.0/kappa) * (1/gammafn(1.0/kappa)) * pow(theta,4.0);
    ans = A4 - 4.0*A3*m + 6.0*pow(m,2.0) + 3.0*pow(m,4.0);
  }
  if(ndis == 3){
    kappa = 2.0;
    double beta1, beta2, beta3, nuAdj, g, rho, theta, delta, M3, M4;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1)/kappa);
    beta3 = beta(3.0/kappa,(nu-2)/kappa);
    nuAdj = (nu+1.0)/kappa;
    g = (1.0 + 3.0*pow(lambda,2.0)) * 1.0/beta1 * pow(nuAdj,(2.0/kappa)) * beta3;
    rho = 2.0 * lambda * 1.0/beta1 * pow(nuAdj,(1.0/kappa)) * beta2;
    theta = 1.0/sqrt(g - pow(rho,2));
    delta = rho * theta;
    M3 = Mr(3,theta,lambda,kappa,nu);
    M4 = Mr(4,theta,lambda,kappa,nu);
    ans = M4 - 4*M3*delta +6*pow(delta,2) + 3*pow(delta,4);
  }
  return(ans);
}
double Pskewness(const double sk, const double ku1, const double ku2, const int ndis)
{
  double lambda = sk;
  double kappa = ku1;
  double nu = ku2;
  double A, S;
  double ans = 0.0;
  if(ndis ==1){
    double beta1, beta2, beta3;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1.0)/kappa);
    beta3 = beta(3.0/kappa,(nu-2.0)/kappa);
    A = beta2 * 1/sqrt(beta1) * 1/sqrt(beta3);
    S = sqrt(1 + 3*pow(lambda,2) - 4 * pow(A,2) * pow(lambda,2));
    ans = 2 * lambda * A * 1/S;
  }
  if(ndis == 2){
    double gamma1, gamma2, gamma3;
    gamma1 = gammafn(1.0/kappa);
    gamma2 = gammafn(2.0/kappa);
    gamma3 = gammafn(3.0/kappa);
    A = gamma2/sqrt(gamma1 * gamma3);
    S = sqrt(1 + 3*pow(lambda,2) - 4 * pow(A,2) * pow(lambda,2));
    ans = 2 * lambda * A/S;
  }
  if(ndis == 3){
    kappa = 2.0;
    double beta1, beta2, beta3;
    beta1 = beta(1.0/kappa,nu/kappa);
    beta2 = beta(2.0/kappa,(nu-1.0)/kappa);
    beta3 = beta(3.0/kappa,(nu-2.0)/kappa);
    A = beta2 * 1/sqrt(beta1) * 1/sqrt(beta3);
    S = sqrt(1 + 3*pow(lambda,2) - 4 * pow(A,2) * pow(lambda,2));
    ans = 2 * lambda * A * 1/S;
  }
  return(ans);
}
