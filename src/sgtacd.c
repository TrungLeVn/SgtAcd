/*################################################################################
##
##   R package rarcd by Alexios Ghalanos Copyright (C) 2012,2013
##   This file is part of the R package rarcd.
##
##   The R package rarcd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rarcd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
# include <R.h>
# include <math.h>
# include "sgtacd.h"
# include "sgtfilters.h"
# include "Sacd.h"
# include "sgtdistributions.h"

void armafilterC(int *model, double *pars, int *idx, double *hm, double *skm, double *psk, double *x, double *res,
		double *zrf, double *constm, double *condm,int *m, int *T)
{
	int i;
  double h = 0;
  double sk = 0;
  double sh = 0;
	for(i=0; i<*T; i++)
	{
	    h = hm[i];
	    sk = skm[i];
	    sh = psk[i];
	    armafilter(model, pars, idx, h, sk, sh, x, res,
                zrf, constm, condm, *m, i, *T);
	}
}
void sacd(int *model, double *pars, int *idx, double *hEst, double *x,
                 double *res, double *e, double *zrf,
                 double *constm, double *condm, int *m, int *T, double *h, double *z,
                 double *tempskew, double *tempshape1,double *tempshape2, double *skhEst, double *tskew,
                 double *tshape1, double *tshape2,double *sbounds,
                 double *llh, double *LHT, double *skew, double *kurt, double *pskew)
{

  int i;
  double lk=0;
// Handle the initial days with lagged. Start the recursion from day m+1 in which m is maximum order
  for(i=0; i<*m; i++)
  {
    if(model[13]>0)
    {
      tempskew[i] = skhEst[0];
      tskew[i] = logmap1dT(sbounds[0], sbounds[1], tempskew[i]);
    } else{
      tskew[i] = pars[idx[10]];
    }
    if(model[18]>0)
    {
      tempshape1[i] = skhEst[1];
       tshape1[i] = expmap1dT(sbounds[2], tempshape1[i], sbounds[6]);
      /* tshape1[i] = logmap1dT(sbounds[2], sbounds[3], tempshape1[i]); */
    } else{
      tshape1[i] = pars[idx[11]];
    }
    if(model[23]>0)
    {
      tempshape2[i] = skhEst[2];
      /* tshape2[i] = expmap2dT(sbounds[4], sbounds[5],tempshape2[i], sbounds[6]);*/
       tshape2[i] = logmap1dT(sbounds[4], sbounds[5],tempshape2[i]);
    } else{
      tshape2[i] = pars[idx[12]];
    }
    skew[i] = skewness( tskew[i], tshape1[i],tshape2[i], model[34]);
    kurt[i] = kurtosis(tskew[i], tshape1[i],tshape2[i], model[34]);
    pskew[i] = Pskewness(tskew[i],tshape1[i],tshape2[i],model[34]);
    h[i] = *hEst;
    armafilter(model, pars, idx, sqrt(fabs(*hEst)),tskew[i],pskew[i], x, res, zrf, constm, condm, *m, i, *T);
    e[i] = res[i] * res[i];
    z[i] = res[i]/sqrt(fabs(h[i]));
    LHT[i] = log(ddist(z[i], sqrt(fabs(h[i])), tskew[i], tshape1[i],tshape2[i], model[34]));
    lk = lk - LHT[i];
  }
  for (i=*m; i<*T; i++)
  {
    sacdfilter(model, pars, idx, e, *T, i, h);
// After this step, we have h is a column vector of conditonal variance
   /* hm = sqrt(fabs(h[i])); */
    if(model[31]==1){
      skewfilter(model, pars, idx, z, tempskew, tskew, sbounds, h, i, *T);
    } else{
      skewfilter(model, pars, idx, res, tempskew, tskew, sbounds, h, i, *T);
    }
    if(model[32]==1){
      acdshape1filter(model, pars, idx, z, tempshape1, tshape1, sbounds, h, i, *T);
    } else{
      acdshape1filter(model, pars, idx, res, tempshape1, tshape1, sbounds, h, i, *T);
    }
    if(model[33]==1){
      acdshape2filter(model, pars, idx, z, tempshape2, tshape2, sbounds, h, i, *T);
    } else{
      acdshape2filter(model, pars, idx, res, tempshape2, tshape2, sbounds, h, i, *T);
    }
    skew[i] = skewness(tskew[i], tshape1[i],tshape2[i], model[34]);
    kurt[i] = kurtosis(tskew[i], tshape1[i],tshape2[i], model[34]);
    pskew[i] = Pskewness(tskew[i],tshape1[i],tshape2[i],model[34]);
    armafilter(model, pars, idx, sqrt(fabs(h[i])),tskew[i],pskew[i],x, res, zrf, constm, condm, *m, i, *T);
    e[i] = res[i] * res[i];
    z[i] = res[i]/sqrt(fabs(h[i]));
    LHT[i] = log(ddist(z[i], sqrt(fabs(h[i])), tskew[i], tshape1[i],tshape2[i], model[34]));
    lk = lk - LHT[i];
  }
  *llh=lk;
}
void gjracd(int *model, double *pars, int *idx, double *hEst, double *x,
                 double *res, double *e, double *nres, double *zrf,
                 double *constm, double *condm, int *m, int *T, double *h, double *z,
                 double *tempskew, double *tempshape1,double *tempshape2, double *skhEst, double *tskew,
                 double *tshape1,double *tshape2, double *sbounds,
                 double *llh, double *LHT, double *skew, double *kurt, double *pskew)
{

  int i;
  double lk=0;
  // Handle the initial days with lagged. Start the recursion from day m+1 in which m is maximum order
  for(i=0; i<*m; i++)
  {
    if(model[13]>0)
    {
      tempskew[i] = skhEst[0];
      tskew[i] = logmap1dT(sbounds[0], sbounds[1], tempskew[i]);
    } else{
      tskew[i] = pars[idx[10]];
    }
    if(model[18]>0)
    {
      tempshape1[i] = skhEst[1];
      tshape1[i] = expmap1dT(sbounds[2], tempshape1[i], sbounds[6]);
      /* tshape1[i] = logmap1dT(sbounds[2], sbounds[3], tempshape1[i]); */
      } else{
      tshape1[i] = pars[idx[11]];
    }
    if(model[23]>0)
    {
      tempshape2[i] = skhEst[2];
      /*  tshape2[i] = expmap2dT(sbounds[4], sbounds[5],tempshape2[i], sbounds[6]);*/
      tshape2[i] = logmap1dT(sbounds[4], sbounds[5],tempshape2[i]);
    } else{
      tshape2[i] = pars[idx[12]];
    }
    skew[i] = skewness( tskew[i], tshape1[i],tshape2[i], model[34]);
    kurt[i] = kurtosis(tskew[i], tshape1[i],tshape2[i], model[34]);
    pskew[i] = Pskewness(tskew[i],tshape1[i],tshape2[i],model[34]);
    h[i] = *hEst;
    armafilter(model, pars, idx,sqrt(fabs(h[i])),tskew[i],pskew[i], x, res, zrf, constm, condm, *m, i, *T);
    e[i] = res[i] * res[i];
    nres[i] = res[i] < 0.0 ? e[i] : 0.0;
    z[i] = res[i]/sqrt(fabs(h[i]));
    LHT[i] = log(ddist(z[i], sqrt(fabs(h[i])), tskew[i], tshape1[i],tshape2[i], model[34]));
    lk = lk - LHT[i];
  }
  for (i=*m; i<*T; i++)
  {
    gjracdfilter(model, pars, idx, nres, e, *T, i, h);
    // After this step, we have h is a column vector of conditonal variance
   /* hm = sqrt(fabs(h[i])); */
    if(model[31]==1){
      skewfilter(model, pars, idx, z, tempskew, tskew, sbounds, h, i, *T);
    } else{
      skewfilter(model, pars, idx, res, tempskew, tskew, sbounds, h, i, *T);
    }
    if(model[32]==1){
      acdshape1filter(model, pars, idx, z, tempshape1, tshape1, sbounds, h, i, *T);
    } else{
      acdshape1filter(model, pars, idx, res, tempshape1, tshape1, sbounds, h, i, *T);
    }
    if(model[33]==1){
      acdshape2filter(model, pars, idx, z, tempshape2, tshape2, sbounds, h, i, *T);
    } else{
      acdshape2filter(model, pars, idx, res, tempshape2, tshape2, sbounds, h, i, *T);
    }
    skew[i] = skewness(tskew[i], tshape1[i],tshape2[i], model[34]);
    kurt[i] = kurtosis(tskew[i], tshape1[i],tshape2[i], model[34]);
    pskew[i] = Pskewness(tskew[i],tshape1[i],tshape2[i],model[34]);
    armafilter(model, pars, idx, sqrt(fabs(h[i])),tskew[i],pskew[i],x, res, zrf, constm, condm, *m, i, *T);
    e[i] = res[i] * res[i];
    nres[i] = (res[i] < 0.0) ? e[i] : 0.0;
    z[i] = res[i]/sqrt(fabs(h[i]));
    LHT[i] = log(ddist(z[i], sqrt(fabs(h[i])), tskew[i], tshape1[i],tshape2[i], model[34]));
    lk = lk - LHT[i];
  }
  *llh=lk;

}
