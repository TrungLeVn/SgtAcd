/*################################################################################
##
##   R package rgarch by Alexios Ghalanos Copyright (C) 2009
##   This file is part of the R package rgarch.
##
##   The R package rgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
# include <R.h>
# include <math.h>
# include "sgtfilters.h"
# include "Sacd.h"
 void sacdfilter(int *model, double *pars, int *idx, double *e, int T, int i, double *h)
 {
   int j;
   h[i] = h[i] + pars[idx[3]];//pars[idx[3]] is starting value of omega
   for( j=0; j<model[4]; j++ )
   {
     h[i] = h[i] + pars[idx[4]+j]*e[i-(j+1)]; //pars[idx[4]] is the starting values of alpha1 in variance process
   }
   for( j=0; j<model[5]; j++ )
   {
     h[i] = h[i] + pars[idx[5]+j]*h[i-(j+1)]; //pars[idx[5]] is the starting values of beta 1 in variance process
   }
 }
 // NOTE the value of nres
 void gjracdfilter(int *model, double *pars, int *idx, double *nres, double *e, int T, int i, double *h)
 {
   int j;
   h[i] = h[i] + pars[idx[3]];
   for( j=0; j<model[4]; j++ )
   {
     h[i] = h[i] + pars[idx[4]+j]*e[i-(j+1)]+pars[idx[6]+j]*nres[i-(j+1)];
   }
   for( j=0; j<model[5]; j++ )
   {
     h[i] = h[i] + pars[idx[5]+j]*h[i-(j+1)];
   }
 }
void armafilter(int* model, double *pars, int *idx, double *x, double *res,
                      double *zrf, double *constm, double *condm, int m, int i, int T)
{
/* --------------------------------------------------------------------------------
 * ARFIMA Process :
 * (1-L)^(-darfima).e[t] = phi(1-L)(y[t] - mu[t]) - psi(L).e[t]
 * where mu[t] includes constant, external and arch-in-mean
 * L is the lag operator
 * --------------------------------------------------------------------------------
 * */
	/*0 constm, 1 condm, 2 res*/
	constm[i] = pars[0];
	condm[i]+=constm[i];
	//ARMA initialization
	if(model[1]>0 || model[2]>0)
	{
		if(i>=model[1])
		{
			if(model[1]>0)
			{
			  condm[i] += pars[idx[1]]*x[i-1];
			}
			if(model[2]>0)
			{
			  condm[i]+=pars[idx[2]]*(x[i-1]-condm[i-1]);
			}
		}
	}
	res[i]=x[i]-condm[i];
}

