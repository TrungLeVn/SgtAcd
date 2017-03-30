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
# include "sgtdistributions.h"
# include "Sacd.h"

void skewfilter(int *model, double *pars, int *idx, double *z, double *tempskew, double *tskew,
		double *sbounds, double *h, int i, int T)
{
  int k = 0;
	if( model[10]>0 )
	{
	  tempskew[i] += skewdynamicsT(model,idx,pars,z[i-1],fabs(h[i-1]),tempskew[i-1],k);
		tskew[i] = logmap1dT(sbounds[0], sbounds[1], tempskew[i]);
	}
	else
	{
		tskew[i] = pars[idx[7]];
	}
}

void acdshape1filter(int *model, double *pars, int *idx, double *z, double *tempshape1, double *tshape1,
		double *sbounds, double *h, int i, int T)
{
  int k = 0;
	if( model[14]>0 )
	{
   tempshape1[i] += shape1dynamics(model,idx,pars,z[i-1],fabs(h[i-1]),tempshape1[i-1],k);
	 tshape1[i] = expmap1dT(sbounds[2], tempshape1[i], sbounds[6]);
/*	tshape1[i] = logmap1dT(sbounds[2],sbounds[3],tempshape1[i]);*/
	}
	else
	{
		tshape1[i] = pars[idx[8]];
	}
}
void acdshape2filter(int *model, double *pars, int *idx, double *z, double *tempshape2, double *tshape2,
                     double *sbounds, double *h, int i, int T)
{
  int k = 0;
  if( model[18]>0 )
  {
    tempshape2[i] += shape2dynamics(model, idx, pars, z[i-1],fabs(h[i-1]),tempshape2[i-1],k);
    /*  tshape2[i] = expmap2dT(sbounds[4],sbounds[5], tempshape2[i], sbounds[6]);*/
  tshape2[i] = logmap1dT(sbounds[4],sbounds[5],tempshape2[i]);
  }
  else
  {
    tshape2[i] = pars[idx[9]];
  }
}
// case 1 is pwl, case 2 is quad with squared, 3 is absolute. Shock = 2 for residuals, shock= 1 for standardized residuals

double skewdynamicsT(const int *model, const int *idx, const double *pars, const double z, const double h, const double lagval, const int lag)
{
	double res=0.0;
	switch(model[22])
	{
		case 0:
		{
			// quadratic model
			res = pars[idx[10]];
			if(model[11]>lag) res+=pars[idx[11]+lag]*z;
			if(model[12]>lag) res+=pars[idx[12]+lag]*(z*z);
			if(model[13]>lag) res+=pars[idx[13]+lag]*lagval;
			break;
		}

		case 1:
		{
			// quadratic model
			res = pars[idx[10]];
			if(model[11]>lag) res+=pars[idx[11]+lag]*z;
			if(model[12]>lag) res+=pars[idx[12]+lag]*fabs(z);
			if(model[13]>lag) res+=pars[idx[13]+lag]*lagval;
			break;
		}
	case 2:
		{
		  // piece-wise linear threshold model (pwl)
		  res = pars[idx[10]];
		  double z1 = (z <  0)? z : 0;
		  double z2 = (z >= 0)? z : 0;
		  if(model[11]>lag) res+=pars[idx[11]+lag]*z1;
		  if(model[12]>lag) res+=pars[idx[12]+lag]*z2;
		  if(model[13]>lag) res+=pars[idx[13]+lag]*lagval;
		  break;
		}
	}
	return res;
}
// case 1 is pwl, case 2 is quad with squared, 3 is absolute. Shock = 2 for residuals, shock= 1 for standardized residuals

double shape1dynamics(const int *model, const int *idx, const double *pars, const double z, const double h, const double lagval, const int lag)
{
	double res=0.0;
	switch(model[23])
	{

	  // quadric model with squared
		case 0:
		{
			res = pars[idx[14]];
			if(model[15]>lag) res+=pars[idx[15]+lag]*z;
			if(model[16]>lag) res+=pars[idx[16]+lag]*(z*z);
			if(model[17]>lag) res+=pars[idx[17]+lag]*lagval;
			break;
		}

		  // quadric model with absolute
		case 1:
		{
			res = pars[idx[14]];
			if(model[15]>lag) res+=pars[idx[15]+lag]*z;
			if(model[16]>lag) res+=pars[idx[16]+lag]*fabs(z);
			if(model[17]>lag) res+=pars[idx[17]+lag]*lagval;
			break;
		}
	case 2:
		{
		  res = pars[idx[14]];
		  double z1 = (z <  0)? (z*z) : 0;
		  double z2 = (z >= 0)? (z*z) : 0;
		  if(model[15]>lag) res+=pars[idx[15]+lag]*z1;
		  if(model[16]>lag) res+=pars[idx[16]+lag]*z2;
		  if(model[17]>lag) res+=pars[idx[17]+lag]*lagval;
		  break;
		}
	case 3:
		  {
		    res = pars[idx[14]];
		    double z1 = (z < 0)? fabs(z) : 0;
		    double z2 = (z >=0)? fabs(z) : 0;
		    if(model[15]>lag) res+=pars[idx[15]+lag]*z1;
		    if(model[16]>lag) res+=pars[idx[16]+lag]*z2;
		    if(model[17]>lag) res+=pars[idx[17]+lag]*lagval;
		    break;
		  }
	}
	return res;
}

double shape2dynamics(const int *model, const int *idx, const double *pars, const double z, const double h, const double lagval, const int lag)
{
  double res=0.0;
  switch(model[24])
  {
    // pwl model with standardized residuals

    // quadric model with squared
  case 0:
  {
    res = pars[idx[18]];
    if(model[19]>lag) res+=pars[idx[19]+lag]*z;
    if(model[20]>lag) res+=pars[idx[20]+lag]*(z*z);
    if(model[21]>lag) res+=pars[idx[21]+lag]*lagval;
    break;
  }

    // quadric model with absolute
  case 1:
  {
    res = pars[idx[18]];
    if(model[19]>lag) res+=pars[idx[19]+lag]*z;
    if(model[20]>lag) res+=pars[idx[20]+lag]*fabs(z);
    if(model[21]>lag) res+=pars[idx[21]+lag]*lagval;
    break;
  }
  case 2:
  {
    res = pars[idx[18]];
    double z1 = (z <  0)? (z*z) : 0;
    double z2 = (z >= 0)? (z*z) : 0;
    if(model[19]>lag) res+=pars[idx[19]+lag]*z1;
    if(model[20]>lag) res+=pars[idx[20]+lag]*z2;
    if(model[21]>lag) res+=pars[idx[21]+lag]*lagval;
    break;
  }
  case 3:
  {
    res = pars[idx[18]];
    double z1 = (z < 0)? fabs(z) : 0;
    double z2 = (z >=0)? fabs(z) : 0;
    if(model[19]>lag) res+=pars[idx[19]+lag]*z1;
    if(model[20]>lag) res+=pars[idx[20]+lag]*z2;
    if(model[21]>lag) res+=pars[idx[21]+lag]*lagval;
    break;
  }
  }
  return res;
}

double logmap1dT(const double LB, const double UB, const double val)
{
	double res=0.0;
	res=LB+(UB-LB)/(1.0 + exp(-1*val));
	return res;
}

double invlogmap1dT(const double LB, const double UB, const double res)
{
	double val=0.0;
	val=-1*log(-(UB-res)/(-res+LB));
	return val;
}

double expmap2dT(const double LB, const double UB, const double val, const double rate)
{
  double res=0.0;
  res = LB + exp(-1.0*rate*val)*UB;
  return res;
}

double invexpmap2dT(const double LB, const double UB, const double res, const double rate)
{
  double val=0.0;
  val=-(1.0/rate)*log((res - LB)/UB);
  return val;
}

double expmap1dT(const double LB, const double val, const double rate)
{
  double res=0.0;
  res = LB + exp(1.0*rate*val);
  return res;
}

double invexpmap1dT(const double LB, const double res, const double rate)
{
  double val=0.0;
  val=(1.0/rate)*log((res - LB));
  return val;
}
