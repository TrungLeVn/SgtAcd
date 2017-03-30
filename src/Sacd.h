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
#ifndef __ARCD_H__
#define __ARCD_H__
void skewfilter(int *model, double *pars, int *idx, double *z, double *tempskew,
		double *tskew, double *sbounds, double *h, int i, int T);
void acdshape1filter(int *model, double *pars, int *idx, double *z, double *tempshape1,
		double *tshape1, double *sbounds, double *h, int i, int T);
void acdshape2filter(int *model, double *pars, int *idx, double *z, double *tempshape2,
                    double *tshape2, double *sbounds, double *h, int i, int T);
double skewdynamicsT(const int *, const int *, const double *, const double ,
		const double , const double , const int );
double shape1dynamics(const int *, const int *, const double *, const double ,
		const double , const double, const int );
double shape2dynamics(const int *, const int *, const double *, const double ,
                     const double , const double, const int );
double logmap1dT(const double , const double , const double );
double invlogmap1dT(const double , const double , const double );
double expmap1dT(const double, const double, const double );
double invexpmap1dT(const double , const double, const double );
double expmap2dT(const double, const double, const double, const double);
double invexpmap2dT(const double, const double, const double, const double);
#endif
