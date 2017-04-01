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
#ifndef __FILTERS_H__
#define __FILTERS_H__
 void sacdfilter(int *model, double *pars, int *idx, double *e, int T, int i, double *h);
 void gjracdfilter(int *model, double *pars, int *idx, double *nres, double *e, int T, int i, double *h);
 void armafilter(int* model, double *pars, int *idx, double *x, double *res, double *zrf,
                       double *constm, double *condm, double hm, double sk, double ku, int m, int i, int T);
#endif
