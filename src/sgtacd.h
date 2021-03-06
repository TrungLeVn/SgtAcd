/*################################################################################
##
##   R package racd by Alexios Ghalanos Copyright (C) 2012, 2013
##   This file is part of the R package racd.
##
##   The R package racd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package racd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
#ifndef SGTACD_H
#define SGTACD_H
 void armafilterC(int *model, double *pars, int *idx, double *hm, double *skm, double *psk, double *x, double *res,
                  double *zrf, double *constm, double *condm,
                  int *m, int *T);
 void sacd(int *model, double *pars, int *idx, double *hEst, double *x,
           double *res, double *e, double *zrf,
           double *constm, double *condm, int *m, int *T, double *h, double *z,
           double *tempskew, double *tempshape1,double *tempshape2, double *skhEst, double *tskew,
           double *tshape1,double *tshape2, double *sbounds,
           double *llh, double *LHT, double *skew, double *kurt, double *pskew);
 void gjracd(int *model, double *pars, int *idx, double *hEst, double *x,
             double *res, double *e, double *nres, double *zrf,
             double *constm, double *condm, int *m, int *T, double *h, double *z,
             double *tempskew, double *tempshape1,double *tempshape2, double *skhEst, double *tskew,
             double *tshape1, double *tshape2, double *sbounds,
             double *llh, double *LHT, double *skew, double *kurt, double *pskew);
 void sacdsimC(int *model, double *pars, int *idx,double *x, double *constm, double *h, double *z,
               double *res, double *e, double *tempskew, double *tempshape1, double *tempshape2,
               double *tskew, double *tshape1, double *tshape2,double *sbounds,double *pskew,
               int *T, int *m);
 void gjracdsimC(int *model, double *pars, int *idx,double *x, double *constm, double *h, double *z,
                 double *res, double *e, double *nres, double *tempskew, double *tempshape1, double *tempshape2,
                 double *tskew, double *tshape1, double *tshape2,double *sbounds,double *pskew,
                 int *T, int *m);

 void meansimC(int *model, double *pars, int *idx, double *h, double *sk, double *ps, double *x, double *res,
               double *constm, int *m, int *T);
#endif /* SGTACD_H */
