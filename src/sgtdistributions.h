#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
double sgn(double x);
double dsged(const double, const double, const double , const double , const double );
double qsged(const double , const double , const double , const double , const double );
double psged(const double , const double ,const double ,const double ,const double );
double rsged(const double ,const double , const double , const double );
double dsgt(const double ,const double ,const double ,const double ,const double ,const double );
double qsgt(const double ,const double ,const double ,const double ,const double ,const double );
double psgt(const double ,const double ,const double ,const double ,const double ,const double );
double rsgt(const double ,const double ,const double ,const double ,const double );
double ddist(const double ,const double , const double , const double , const double , const int );
double rdist(const double , const double , const double , const int );
double pdist(const double , const double , const double , const double , const double , const double , const int );
double qdist(const double ,const double , const double , const double , const double , const double , const int );
double dsgtB(const double,const double,const double,const double);
double skewness(const double, const double, const double, const int);
double kurtosis(const double, const double, const double, const int);
#endif /* DISTRIBUTIONS_H */

