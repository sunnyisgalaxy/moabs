#include <map>

//CI singleCI(int  k, int  n, double alpha, int method);
CI singleCI(int  k, int  n, double alpha, int method); //right now, method = symmetric width
double pdiff(int  k1, int  n1, int  k2, int  n2, double  d, double tolerance=TOLER);
CI pdiffCI(int  k1, int  n1, int  k2, int  n2, double alpha, int method);
//float AdaptLob(float (*f)(float, const void*), float a, float b, float Tol, const void *Param);
//int fet_1k( double & pvalue, int contingency_table[], int nrow, int ncol);

extern std::map < MultiKey, double > lut_fet_1k;

double pdiffInRegion(int m1, int t1, int m2, int t2, double d);
double pdiffInRegionOpt(int m1, int t1, int m2, int t2, double d);
CI pdiffCIOpt(int  k1, int  n1, int  k2, int  n2, double alpha, int method);
