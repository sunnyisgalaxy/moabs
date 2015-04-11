#include <map>

extern std::map < MultiKey, double > lut_fet_1k;

int fet_1k( double & pvalue, int contingency_table[], int nrow, int ncol);
