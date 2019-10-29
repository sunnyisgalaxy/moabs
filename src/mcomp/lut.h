#include <map>
#include <utility>

extern std::map < MultiKey, double > lut_pdiffInRegion;
void build_lut_pdiffInRegion(int  tableMax, int  numThreads);
void build_lut_pdiffInRegion_singleThread(int tableMax);

extern std::map < MultiKey, CI > lut_pdiffCI;
void build_lut_pdiffCI_singleThread(int tableMax);
void build_lut_pdiffCI(int  tableMax, int  numThreads);

extern std::map < std::pair <int, int>, CI > lut_singleCI;
void build_lut_singleCI_singleThread(int tableMax);
void build_lut_singleCI(int  tableMax, int  numThreads);


// `lut_fet` is not used
// extern std::map < MultiKey, double > lut_fet;
// void build_lut_fet(int  tableMax, int  numThreads);
// void build_lut_fet_singleThread(int tableMax);
