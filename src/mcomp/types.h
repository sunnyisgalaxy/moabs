#pragma once

#include <vector>

#define TOLER 1e-10
#define TENEPS 1e-16
//you may change this to 1e-18 by change double to long double for toler in Adapt.h
#define MAXITER 100
#define DIFFDELTA 0.01
//make it option?



struct CI
{
	double a;
	double b;
};

class MultiKey {
  public:
    int  n1;
    int  k1;
    int  n2;
    int  k2;

    MultiKey(int N1, int K1, int N2, int K2)
      : n1(N1), k1(K1), n2(N2), k2(K2) {}

    bool operator<(const MultiKey &right) const
    {
        if ( n1 == right.n1 ) {
            if ( k1 == right.k1 ) {
                if ( n2 == right.n2 ) {
                    return k2 < right.k2;
                }
                else {
                    return n2 < right.n2;
                }
            }
            else {
                return k1 < right.k1;
            }
        }
        else {
            return n1 < right.n1;
        }
    }
};

class BiKey {
  public:
    int  n1;
    int  k1;

    BiKey()
      : n1(-1), k1(-1) {}

    BiKey(int N1, int K1)
      : n1(N1), k1(K1) {}

    bool operator<(const BiKey &right) const
    {
        if ( n1 == right.n1 ) {
        	return k1 < right.k1;
        }
        else {
            return n1 < right.n1;
        }
    }
};

struct MergeLaneElement {
	int start;
	int end;
	int k;
	int n;
	int tcp;
	int mcp;
	int tcm;
	int mcm;
	char strand;
	char next;
	std::string chr;
	MergeLaneElement(): start(-1), end(-1), k(-1), n(-1), tcp(-1), mcp(-1), tcm(-1), mcm(-1), strand('x'), next('X') { }
};


////hmm dmr method parameters
//typedef struct
//{
//	int minDepth;			//sequencing depth for initial DMR calling
//	double minDiff;  		//ratio difference for initial DMR calling
//	int maxIterations;   	//maximum number of iterations for training
//	int minRegionDist;     	//minimum distance between two DMRs. DMCs with distance smaller than this threshold will be merged into a region.
//	double minP;           	//threshold for confidence
//	int trainPct;      		// 1/trainPct of initial regions for training.
//	int samplingSaturation; // for hmmDmr calling, precision is not that strict. 100 is good enough to set as max sequencing depth.
//	int minDmcs;			//min number of Dmcs to form a Dmr
//	double relaxedCondition;	//relaxedCondition > 0 : allow Dmr to contain relaxedCondtion * Num_DMCs of sites at middle hmm state.
//}HMMCONFIG;

//	config.maxIterations = 100;
//	config.minP = 0.95;
//	config.trainPct = 10; //25 optimal; training is about 1/25 of total initial DMRs. ~ 1 chrom ~.
//	config.minDiff = 0.33; //0.3333 optimal
//	config.minDepth = 5; //3 optimal;
//	config.minRegionDist = 300; //500 optimal
//	config.samplingSaturation = 100;
//	config.minDmcs = 3;
//	config.relaxedCondition = 0.2;

//hmm dmr method parameters
class HMMCONFIG
{
public:
	int trainPct;      		// 1/trainPct of initial regions for training.
	int maxIterations;   	//maximum number of iterations for training
	int minDepth;			//sequencing depth for initial DMR calling
	double minDiff;  		//ratio difference for initial DMR calling
	int minRegionDist;     	//minimum distance between two DMRs. DMCs with distance smaller than this threshold will be merged into a region.
	int minDmcs;			//min number of Dmcs to form a Dmr
	double relaxedCondition;	//relaxedCondition > 0 : allow Dmr to contain relaxedCondtion * Num_DMCs of sites at middle hmm state.
	double minP;           	//threshold for confidence
	int samplingSaturation; // for hmmDmr calling, precision is not that strict. 100 is good enough to set as max sequencing depth.
	HMMCONFIG(){
		trainPct = 25;
		maxIterations = 100;
		minDepth = 3;
		minDiff = 0.3333;
		minRegionDist = 500;
		minDmcs = 3;
		relaxedCondition = 0.8;
		minP = 0.95;
		samplingSaturation = 100;
	}
};


std::string itos(int i);
int string_to_int( std::string s);
std::string get_exepath();
std::string do_readlink(std::string const& path);
bool isconsensus(std::vector< int > & tcs, std::vector< int > & mcs, int depth=10); // depth>=10, max-min<0.1
