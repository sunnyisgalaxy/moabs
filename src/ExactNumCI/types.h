#pragma once

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

    MultiKey(long N1, int K1, long N2, long K2)
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




//std::string itos(int i);
//int string_to_int( std::string s);
