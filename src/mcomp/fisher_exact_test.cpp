#include<iostream>
#include "types.h"
#include <map>

extern "C" {
  #include "fexact.h" //The C code, modified from R, for Fisher's test.
}

using namespace std;

std::map < MultiKey, double > lut_fet_1k;

void fisher_exact_test( double & pvalue, int contingency_table[], int nrow, int ncol) //preset workspace=1e5;mult=30;hybrid setting; for the R.fisher.test() function.
{
	//output paramters
	double prt = 1.0;
	double *prtp = &prt;
	//double pvalue = 1.0;
	double *pvaluep = &pvalue;

	//preset parameters
	int workspace = 200000;
	int mult = 30;
	int hybrid = (nrow > 2 || ncol >2) ? 1 : 0;

	//input parameters to use non-hybrid method
	double expected = -1.0;
	double percnt = 100.0;
	double emin = 0.0;
	
	if (hybrid) {
		expected = 5.0;
		percnt = 80.0;
		emin = 1.0;
        }
	//use hybrid if table is not 2x2 and if 80% of cells >= expected and if remaining >=emin;

	int *workspacep = &workspace;
	int *multp = &mult;
	double *expectedp = &expected;
	double *percntp = &percnt;
	double *eminp = &emin;
	int *contingency_table_p = contingency_table;
	int *nrowp = &nrow;
	int *ncolp = &ncol;
	
	fexact(nrowp, ncolp, contingency_table_p, nrowp, expectedp, percntp,eminp, prtp, pvaluep, workspacep, multp);

}


int fet_1k( double & pvalue, int table[], int nrow, int ncol ) //table might be changed
{	
	try{
		//if it's 2x2 table, check lookup table first
		if(nrow==2 && ncol==2){
			//table format : {m1, t1 - m1, m2, t2 - m2};
			int k1 = table[0];
			int n1 = table[0]+table[1];
			int k2 = table[2];
			int n2 = table[2]+table[3];
			MultiKey combi(n1,k1,n2,k2);
			if(lut_fet_1k.count(combi) == 1){
				pvalue = lut_fet_1k[combi];
				return 0;
			}
		}
		int sum = 0;
		for(int i = 0; i < nrow * ncol; i += 2){
			sum += table[i] + table[i+1];
			if(table[i] + table[i+1] > 1000){
				table[i]   = int( double(table[i])   / double(table[i] + table[i+1]) * 1000 + 0.5 );
				//table[i+1] = int( double(table[i+1]) / double(table[i] + table[i+1]) * 1000 + 0.5 ); if for 0.5 vs 999.5 this will generate 1 vs 1000;
				table[i+1] = 1000 - table[i];
			}
		}

		fisher_exact_test( pvalue, table, nrow, ncol );

		//if it's 2x2 table and sum <=100, write to lookup table.
		//todo
		//possible bug: this works for single thread.
		//if I can run fet_1k under multi threads, this function need be rewritten.
		if(sum < 101 && nrow==2 && ncol==2){
			int k1 = table[0];
			int n1 = table[0]+table[1];
			int k2 = table[2];
			int n2 = table[2]+table[3];
			MultiKey combi(n1,k1,n2,k2);
			lut_fet_1k[combi] = pvalue;
		}
		return 0;
	}
	catch (exception& e){
	  cout << e.what() << endl;
	}
	catch(...){
		std::cerr <<"unknown error" <<std::endl;
	}
	
}
