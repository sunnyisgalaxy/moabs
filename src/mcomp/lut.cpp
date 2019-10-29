//lookup tables: lut.cpp

#define ALPHA 0.05
#define BATCHMAX 2000
//#include <RInside.h>                            // for the embedded R via RInside
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>

#include "types.h"
#include "lut.h"
#include "stats.h"
#include "hmm.h"
#include "fisher_exact_test.h"
#include "fet2x2.h"

#include <iostream>
#include <iomanip>
#include <time.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>

#include <math.h> //round() function

namespace po = boost::program_options;
using namespace std;


////Table for pdiffInRegion
map<MultiKey, double> lut_pdiffInRegion;


void build_lut_pdiffInRegion_singleThread(int tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){

			for(int n2 = 1; n2 <= tableMax; n2 ++){

				for(int k2 = 0; k2 <= n2; k2 ++){

					MultiKey combi(n1,k1,n2,k2);
					lut_pdiffInRegion[combi] = pdiffInRegion(k1,n1,k2,n2,DIFFDELTA);
					//cout << lut_pdiffInRegion[combi] << endl;
				}
			}
		}
	}
}


void build_lut_pdiffInRegion_init(int & tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){

			for(int n2 = 1; n2 <= tableMax; n2 ++){

				for(int k2 = 0; k2 <= n2; k2 ++){

					MultiKey combi(n1,k1,n2,k2);
					//lut_pdiffInRegion[combi] = -1;//pdiffInRegion(k1,n1,k2,n2,DIFFDELTA);
					lut_pdiffInRegion[combi] = -1;
				}
			}
		}
	}
}

void build_lut_pdiffInRegion_step(int & offset, int & step)
{
	map<MultiKey, double>::iterator it = lut_pdiffInRegion.begin();
	//for(map<MultiKey, double>::iterator it = lut_pdiffInRegion.begin() + offset;it != lut_pdiffInRegion.end(); it += step ){ //it + x not allowd?
	for(int i = 0; i < offset && it != lut_pdiffInRegion.end(); i++){
				it++;
	}
	for( ;it != lut_pdiffInRegion.end();){
		it->second = pdiffInRegion(it->first.k1, it->first.n1, it->first.k2, it->first.n2, DIFFDELTA);
		for(int i = 0; i < step && it != lut_pdiffInRegion.end(); i++){
			it++;
		}
	}
}


void build_lut_pdiffInRegion(int  tableMax, int  numThreads)
{
	build_lut_pdiffInRegion_init(tableMax);
	//int tableSize = lut_pdiffInRegion.size();

	boost::thread_group g;

	for(int i = 0; i < numThreads; i++){
		boost::thread *tp = new boost::thread( build_lut_pdiffInRegion_step, i, numThreads);
		g.add_thread(tp);
	}
	g.join_all();

//	for(map<MultiKey, double>::iterator it = lut_pdiffInRegion.begin() ;it != lut_pdiffInRegion.end(); it ++ ){
//		cout << it->second << endl;
//	}
}


////Table for pdiffCI
map<MultiKey, CI> lut_pdiffCI;


void build_lut_pdiffCI_singleThread(int tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){

			for(int n2 = 1; n2 <= tableMax; n2 ++){

				for(int k2 = 0; k2 <= n2; k2 ++){

					MultiKey combi(n1,k1,n2,k2);
					lut_pdiffCI[combi] = pdiffCI(k1,n1,k2,n2,ALPHA, 1);
					//cout << lut_pdiffCI[combi] << endl;
				}
			}
		}
	}
}


void build_lut_pdiffCI_init(int & tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){

			for(int n2 = 1; n2 <= tableMax; n2 ++){

				for(int k2 = 0; k2 <= n2; k2 ++){

					MultiKey combi(n1,k1,n2,k2);
					CI ci; ci.a = -1; ci.b = -1;
					lut_pdiffCI[combi] = ci;
				}
			}
		}
	}
}

void build_lut_pdiffCI_step(int & offset, int & step)
{
	map<MultiKey, CI>::iterator it = lut_pdiffCI.begin();
	//for(map<MultiKey, double>::iterator it = lut_pdiffCI.begin() + offset;it != lut_pdiffCI.end(); it += step ){ //it + x not allowd?
	for(int i = 0; i < offset && it != lut_pdiffCI.end(); i++){
				it++;
	}
	for( ;it != lut_pdiffCI.end();){
		it->second = pdiffCI(it->first.k1, it->first.n1, it->first.k2, it->first.n2, ALPHA, 1);
		for(int i = 0; i < step && it != lut_pdiffCI.end(); i++){
			it++;
		}
	}
}


void build_lut_pdiffCI(int  tableMax, int  numThreads)
{
	build_lut_pdiffCI_init(tableMax);
	//int tableSize = lut_pdiffCI.size();

	boost::thread_group g;

	for(int i = 0; i < numThreads; i++){
		boost::thread *tp = new boost::thread( build_lut_pdiffCI_step, i, numThreads);
		g.add_thread(tp);
	}
	g.join_all();

//	for(map<MultiKey, CI>::iterator it = lut_pdiffCI.begin() ;it != lut_pdiffCI.end(); it ++ ){
//		cout << it->first.n1 << "\t" << it->first.k1 << "\t" << it->first.n2 << "\t" << it->first.n2 << "\t=";
//		cout << it->second.a << "\t" << it->second.b << endl;
//	}
}






//result from R version:
//1	0	3	3	0.24999999999999994449
//1	0	4	1	1.0000000000000004441
//acutally the result should be:
//1	0	3	3	0.25
//1	0	4	1	1

////Table for singleCI
map<pair <int, int>, CI> lut_singleCI;

void build_lut_singleCI_singleThread(int tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){
			pair <int, int> combi(n1,k1);
			lut_singleCI[combi] = singleCI(k1,n1,ALPHA, 1);
		}
	}
}


void build_lut_singleCI_init(int & tableMax)
{
	for(int n1 = 1; n1 <= tableMax; n1 ++){
		//cout << n1 << endl;
		for(int k1 = 0; k1 <= n1; k1 ++){
			pair <int, int> combi(n1,k1);
			CI ci; ci.a = -1; ci.b = -1;
			lut_singleCI[combi] = ci;
		}
	}
}

void build_lut_singleCI_step(int & offset, int & step)
{
	map<pair<int, int>, CI>::iterator it = lut_singleCI.begin();
	//for(map<pair<int, int>, CI>::iterator it = lut_singleCI.begin() + offset;it != lut_singleCI.end(); it += step ){ //it + x not allowd?
	for(int i = 0; i < offset && it != lut_singleCI.end(); i++){
				it++;
	}
	for( ;it != lut_singleCI.end();){
		pair <int, int> combi = it->first;
		it->second = singleCI(combi.second, combi.first, ALPHA, 1);
		for(int i = 0; i < step && it != lut_singleCI.end(); i++){
			it++;
		}
	}
}


void build_lut_singleCI(int  tableMax, int  numThreads)
{
	build_lut_singleCI_init(tableMax);
	//int tableSize = lut_singleCI.size();

	boost::thread_group g;

	for(int i = 0; i < numThreads; i++){
		boost::thread *tp = new boost::thread( build_lut_singleCI_step, i, numThreads);
		g.add_thread(tp);
	}
	g.join_all();

//	for(map<MultiKey, double>::iterator it = lut_singleCI.begin() ;it != lut_singleCI.end(); it ++ ){
//		cout << it->second << endl;
//	}
}





// `lut_fet` is not used, remove it for now
//
//      ////because the c code from R for fisher.test ('fexact.c') is not thread safe. I replaced it fet_1k with functions from 'fisher.c'
//      ////Table for fisher exact test
//      map<MultiKey, double> lut_fet;
//      
//      void build_lut_fet_singleThread(int tableMax)
//      {
//      	for(int n1 = 1; n1 <= tableMax; n1 ++){
//      		//cout << n1 << endl;
//      		for(int k1 = 0; k1 <= n1; k1 ++){
//      
//      			for(int n2 = 1; n2 <= tableMax; n2 ++){
//      
//      				for(int k2 = 0; k2 <= n2; k2 ++){
//      
//      					MultiKey combi(n1,k1,n2,k2);
//      					double value = -1;
//      					int varray[] = {k1, n1-k1, k2, n2-k2};
//      					fet_1k( value, varray, 2, 2 );
//      					lut_fet[combi] = value;
//      					//cout << lut_pdiffInRegion[combi] << endl;
//      				}
//      			}
//      		}
//      	}
//      }
//      
//      void build_lut_fet_init(int & tableMax)
//      {
//      	for(int n1 = 1; n1 <= tableMax; n1 ++){
//      		//cout << n1 << endl;
//      		for(int k1 = 0; k1 <= n1; k1 ++){
//      
//      			for(int n2 = 1; n2 <= tableMax; n2 ++){
//      
//      				for(int k2 = 0; k2 <= n2; k2 ++){
//      
//      					MultiKey combi(n1,k1,n2,k2);
//      					//lut_pdiffInRegion[combi] = -1;//pdiffInRegion(k1,n1,k2,n2,DIFFDELTA);
//      					lut_fet[combi] = -1;
//      				}
//      			}
//      		}
//      	}
//      }
//      
//      void build_lut_fet_step(int & offset, int & step)
//      {
//      	map<MultiKey, double>::iterator it = lut_fet.begin();
//      	//for(map<MultiKey, double>::iterator it = lut_fet.begin() + offset;it != lut_fet.end(); it += step ){ //it + x not allowd?
//      	for(int i = 0; i < offset && it != lut_fet.end(); i++){
//      				it++;
//      	}
//      	for( ;it != lut_fet.end();){
//      		//double value = -1;
//      		//int varray[] = {it->first.k1, it->first.n1 - it->first.k1, it->first.k2, it->first.n2 - it->first.k2};
//      		//int varray[] = {k1, n1-k1, k2, n2-k2};
//      		//fet_1k( value, varray, 2, 2 );
//      		//cout << it->first.k1 << "\t" << it->first.n1 << "\t" << it->first.k2 << "\t" << it->first.n2 << "\t" << value << endl;
//      		it->second = fet2x2(it->first.k1, it->first.n1 - it->first.k1, it->first.k2, it->first.n2 - it->first.k2); //value;
//      		for(int i = 0; i < step && it != lut_fet.end(); i++){
//      			it++;
//      		}
//      	}
//      }
//      
//      
//      void build_lut_fet(int  tableMax, int  numThreads)
//      {
//      	build_lut_fet_init(tableMax);
//      	//int tableSize = lut_fet.size();
//      	//cout << "finish init" << endl;
//      	boost::thread_group g;
//      
//      	for(int i = 0; i < numThreads; i++){
//      		boost::thread *tp = new boost::thread( build_lut_fet_step, i, numThreads);
//      		g.add_thread(tp);
//      	}
//      	g.join_all();
//      
//      //	for(map<MultiKey, double>::iterator it = lut_fet.begin() ;it != lut_fet.end(); it ++ ){
//      //		cout << it->second << endl;
//      //	}
//      }
//      
