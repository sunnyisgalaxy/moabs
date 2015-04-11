#include "types.h"
#include "stats.h"
#include <iostream>
//#include <math.h> //round() function
using namespace std;

int main(){
	//double pdiff(int  k1, int  n1, int  k2, int  n2, double  D)
	int  k1; int  n1; int  k2; int  n2; double  D;
	cout << "Function in request is pdiff(int k1, int n1, int k2, int n2, double D)" << endl;
	cout << "Enter these parameters separated by space/tab/return: \n";
	cin>>k1>>n1>>k2>>n2>>D;
	// the >> operator skips whitespace characters such as tabs, blank space and newline.  When eof is encountered, zero (false) is returned.
	cout << "k1,n1,k2,n2,D=" << k1 << "," << n1 << "," << k2 << "," << n2 << "," << D << endl; 
	cout << pdiff(k1,n1,k2,n2,D) << endl;	
	cout << pdiff(172,172,115,125, 0.01) << endl;	
	cout << pdiff(172,172,115,125,-0.01) << endl;	
	cout << "good" <<endl;
	//The other functions are pdiffCI() and singleCI(). 
}
