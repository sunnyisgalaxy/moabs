//#include <boost/function.hpp>//dsun: must be put before nr3.h gamma.h and incgammabeta.h, otherwise can not compile.!!!!
//#include <boost/format.hpp>
#include <iostream>
#include <boost/math/distributions.hpp>
//#include "nr3.h"
#include "types.h"
#include "lut.h"

//#include "adapt.h"

using namespace std;
using namespace boost;
using namespace boost::math;



////////////////////adapt.h##################
struct Adapt {
	long double TOL,toler;
	static const double alpha,beta,x1,x2,x3,x[12];
	bool terminate,out_of_tolerance;
	Adapt(long double tol);
	template <class T>
	double integrate(T &func, const double a, const double b);
	template <class T>
	double adaptlob(T &func, const double a, const double b, const double fa,
		const double fb, const double is);
};

Adapt::Adapt(long double tol) : TOL(tol),terminate(true),out_of_tolerance(false)
{
	const long double EPS=numeric_limits<long double>::epsilon();
	if (TOL < 10.0*EPS)
		TOL=10.0*EPS;
}

template <class T>
double Adapt::integrate(T &func, const double a, const double b)
{
	double m,h,fa,fb,i1,i2,is,erri1,erri2,r,y[13];
	m=0.5*(a+b);
	h=0.5*(b-a);
	fa=y[0]=func(a);
	fb=y[12]=func(b);
	for (int i=1;i<12;i++)
		y[i]=func(m+x[i]*h);
	i2=(h/6.0)*(y[0]+y[12]+5.0*(y[4]+y[8]));
	i1=(h/1470.0)*(77.0*(y[0]+y[12])+432.0*(y[2]+y[10])+
		625.0*(y[4]+y[8])+672.0*y[6]);
	is=h*(0.0158271919734802*(y[0]+y[12])+0.0942738402188500*
		(y[1]+y[11])+0.155071987336585*(y[2]+y[10])+
		0.188821573960182*(y[3]+y[9])+0.199773405226859*
		(y[4]+y[8])+0.224926465333340*(y[5]+y[7])+
		0.242611071901408*y[6]);
	erri1=abs(i1-is);
	erri2=abs(i2-is);
	r=(erri2 != 0.0) ? erri1/erri2 : 1.0;
	toler=(r > 0.0 && r < 1.0) ? TOL/r : TOL;
	if (is == 0.0)
		is=b-a;
	is=abs(is);
	return adaptlob(func,a,b,fa,fb,is);
}

template <class T>
double Adapt::adaptlob(T &func, const double a, const double b, const double fa,
		const double fb, const double is)
{
	double m,h,mll,ml,mr,mrr,fmll,fml,fm,fmrr,fmr,i1,i2;
	m=0.5*(a+b);
	h=0.5*(b-a);
	mll=m-alpha*h;
	ml=m-beta*h;
	mr=m+beta*h;
	mrr=m+alpha*h;
	fmll=func(mll);
	fml=func(ml);
	fm=func(m);
	fmr=func(mr);
	fmrr=func(mrr);
	i2=h/6.0*(fa+fb+5.0*(fml+fmr));
	i1=h/1470.0*(77.0*(fa+fb)+432.0*(fmll+fmrr)+625.0*(fml+fmr)+672.0*fm);
	if (abs((long double)i1-(long double)i2) <= toler*(long double)is || mll <= a || b <= mrr) {
		if ((mll <= a || b <= mrr) && terminate) {
			out_of_tolerance=true;
			terminate=false;
		}
		return i1;
	}
	else
		return adaptlob(func,a,mll,fa,fmll,is)+
			adaptlob(func,mll,ml,fmll,fml,is)+
			adaptlob(func,ml,m,fml,fm,is)+
			adaptlob(func,m,mr,fm,fmr,is)+
			adaptlob(func,mr,mrr,fmr,fmrr,is)+
			adaptlob(func,mrr,b,fmrr,fb,is);
}
const double Adapt::alpha=sqrt(2.0/3.0);
const double Adapt::beta=1.0/sqrt(5.0);
const double Adapt::x1=0.942882415695480;
const double Adapt::x2=0.641853342345781;
const double Adapt::x3=0.236383199662150;
const double Adapt::x[12]={0,-x1,-alpha,-x2,-beta,-x3,0.0,x3,beta,x2,alpha,x1};

////////////////////////////adapt.h####################################


/* * singleCI.cpp
 * Given k positive results from n observations, what is the confidence interval at confidence level 1-alpha ?
 * 	TODO: need immplement other methods;
 * Created on: Apr 28, 2011
 * Author: deqiangs
 */
CI singleCI(int  k, int  n, double alpha, int method) //right now, method = symmetric width
{
	double pNormal = double(k)/n;
	double d = 0.5 * ( pNormal - (-1) ); //initial guess of symmetric width	
	double x1 = 0;
	double x2 = 1.0;
	double dif = 1;

	int iter = 0;

	int a = k + 1;
	int b = n - k + 1;
	beta_distribution <> f(a,b);

	double h = 1.0;
	double l = 0.0;

	while( abs(dif) > TOLER && iter < MAXITER )
	{
		x1 = (pNormal - d < 0) ? 0 : (pNormal - d);
		x2 = (pNormal + d > 1) ? 1 : (pNormal + d);

		dif = ( cdf(f, x2) - cdf(f, x1) ) - ( 1 - alpha );
		
		if(dif > 0){ //d is too big
			h = d;
			d = ( l + h ) * 0.5;
			
		}
		else if(dif < 0){ //d is too small
			l = d;
			d = ( l + h ) * 0.5;
		}
		
		iter ++;
	}
	CI ci;
	ci.a = x1;
	ci.b = x2;
	return ci;
};

/* * pdiff.cpp

 * Given k1(k2) positive results from n1(n2) observations, what is the probability of 
 *  r1,i.e., ratio 1 greater than r2 by d? P(r1-r2>d)=?
 * Created on: Apr 28, 2011
 * Author: deqiangs
 */

class JointProb {
	int a1, b1, a2, b2;
	double D;
public:
	JointProb(int A1, int B1, int A2, int B2, double d) : a1(A1), b1(B1), a2(A2), b2(B2), D(d) { };
	double operator()(double x) const {
		beta_distribution <> f1(a1, b1);
		beta_distribution <> f2(a2, b2);
		
		//extend pdf(f, x)=0 for x=[-Infinity, 0) , (1, +Infinity];
		//extend cdf(f, x)=0 for x=[-Infinity, 0), cdf(f, x)=1 for x=(1, +Infinity];
		if( D >= 0.0 ) {					//D>=0;
			if( x < D ) 		{ return 0.0; } 	//cdf(f2, x-D)=0 for x<D;
			else if ( x > 1.0 ) 	{ return 0.0; } 	//pdf(f1, x)=0 for x>1;
			else 			{ return pdf(f1, x) * cdf(f2, x-D); }
		} else { 						//D<0;
			if( x < 0.0 ) 		{ return 0.0; } 	//pdf(f1, x)=0 for x<0;
			else if ( x <= 1.0+D ) 	{ return pdf(f1, x) * cdf(f2, x-D); } //for 0 <= x <= 1+D;
			else if ( x <= 1.0 ) 	{ return pdf(f1, x); } 	//for 1+D < x <= 1.0, cdf(f2, x-D)=1; 
			else 			{ return 0.0; }		//pdf(f1,x)=0 for x>1.
		}
		//return pdf(f1, x) * cdf(f2, x-D);
	}
};


double pdiff(int  k1, int  n1, int  k2, int  n2, double  d, double tolerance=TOLER)
{
	int a1 = k1 + 1; 
	int b1 = n1 - k1 + 1; 
	int a2 = k2 + 1; 
	int b2 = n2 - k2 + 1;
	JointProb jointProb(a1, b1, a2, b2, d);
	//cout << "IN PDIFF " << k1 << "\t" << n1 << "\t" << k2 << "\t" << n2 << "\t" << d << endl;
 	int no_std_dev = 25;
 	double mean1 = a1/double(a1+b1); 
	double mean2 = min(a2/double(a2+b2) + d, 1.0);
	//for this dataset (34,34,2551,2551,0.01), mean2 is 1.01! 
	double sd1 = sqrt(double(a1)*b1/(double(a1+b1)*(a1+b1)*(a1+b1+1)));
 	double sd2 = sqrt(double(a2)*b2/(double(a2+b2)*(a2+b2)*(a2+b2+1)));
	double lower = max( max(mean1-no_std_dev*sd1, mean2-no_std_dev*sd2), max(0.0, d) );
	double upper = min( mean1 + no_std_dev*sd1, 1.0 );
	
	Adapt Ad(tolerance); 

	double mainpart = Ad.integrate(jointProb, lower, upper);
	if(lower > upper){ double tmp = lower; lower = upper; upper = tmp; }
	
	double testLower = Ad.integrate( jointProb, max(0.0, lower - no_std_dev * sd1), lower );
	double testUpper = Ad.integrate( jointProb, upper, min(1.0, upper + no_std_dev * sd1) );
	if( testLower + testUpper > tolerance && testLower + testUpper > 0.001*mainpart ){ 
		std::cerr<<"please increase no_std_dev at this tolerance for pdiff:" << k1 <<","<< n1 <<","<< k2 <<","<< n2 <<","<< d << " AT "<< testLower <<","<< testUpper <<","<< mainpart <<"\t" << tolerance << endl;
		std::cerr << "For " << mean1 << "\t" << sd1 << "\t" << mean2 << "\t" << sd2 << "\t" << lower << "\t" << upper << endl;
		//exit(1);
	}
	//cout << "mainpart, testLower, testUpper=" << mainpart << "\t" << testLower <<"\t" << testUpper << "\tAx1.0" << mainpart*1.0 << "\t" << testLower*1.0<< "\t" <<testUpper*1.0<<endl;
	
	//cout << "lower - upper A:\t" << mainpart << endl; 
	//double remainingpart = Ad.integrate(jointProb, 0, lower) + Ad.integrate(jointProb, upper, 1);
	
	//cout << "0, lower, upper, 1 A:\t" << final << endl;
	//cout << "DIFFERENCE=" << final - mainpart << endl;
	//double final = 
	return mainpart + testLower + testUpper;
	//return Ad.integrate(jointProb, lower, upper);
}


//probability of P(|r1-r2|<d)
double pdiffInRegion(int m1, int t1, int m2, int t2, double d)
{
	//double tolerance = 10.0 * numeric_limits<double>::epsilon();
	double nominalDif = double(m2)/t2-double(m1)/t1;
	double p_acc = 0;

	if( nominalDif > 0 ) //P(-d<p1-p2<d)=P(-d<p2-p1<d)= P(p1-p2>-d)-P(p1-p2>d) = P(p2-p1>-d)-P(p2-p1>d) for accruacy reason
	{
		p_acc = pdiff(m2,t2,m1,t1,-d, TENEPS)-pdiff(m2,t2,m1,t1,d, TENEPS);
	} else 
	{
		p_acc = pdiff(m1,t1,m2,t2,-d, TENEPS)-pdiff(m1,t1,m2,t2,d, TENEPS);
	}
	if(p_acc + TOLER <0 || p_acc -TOLER >1){
		std::cerr << "check the function pdiffInRegion(" << m1 <<","<< t1 <<","<< m2 <<","<< t2 <<","<< d << ")" << endl;
		exit(0);
	}
	p_acc = min(p_acc, 1.0);
	p_acc = max(p_acc, 0.0);
	return p_acc;
}

double pdiffInRegionOpt(int m1, int t1, int m2, int t2, double d)
{
	MultiKey combi(t1,m1,t2,m2);
	//MultiKey combi(n1,k1,n2,k2);
	if(lut_pdiffInRegion.count(combi) == 1){
		return lut_pdiffInRegion[combi];
	} else {
		return pdiffInRegion( m1, t1, m2, t2, d);
	}
}

/* * pdiffCI.cpp
 * Given k1/n1, k2/n2, and alpha, what is the confidence interval of r1-r2 at confidence level of 1-alpha ?
 *
 * TODO: method is not currently immplimented. As of now 20110708, method is based on one side.
 *
 * For Adapt integration method, for 2000 calculations of pdiffCI(5,10,0,10,0.05,1) at tol=1e-6, sandbox takes 27 seconds; So for 2,000,000 such calculations, it will take 7.5 hours. By setting tol = 1e-5, it is 15 seconds;
 * For GaussLobattoFor   method, for 2000 calculations of pdiffCI(5,10,0,10,0.05,1) at tol=1e-5, sandbox takes 28 seconds; So for 2,000,000 such calculations, it will take 7.5 hours.
 * Created on: Apr 28, 2011
 * Author: deqiangs
 */

//P(t>d) = 0.95, where is alpha = 0.05
double pdiffCI_low(int  k1, int  n1, int  k2, int  n2, double alpha, int method)
{
	double dNominal = double(k1)/n1 - double(k2)/n2;
	double alphaComp = 1 - alpha;

	double l = -1.0;
	double h = dNominal;
	double guess = ( l + h ) /2.0;

	int iter = 0;

	double p = pdiff(k1,n1,k2,n2, guess);
	while(abs(p - alphaComp) > TOLER && iter < MAXITER){ //cout <<"THIS LOW: "<< iter << "\t" << h << "\t" << l << "\t" << p << "\t" << p-alphaComp << endl;
		if( p < alphaComp - TOLER ){
			h = guess;
			guess = (l + guess) / 2.0;
		}
		else if( p > alphaComp + TOLER )
		{
			l = guess;
			guess = (guess + h ) / 2.0;
		}
		iter ++;
		
		p = pdiff(k1,n1,k2,n2, guess);
	}
	//cout << "LOW" << endl;
	return guess;
};

//P(t<d) = 0.95 or P(t>d) = 0.05, where alpha = 0.05
double pdiffCI_high(int  k1, int  n1, int  k2, int  n2, double alpha, int method)
{
	double dNominal = double(k1)/n1 - double(k2)/n2;

	double l = dNominal;
	double h = 1.0;
	double guess = ( l + h ) /2.0;

	int iter = 0;

	double p = pdiff(k1,n1,k2,n2, guess);
	while(abs(p - alpha) > TOLER && iter < MAXITER){ //cout <<"THIS HIGH: "<< iter << "\t" << h << "\t" << l << "\t" << p << "\t" << p-alpha << endl;
		if( p < alpha - TOLER ){
			h = guess;
			guess = (l + guess) / 2.0;
		}
		else if( p > alpha + TOLER )
		{
			l = guess;
			guess = (guess + h ) / 2.0;
		}
		iter ++;
		//cout << "INSIDE: " << k1 <<"," << n1 <<"," << k2 <<"," << n2 <<"," << guess << endl;
		p = pdiff(k1,n1,k2,n2, guess);
	}
	//cout << "HIGH"<<endl;
	return guess;
};


CI pdiffCI(int  k1, int  n1, int  k2, int  n2, double alpha, int method)
{
	CI ci;
	ci.a = pdiffCI_low(k1, n1, k2,  n2, alpha, method);
	ci.b = pdiffCI_high(k1, n1, k2,  n2, alpha, method);
	return ci;
};

CI pdiffCIOpt(int  k1, int  n1, int  k2, int  n2, double alpha, int method)
{
	MultiKey combi(n1,k1,n2,k2);

	//added for fitting model where the fit could be very broad/divergent
	if(n1 == 0 || n2 == 0){
		CI ci;
		ci.a = 0;
		ci.b = 0;
	}

	if(lut_pdiffCI.count(combi) == 1){
		return lut_pdiffCI[combi];
	} else {
		return pdiffCI( k1, n1, k2, n2, alpha, method);
	}
}
