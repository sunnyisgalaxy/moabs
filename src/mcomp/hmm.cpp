#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iterator>
#include "float.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h> //round() function
#include <sstream>
//#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
//#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
//#include <boost/thread.hpp>
//#include <boost/program_options.hpp>
//#include <boost/date_time.hpp>

#include <math.h> //round() function
#include "types.h"
#include "lut.h"

//namespace po = boost::program_options;
using namespace std;

#define NUM_CHROM 100      //maximum chrom numbers
#define NSTATE     3       //Number of hidden states; 0|1|2 = L1 small | no diff | L2 small ?

#define TRANS_TABLE_ELEM_PREC 0.00001


//struct for the initial differential methylation regions
typedef struct
{
	int chromIndex;
	int start;				//index of genomic coordinate of start of region
	int end;				//index of genomic coordinate of end of region
	int cxSites;			//number of cpg sites in this region.
	int *siteCountsM1;      //mC counts in L1
	int *siteCountsT1;      //tC counts in L1
	int *siteCountsM2;	   	//mC counts in L2
	int *siteCountsT2;	   	//tC counts in L2
	double *prob;          	//posterior probabilities of states
	int *isSignificant;    	//if the cpg is a dmc: 1 if L1 small, -1 if L2 small, 0 if non-differential. ?
}DMR_REGION;




map<int, map<int, int> > State;
map<int, map<int, int> > loci; //chromIndex->startIndex => start coordinate
map<int, string> chromNameByIndex; // index -> name

DMR_REGION *regions;


//prior prability of the first bin in HMM.
double priorProb[NSTATE];

//the transition probability table
double transitionTable[NSTATE*NSTATE];

//the emmission probability table
map<int, map<MultiKey, double> > emmTable;


int regionNum;

HMMCONFIG config;



int getSiteCount(string & inputFile);
void getSignificantSites();
void getOptimalStates();
int buildLookupTable();
void initHMM();
void forwardBackward(double *newTransition, int step);
double stopClimb(double *newTransition);
void outputFile(string & outputName);
void freeMem();

//int string_to_int( string s){return atoi(s.c_str());}
//string itos(int i)
//{
//    stringstream s;
//    s << i;
//    return s.str();
//}


class cpgMeth
{
public:
	int totalC;
	int methC;
	cpgMeth()
	{
		totalC = 0;
		methC = 0;
	};
	cpgMeth(int t, int m) : totalC(t), methC(m) {};
	void out(){
		cout << totalC << "\t" << methC << endl;
	}

};

void readCompToHashSimple(string fileName, int i, int j, map <int, map <string, map<int, cpgMeth> > > & lane ){

	std::stringstream fakeLaneId;
	fakeLaneId << "100" << i << "100" << j;
	int laneIdFake = string_to_int(fakeLaneId.str());
	int chrIndex = -1;
	int startIndex = -1;
	int t1Index = -1;
	int r1Index = -1;
	int t2Index = -1;
	int r2Index = -1;

	string line = "";
	ifstream inputf(fileName.c_str(), ios::in);

	while (inputf.good()) {
		getline(inputf, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));

		if( find(line.begin(), line.begin()+1, '#') == line.begin() )
		{
			chrIndex = find (fields.begin(), fields.end(), "#chrom") - fields.begin();
			startIndex = find (fields.begin(), fields.end(), "start") - fields.begin();
			t1Index = find (fields.begin(), fields.end(), "totalC_" + itos(i-1)) - fields.begin();
			r1Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(i-1)) - fields.begin();
			t2Index = find (fields.begin(), fields.end(), "totalC_" + itos(j-1)) - fields.begin();
			r2Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(j-1)) - fields.begin();
			//cout << i << "\t" << j << endl;
			//cout << t1Index << "\t" << r1Index << endl;
			//cout << t2Index << "\t" << r2Index << endl;
		}
		else{
			//if(fields[fetPIndex] == "NA"){continue;} //if fetP field is "NA" then next;

			string chr = fields[chrIndex];
			int start = string_to_int(fields[startIndex]);
			int t1 = string_to_int(fields[t1Index]);
			double r1 = atof((fields[r1Index]).c_str());
			int t2 = string_to_int(fields[t2Index]);
			double r2 = atof((fields[r2Index]).c_str());

			if(t1 > config.samplingSaturation){t1=config.samplingSaturation;}
			if(t2 > config.samplingSaturation){t2=config.samplingSaturation;}

			lane[i][chr][start] = cpgMeth(t1, int(t1*r1+0.5));
			lane[j][chr][start] = cpgMeth(t2, int(t2*r2+0.5));
			//cout << chr << "\t" << start << "\t" << t1 << "\t" << int(t1*r1+0.5) << "\t" << t2 << "\t" << int(t2*r2+0.5) << endl;
			//lane[i][chr][start].out();
			//lane[j][chr][start].out();

		}
	}

	inputf.close();

}


int hmm(string & inputFile, string & outputName, HMMCONFIG & conf)
{
	double newTransitionMatrix[NSTATE*NSTATE];
	double transMatrixElemPrec;
	int count=0;

	//Initialize the configurations
//	config.maxIterations = 100;
//	config.minP = 0.95;
//	config.trainPct = 10; //25 optimal; training is about 1/25 of total initial DMRs. ~ 1 chrom ~.
//	config.minDiff = 0.33; //0.3333 optimal
//	config.minDepth = 5; //3 optimal;
//	config.minRegionDist = 300; //500 optimal
//	config.samplingSaturation = 100;
//	config.minDmcs = 3;
//	config.relaxedCondition = 0.2;

	getSiteCount(inputFile);

	printf("reading file completed. regionNum=%d\n", regionNum);

	printf("building lookup table...\n");

	//calculate the emission probability table
	buildLookupTable();

	cout << "Initializing HMM..." << endl;

	//initialize the state of the HMM
	initHMM();


	//train the transition probability table
	//todo: use supplied dmr for training; use specified transition table to calculate hidden states
	while (1)
	{
		cout << "Iteration " << count << endl;
		count++;

		if (count>config.maxIterations)
		{
			break;
		}

		forwardBackward(newTransitionMatrix, config.trainPct);

		transMatrixElemPrec = stopClimb(newTransitionMatrix);

		cout << "Iteration " << count << endl;

		if (transMatrixElemPrec<TRANS_TABLE_ELEM_PREC)
		{
			break;
		}

		memcpy(transitionTable, newTransitionMatrix, NSTATE*NSTATE*sizeof(double));
	}

	//compute the posterior probability of states for all sites
	forwardBackward(newTransitionMatrix, 1);

	//get the DMCs if prob is > 0.95 and get optimal hidden states
	getSignificantSites();
	getOptimalStates();

	outputFile(outputName);

	freeMem();

	return 0;
}


//int main()
//{
//	string inf = "xxx";
//	string ouf = "testx";
//	CONFIG conf;
//	hmm(inf, ouf, conf);
//	return 0;
//}

//read methcall file or dmc file to profile each cytosine for two samples
int getSiteCount(string & inputFile)
{
	int tmpStart, tmpEnd;
	int i,j,k;

	map<int, map<int, int> > selected; //chromIndex->startIndex => 1|0 : selected | not selected for consideration

	cout << "start reading" << endl;
	string fileName = inputFile;
	map <int, map <string, map<int, cpgMeth> > > lane; //lane = laneId -> chrom -> position -> (n,k)
	readCompToHashSimple(fileName, 1, 2, lane ); //read fileName into lane
	cout << "read done" << endl;

	//select sites for consideration by looping through chromIndex->startIndex

	regionNum = 0;
	i = 0;
	for( map<string, map<int, cpgMeth> >::iterator pchr = lane[1].begin(); pchr != lane[1].end(); pchr++, i++){
	//for (i=0;i<chromNum;i++){
		string chr = pchr->first;
		//i = pchr - lane[1].begin(); //error. why?
		chromNameByIndex[i] = chr;

		//read coordinates into loci by chromIndex->startIndex => start
		j = 0;
		for(map<int, cpgMeth>::iterator it = lane[1][chr].begin(); it != lane[1][chr].end(); ++it) {
			//j = it - lane[1][chr].begin(); //why compile error?
			loci[i][j] = it->first;
			j ++;
		}

		int locusSize = lane[1][chr].size();

		//set selected for each site
		//for(map<int, cpgMeth>::iterator it = lane[1][chr].begin(); it != lane[1][chr].end(); ++it) {
		for (j=0;j<locusSize;j++)
		{
			int n1 = lane[1][chr][loci[i][j]].totalC;
			int k1 = lane[1][chr][loci[i][j]].methC;
			int n2 = lane[2][chr][loci[i][j]].totalC;
			int k2 = lane[2][chr][loci[i][j]].methC;

			//worry about direction later //todo
			//determine the initial DMRs based on rough condition
			//(selected = 1) is selected for training, and call DMC/DMR. Neglect sites with selected = 0;
			if ( (abs(double(k1)/n1 - double(k2)/n2) > config.minDiff) && n1 >= config.minDepth && n2 >= config.minDepth )
			{
				selected[i][j] = 1;
			}
			else
			{
				selected[i][j] = 0;
			}
			//cout << chr << "\t" << i << "\t" << j << "\t" << locus[j] << selected[i][j] << endl;
		}

		//merge close regions
		//if distance between two consective Dmcs is larger than minRegionDist, these two are still in one initial region.
		//That's because the CpG before a region is used as null.
		//Further sceening will be applied later for DMR output.
		for (j=1;j<locusSize;j++)
		{
			if ((selected[i][j-1]==1)&&(selected[i][j]==0))
			{
				tmpStart=j;
				for (tmpEnd=tmpStart+1;tmpEnd<locusSize;tmpEnd++)
				{
					if (selected[i][tmpEnd]==1)
					{
						break;
					}
				}
				if ((tmpEnd!=locusSize)&&(loci[i][tmpEnd]-loci[i][tmpStart]<=config.minRegionDist))
				{
					for (k=tmpStart;k<=tmpEnd;k++)
					{
						selected[i][k] = 1;
					}
				}
				j=tmpEnd;
			}
		}

		for (j=0;j<locusSize;j++)
		{
			if (selected[i][j]==1)
			{
				tmpStart = j;

				for (k=tmpStart+1;k<locusSize;k++)
				{
					if (selected[i][k]==0)
					{
						break;
					}
				}

				j=k;

				regionNum++;
			}
		}
	}

	//get initially identified DMRs into regions array by looping through chromIndex->startIndex
	regions = (DMR_REGION *)malloc(regionNum*sizeof(DMR_REGION));
	regionNum = 0;
	i = 0;
	for( map<string, map<int, cpgMeth> >::iterator pchr = lane[1].begin(); pchr != lane[1].end(); pchr++, i++){
	//for (i=0;i<chromNum;i++){
		string chr = pchr->first;
		int locusSize = lane[1][chr].size();

		//merge all consective sites with selected == 1
		for (j=0;j<locusSize;j++)
		{
			if (selected[i][j]==1)
			{
				tmpStart = j;

				for (k=tmpStart+1;k<locusSize;k++)
				{
					if (selected[i][k]==0)
					{
						break;
					}
				}

				j=k; //j is reset here, in addition to the for() loop

				tmpEnd = k-1;

				if (tmpStart>0)
				{
					tmpStart--;
				}

				regions[regionNum].chromIndex=i;
				regions[regionNum].start = tmpStart;
				regions[regionNum].end = tmpEnd;
				regions[regionNum].cxSites = tmpEnd-tmpStart+1;
				regions[regionNum].siteCountsM1 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regions[regionNum].siteCountsM2 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regions[regionNum].siteCountsT1 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regions[regionNum].siteCountsT2 = (int *)malloc(regions[regionNum].cxSites*sizeof(int));

				for (k=0;k<regions[regionNum].cxSites;k++)
				{
					regions[regionNum].siteCountsM1[k] = lane[1][chr][loci[i][tmpStart+k]].methC;
					regions[regionNum].siteCountsM2[k] = lane[2][chr][loci[i][tmpStart+k]].methC;
					regions[regionNum].siteCountsT1[k] = lane[1][chr][loci[i][tmpStart+k]].totalC;
					regions[regionNum].siteCountsT2[k] = lane[2][chr][loci[i][tmpStart+k]].totalC;
				}
				regions[regionNum].prob = (double *)malloc(regions[regionNum].cxSites*NSTATE*sizeof(double));
				regions[regionNum].isSignificant = (int *)malloc(regions[regionNum].cxSites*sizeof(int));
				regionNum++;
			}
		}
	}

	if(0){ //this is for debug

		for (i=0;i<regionNum;i++){
			printf( "%i\t%i\t%i\t%i\t%i\t%i\t%f\t%i\n",regions[i].chromIndex, regions[i].start, regions[i].end, regions[i].cxSites, *(regions[i].siteCountsM1), *(regions[i].siteCountsM2), *(regions[i].prob), *(regions[i].isSignificant));
		}

		ofstream locFile;
		locFile.open("binsInConsideration", ios_base::out);
		//output bin file
		for (i=0;i<regionNum;i++)
		{
			string tmpChrom = chromNameByIndex[regions[i].chromIndex];

			for (j=0;j<regions[i].cxSites;j++)
			{
				int st = loci[regions[i].chromIndex][regions[i].start];
				int en = loci[regions[i].chromIndex][regions[i].end];
				locFile << tmpChrom << "\t";
				locFile << i << "\t" << j << "\t";
				locFile << selected[regions[i].chromIndex][regions[i].start+j] << "\t";
				locFile << st << "\t";
				locFile << en << "\t";
				locFile << loci[regions[i].chromIndex][regions[i].start+j] << "\t";
				locFile << regions[i].siteCountsT1[j] << "\t";
				locFile << regions[i].siteCountsM1[j] << "\t";
				locFile << regions[i].siteCountsT2[j] << "\t";
				locFile << regions[i].siteCountsM2[j] << "\t";


				for (k=0;k<NSTATE;k++)
				{
					//fprintf(fh, "%1.3f\t",regions[i].prob[j*NSTATE+k]);
					locFile << regions[i].prob[j*NSTATE+k] << "\t";

				}
				//fprintf(fh, "%d\n", regions[i].isSignificant[j]);
				locFile << regions[i].isSignificant[j] << endl;
			}
		}

		locFile.close();

	}

	return 0;
}

void freeMem()
{
	int i;

	for (i=0;i<regionNum;i++)
	{
		free(regions[i].siteCountsM1);
		free(regions[i].siteCountsM2);
		free(regions[i].prob);
		free(regions[i].isSignificant);
	}

	free(regions);
}

//calculate the emission probability for each possible Multikey
int buildLookupTable()
{
	int i,j,k,l;
	map <MultiKey, int> marked;

	double sum, sum1,p;
	double alpha1,beta1, alpha2,beta2;
	double p1[100], p2[100], p1Prob[100],p2Prob[100], upperBound, lowerBound, mean, stdev, max;
	double denominator[3];

	//parameters for prior beta distribution
	alpha1 = 1.0;
	alpha2 = 1.0;
	beta1 = 1.0;
	beta2 = 1.0;


	for (i=0;i<regionNum;i++)
	{
		for (j=0;j<regions[i].cxSites;j++)
		{
			MultiKey obs(regions[i].siteCountsT1[j], regions[i].siteCountsM1[j], regions[i].siteCountsT2[j], regions[i].siteCountsM2[j]);
			marked[obs] = 1;
		}
	}

	//computing the integral in denominator
	mean = (alpha1)/(alpha1+beta1);
	stdev = sqrt(alpha1*beta1/(alpha1+beta1)/(alpha1+beta1)/(alpha1+beta1+1));

	lowerBound = mean-10*stdev<0?0:mean-10*stdev;
	upperBound = mean+10*stdev>1?1:mean+10*stdev;

	p1[0] = lowerBound+(upperBound-lowerBound)/200;

	for (i=1;i<100;i++)
	{
		p1[i] = p1[i-1]+(upperBound-lowerBound)/100;
	}

	max = -DBL_MAX;

	for (i=0;i<100;i++)
	{
		p1Prob[i] = (alpha1-1)*log(p1[i])+(beta1-1)*log(1-p1[i]);

		if (p1Prob[i]>max)
		{
			max = p1Prob[i];
		}
	}

	sum = 0;

	for (i=0;i<100;i++)
	{
		p1Prob[i] = exp(p1Prob[i]-max);

		sum+=p1Prob[i];
	}

	for (i=0;i<100;i++)
	{
		p1Prob[i] /= sum;
	}

	for (i=0;i<3;i++)
	{
		denominator[i] = 0;
	}

	for (i=0;i<100;i++)
	{
		for (l=0;l<100;l++)
		{
			if ((p1[i]-p1[l]>config.minDiff))
			{
				denominator[2]+=p1Prob[i]*p1Prob[l];
			}
			else if ((p1[l]-p1[i]>config.minDiff))
			{
				denominator[0]+=p1Prob[i]*p1Prob[l];
			}
			else
			{
				denominator[1]+=p1Prob[i]*p1Prob[l];
			}
		}
	}
	//cout << "denom: " << denominator[0] << "\t" << denominator[1] << "\t" << denominator[2] << endl;

	//calculate the emission probability
	for (int JJ=1; JJ<config.samplingSaturation+1; JJ++){
	for (j=0;j<=JJ;j++)
	{
		for (int KK=1; KK<config.samplingSaturation+1; KK++){
		for (k=0;k<=KK;k++)
		{
			MultiKey combi(JJ, j, KK, k);
			if (!marked[combi])
			{
				continue;
			}

			mean = (j+alpha1)/(JJ+alpha1+beta1);
			stdev = sqrt((j+alpha1)*(JJ-j+beta1)/(JJ+alpha1+beta1)/(JJ+alpha1+beta1)/(JJ+alpha1+beta1+1));

			lowerBound = mean-10*stdev<0?0:mean-10*stdev;
			upperBound = mean+10*stdev>1?1:mean+10*stdev;

			p1[0] = lowerBound+(upperBound-lowerBound)/200;

			for (i=1;i<100;i++)
			{
				p1[i] = p1[i-1]+(upperBound-lowerBound)/100;
			}

			max = -DBL_MAX;

			for (i=0;i<100;i++)
			{
				p1Prob[i] = (j+alpha1-1)*log(p1[i])+(JJ-j+beta1-1)*log(1-p1[i]);

				if (p1Prob[i]>max)
				{
					max = p1Prob[i];
				}
			}

			sum = 0;

			for (i=0;i<100;i++)
			{
				p1Prob[i] = exp(p1Prob[i]-max);

				sum+=p1Prob[i];
			}

			for (i=0;i<100;i++)
			{
				p1Prob[i] /= sum;
			}

			mean = (k+alpha2)/(KK+alpha2+beta2);
			stdev = sqrt((k+alpha2)*(KK-k+beta2)/(KK+alpha2+beta2)/(KK+alpha2+beta2+1)/(KK+alpha2+beta2+1));

			lowerBound = mean-10*stdev<0 ? 0:mean-10*stdev;
			upperBound = mean+10*stdev>1 ? 1:mean+10*stdev;

			p2[0] = lowerBound+(upperBound-lowerBound)/200;

			for (i=1;i<100;i++)
			{
				p2[i] = p2[i-1]+(upperBound-lowerBound)/100;
			}

			max = -DBL_MAX;

			for (i=0;i<100;i++)
			{
				p2Prob[i] = (k+alpha2-1)*log(p2[i])+(KK-k+beta2-1)*log(1-p2[i]);

				if (p2Prob[i]>max)
				{
					max = p2Prob[i];
				}
			}

			sum = 0;

			for (i=0;i<100;i++)
			{
				p2Prob[i] = exp(p2Prob[i]-max);

				sum+=p2Prob[i];
			}

			for (i=0;i<100;i++)
			{
				p2Prob[i] /= sum;
			}

			for (i=0;i<3;i++)
			{
				emmTable[i][combi] = 0;

			}

			for (i=0;i<100;i++)
			{
				for (l=0;l<100;l++)
				{
					if   ((p1[i]-p2[l]>config.minDiff))
					{
						//emmTable[2][combi] += p1Prob[i]*p2Prob[l]/denominator[2];
						emmTable[2][combi] += p1Prob[i]*p2Prob[l];
					}
					else if   ((p2[l]-p1[i]>config.minDiff))
					{
						//emmTable[0][combi] += p1Prob[i]*p2Prob[l]/denominator[0];
						emmTable[0][combi] += p1Prob[i]*p2Prob[l];
					}
					else
					{
						//emmTable[1][combi] += p1Prob[i]*p2Prob[l]/denominator[1];
						emmTable[1][combi] += p1Prob[i]*p2Prob[l];
					}
				}
			}

			if(0){ //for debug
				printf("%d\t%d\t%d\t%d\t", JJ, j, KK, k);
				for (i=0;i<NSTATE;i++){
					printf( "%f\t",emmTable[i][combi]);
				}
				printf("\n");
			}
		}
		}
	}
	}
	//cout << "normalization" << endl;
	//normalization
	for (int JJ=1; JJ<config.samplingSaturation+1; JJ++){
	for (j=0;j<=JJ;j++)
	{

		for (int KK=1; KK<config.samplingSaturation+1; KK++){
		for (k=0;k<=KK;k++)
		{

			MultiKey combi(JJ, j, KK, k);
			if (!marked[combi])
			{
				continue;
			}

			sum = 0;

			for (i=0;i<NSTATE;i++)
			{
				sum+=emmTable[i][combi];
			}

			for (i=0;i<NSTATE;i++)
			{
				emmTable[i][combi]/=sum;
			}

			if(0){
				printf("%d\t%d\t%d\t%d\t", JJ, j, KK, k);
				for (i=0;i<NSTATE;i++){
					printf( "%f\t",emmTable[i][combi]);
				}
				printf("\n");
			}

		}
		}
	}
	}

	//cout << "return.." << endl;
	//exit;
	return 1;
}

//initialize HMM
void initHMM()
{
	int i,j;

	priorProb[0] = 0.0;
	priorProb[1] = 1.0;
	priorProb[2] = 0.0;

	for (i=0;i<NSTATE;i++)
	{
		for (j=0;j<NSTATE;j++)
		{
			transitionTable[i*NSTATE+j]=1.0/NSTATE;
		}
	}
}

//Forward-Backward algorithm for training transition table. step is used for selecting a subset of regions for training
void forwardBackward(double *newTransition, int step)
{
	double *alpha, *beta, *eps;
	double sum, sum1;
	double priorSum[NSTATE];
	double transitionSum[NSTATE*NSTATE],postProbSum[NSTATE*NSTATE];
	int maxLen;
	int i,j,k,l;

	//allocate memory and initializing arrays
	memset(priorSum,0,NSTATE*sizeof(double));
	memset(transitionSum, 0, NSTATE*NSTATE*sizeof(double));
	memset(postProbSum, 0 , NSTATE*NSTATE*sizeof(double));

	maxLen = 0;

	for (i=0;i<regionNum;i+=step)
	{
		if (regions[i].cxSites>maxLen)
		{
			maxLen = regions[i].cxSites;
		}
	}

	alpha = (double *)malloc(NSTATE*maxLen*sizeof(double));
	beta = (double *)malloc(NSTATE*maxLen*sizeof(double));
	eps = (double *)malloc(NSTATE*NSTATE*maxLen*sizeof(double));

	for (l=0;l<regionNum;l+=step)
	{
		//Forward algorithm

		MultiKey obs(regions[l].siteCountsT1[0], regions[l].siteCountsM1[0], regions[l].siteCountsT2[0], regions[l].siteCountsM2[0]);
		sum=0;
		for (k=0;k<NSTATE;k++)
		{
			//alpha[k] = emmissionProbTable[k][ regions[l].siteCountsM1[0] ][ regions[l].siteCountsM2[0] ]*priorProb[k];
			alpha[k] = emmTable[k][ obs ]*priorProb[k];
			//cout << "region=" << l << " state=" << k << " alpha[k]=" << alpha[k] ;
			//cout << " emmTable[k][ obs ]*priorProb[k]=" << emmTable[k][ obs ] << " * " << priorProb[k] << endl;
		}

		for (i=1;i<regions[l].cxSites;i++)
		{
			sum1 = 0;

			for (j=0;j<NSTATE;j++)
			{
				sum = 0;

				for (k=0;k<NSTATE;k++)
				{
					sum+=alpha[(i-1)*NSTATE+k]*transitionTable[k*NSTATE+j];
				}

				//alpha[i*NSTATE+j] = sum*emmissionProbTable[j][ regions[l].siteCountsM1[i] ][ regions[l].siteCountsM2[i] ];
				MultiKey obs1(regions[l].siteCountsT1[i], regions[l].siteCountsM1[i], regions[l].siteCountsT2[i], regions[l].siteCountsM2[i]);
				alpha[i*NSTATE+j] = sum*emmTable[j][obs1];
				sum1+=alpha[i*NSTATE+j];
				//cout << "region=" << l << " state=" << i << "," << j << " alpha[i*NSTATE+j]=" << alpha[i*NSTATE+j] << endl;
			}

			//cout << "scale:" << endl;
			for (j=0;j<NSTATE;j++)
			{
				alpha[i*NSTATE+j]/=sum1;
				//cout << "region=" << l << " state=" << i << "," << j << " alpha[i*NSTATE+j]=" << alpha[i*NSTATE+j] << endl;
			}
		}

		//Backward algorithm

		for (j=0;j<NSTATE;j++)
		{
			beta[(regions[l].cxSites-1)*NSTATE+j] = 1.0/NSTATE;
		}

		for (i=regions[l].cxSites-2;i>=0;i--)
		{
			sum1 = 0;

			for (j=0;j<NSTATE;j++)
			{
				sum = 0;

				for (k=0;k<NSTATE;k++)
				{
					//sum+=beta[(i+1)*NSTATE+k]*transitionTable[j*NSTATE+k]*emmissionProbTable[k][ regions[l].siteCountsM1[i+1] ][ regions[l].siteCountsM2[i+1] ];
					MultiKey obs1(regions[l].siteCountsT1[i+1], regions[l].siteCountsM1[i+1], regions[l].siteCountsT2[i+1], regions[l].siteCountsM2[i+1] );
					sum+=beta[(i+1)*NSTATE+k]*transitionTable[j*NSTATE+k]*emmTable[k][obs1];
				}

				beta[i*NSTATE+j]=sum;
				sum1+=beta[i*NSTATE+j];
			}

			for (j=0;j<NSTATE;j++)
			{
				beta[i*NSTATE+j]/=sum1;
			}
		}

		//compute posterior probability

		for (i=0;i<regions[l].cxSites;i++)
		{
			sum1=0;
			for (j=0;j<NSTATE;j++)
			{
				regions[l].prob[i*NSTATE+j]=alpha[i*NSTATE+j]*beta[i*NSTATE+j];
				sum1+=regions[l].prob[i*NSTATE+j];
			}

			for (j=0;j<NSTATE;j++)
			{
				regions[l].prob[i*NSTATE+j]/=sum1;
			}
		}

		//compute eps

		for (i=0;i<regions[l].cxSites-1;i++)
		{
			sum1=0;
			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					//eps[i*NSTATE*NSTATE+j*NSTATE+k]=alpha[i*NSTATE+j]*beta[(i+1)*NSTATE+k]*emmissionProbTable[k][ regions[l].siteCountsM1[i+1] ][ regions[l].siteCountsM2[i+1] ]*transitionTable[j*NSTATE+k];
					MultiKey obs1(regions[l].siteCountsT1[i+1], regions[l].siteCountsM1[i+1], regions[l].siteCountsT2[i+1], regions[l].siteCountsM2[i+1] );
					eps[i*NSTATE*NSTATE+j*NSTATE+k]=alpha[i*NSTATE+j]*beta[(i+1)*NSTATE+k]*emmTable[k][obs1]*transitionTable[j*NSTATE+k];
					sum1+=eps[i*NSTATE*NSTATE+j*NSTATE+k];
				}
			}

			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					eps[i*NSTATE*NSTATE+j*NSTATE+k]/=sum1;

				}
			}
		}

		//update the cumulation for caluculating transition probability table

		for (i=0;i<regions[l].cxSites-1;i++)
		{
			for (j=0;j<NSTATE;j++)
			{
				for (k=0;k<NSTATE;k++)
				{
					transitionSum[j*NSTATE+k]+=eps[i*NSTATE*NSTATE+j*NSTATE+k];
					postProbSum[j*NSTATE+k]+=regions[l].prob[i*NSTATE+j];
				}
			}
		}

	}

	//update transition table
	for (j=0;j<NSTATE;j++)
	{
		for (k=0;k<NSTATE;k++)
		{
			newTransition[j*NSTATE+k] = transitionSum[j*NSTATE+k]/postProbSum[j*NSTATE+k];
		}
	}

	free(alpha);
	free(beta);
	free(eps);
}

//calculate the biggest difference for transition matrix elements
double stopClimb(double *newTransition)
{
	int i,j;
	double maxDiff = 0.0;

	for (i=0;i<NSTATE;i++)
	{

		for (j=0;j<NSTATE;j++)
		{
			double diff = abs(newTransition[i*NSTATE+j] - transitionTable[i*NSTATE+j]);
			if(diff > maxDiff){
				maxDiff = diff;
			}
		}
	}

	return maxDiff;
}

//output files
void outputFile(string & outputName)
{

	int i,j,k;
	string tmpChrom;
	int start,end;

	//get file names
	string locFileName = outputName + ".loc";
	string hmmFileName = outputName + ".hmm";
	string regionFileName = outputName + ".region";
	string dmrFileName = outputName + ".dmr";

	ofstream locFile;
	locFile.open(locFileName.c_str(), ios_base::out);


	//output loc file
	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{
			int st = loci[regions[i].chromIndex][regions[i].start];
			int en = loci[regions[i].chromIndex][regions[i].end];
			locFile << tmpChrom << "\t";
			locFile << st << "\t";
			locFile << en << "\t";
			locFile << loci[regions[i].chromIndex][regions[i].start+j] << "\t";
			locFile << regions[i].siteCountsT1[j] << "\t";
			locFile << regions[i].siteCountsM1[j] << "\t";
			locFile << regions[i].siteCountsT2[j] << "\t";
			locFile << regions[i].siteCountsM2[j] << "\t";


			for (k=0;k<NSTATE;k++)
			{
				locFile << regions[i].prob[j*NSTATE+k] << "\t";
			}

			MultiKey combi(regions[i].siteCountsT1[j], regions[i].siteCountsM1[j], regions[i].siteCountsT2[j], regions[i].siteCountsM2[j]);
			for (k=0;k<NSTATE;k++)
			{
				locFile << emmTable[k][combi] << "\t";
			}
			locFile << State[i][j] << "\t";
			locFile << regions[i].isSignificant[j] << endl;
		}
	}
	locFile.close();

	//output regions formed by DMCs only
	ofstream regionFile;
	regionFile.open(regionFileName.c_str(), ios_base::out);

	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{
			if (regions[i].isSignificant[j])
			{
				start = j;
			}
			else
			{
				continue;
			}

			for (k=start+1;k<regions[i].cxSites;k++)
			{
				if (!(regions[i].isSignificant[k]==regions[i].isSignificant[k-1]))
				{
					break;
				}
			}

			end = k;

			regionFile << tmpChrom << "\t";
			regionFile << loci[regions[i].chromIndex][regions[i].start+start] << "\t";
			regionFile << loci[regions[i].chromIndex][regions[i].start+end] << "\t";
			regionFile << (regions[i].isSignificant[j]==1?"+":"-") << endl;


			j=k-1;
		}
	}

	regionFile.close();


	//output DMRs by joining consective differential states, allowing 20% sites (around 1 or 2 cpgs) as middle states
	ofstream dmrFile;
	dmrFile.open(dmrFileName.c_str(), ios_base::out);

	for (i=0;i<regionNum;i++)
	{
		tmpChrom = chromNameByIndex[regions[i].chromIndex];

		for (j=0;j<regions[i].cxSites;j++)
		{

			if (State[i][j] != 0)
			{
				start = j;
			}
			else
			{
				continue;
			}
			int dmcCount = 1;
			int cpgCount = 0;
			int type = State[i][j];
			int lastDmc = j;
			for (k=start+1;k<regions[i].cxSites;k++)
			{

				//int relaxedCondition = 1;//make it an option todo
				if(config.relaxedCondition > 0){
					if( (regions[i].prob[(k-1)*NSTATE+0] - regions[i].prob[(k-1)*NSTATE+2]) * (regions[i].prob[k*NSTATE+0] - regions[i].prob[k*NSTATE+2]) < 0  ){
					//if( State[i][k] == -1 * type ){ //if State[i][k] == 0, it's fine; Stops only if a reversed dmc is met.
						cpgCount -= k - lastDmc - 1;
						k = lastDmc + 1;
						//dmrFile << k << endl;
						break;
					}
				} else {
					if (! (State[i][k] == State[i][k-1]) )
					{
						break;
					}
				}

				if(State[i][k]==type){
					if(loci[regions[i].chromIndex][regions[i].start + k] - loci[regions[i].chromIndex][regions[i].start + lastDmc] > config.minRegionDist){
						cpgCount -= k - lastDmc - 1;
						k = lastDmc + 1;
						break;
					}
					lastDmc = k;
					dmcCount ++;
				} else if(State[i][k]==0){
					cpgCount ++;
					if( config.relaxedCondition > 0 && cpgCount >int( dmcCount * config.relaxedCondition + 1 ) ){
						cpgCount -= k - lastDmc - 1;
						k = lastDmc + 1;
						break;
					}
				}

			}
			if(k== regions[i].cxSites ){
				cpgCount -= k - lastDmc - 1;
				k = lastDmc + 1;
			}

			end = k - 1;
			//dmrFile << end << endl;
			if(dmcCount >=config.minDmcs){ //3 dmc and 500bp seem to be good parameters
				dmrFile << tmpChrom << "\t";
				dmrFile << loci[regions[i].chromIndex][regions[i].start+start] << "\t";
				dmrFile << loci[regions[i].chromIndex][regions[i].start+end] + 2 << "\t";
				dmrFile << (regions[i].isSignificant[j]==1?"+":"-") << "\t";
				dmrFile << end -start + 1 << "\t" << dmcCount << "\t" << cpgCount << "\t";
				dmrFile << State[i][j] << endl;
			}

			j=k-1;
		}
	}

	dmrFile.close();


	//output the transition probability table
	ofstream hmmFile;
	hmmFile.open(hmmFileName.c_str(), ios_base::out);
	hmmFile << "transition matrix:" << endl;
	for (i=0;i<NSTATE;i++)
	{
		for (j=0;j<NSTATE;j++)
		{
			hmmFile << transitionTable[i*NSTATE+j] << "\t";
		}
		hmmFile << endl;
	}
	hmmFile.close();
}

// get DMCs if prob > minP(0.95)
void getSignificantSites()
{
	int i,j;
	double sum;

	for (i=0;i<regionNum;i++)
	{
		for (j=0;j<regions[i].cxSites;j++)
		{
			sum = 0;

			if (regions[i].prob[j*NSTATE]>config.minP)
			{
				regions[i].isSignificant[j] = -1;
			}
			else if (regions[i].prob[j*NSTATE+2]>config.minP)
			{
				regions[i].isSignificant[j] = 1;
			}
			else
			{
				regions[i].isSignificant[j] = 0;
			}
		}
	}
}

//get optimal state sequence
void getOptimalStates()
{
	int i,j;
	double sum;

	for (i=0;i<regionNum;i++)
	{
		for (j=0;j<regions[i].cxSites;j++)
		{
			sum = 0;

			double s0 = regions[i].prob[j*NSTATE];
			double s1 = regions[i].prob[j*NSTATE+1];
			double s2 = regions[i].prob[j*NSTATE+2];
			if(s0 > s1 && s0 > s2){
				State[i][j] = -1;
			} else if (s2 > s0 && s2 > s1){
				State[i][j] = 1;
			} else {
				State[i][j] = 0;
			}

			if (regions[i].prob[j*NSTATE]>config.minP)
			{
				regions[i].isSignificant[j] = -1;
			}
			else if (regions[i].prob[j*NSTATE+2]>config.minP)
			{
				regions[i].isSignificant[j] = 1;
			}
			else
			{
				regions[i].isSignificant[j] = 0;
			}
		}
	}
}
