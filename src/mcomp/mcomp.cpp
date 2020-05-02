
//need complete all the methods for pdiff, pdiffCI and etc.
//bug: when threads >1, \t\t shows up.
//TODO: times options
//TODO: qvalue
//TODO: reproducibility check for replicates
//TODO: construct tag/header data so that the tags are not hard coded or created/recreaded again.//replace i/j with lane names
//TODO: This is checking epimutation on C and you could make it universal for use of checking mutation on other nucleotides.
//TODO: statistics when merging ratio files
#define ALPHA 0.05
#define BATCHMAX 2000
#define BUFFSIZE 2048000                          // 2000 batch size, each with 1024 bytes for current replicates
//#include <RInside.h>                            // for the embedded R via RInside
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iterator>
//#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
//#include <sys/auxv.h>

//#include "nr.h"
//#include "nr3.h"
//#include "adapt.h"
//#include "qgaus.h"
#include "types.h"
#include "stats.h"
#include "hmm.h"
#include "lut.h"
#include "fisher_exact_test.h"
#include "fet2x2.h"
#include "bbf.h"

//#include "pdiffCI.h"
//#include "sci.h" //do not understand why it generates error undefined reference to `singleCI(int, int, double, int)'
//#include "gamma.h"
//#include "incgammabeta.h"

//#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
////#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
//#include <boost/date_time.hpp>

// For boost::process, not working in CentOS due to no automatic pipe closing
// #include <boost/process.hpp>
// #include <boost/regex.hpp>
// #include <boost/asio/io_service.hpp>
// #include <future>
// #include <sstream>

#include <math.h> //round() function

// For two-way pipe
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <sys/wait.h>

namespace po = boost::program_options;
using namespace std;


// Moved to `types.h` and `types.cpp`
// std::string do_readlink(std::string const& path) {
//     char buff[1024];
//     ssize_t len = ::readlink(path.c_str(), buff, sizeof(buff)-1);
//     if (len != -1) {
//       buff[len] = '\0';
//       return std::string(buff);
//     } else {
//      /* handle error condition */
//     }
// }


//RInside R;
//SEXP ans;

//global variable opts
class Opts //Options
{
public:
	string email;
	vector <string> ratiosFiles;
	vector <string> mergedRatiosFiles;
	vector <string> labels;
	string outputDir;
	string webOutputDir;
	string compFile;
	string inGenome;
	string reference;
	int precision;
	int threads;
	int lmFit;
	int doMergeRatioFiles;
	int doStrandSpecifiMeth;
	int doComp;
//	vector <double> xVector;
	int mergeNotIntersect;
	int withVariance;
	int minDepthForComp;
	int doDmcScan;
	int doDmrScan;
	int dmrMethods;
	int dopredefinedFeature;
	double filterCredibleDif;
	double pFetDmc;
	double pFetDmr;
	double minNominalDif;
	double pSimDmc;
	double pSimDmr;
	double minCredibleDif;
	double topRankByCDif;
	double topRankByPSim;
	double minDmcsInDmr;
	double maxDistConsDmcs;
	vector <string> predefinedFeature;

	Opts(){
		email = "";
		outputDir = "./";
		webOutputDir = "./";
		compFile = "mSuite";
		inGenome = "";
		reference = "";
		precision = 3;
		threads = 1;
		lmFit = 1;
		mergeNotIntersect = 1; //for 2 samples, this option defaults at 0(intersect); for 3 samples it defaults to 1(merge)//TODO this option
		doMergeRatioFiles = 0;
		doStrandSpecifiMeth = 0;
		doComp = 1;
		minDepthForComp = 3;
		doDmcScan = 1;
		doDmrScan = 1;
		dmrMethods = 7;
		filterCredibleDif = -10.0;
		pFetDmc = 0.05;
		pFetDmr = 0.05;
		minNominalDif = 0.3333;
		pSimDmc = 1;
		pSimDmr = 1;
		minCredibleDif = 0.2;
		topRankByCDif = 1;//not yet implemented
		topRankByPSim = 1;//not yet implemented
		minDmcsInDmr = 3;
		maxDistConsDmcs = 300;
	};
	void out(){
		cout << endl;
		cout << "Printing options for compFile, threads, lmFit, mergeNotIntersect " << endl;
		cout << "labels, size of labels" << endl;
		cout << "ratioFiles, size of ratioFiles" << endl;
		cout << "doMergeRatioFiles, doStrandSpecificMeth, doComp" << endl;
		cout << compFile << "\t" << threads << "\t" << lmFit << "\t" << mergeNotIntersect << endl;
//		copy(xVector.begin(), xVector.end(), ostream_iterator<double>(cout, ",")); cout << "size of xVector=" << xVector.size() << endl;
		copy(labels.begin(), labels.end(), ostream_iterator<string>(cout, "\t")); cout << "size of labels=" << labels.size() << endl;
		copy(ratiosFiles.begin(), ratiosFiles.end(), ostream_iterator<string>(cout, "\t")); cout << "size of ratiosFiles=" << ratiosFiles.size() << endl;
		cout << doMergeRatioFiles << "\t" << doStrandSpecifiMeth << "\t" << doComp << "\t" << endl;
	}
} opts;


int parse_options(int ac, char * av[]){

	po::variables_map options;

	po::options_description desc("Allowed options for methComp");
	desc.add_options()
	("help,h", 								"produce help message;")
	("email",								po::value<string>(),"Specify email;")
	("ratiosFiles,r",						po::value< vector<string> >()->multitoken(), "Specify the names of ratio files from methCall. Multiple lane files can be separated by , to be combined into a single track; example: -r sample1 -r sample2 -r s_r1,s_r2,s_r3;")
	("mergedRatiosFiles,m",					po::value< vector<string> >()->multitoken(), "If --ratiosFiles is ',' separated, then this option must be set;")
	("labels,l",							po::value< vector<string> >()->multitoken(), "Name labels for samples, defaut 0, 1, ...;")
	("outputDir",							po::value<string>(), "Specify the name of the output directory;")
	("webOutputDir",						po::value<string>(), "Specify the name of the web-accessible output directory for UCSC Genome Browser tracks;")
	("compFile,c", 							po::value<string>(), "Name of the comparison file resulted from statistical tests;")
	("inGenome", 							po::value<string>(), "Specify the UCSC Genome Browser identifier of source genome assembly;")
	("reference", 							po::value<string>(), "Specify the path to the reference genome for example mm9.fa; mm9.chrom.sizes must be in the same dir;")
//	("xVector", 							po::value<string>(), "Specify the x vector for R.lm() function;x is comma(,) separated float numbers; default, 1.0,2.0,...,;")
	("precision", 							po::value<int>()->default_value(3), "Specify the precision of float numbers in output files (default: 3);")
	("threads,p", 							po::value<int>()->default_value(6), "Specify number of threads; suggest number 6-12; default 6;")
	("lmFit", 								po::value<int>()->default_value(1), "Specify if lenear model fitting is performed; default true; Note that 'na' is generated if slope is 0;")
	("mergeNotIntersect", 					po::value<int>()->default_value(1), "Specify if genomic locations are merged or intersected among samples; 1 for merge(default) and 0 for intersect;")
	("withVariance", 						po::value<int>()->default_value(0), "Specify if there's individual biological variance among the same condition; default 0; Should be 0 for most animal models 1 for most patient studies; WithVariance=1 is not effective if only 1 or 2 replicates.")
	("doMergeRatioFiles", 					po::value<int>()->default_value(0), "Internal parameter. Is true when -m parameter is ',' separated and program will merge ratio Files that are separated by ',' and the output files are named according to option -x;")
	("doStrandSpecifiMeth", 				po::value<int>()->default_value(0), "whether strand specific methylation analysis will be performed;")
	("doComp", 								po::value<int>()->default_value(1), "doComp;")
	("minDepthForComp,d", 					po::value<int>()->default_value(3), "If a site has depth < d then this site is ignored for statistical tests; This option affects much of nominal ratios but none of credible ratios; Suggest 10 for method 2 and 3 for method 2; You may also reset this option during later DMC/DMR rescan to filter sites with depth < d;")
	("----------",							"Below are options for Dmc and Dmr scan;")
	("doDmcScan", 							po::value<int>()->default_value(1), "doDmcScan;")
	("doDmrScan", 							po::value<int>()->default_value(1), "doDmrScan;")

	("filterCredibleDif", 					po::value<double>()->default_value(-10), "if absolute value of cDif for a site < filterCredibleDif, then this site is ignored for regional calculation. use 0.01(for example) to filter all sites with no difference; use 0.20(for example) to select DMCs; Any negative number = no filter;")
	("dmrMethods", 							po::value<int>()->default_value(7), "dmrMethods: add 2^x  method x; examples: 7 for three methods, 4 for method 3 only;")
	("pFetDmc", 							po::value<double>()->default_value(0.05, "0.05"), "Cutoff of P value from Fisher Exact Test for Dmc scan;")
	("pFetDmr", 							po::value<double>()->default_value(0.05, "0.05"), "Cutoff of P value from Fisher Exact Test for Dmr scan;")
	("minNominalDif", 						po::value<double>()->default_value(0.3333, "0.3333"), "min nominal meth diff for Dmc Dmr;")
	("pSimDmc", 							po::value<double>()->default_value(1.0, "1.0"), "Cutoff P value from Similarity Test for Dmc scan; Since p is alwasys less than 1, default 1 means not a criteria;")
	("pSimDmr", 							po::value<double>()->default_value(1.0, "1.0"), "Cutoff P value from Simlarity Test for Dmr scan;")
	("minCredibleDif", 						po::value<double>()->default_value(0.2, "0.2"), "min credible meth diff for Dmc calling, used in M2 or predefined regions;")
	("topRankByCDif", 						po::value<double>()->default_value(1.0, "1.0"), "filter Dmc by asking it to be in top (default 100%) percent by ranking absolute value of credibleDif; suggest 0.05 as the only condition to call Dmc if cDif condition is not prefered; The cutoff cDif will be used as Dmr criteria;")
	("topRankByPSim", 						po::value<double>()->default_value(1.0, "1.0"), "filter Dmc by asking it to be in top (default 100%) percent by ranking P value from Similarity Test;")
	("minDmcsInDmr", 						po::value<int>()->default_value(3), "minimum number of Dmcs in a Dmr;")
	("maxDistConsDmcs", 					po::value<int>()->default_value(300), "max distance between two consective Dmcs for them to be considered in a Dmr;")
	("predefinedFeature,f", 				po::value< vector<string> >()->multitoken(), "supply bed files as predefined feature; -f promoter.bed -f CpgIsland.bed -f LINE.bed is same as -f promoter.bed,CpgIsland.bed,Line.bed")
	;

	po::store(po::parse_command_line(ac, av, desc), options);
	po::notify(options);

	//////////////////////////////////
	if (options.count("help") || ac == 1) {
		cout << desc << endl;
		exit(1);
	}

	//////////////////////////////////
	cout<<"Options are saved in file run.config and printed here:"<<endl;
	ofstream configFile;
	//configFile.open("run.config", ios_base::app);
	stringstream configFileName;
	configFileName << "run.config." << getpid();
	configFile.open(configFileName.str().c_str(), ios_base::app);

	for(std::map<string,po::variable_value>::iterator iter = options.begin(); iter != options.end(); ++iter)
	{
		string k =  (*iter).first;

		cout		<<k<<"=";
		configFile 	<<k<<"=";


		if( k == "email"){
			opts.email 	= 	options[k].as<string>();
			configFile 	<<	options[k].as<string>();
			cout 		<<	options[k].as<string>();
		}
		else if( k == "ratiosFiles"){
			opts.ratiosFiles = options[k].as< vector<string> >();
			for(unsigned int i = 0 ; i < opts.ratiosFiles.size(); i++)
			{
				cout 		<< opts.ratiosFiles[i]<<" ";//copy(myvector.begin(), myvector.end(), ostream_iterator<int>(cout, "\n"));
				configFile 	<< opts.ratiosFiles[i]<<" ";

				if( find(opts.ratiosFiles[i].begin(),opts.ratiosFiles[i].end(), ',') != opts.ratiosFiles[i].end() ){
					opts.doMergeRatioFiles = 1;
					vector <string> toMergeFiles;
					boost::split(toMergeFiles, opts.ratiosFiles[i], boost::is_any_of(","));
					for(unsigned int j = 0; j < toMergeFiles.size(); j++ ){
						if ( !boost::filesystem::exists( toMergeFiles[j] ) )
						{
						  cout << endl << "Can't find file " << toMergeFiles[j] << endl;
						  exit(1);
						}
					}
				} else {
					if ( !boost::filesystem::exists( opts.ratiosFiles[i] ) )
					{
					  cout << endl << "Can't find file " << opts.ratiosFiles[i] << endl;
					  exit(1);
					}
				}
			}
		}
		else if( k == "mergedRatiosFiles"){
			opts.mergedRatiosFiles = options[k].as< vector<string> >();
			for(unsigned int i = 0 ; i < opts.mergedRatiosFiles.size(); i++)
			{
				cout 		<< opts.mergedRatiosFiles[i]<<" ";//copy(myvector.begin(), myvector.end(), ostream_iterator<int>(cout, "\n"));
				configFile 	<< opts.mergedRatiosFiles[i]<<" ";

				if ( boost::filesystem::exists( opts.mergedRatiosFiles[i] ) )
				{
				  cout << endl << "Already existing: " << opts.mergedRatiosFiles[i] << endl;
				  exit(1);
				}

			}
		}
		else if( k == "labels"){
			opts.labels = options[k].as< vector<string> >();
			for(unsigned int i = 0 ; i < opts.labels.size(); i++)
			{
				cout 		<< opts.labels[i]<<" ";//copy(myvector.begin(), myvector.end(), ostream_iterator<int>(cout, "\n"));
				configFile 	<< opts.labels[i]<<" ";
			}
		}
//		else if( k == "xVector"){
//			string xVector = options[k].as<string>();
//			vector <string> fields;
//			boost::split(fields, xVector, boost::is_any_of(","));
//			vector<double>v;
//			for(unsigned int i = 0 ; i < fields.size(); i++)
//			{
//				v.push_back( atof( (fields[i]).c_str() ) );
//			}
//			opts.xVector = v;
//		}
		else if( k == "outputDir"){
			opts.outputDir 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if(k== "compFile"){
			opts.compFile 		= 	options[k].as<string>();
			configFile 			<<	options[k].as<string>();
			cout 				<<	options[k].as<string>();
		}
		else if( k == "inGenome"){
			opts.inGenome 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k == "reference"){
			opts.reference 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k == "precision"){
			opts.precision 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "threads"){
			opts.threads 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}

		else if( k == "lmFit"){
			opts.lmFit 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "withVariance"){
			opts.withVariance 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "mergeNotIntersect"){
			opts.mergeNotIntersect 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "doMergeRatioFiles"){
			opts.doMergeRatioFiles 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "doStrandSpecifiMeth"){
			opts.doStrandSpecifiMeth 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "doComp"){
			opts.doComp 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "dmrMethods"){
			opts.dmrMethods 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}

		else if( k == "minDepthForComp"){
			opts.minDepthForComp 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "doDmcScan"){
			opts.doDmcScan 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "doDmrScan"){
			opts.doDmrScan 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "filterCredibleDif"){
			opts.filterCredibleDif 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "pFetDmc"){
			opts.pFetDmc 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "pFetDmr"){
			opts.pFetDmr 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "minNominalDif"){
			opts.minNominalDif 		= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}

		else if( k == "pSimDmc"){
			opts.pSimDmc 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "pSimDmr"){
			opts.pSimDmr 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "minCredibleDif"){
			opts.minCredibleDif 	= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "topRankByCDif"){
			opts.topRankByCDif 			= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}
		else if( k == "topRankByPSim"){
			opts.topRankByPSim		= 	options[k].as<double>();
			configFile 				<<	options[k].as<double>();
			cout 					<<	options[k].as<double>();
		}

		else if( k == "minDmcsInDmr"){
			opts.minDmcsInDmr 		= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "maxDistConsDmcs"){
			opts.maxDistConsDmcs 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "predefinedFeature"){
			opts.predefinedFeature 	= options[k].as< vector<string> >();
			for(unsigned int i = 0 ; i < opts.predefinedFeature.size(); i++)
			{
				cout 		<< opts.predefinedFeature[i]<<" ";//copy(myvector.begin(), myvector.end(), ostream_iterator<int>(cout, "\n"));
				configFile 	<< opts.predefinedFeature[i]<<" ";
			}
		}
		else{
			cerr << "Please modify code for this new type: " << k <<endl;
		}

		cout			<<endl;
		configFile 		<<endl;
	}
	configFile.close();
	configFile.clear();

	if( options["ratiosFiles"].empty() )
	{
		//cout << "Mandatory parameters missing. Program will terminate now."<<endl;
		//exit(1);
		if( ( ! options["compFile"].empty() ) && boost::filesystem::exists( opts.compFile ) ){
			cout << "Since ratiosFiles are not specified but compFile " << opts.compFile << " is specified and exists, program will go to DMR calling." << endl;
			opts.out();
			opts.doMergeRatioFiles = 0;
			opts.doStrandSpecifiMeth = 0;
			opts.doComp = 0;
			opts.out();
		} else {
			cout << "Mandatory parameters either ratiosFiles or compFile missing. Program will terminate now."<<endl;
			exit(1);
		}
	}

	if(opts.doMergeRatioFiles == 0){
		opts.mergedRatiosFiles = opts.ratiosFiles;
	} else {
		if( options["mergedRatiosFiles"].empty() ){
			cout << "You request to merge files before comparison so please specify the filename for the merged file" << endl;
			exit(1);
		}
	}
	if(opts.labels.empty()){
		if( ! opts.mergedRatiosFiles.empty() ){
			opts.labels = opts.mergedRatiosFiles;
		} else {
		    string line = "";
			ifstream inputf(opts.compFile.c_str(), ios::in);
			getline(inputf, line);
			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			for( unsigned int i = 0; i < fields.size(); i ++ ){
				if( fields[i].find("totalC") != string::npos ){
					opts.labels.push_back(fields[i].erase(0, 7));
				}
			}
			inputf.close();
		}
	}

//	if( opts.xVector.empty() ){
//		vector <double> v;
//		double seed = 1.0;
//		for(unsigned int i = 0; i < opts.ratiosFiles.size(); i++){
//			v.push_back(seed+i);
//		}
//		opts.xVector = v;
//	}
//	if( opts.ratiosFiles.size() != opts.xVector.size() ){
//		cout << "xVector size is not same as number of ratiosFiles." << endl;
//		exit(1);
//	}
	//opts.out();
	return 0;
};

//int string_to_int( string s){return atoi(s.c_str());}
//string itos(int i)
//{
//    stringstream s;
//    s << i;
//    return s.str();
//}

class cMeth
{
public:
	int totalC;
	int methC;
	char strand;
	char nextN;
	cMeth()
	{
		totalC = 0;
		methC = 0;
		strand = '*';
		nextN = 'N';
	};
	cMeth(int t, int m, char s, char n) : totalC(t), methC(m), strand(s), nextN(n) {};

};

class statSim
{
public:
	long double genomeSimP;
	int cytosinesCountGenome;
	statSim(){
		genomeSimP = 0.0;
		cytosinesCountGenome = 0;
	}
	statSim(long double g, int c) : genomeSimP(g), cytosinesCountGenome(c) {};
};

void readLaneToHash( string file, Opts & option, map <string, map<int, cMeth> > & meth ) //read file back to hash //for around 6,000,000 records, this consumes 500M ram and takes 40 seconds.
{
	int colIdForChr = 0;
	int colIdForTotalC = 4;
	int colIdForMethC = 5;
	int colIdForStart = 1;
	int colIdForStrand = 6;
	int colIdForNext = 7;
	cout << "if no header #chom in file, default index order is #chrom, start, end, ratio, totalC, methC, strand, next" << endl;
	string chrom;
	int start;
	int totalC;
	int methC;
	char strand;
	char next;

	//ofstream testFile;
	//testFile.open("test", ios_base::app);

    int count = 0;
    string line = "";
	ifstream inputf(file.c_str(), ios::in);
	while (inputf.good()) {
		getline(inputf, line);
		if(line == ""){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		if( fields[0] == "#chrom" )		//if there's header line, reset the index values; otherwise just use defalut index values
		{
			vector<string>::iterator it;

			it = find (fields.begin(), fields.end(), "totalC");
			if( it != fields.end() ){ colIdForTotalC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "methC");
			if( it != fields.end() ){ colIdForMethC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "#chrom");
			if( it != fields.end() ){ colIdForChr = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "start");
			if( it != fields.end() ){ colIdForStart = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "strand");
			if( it != fields.end() ){ colIdForStrand = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

			it = find (fields.begin(), fields.end(), "next");
			if( it != fields.end() ){ colIdForNext = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

		} else {				//content
			chrom = fields[colIdForChr];
			start = string_to_int(fields[colIdForStart]);
			totalC = string_to_int(fields[colIdForTotalC]);
			methC = string_to_int(fields[colIdForMethC]);
			strand = fields[colIdForStrand][0];
			next = fields[colIdForNext][0];
			if(totalC >= option.minDepthForComp){
				meth[chrom][start] = cMeth(totalC, methC, strand, next);
			}
			//testFile << totalC <<"\t"<< methC <<"\t"<< strand <<"\t" << next << "" << endl;
			//cout << "totalC\tmethC\tstrand\tnext" << endl;
		}
		count += 1 ;
	}
	inputf.close();
}


//R cannot be threaded and hence the R part has be in one thread. In other words, anything in one subruitine with R must be processed line by line
string compOneLocSingleThread(map <int, map <string, map<int, cMeth> > > & lane, string chr, int start)
{
	char next = 'X';
	for( unsigned int i = 0; i < lane.size(); i++)
	{
		if(lane[i].count(chr) && lane[i][chr].count(start))
		{
			next = lane[i][chr][start].nextN;
			break;
		}
	}
	int end = (next == 'G') ? (start + 2) : (start + 1);

	stringstream toPrint ;
	toPrint << setprecision(3) << chr << "\t" << start << "\t" << end;

	//single comp to get CI
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 ){//chr or start does not exist;
			toPrint << "\tsingle:\t" << "NA" << "\t" << "NA" << "\t" << setprecision(3) << "NA" << "," << "NA" ;
			continue;
		}
		int methC = lane[i][chr][start].methC;
		int totalC = lane[i][chr][start].totalC;
		double nominal = double(methC) / totalC;

		//time_t tstart,tend;
		//time (&tstart);
		CI sCI = singleCI(methC, totalC, ALPHA, 1);
		//time (&tend);
		//double tdif = difftime (tend,tstart);
		//if(tdif>1)cout << "singleCI: " << start << "\t" << i << "\t" << tdif << endl;
		toPrint << "\tsingle:\t" << totalC << "\t" << nominal << "\t" << setprecision(3) << sCI.a << "," << sCI.b ;
	}

	//pair comp to get diffCI, and pvalue.
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tpair:\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "," << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

			CI ci = pdiffCIOpt(m2, t2, m1, t1,ALPHA, 1); //pdiffCI returns ci for r2-r1
			double nominalDif = double(m2)/t2 - double(m1)/t1;
			//credibleDif is minimum value of abs(p2-p1)
			double credibleDif = abs(ci.a) < abs(ci.b) ? abs(ci.a) : abs(ci.b);
			if( ci.a * ci.b <= 0 ){credibleDif = 0;}
			toPrint << "\tpair:\t" << nominalDif << "\t" << credibleDif << "\t" << ci.a << "," << ci.b ;

		}
	}
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tp:\t" << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

			//fisher exact test p value//commentted out method is Rcpp call.
			int varray[] = {m1, t1 - m1, m2, t2 - m2};
			//vector<int> v (varray, varray + sizeof(varray) / sizeof(int) );
			//R.assign(v,"v");
			//string txt = "test=fisher.test(matrix(c( v[1],v[2], v[3], v[4] ), nrow=2)); test$p.value";
			//ans = R.parseEval(txt);
			//double p = Rcpp::as< double >(ans);
			double p = 1.0;
			fet_22( p, varray, 2, sizeof(varray) / sizeof(int) / 2 );
			toPrint << "\tp:\t" << p ;

			//TODO: my p value;
		}
	}


	//group comp to get model fitting
//	unsigned int lanesCoveredCount  = 0;
//	for( unsigned int i = 0; i < lane.size(); i++ )
//	{
//		//fisher exact p value
//		if( lane[i].count(chr) != 0 && lane[i][chr].count(start) != 0 ){
//			lanesCoveredCount += 1;
//		}
//	}

//	if( lanesCoveredCount > 2 ){
//		if( lanesCoveredCount == lane.size() )// all lanes covering this location to do group comp.
//		{
//			vector <int> pairs;
//			vector <double> ratios;
//			for( unsigned int i = 0; i < lane.size(); i++ )
//			{
//				int mc = lane[i][chr][start].methC;
//				int tc = lane[i][chr][start].totalC;
//				int uc = tc - mc;
//				double ratio = double( mc ) / tc;
//				pairs.push_back(mc);
//				pairs.push_back(uc);
//				ratios.push_back(ratio);
//			}
//
//			//fisher's exact test p-value for all lanes //commentted out method is Rcpp call.
//			//R.assign(pairs, "pairs");
//			//string txt = "test=fisher.test(matrix(pairs, nrow=2),hybrid=TRUE, workspace = 2e8); test$p.value";
//			//ans = R.parseEval(txt);
//			//double p_all = Rcpp::as< double >(ans);
//			int *varray = & pairs[0];
//			int ncol = pairs.size() / 2;
//			double p_all = 1.0;
//			fet_1k( p_all, varray, 2, ncol );
//			toPrint << "\tfet_p_all:\t" << p_all ;
//			//cout << "xxx: " << pairs[0] << "\t" << pairs[1] << endl;
	
			//TODO: my method
	
//			//Linear model to get slope, correlation, R^2, and anova p-value
//			R.assign(ratios, "ratios");
//			//double v[] = {4.0, 12.0, 24.0};
//			//vector<double> times (v, v + sizeof(v) / sizeof(double) );
//			vector <double> times(opts.xVector);
//			R.assign(times, "times");
//			string txt = "fit=lm(ratios~times); coef(fit)[2]";
//			ans = R.parseEval(txt);
//			double rate = Rcpp::as< double >(ans);
//			txt = "cor(ratios,times)";
//			ans = R.parseEval(txt);
//			double cor = Rcpp::as< double >(ans);
//			txt = "summary(fit)$adj.r.squared";
//			ans = R.parseEval(txt);
//			double ars = Rcpp::as< double >(ans);
//			txt = "anova(fit)$P[1]";
//			ans = R.parseEval(txt);
//			double anovaP = Rcpp::as< double >(ans);
//
//			toPrint << "\tlm_rate_cor_ars_anovaP:\t" << rate << "\t" << cor << "\t" << ars << "\t" << anovaP ;
	
//		}
//		else {
//			toPrint << "\tfet_p_all:\t" << "NA" << "\tlm_rate_cor_ars_anovaP:\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" ;
//		}
//	}
	//toPrint << endl;
	return toPrint.str();
};


string compOneLocNonR(map <int, map <string, map<int, cMeth> > > & lane, string chr, int start)
{
	char next = 'X';
	for( unsigned int i = 0; i < lane.size(); i++)
	{
		if(lane[i].count(chr) && lane[i][chr].count(start))
		{
			next = lane[i][chr][start].nextN;
			break;
		}
	}
	int end = (next == 'G') ? (start + 2) : (start + 1);

	stringstream toPrint ;
	toPrint << setprecision(3) << chr << "\t" << start << "\t" << end;
	
	//single comp to get CI
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 ){//chr or start does not exist;
			toPrint << "\tsingle:\t" << "NA" << "\t" << "NA" << "\t" << setprecision(3) << "NA" << "," << "NA" ;
			continue;
		}
		int methC = lane[i][chr][start].methC;
		int totalC = lane[i][chr][start].totalC;
		double nominal = double(methC) / totalC;

		//toPrint << "\tsingle:\t" << "NA" << "\t" << "NA" << "\t" << setprecision(3) << "NA" << "," << "NA" ;
		//cout << "START singleCI " << totalC << "\t" << methC << endl;

		CI sCI = singleCI(methC, totalC, ALPHA, 1);

		toPrint << "\tsingle:\t" << totalC << "\t" << nominal << "\t" << setprecision(3) << sCI.a << "," << sCI.b ;
		//cout << "END singleCI " << totalC << "\t" << methC << endl;

	};



	//pair comp to get diffCI, and pvalue.
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tpair:\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "," << "NA";
				toPrint << "\tp_acc:\t" << "NA" ;
				//toPrint << "\tp_fet:\t" << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

			//toPrint << "\tpair:\t" << 0.55 << "\t" << 0.45 << "\t" << 0.3 << "," << 0.7 ;
			//cout << "START pdiffCI " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;

			CI ci = pdiffCIOpt(m2, t2, m1, t1,ALPHA, 1); //pdiffCI returns ci for r2-r1
			double nominalDif = double(m2)/t2 - double(m1)/t1;
			//credibleDif is minimum value of abs(p2-p1)
			double credibleDif = abs(ci.a) < abs(ci.b) ? abs(ci.a) : abs(ci.b);
			if( ci.a * ci.b <= 0 ){credibleDif = 0;}
			if( nominalDif < 0 ) { credibleDif = -1 * credibleDif; }
			toPrint << "\tpair:\t" << nominalDif << "\t" << credibleDif << "\t" << ci.a << "," << ci.b ;
			//cout << "END pdiffCI " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;


/*
			//the pvalue calculation
			//fisher exact test p value//commentted out method is Rcpp call.
			//cout << "pair test for " << m1 << "\t" << t1-m1 << "\t" << m2 << "\t" << t2-m2 << endl;
			//saturate reads at 1000; otherwise fet will fail for those sites covered like 30k+ times
			int varray[] = {m1, t1 - m1, m2, t2 - m2};
			double p = 1.0;
			fet_1k( p, varray, 2, sizeof(varray) / sizeof(int) / 2 );
			toPrint << "\tp_fet:\t" << p ;
*/
			//TODO: my p value;
			double p_acc = pdiffInRegionOpt(m1,t1,m2,t2,0.01);
//			This selection works for GaussLobatto method but Adapt methods hangs on for dataset 2839, 3028, 1548, 1653.
// 			if( nominalDif < 0 ) //P(-d<p1-p2<d)=P(-d<p2-p1<d)= P(p1-p2>-d)-P(p1-p2>d) = P(p2-p1>-d)-P(p2-p1>d) for accruacy reason
// 			{
// 				p_acc = pdiff(m2,t2,m1,t1,-0.01)-pdiff(m2,t2,m1,t1,0.01);
// 			} else 
// 			{
// 				p_acc = pdiff(m1,t1,m2,t2,-0.01)-pdiff(m1,t1,m2,t2,0.01);
// 			}
			//cout << "starting to do " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;

			//cout << "finished to do " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;
			toPrint << "\tp_acc:\t" << p_acc;

		}
	}

//	stringstream x;
//	x << "..." << chr << "..." << start << "..." <<  "ZXCVASDFqwer";
//	return toPrint.str() + x.str();
/*
	//group comp to get model fitting
	unsigned int lanesCoveredCount  = 0;
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		//fisher exact p value
		if( lane[i].count(chr) != 0 && lane[i][chr].count(start) != 0 ){
			lanesCoveredCount += 1;
		}
	}
	if( lanesCoveredCount > 2 ){
		if( lanesCoveredCount == lane.size() )// all lanes covering this location to do group comp.
		{
			string txt = "";
	
			vector <int> pairs;
			vector <double> ratios;
			for( unsigned int i = 0; i < lane.size(); i++ )
			{
				int mc = lane[i][chr][start].methC;
				int tc = lane[i][chr][start].totalC;
				int uc = tc - mc;
				double ratio = double( mc ) / tc;
				pairs.push_back(mc);
				pairs.push_back(uc);
				ratios.push_back(ratio);
			}
	
			//saturate reads at 1000; otherwise fet will fail for those sites covered like 30k+ times
			//fisher's exact test p-value for all lanes //commentted out method is Rcpp call.
			int *varray = & pairs[0];
			int ncol = pairs.size() / 2;
			double p_all = 1.0;
			//cout << "group test for " << varray[0] << "\t" << varray[1] << "\t" << varray[2] << "\t" << varray[3] << endl;
			
			fet_1k( p_all, varray, 2, ncol );
			//cout << "\t....done" << endl;
			toPrint << "\tfet_p_all:\t" << p_all ;
		} 
		else {
			toPrint << "\tfet_p_all:\t" << "NA";
		}
	}
*/
	return toPrint.str();
}

void compBatchLocNonR(map <int, map <string, map<int, cMeth> > > & lane, string chr, set <int> & starts, int pstart, int batchSize, vector<string> & out)
{
	set <int>::iterator it = starts.begin();
	for(int j = 0; j <  pstart; j++) it++; //now it points to starts.begin() + pstart// do not understand why "it=starts.begin() + pstart" does not compile.

	for(int i = 0; i < batchSize; i++)
	{
		int start = *it;
		//
//		stringstream x;
//		x << chr << "..." << pstart << "..." << i << "..." << start << "..." <<  "ZXCVASDFqwer";
//		out.push_back(x.str());
		string outOne = compOneLocNonR(lane, chr, start);
		out.push_back(outOne);
		it++;
	}
	//cout << "done with subthread with pstart: " << pstart <<endl;
}


string compOneLocR(map <int, map <string, map<int, cMeth> > > & lane, string chr, int start)
{
	stringstream toPrint ;
	toPrint << setprecision(3);

	//now this is nonR calculation
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tp_fet:\t" << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

//			int useRcpp = 0;
//			string txt = "";
//			//fisher exact test p value//commentted out method is Rcpp call.
			int varray[] = {m1, t1 - m1, m2, t2 - m2};
//			if(useRcpp){
//				vector<int> v (varray, varray + sizeof(varray) / sizeof(int) );
//				R.assign(v,"v");
//				txt = "test=fisher.test(matrix(c( v[1],v[2], v[3], v[4] ), nrow=2)); test$p.value";
//				ans = R.parseEval(txt);
//				double p = Rcpp::as< double >(ans);
//				toPrint << "\tp_fet:\t" << p ;
//			}
//			else{
				double p = 1.0;
				fet_22( p, varray, 2, sizeof(varray) / sizeof(int) / 2 );
				toPrint << "\tp_fet:\t" << p ;
//			}
			//TODO: my p value;
		}
	}


//	//group comp to get model fitting
//	unsigned int lanesCoveredCount  = 0;
//	for( unsigned int i = 0; i < lane.size(); i++ )
//	{
//		//fisher exact p value
//		if( lane[i].count(chr) != 0 && lane[i][chr].count(start) != 0 ){
//			lanesCoveredCount += 1;
//		}
//	}

//	if( lanesCoveredCount > 2 ) {
//		if( lanesCoveredCount == lane.size() )// all lanes covering this location to do group comp.
//		{
//			string txt = "";
//
//			vector <int> pairs;
//			vector <double> ratios;
//			for( unsigned int i = 0; i < lane.size(); i++ )
//			{
//				int mc = lane[i][chr][start].methC;
//				int tc = lane[i][chr][start].totalC;
//				int uc = tc - mc;
//				double ratio = double( mc ) / tc;
//				pairs.push_back(mc);
//				pairs.push_back(uc);
//				ratios.push_back(ratio);
//				//cout << start << "\t" << tc << "\t==" << mc << endl;
//			}
//
//			int useRcpp = 0;
//
//			//fisher's exact test p-value for all lanes //commentted out method is Rcpp call.
//			if( useRcpp ){
//				R.assign(pairs, "pairs");
//				txt = "test=fisher.test(matrix(pairs, nrow=2),hybrid=TRUE, workspace = 2e8); test$p.value";
//				ans = R.parseEval(txt);
//				double p_all = Rcpp::as< double >(ans);
//				toPrint << "\tfet_p_all:\t" << p_all ;
//			}
//			else{
//				int *varray = & pairs[0];
//				int ncol = pairs.size() / 2;
//				double p_all = 1.0;
//				fet_1k( p_all, varray, 2, ncol );
//				toPrint << "\tfet_p_all:\t" << p_all ;
//			}

			//TODO: my method

//			copy(ratios.begin(), ratios.end(), ostream_iterator<double>(cout, ","));
//			//Linear model to get slope, correlation, R^2, and anova p-value
//			R.assign(ratios, "ratios");
//			//double v[] = {4.0, 12.0, 24.0};
//			//vector<double> times (v, v + sizeof(v) / sizeof(double) );
//			vector <double> times(opts.xVector);
//			copy(times.begin(), times.end(), ostream_iterator<double>(cout, ","));
//			R.assign(times, "times");
//			txt = "fit=lm(ratios~times); coef(fit)[2]";
//			ans = R.parseEval(txt);
//			double rate = Rcpp::as< double >(ans);
//			//cout << rate <<"rrr" << endl;
//			txt = "cor(ratios,times)";
//			ans = R.parseEval(txt);
//			double cor = Rcpp::as< double >(ans);
//			txt = "summary(fit)$adj.r.squared";
//			ans = R.parseEval(txt);
//			double ars = Rcpp::as< double >(ans);
//			txt = "anova(fit)$P[1]";
//			ans = R.parseEval(txt);
//			double anovaP = Rcpp::as< double >(ans);
//
//			toPrint << "\tlm_rate_cor_ars_anovaP:\t" << rate << "\t" << cor << "\t" << ars << "\t" << anovaP ;
//		}
//		else {
//			toPrint << "\tlm_rate_cor_ars_anovaP:\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" ;
//		}
//	}
	//toPrint << endl;
	//cout << toPrint.str() << endl;;
	return toPrint.str();
}

void compBatchLocR(map <int, map <string, map<int, cMeth> > > & lane, string chr, set <int> & starts, int thisStart, int length, vector<string> & rout)
{
	set <int>::iterator it = starts.begin();
	for(int j = 0; j <  thisStart; j++) it++; //now it points to starts.begin() + pstart// do not understand why "it=starts.begin() + pstart" does not compile.

	for(int i = 0; i < length; i++)
	{
		int start = *it;
		string outOne = compOneLocR(lane, chr, start);
		rout.push_back(outOne);
		//cout << i << "=i and rout.size()=" << rout.size() << endl;
		it++;
		//cout << "HERE: " << i << endl;
	}
}






string compOneLoc(map <int, map <string, map<int, cMeth> > > & lane, string chr, int start)
{
	char next = 'X';
	for( unsigned int i = 0; i < lane.size(); i++)
	{
		if(lane[i].count(chr) && lane[i][chr].count(start))
		{
			next = lane[i][chr][start].nextN;
			break;
		}
	}
	int end = (next == 'G') ? (start + 2) : (start + 1);

	stringstream toPrint ;
	toPrint << setprecision(3) << chr << "\t" << start << "\t" << end;

	//single comp to get CI
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 ){//chr or start does not exist;
			toPrint << "\tsingle:\t" << "NA" << "\t" << "NA" << "\t" << setprecision(3) << "NA" << "," << "NA" ;
			continue;
		}
		int methC = lane[i][chr][start].methC;
		int totalC = lane[i][chr][start].totalC;
		double nominal = double(methC) / totalC;

		//toPrint << "\tsingle:\t" << "NA" << "\t" << "NA" << "\t" << setprecision(3) << "NA" << "," << "NA" ;
		//cout << "START singleCI " << totalC << "\t" << methC << endl;

		CI sCI = singleCI(methC, totalC, ALPHA, 1);

		toPrint << "\tsingle:\t" << totalC << "\t" << nominal << "\t" << setprecision(3) << sCI.a << "," << sCI.b ;
		//cout << "END singleCI " << totalC << "\t" << methC << endl;

	};



	//pair comp to get diffCI, and pvalue.
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tpair:\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "," << "NA";
				toPrint << "\tp_acc:\t" << "NA" ;
				//toPrint << "\tp_fet:\t" << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

			//toPrint << "\tpair:\t" << 0.55 << "\t" << 0.45 << "\t" << 0.3 << "," << 0.7 ;
			//cout << "START pdiffCI " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;

			CI ci = pdiffCIOpt(m2, t2, m1, t1,ALPHA, 1); //pdiffCI returns ci for r2-r1
			double nominalDif = double(m2)/t2 - double(m1)/t1;
			//credibleDif is minimum value of abs(p2-p1)
			double credibleDif = abs(ci.a) < abs(ci.b) ? abs(ci.a) : abs(ci.b);
			if( ci.a * ci.b <= 0 ){credibleDif = 0;}
			if( nominalDif < 0 ) { credibleDif = -1 * credibleDif; }
			toPrint << "\tpair:\t" << nominalDif << "\t" << credibleDif << "\t" << ci.a << "," << ci.b ;
			//cout << "END pdiffCI " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;


/*
			//the pvalue calculation
			//fisher exact test p value//commentted out method is Rcpp call.
			//cout << "pair test for " << m1 << "\t" << t1-m1 << "\t" << m2 << "\t" << t2-m2 << endl;
			//saturate reads at 1000; otherwise fet will fail for those sites covered like 30k+ times
			int varray[] = {m1, t1 - m1, m2, t2 - m2};
			double p = 1.0;
			fet_1k( p, varray, 2, sizeof(varray) / sizeof(int) / 2 );
			toPrint << "\tp_fet:\t" << p ;
*/
			//TODO: my p value;
			double p_acc = pdiffInRegionOpt(m1,t1,m2,t2,0.01);
//			This selection works for GaussLobatto method but Adapt methods hangs on for dataset 2839, 3028, 1548, 1653.
// 			if( nominalDif < 0 ) //P(-d<p1-p2<d)=P(-d<p2-p1<d)= P(p1-p2>-d)-P(p1-p2>d) = P(p2-p1>-d)-P(p2-p1>d) for accruacy reason
// 			{
// 				p_acc = pdiff(m2,t2,m1,t1,-0.01)-pdiff(m2,t2,m1,t1,0.01);
// 			} else
// 			{
// 				p_acc = pdiff(m1,t1,m2,t2,-0.01)-pdiff(m1,t1,m2,t2,0.01);
// 			}
			//cout << "starting to do " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;

			//cout << "finished to do " << m2 << "\t" << t2 << "\t" << m1 << "\t" << t1 << endl;
			toPrint << "\tp_acc:\t" << p_acc;

		}
	}


	//now this is nonR calculation
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			if( lane[i].count(chr) == 0 || lane[i][chr].count(start) == 0 || lane[j].count(chr) == 0 || lane[j][chr].count(start) == 0 ){//chr or start does not exist;
				toPrint << "\tp_fet:\t" << "NA" ;
				continue;
			}

			int m1 = lane[i][chr][start].methC;
			int t1 = lane[i][chr][start].totalC;
			int m2 = lane[j][chr][start].methC;
			int t2 = lane[j][chr][start].totalC;

//			int useRcpp = 0;
//			string txt = "";
//			//fisher exact test p value//commentted out method is Rcpp call.
			int varray[] = {m1, t1 - m1, m2, t2 - m2};
//			if(useRcpp){
//				vector<int> v (varray, varray + sizeof(varray) / sizeof(int) );
//				R.assign(v,"v");
//				txt = "test=fisher.test(matrix(c( v[1],v[2], v[3], v[4] ), nrow=2)); test$p.value";
//				ans = R.parseEval(txt);
//				double p = Rcpp::as< double >(ans);
//				toPrint << "\tp_fet:\t" << p ;
//			}
//			else{
				double p = 1.0;
				fet_22( p, varray, 2, sizeof(varray) / sizeof(int) / 2 );
				toPrint << "\tp_fet:\t" << p ;
//			}
			//TODO: my p value;
		}
	}

//	stringstream x;
//	x << "..." << chr << "..." << start << "..." <<  "ZXCVASDFqwer";
//	return toPrint.str() + x.str();
/*
	//group comp to get model fitting
	unsigned int lanesCoveredCount  = 0;
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		//fisher exact p value
		if( lane[i].count(chr) != 0 && lane[i][chr].count(start) != 0 ){
			lanesCoveredCount += 1;
		}
	}
	if( lanesCoveredCount > 2 ){
		if( lanesCoveredCount == lane.size() )// all lanes covering this location to do group comp.
		{
			string txt = "";

			vector <int> pairs;
			vector <double> ratios;
			for( unsigned int i = 0; i < lane.size(); i++ )
			{
				int mc = lane[i][chr][start].methC;
				int tc = lane[i][chr][start].totalC;
				int uc = tc - mc;
				double ratio = double( mc ) / tc;
				pairs.push_back(mc);
				pairs.push_back(uc);
				ratios.push_back(ratio);
			}

			//saturate reads at 1000; otherwise fet will fail for those sites covered like 30k+ times
			//fisher's exact test p-value for all lanes //commentted out method is Rcpp call.
			int *varray = & pairs[0];
			int ncol = pairs.size() / 2;
			double p_all = 1.0;
			//cout << "group test for " << varray[0] << "\t" << varray[1] << "\t" << varray[2] << "\t" << varray[3] << endl;

			fet_1k( p_all, varray, 2, ncol );
			//cout << "\t....done" << endl;
			toPrint << "\tfet_p_all:\t" << p_all ;
		}
		else {
			toPrint << "\tfet_p_all:\t" << "NA";
		}
	}
*/
	return toPrint.str();
}

void compBatchLoc(map <int, map <string, map<int, cMeth> > > & lane, string chr, set <int> & starts, int pstart, int batchSize, vector<string> & out)
{
	set <int>::iterator it = starts.begin();
	for(int j = 0; j <  pstart; j++) it++; //now it points to starts.begin() + pstart// do not understand why "it=starts.begin() + pstart" does not compile.

	for(int i = 0; i < batchSize; i++)
	{
		int start = *it;
		//
//		stringstream x;
//		x << chr << "..." << pstart << "..." << i << "..." << start << "..." <<  "ZXCVASDFqwer";
//		out.push_back(x.str());
		string outOne = compOneLoc(lane, chr, start);
		out.push_back(outOne);
		it++;
	}
	//cout << "done with subthread with pstart: " << pstart <<endl;
}



//can not use reference for out or rout because they might be destroyed in a new batch.
void appendToFile( ofstream & allCompFile, vector<vector<string> >  out, vector<string>  rout )
{
	unsigned int count  = 0;
	for(unsigned int i = 0; i < out.size(); i++)
	{
		vector <string> oneBatch = out[i];
		for(unsigned int j = 0; j < oneBatch.size(); j++)
		{
			if( rout.empty() ){
				allCompFile << oneBatch[j] << endl;
			} else {
				allCompFile << oneBatch[j] << rout[count] << endl;
			}
			count ++;
		}
	}
//	if(count != rout.size()){cerr << "unknown error here..." << count << " and " << rout.size() << endl; exit(1);}
};

//can not use reference for out or rout because they might be destroyed in a new batch.
void appendToFile2( ofstream & allCompFile, vector<vector<string> >  out)
{
	for(unsigned int i = 0; i < out.size(); i++){
		vector <string> oneBatch = out[i];
		for(unsigned int j = 0; j < oneBatch.size(); j++){
				allCompFile << oneBatch[j] << endl;
		}
	}
};

class PairCompOut
{
public:
	//string loc; //"chr\tstart\tend"
	string chr;
	int start;
	int end;
	int ta;		//totalC_a;
	double ra; 	//nominalRatio_a;
	string cia; 	//ratioCI_a;
	int	tb;		//totalC_b;
	double rb;	//nominalRatio_b;
	string cib;	//ratioCI_b;
	double simP;
	double fetP;
	double nDif;	//nominalDif;
	string difCI;
	double cDif;	//credibleDif;
	string cls;		//classification of DMC
	PairCompOut(){};
	PairCompOut(string chrom, int st, int en, int t1, double r1, string ci1, int t2, double r2, string ci2, double sp, double fp, double NDif, string DifCI, double CDif, string c) : chr(chrom), start(st), end(en), ta(t1), ra(r1), cia(ci1),tb(t2), rb(r2), cib(ci2), simP(sp), fetP(fp), nDif(NDif), difCI(DifCI), cDif(CDif), cls(c) {};
	string ToString(){
		stringstream ss;
		ss 	<< chr << "\t" << start << "\t" << end << "\t" << ta << "\t" << ra << "\t" << cia << "\t" << tb << "\t" << rb << "\t" << cib
			<< "\t" << nDif << "\t" << cDif << "\t" << difCI << "\t" << simP << "\t" << fetP << "\t"  << cls;
		return ss.str();
	}
	string ToBedTrackString(){
		stringstream ss;
		ss 	<< setprecision(3) << chr << "\t" << start << "\t" << end
			<< "\t" << "meth" << rb*100 << "vs" << ra*100 << "\t" << 1000 << "\t" << "+"
			<< "\t" << 0 << "\t" << 0 << "\t" << (nDif > 0 ? "160,0,0" : "0,0,160");
		return ss.str();
	}
};

class Dmr
{
public:
	string chr;
	int start;
	int end;
	int allSites;//# of all cytosine sites
	int dmcSites;
	int mC1;
	int uC1;
	int mC2;
	int uC2;
	double meanRatio1;
	double meanRatio2;
	double methDif;
	double fetP;
	string cls;
	Dmr(){};
	Dmr(string Chr, int Start, int End, int AllSites, int DmcSites, int MC1, int UC1, int MC2, int UC2, double mr1, double mr2, double MethDif, double FetP, string c) : chr(Chr), start(Start), end(End), allSites(AllSites), dmcSites(DmcSites), mC1(MC1), uC1(UC1), mC2(MC2), uC2(UC2), meanRatio1(mr1), meanRatio2(mr2), methDif(MethDif), fetP(FetP), cls(c) {};
	//replace end + 2 by end in function ToString() and ToBedTrackString() for consideration of predefined region.
	string ToString(){
		stringstream ss;
		ss 	<< setprecision(3) << chr << "\t" << start << "\t" << end + 2
			<< "\t" << meanRatio1 << "\t" << mC1 + uC1 << "\t" << allSites
			<< "\t" << meanRatio2 << "\t" << mC2 + uC2 << "\t" << allSites
			<< "\t" << methDif << "\t" << fetP << "\t" << cls;
		return ss.str();
	}
	string ToBedTrackString(){
		stringstream ss;
		ss 	<< setprecision(3) << chr << "\t" << start << "\t" << end + 2
			<< "\t" << "meth" << meanRatio2*100 << "vs" << meanRatio1*100 << "\t" << 1000 << "\t" << "+"
			<< "\t" << 0 << "\t" << 0 << "\t" << (methDif > 0 ? "160,0,0" : "0,0,160");
		return ss.str();
	}

};

bool sortDmrByMethDif(const Dmr& d1, const Dmr& d2)
{
  return abs(d1.methDif) > abs(d2.methDif);
}

void getDmr(string chr, vector<int> & allLocs, map<int, PairCompOut> & thisChr, vector<Dmr> & dmrs )
{
	double pCut = 0.005;//TODO parameter input;
	int minDepth = 0;
	double minNominalDiff = 0;
	double minConfidentRatio = 0;

	Dmr dmr;
	for(unsigned int i = 0; i < allLocs.size(); i ++)
	{


	}

	//PairCompOut dimer = thisChr[x];

};

// copied from void readLaneToStrandSpecificHash( string file, map <string, map<int, cMeth> > & methPs, map <string, map<int, cMeth> > & methMs)
void readLaneToStrandSpecificHash( string file, map <string, map<int, cMeth> > & methPs, map <string, map<int, cMeth> > & methMs, string chr )
{
	int colIdForChr = 0;
	int colIdForTotalCPs = 9;
	int colIdForMethCPs = 10;
	int colIdForTotalCMs = 12;
	int colIdForMethCMs = 13;
	int colIdForStart = 1;
	//int colIdForStrand = -1;
	int colIdForNext = 7;

	string chrom;
	int start;
	int totalCPs;
	int methCPs;
	int totalCMs;
	int methCMs;
	char strand;
	char next;

	int count = 0;
	string line = "";
	ifstream inputf(file.c_str(), ios::in);
	while (inputf.good()) {
		getline(inputf, line);
		if(line == ""){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		if( count == 0 )		//header line
		{
		} else {				//content
			chrom = fields[colIdForChr];
			if (chrom != chr) continue;
			start = string_to_int(fields[colIdForStart]);
			totalCPs = string_to_int(fields[colIdForTotalCPs]);
			methCPs = string_to_int(fields[colIdForMethCPs]);
			totalCMs = string_to_int(fields[colIdForTotalCMs]);
			methCMs = string_to_int(fields[colIdForMethCMs]);
			//strand = fields[colIdForStrand][0];
			next = fields[colIdForNext][0];
			if(totalCPs > 0){ methPs[chrom][start] = cMeth(totalCPs, methCPs, '+', next);}
			if(totalCMs > 0){ methMs[chrom][start] = cMeth(totalCMs, methCMs, '-', next);}
		}
		count += 1 ;
	}
	inputf.close();
}

void readLaneToStrandSpecificHash( string file, map <string, map<int, cMeth> > & methPs, map <string, map<int, cMeth> > & methMs ) //read file back to hash //for around 6,000,000 records, this consumes 500M ram and takes 40 seconds.
{
	int colIdForChr = 0;
	int colIdForTotalCPs = 9;
	int colIdForMethCPs = 10;
	int colIdForTotalCMs = 12;
	int colIdForMethCMs = 13;
	int colIdForStart = 1;
	//int colIdForStrand = -1;
	int colIdForNext = 7;

	string chrom;
	int start;
	int totalCPs;
	int methCPs;
	int totalCMs;
	int methCMs;
	char strand;
	char next;

	//ofstream testFile;
	//testFile.open("test", ios_base::app);

    int count = 0;
    string line = "";
	ifstream inputf(file.c_str(), ios::in);
	while (inputf.good()) {
		getline(inputf, line);
		if(line == ""){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		if( count == 0 )		//header line
		{
			//modify this after header is changed to unique header fields
//			vector<string>::iterator it;
//
//			it = find (fields.begin(), fields.end(), "totalC");
//			if( it != fields.end() ){ colIdForTotalC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "methC");
//			if( it != fields.end() ){ colIdForMethC = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "#chrom");
//			if( it != fields.end() ){ colIdForChr = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "start");
//			if( it != fields.end() ){ colIdForStart = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "strand");
//			if( it != fields.end() ){ colIdForStrand = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }
//
//			it = find (fields.begin(), fields.end(), "next");
//			if( it != fields.end() ){ colIdForNext = it - fields.begin(); }else{ cerr << "header element not found" << endl; exit(1); }

		} else {				//content
			chrom = fields[colIdForChr];
			start = string_to_int(fields[colIdForStart]);
			totalCPs = string_to_int(fields[colIdForTotalCPs]);
			methCPs = string_to_int(fields[colIdForMethCPs]);
			totalCMs = string_to_int(fields[colIdForTotalCMs]);
			methCMs = string_to_int(fields[colIdForMethCMs]);
			//strand = fields[colIdForStrand][0];
			next = fields[colIdForNext][0];
			if(totalCPs > 0){ methPs[chrom][start] = cMeth(totalCPs, methCPs, '+', next);}
			if(totalCMs > 0){ methMs[chrom][start] = cMeth(totalCMs, methCMs, '-', next);}
			//cout << methPs[chrom][start].totalC << endl;
			//testFile << totalC <<"\t"<< methC <<"\t"<< strand <<"\t" << next << "" << endl;
			//cout << "totalC\tmethC\tstrand\tnext" << endl;
		}
		count += 1 ;
	}
	inputf.close();
}




//compare lanes and output result to fileName
void doComp(map <int, map <string, map<int, cMeth> > > & lane, string fileName, Opts & option)
{

///* //comment out this region for dmr only

	ofstream allCompFile;
	allCompFile.open(fileName.c_str(), ios_base::out);
	allCompFile << setprecision(3);


	//prepare the genomic locations that exist in at least one lane
	set <string> chroms;

	for(unsigned int i = 0; i < lane.size(); i++ ){
		for(map<string, map<int, cMeth> >::iterator it = lane[i].begin(); it != lane[i].end(); ++it) {
		  chroms.insert(it->first);
		  //cout << "CCC:" << it->first << endl;
		}
	}
	//cout << "finish merging" << endl;

	//start of header
	allCompFile << "#chrom" << "\t" << "start" << "\t" << "end";
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		allCompFile << "\tsingle:\t" << "totalC_" << i << "\t" << "nominalRatio_" << i << "\t" << "ratioCI_" << i ;
	}
	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			allCompFile << "\tpair:\t" << "nominalDif_" << j << "-" << i << "\t" << "credibleDif_" << j << "-" << i << "\t" << "difCI_" << j << "-" << i ;
			//allCompFile	<< "\tp_fet:\t" << "p_" << j << "_v_" << i ;
		}
	}

	//allCompFile << "\t"; //to make up the bug//tobe fixed

	for( unsigned int i = 0; i < lane.size(); i++ )
	{
		for( unsigned int j = i + 1; j < lane.size(); j++ )	//compare i with j;
		{
			allCompFile << "\tp_sim:\t" << "p_sim_" << j << "_v_" << i ;
			allCompFile << "\tp_fet:\t" << "p_fet_" << j << "_v_" << i ;
		}
	}
//	if(lane.size() > 2 ){
//		allCompFile << "\tfet_p_all:\t" << "fet_p_all" << "\tlm_rate_cor_ars_anovaP:\t" << "rate" << "\t" << "cor" << "\t" << "ars" << "\t" << "anovaP" ;
//	}
	allCompFile << endl;
	//end of header

	//cout << "pp=" << option.threads << " size=" << lane.size() << endl;
	for (set<string>::iterator pchr=chroms.begin(); pchr!=chroms.end(); pchr++)
	{
		string chr = *pchr;
		//cout << chr << endl;
	}
	//if(option.threads == 0)//dsun
	/*
	 * @ 20190804 by Jin Li
	 * threads=1 will generate different comparison table with threads>1
	 * only FET p-value is calculated when threads=1, which will cause segmentation fault in downstream DRM calling
	 * Because both psim and pfet will be read in DMC and DMR calling
	 */
	// if(option.threads == 1)
	// {
	// 	//loop through all the genomic locations by chrom and then start.
	// 	for (set<string>::iterator pchr=chroms.begin(); pchr!=chroms.end(); pchr++)
	// 	{
	// 		string chr = *pchr;
	// 		//cout << "start chr" << endl;

	// 		//get all locations for current chrom
	// 		set <int> starts;
	// 		for(unsigned int i = 0; i < lane.size(); i++ ){
	// 			for(map<int, cMeth>::iterator it = lane[i][chr].begin(); it != lane[i][chr].end(); ++it) {
	// 			  starts.insert(it->first);
	// 			}
	// 		}

	// 		//do a single comparison for chr->start->cMeth
	// 		for(set<int>::iterator pstart = starts.begin(); pstart != starts.end(); pstart++ ){
	// 			int start = *pstart;
	// 			//cout << "this :" << chr << "..." << start << endl;
	// 			allCompFile << compOneLocSingleThread(lane, chr, start) << endl;
	// 		}
	// 	}
	// }
	// else 				// do a comparison of a batch of genomic sites
	{
		//vector < vector<string> > outToPrint;
		//vector<string> routToPrint;

		for (set<string>::iterator pchr=chroms.begin(); pchr!=chroms.end(); pchr++)
		{
			string chr = *pchr;

			//get all locations for current chrom
			set <int> starts;
			for(unsigned int i = 0; i < lane.size(); i++ ){
				for(map<int, cMeth>::iterator it = lane[i][chr].begin(); it != lane[i][chr].end(); ++it) {
				  starts.insert(it->first);
				}
			}

			//do calculation and output to file
			bool endOfChr = false;
			int pstart = 0;
			int batchSize = BATCHMAX;
			boost::thread_group rio;
			//this way of threading makes sure the output is still sorted
			while( ! endOfChr )
			{
				boost::thread_group g;
				vector < vector<string> > out;
				out.resize(option.threads);
				int thisStart = pstart;

				//start nonR calculation on threads
				//cout << "starting subthreads" << endl;
				int i = 0;
				for( ; i < option.threads;  )
				{
					if(starts.size() - pstart < BATCHMAX)
					{
						//cout << i << endl;
						batchSize = starts.size() - pstart;
						boost::thread *tp = new boost::thread( compBatchLoc, boost::ref(lane), chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]) );
						g.add_thread(tp);
						//compBatchLoc(lane, chr, starts, pstart, batchSize, outStr[i]);
						pstart = pstart + batchSize;
						i++;
						//cout << i << endl;
						cout << g.size() << endl;
						endOfChr = true;
						break;//no more location in current chrom. breaks out of inner most for loop
					} else
					{
						//cout << i << endl;
						boost::thread *tp = new boost::thread( compBatchLoc, boost::ref(lane), chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]) );
						g.add_thread(tp);
						//compBatchLoc(lane, chr, starts, pstart, batchSize, outStr[i]);
						pstart = pstart + batchSize;
						i++;
						//cout << i << endl;
						cout << g.size() << endl;
					}
				}

				//start R calculation on main thread while nonR calculation is on backgroup thread
				vector<string> rout;
//				int length = (i-1)*BATCHMAX+batchSize;
				//cout << "starting main thread for length: " << length << endl;
				//rout.resize(length);//this will initialize the vector rout to length and push_back() will make it to length+1.
//				compBatchLocR(lane, chr, starts, thisStart, length, rout);
				//cout << "done with main thread: " << thisStart << "\t" << length << endl;
//				if(option.lmFit){
//					int length = (i-1)*BATCHMAX+batchSize;
//					cout << "starting main thread ." << endl;
//					rout.resize(length);
//					compBatchLocR(lane, chr, starts, thisStart, length, rout);
//					cout << "done with main thread: " << thisStart << "\t" << length <<endl;
//				}
				//cout << "waiting to join sub threads: "<< endl;

				g.join_all();
				//cout << "sub trheads joined......waiting to join output thread" << endl;


				//put both nonR and R calculation to file in a backgroud thread
				rio.join_all();
				//cout << "joined output thread" << endl;

				//out and rout are safe with local variable?
//				if(option.lmFit){
//					outToPrint = out;
//					routToPrint = rout;
//					boost::thread *tprout = new boost::thread( appendToFile, boost::ref(allCompFile), outToPrint, routToPrint);
//					rio.add_thread(tprout);
//				}
//				else{
//					outToPrint = out;
//					boost::thread *tprout = new boost::thread( appendToFile2, boost::ref(allCompFile), outToPrint );
//					rio.add_thread(tprout);
//				}
				//bug here: length should be included because at the end of chrom, out[i] might be empty.
				boost::thread *tprout = new boost::thread( appendToFile, boost::ref(allCompFile), out, rout);
				rio.add_thread(tprout);
				//cout << "added output thread" << endl;
			}
			rio.join_all();
		}
	}
	allCompFile.close();
	//finish of compaire all sites.
//*/


}


//assymmetrically methylated cpg (AMC)
void strandSpecificMethForLane( string inputFile )
{
	map<int, map <string, map<int, cMeth> > > strandLanes ; //strand->chrom->loc->(totalC, methC, strand, nextN);
	map <string, map<int, cMeth> > methPs;
	map <string, map<int, cMeth> > methMs;
	readLaneToStrandSpecificHash( inputFile, methPs, methMs );
	strandLanes[0] = methPs;
	strandLanes[1] = methMs;
	string outF = inputFile + "_ss.bed";
	Opts option(opts);
	option.mergeNotIntersect = 0;
	option.lmFit = 0;//no need too
	//option.xVector.clear();		//no need
	//option.xVector.push_back(1.0);	//no need
	//option.xVector.push_back(1.0);	//no need
	doComp(strandLanes, outF, option);
	string inFile = outF;
	outF = inputFile + "_ssCredible.bed";
	double minCredibleDif = 0.5;
	int minTwoStrands = 14;
	int minOneStrand = 1;
	string cmd = "";
	cmd = "perl -F'\\t' -lane 'print if($_ !~ qq(NA));' " + inFile + " > " + inputFile + "_ssBoth.bed";
	cout << " executing command: " << cmd << endl;
	system(cmd.c_str());
	cmd = "perl -F'\\t' -lane 'print if($F[13] > 0.5);' " + inFile + " > " + inputFile + "_ssCredible.bed";
	cout << " executing command: " << cmd << endl;
	system(cmd.c_str());	
	//selectStrandSpecific(inFile, outF, minCredibleDif);
	//getStats(inFile);
}

string join_vec(vector<int> & vec){
	std::ostringstream oss;
	if (!vec.empty()){
		std::copy(vec.begin(), vec.end()-1, std::ostream_iterator<int>(oss, " "));
		oss << vec.back();
	}
	return oss.str();
}

string join_vec(vector<string> & vec){
	std::ostringstream oss;
	if (!vec.empty()){
		std::copy(vec.begin(), vec.end()-1, std::ostream_iterator<string>(oss, " "));
		oss << vec.back();
	}
	return oss.str();
}

int CallBetaBinomialFit(map< int, vector< vector <int> > > & tcmcs, map< int, BiKey > & fits) {
	string input;
	for (map< int, vector< vector <int> > > :: iterator it=tcmcs.begin(); it!=tcmcs.end(); ++it) {
		input += to_string(it->first) + " " + join_vec(it->second[0]) + " " + join_vec(it->second[1]) + "\n";
	}
	std::cout << '.';
	int fd1[2];
	int fd2[2];
	pid_t cid;
	if(-1==pipe(fd1)) {
		fprintf(stderr, "Pipe failed");
		return 1;
	}
	if(-1==pipe(fd2)) {
		fprintf(stderr, "Pipe failed");
		return 1;
	}

	cid=fork();
	if (cid==-1) {
		fprintf(stderr, "Fork error");
		return 1;
	} else if (cid==0) {
		char buff[BUFFSIZE]={'\0'};
		string result;

		close(fd1[1]);
		read(fd1[0], buff, BUFFSIZE);
		close(fd1[0]);

		std::vector<string> lines;
		boost::split(lines, buff, boost::is_any_of("\n"));
		for (int lno=0; lno<lines.size(); lno++) {
			if (lines[lno].empty()) continue;
			std::vector<string> fields;
			boost::split(fields, lines[lno], boost::is_any_of(" "));
			std::vector<int> n;
			std::vector<int> k;
			int nums=fields.size();
			for(int i=1; i<nums; i++){
				if(i<=nums/2){
					n.push_back(string_to_int(fields[i]));
				} else {
					k.push_back(string_to_int(fields[i]));
				}
			}

			BiKey bk(-1,-1);
			BetaBinomialFit(n, k, bk);
			result += fields[0] + " " + to_string(bk.n1) + " " + to_string(bk.k1) + "\n";
		}

		close(fd2[0]);
		write(fd2[1], result.c_str(), result.size()+1);
		close(fd2[1]);
		exit(0);
	} else {
		char buff[BUFFSIZE]={'\0'};

		close(fd1[0]);
		write(fd1[1], input.c_str(), input.size()+1);
		close(fd1[1]);
		wait(NULL);
		close(fd2[1]);
		read(fd2[0], buff, BUFFSIZE);
		close(fd2[0]);

		std::vector<string> lines;
		boost::split(lines, buff, boost::is_any_of("\n"));
		for (int lno=0; lno<lines.size(); lno++) {
			if (lines[lno].empty()) continue;
			std::vector<string> fields;
			boost::split(fields, lines[lno], boost::is_any_of(" "));
			if (fields.size()==3) {
				int start=string_to_int(fields[0]);
				fits[start].n1=string_to_int(fields[1]);
				fits[start].k1=string_to_int(fields[2]);
			}
		}
	}
	return 0;
}

void mergeRatioFilesWorker(	map<int, map <string, map<int, cMeth> > >  & lanesPlus, map<int, map <string, map<int, cMeth> > >  & lanesMinus, string chr, set <int> & starts, int pstart, int batchSize, vector<string> & out, Opts & option)
{

		set <int>::iterator it = starts.begin();
		for(int j = 0; j <  pstart; j++) it++; //now it points to starts.begin() + pstart// do not understand why "it=starts.begin() + pstart" does not compile.

		map<int, MergeLaneElement> mergedstarts; // start -> merged lane elements
		map<int, vector< vector< int > > > tcmcs; // start -> <tcs, mcs>

		for(int k = 0; k < batchSize; k++){
			int start = *it;
			
		//Todo:It's redundant to do summations for both options. Please detach into two "for" statements.
		//for(set<int>::iterator pstart = starts.begin(); pstart != starts.end(); pstart++ ){
		//	int start = *pstart;
//			int tc;
//			int mc;
			char next = 'X';

			//data structure for methC and totalC count values
			//xx	total-1	plus-2	minus-3
			//r1	x		x		x
			//r2	x		x		x
			//r3	x		x		x
			map<int, map<char, int> > mcount; //laneId, strand ->mcount
			map<int, map<char, int> > tcount; //laneId, strand ->tcount
			for(int i = 0; i < lanesPlus.size(); i++){
				mcount[i]['B'] = 0;
				mcount[i]['+'] = 0;
				mcount[i]['-'] = 0;
				tcount[i]['B'] = 0;
				tcount[i]['+'] = 0;
				tcount[i]['-'] = 0;
			}

			int tcp = 0;
			int mcp = 0;

			for(unsigned int i = 0; i < lanesPlus.size(); i++ ){
				if(lanesPlus[i].count(chr) != 0 && lanesPlus[i][chr].count(start) != 0 && lanesPlus[i][chr][start].totalC > 0 )
				{
					tcp += lanesPlus[i][chr][start].totalC;
					mcp += lanesPlus[i][chr][start].methC;
					tcount[i]['+'] = lanesPlus[i][chr][start].totalC;
					mcount[i]['+'] = lanesPlus[i][chr][start].methC;
					next = lanesPlus[i][chr][start].nextN;
				}
			}

			int tcm = 0;
			int mcm = 0;

			for(unsigned int i = 0; i < lanesMinus.size(); i++ ){
				if(lanesPlus[i].count(chr) != 0 && lanesMinus[i][chr].count(start) != 0 && lanesMinus[i][chr][start].totalC > 0 )
				{
					tcm += lanesMinus[i][chr][start].totalC;
					mcm += lanesMinus[i][chr][start].methC;
					tcount[i]['-'] = lanesMinus[i][chr][start].totalC;
					mcount[i]['-'] = lanesMinus[i][chr][start].methC;
					next = lanesMinus[i][chr][start].nextN;
				}
			}
			
			int end = (next == 'G') ? (start + 2) : (start + 1);
			char strand = 'B';
			if(tcp == 0){ strand = '-'; }
			if(tcm == 0){ strand = '+'; }


			vector <int> tci;
			vector <int> mci;

			//wrong. because tci and mci may have different sizes
//			for(unsigned int i = 0; i < tcpi.size(); i++ ){
//				tci.push_back(tcpi[i]+tcmi[i]);
//				mci.push_back(mcpi[i]+mcmi[i]);
//			}

			//assign value to total column
			int tc = 0;
			for(int i = 0; i < lanesPlus.size(); i++){
				mcount[i]['B'] = mcount[i]['+'] + mcount[i]['-'];
				tcount[i]['B'] = tcount[i]['+'] + tcount[i]['-'];
				if(tcount[i]['B']>0){
					mci.push_back(mcount[i]['B']);
					tci.push_back(tcount[i]['B']);
					tc += tcount[i]['B'];
				}
			}
//			cout << "nRep = " << tci.size() << endl;
//
//			std::cout	 	<< setprecision(3) << chr << "\t" << start << "\t" << end << "\t" << double(mcp+mcm)/(tcp+tcm)
//							<< "\t" << tcp+tcm << "\t" << mcp+mcm << "\t" << strand << "\t" << next
//							<< "\t+\t" << tcp << "\t" << mcp << "\t-\t" << tcm << "\t" << mcm << endl;
			//for debug purpose
			//std::cerr << chr << "\t" << start << "\t" << end << endl;
			//int option_wi_variance = 1;

			if(option.withVariance && tc >= option.minDepthForComp){
				MergeLaneElement element;
				element.start=start;
				element.end=end;
				element.tcp=tcp;
				element.mcp=mcp;
				element.tcm=tcm;
				element.mcm=mcm;
				element.strand=strand;
				element.next=next;
				element.chr=chr;

				if(tci.size()>1){
					if (isconsensus(tci, mci)) { // consistent replicates, no need to fit BBF
						element.k=mcp+mcm;
						element.n=tcp+tcm;
					} else {
						tcmcs[start].push_back(tci);
						tcmcs[start].push_back(mci);
					}
				} else {
					element.k=mci[0];
					element.n=tci[0];
				}

				mergedstarts[start] = element;
			} else if(tc >= option.minDepthForComp) {
				MergeLaneElement element;
				element.start=start;
				element.end=end;
				element.k=mcp+mcm;
				element.n=tcp+tcm;
				element.tcp=tcp;
				element.mcp=mcp;
				element.tcm=tcm;
				element.mcm=mcm;
				element.strand=strand;
				element.next=next;
				element.chr=chr;

				mergedstarts[start] = element;
			} else {
				//less than minDepthForComp
			}
			it++;
		}

		map<int, BiKey > fits; // start -> BiKey
		CallBetaBinomialFit(tcmcs, fits); // Calling Beta-binomial fitting sequentially

		for(map<int, BiKey>::iterator it=fits.begin(); it!=fits.end(); ++it) {
			mergedstarts[it->first].k=it->second.k1;
			mergedstarts[it->first].n=it->second.n1;
		}

		for (map<int, MergeLaneElement> :: iterator it=mergedstarts.begin(); it!=mergedstarts.end(); ++it) {
			MergeLaneElement e=it->second;
			if (e.n>0) {
				stringstream mergedLaneFile;
				mergedLaneFile << setprecision(3) << e.chr << "\t" << e.start << "\t" << e.end << "\t" << double(e.k)/(e.n)
					<< "\t" << e.n << "\t" << e.k << "\t" << e.strand << "\t" << e.next
					<< "\t+\t" << e.tcp << "\t" << e.mcp << "\t-\t" << e.tcm << "\t" << e.mcm;
				string moreThanDepthStr = mergedLaneFile.str();
				if(!   moreThanDepthStr.empty() ) out.push_back(moreThanDepthStr);
			}
		}
}

// From readLaneToStrandSpecificHash
void getChroms(vector<string> & filesToMerge, set<string> & chroms) {
	for(unsigned int i=0; i<filesToMerge.size(); i++) {
		int count = 0;
		string line = "";
		ifstream inputf(filesToMerge[i].c_str(), ios::in);
		while (inputf.good()) {
			getline(inputf, line);
			if(line == ""){continue;}

			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if( count == 0 ) { //header line
			} else {				//content
				chroms.insert(fields[0]);
			}
			count += 1;
		}
		inputf.close();
	}
}

// copied from void mergeRatioFiles(vector<string> & filesToMerge, string outFileName, Opts & option)
void mergeRatioFiles(vector<string> & filesToMerge, string outFileName, Opts & option, string & chr)
{
	map<int, map <string, map<int, cMeth> > > lanesPlus ; //laneId->chrom->loc->(totalC, methC, strand, nextN); for + strand;
	map<int, map <string, map<int, cMeth> > > lanesMinus ; //laneId->chrom->loc->(totalC, methC, strand, nextN); for - strand;

	for(unsigned int i = 0; i < filesToMerge.size(); i++ ){
		cout << " start reading file " << filesToMerge[i] << " " << chr << endl;
		map <string, map<int, cMeth> > methPs;
		map <string, map<int, cMeth> > methMs;
		readLaneToStrandSpecificHash( filesToMerge[i], methPs, methMs, chr );
		lanesPlus[i] = methPs;
		lanesMinus[i] = methMs;
	}
	cout << " finish reading" << endl;
	
	ofstream mergedLaneFile;
	mergedLaneFile.open(outFileName.c_str(), ios_base::out);

	//get all locations for current chrom
	set <int> starts;
	for(unsigned int i = 0; i < lanesPlus.size(); i++ ){
		for(map<int, cMeth>::iterator it = lanesPlus[i][chr].begin(); it != lanesPlus[i][chr].end(); ++it) {
			starts.insert(it->first);
		}
	}
	for(unsigned int i = 0; i < lanesMinus.size(); i++ ){
		for(map<int, cMeth>::iterator it = lanesMinus[i][chr].begin(); it != lanesMinus[i][chr].end(); ++it) {
			starts.insert(it->first);
		}
	}

	bool endOfChr = false;
	int pstart = 0;
	int batchSize = BATCHMAX;
	boost::thread_group mergeio;
	//this way of threading makes sure the output is still sorted
	while( ! endOfChr )
	{
		boost::thread_group g;
		vector < vector<string> > out;
		out.resize(option.threads);
		int thisStart = pstart;

		int i = 0;
		for( ; i < option.threads;  )
		{
			if(starts.size() - pstart < BATCHMAX)
			{
				batchSize = starts.size() - pstart;
				boost::thread *tp = new boost::thread( mergeRatioFilesWorker, boost::ref(lanesPlus), boost::ref(lanesMinus),chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]), option );
				g.add_thread(tp);
				pstart = pstart + batchSize;
				i++;
				endOfChr = true;
				break;//no more location in current chrom. breaks out of inner most for loop
			} else
			{
				boost::thread *tp = new boost::thread( mergeRatioFilesWorker, boost::ref(lanesPlus), boost::ref(lanesMinus),chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]), option );
				g.add_thread(tp);
				pstart = pstart + batchSize;
				i++;
			}
		}

		vector<string> rout;

		g.join_all();

		//put both nonR and R calculation to file in a backgroud thread
		mergeio.join_all();
		//cout << "joined output thread" << endl;

		boost::thread *tprout = new boost::thread( appendToFile, boost::ref(mergedLaneFile), out, rout);
		mergeio.add_thread(tprout);
		//cout << "added output thread" << endl;
	}
	///////////////////////////////////////////
	mergeio.join_all();
	mergedLaneFile.close();
}

void mergeRatioFilesByChrom(vector<string> & filesToMerge, string outFileName, Opts & option)
{
	set<string> chroms;
	getChroms(filesToMerge, chroms);
	vector< string > outfiles;
	for (set<string>::iterator pchr=chroms.begin(); pchr!=chroms.end(); pchr++)
	{
		string chr = *pchr;
		string outfile=outFileName+"_"+chr+".bed";
		if (! boost::filesystem::exists(outfile)) {
			mergeRatioFiles(filesToMerge, outfile, option, chr);
		}
		outfiles.push_back(outfile);
	}

	ofstream mergedLaneFile;
	mergedLaneFile.open(outFileName.c_str(), ios_base::out);
	mergedLaneFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttcP\tmcP\tMinus\ttcM\tmcM" << endl;
	mergedLaneFile.close();

	string sysCmd;
	sysCmd = "cat " + join_vec(outfiles) + " >> " + outFileName;
	system(sysCmd.c_str());
	sysCmd = "rm -f " + join_vec(outfiles);
	system(sysCmd.c_str());
}


void mergeRatioFiles(vector<string> & filesToMerge, string outFileName, Opts & option)
{
	map<int, map <string, map<int, cMeth> > > lanesPlus ; //laneId->chrom->loc->(totalC, methC, strand, nextN); for + strand;
	map<int, map <string, map<int, cMeth> > > lanesMinus ; //laneId->chrom->loc->(totalC, methC, strand, nextN); for - strand;

	//not used as of now
	//map<int, map <string, map<int, cMeth> > > mergedLane ; //strand->chrom->loc->(totalC, methC, strand, nextN); for + strand;
	ofstream mergedLaneFile;
	//configFile.open("run.config", ios_base::app);
	mergedLaneFile.open(outFileName.c_str(), ios_base::out);
	mergedLaneFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttcP\tmcP\tMinus\ttcM\tmcM" << endl;

	for(unsigned int i = 0; i < filesToMerge.size(); i++ ){
		cout << " start reading file " << filesToMerge[i] << endl;
		map <string, map<int, cMeth> > methPs;
		map <string, map<int, cMeth> > methMs;
		readLaneToStrandSpecificHash( filesToMerge[i], methPs, methMs );
		lanesPlus[i] = methPs;
		lanesMinus[i] = methMs;
	}
	cout << " finish reading" << endl;
	
	
	//prepare the genomic locations that exist in at least one lane
	set <string> chroms;

	for(unsigned int i = 0; i < lanesPlus.size(); i++ ){
		for(map<string, map<int, cMeth> >::iterator it = lanesPlus[i].begin(); it != lanesPlus[i].end(); ++it) {
		  chroms.insert(it->first);
		}
	}
	for(unsigned int i = 0; i < lanesMinus.size(); i++ ){
		for(map<string, map<int, cMeth> >::iterator it = lanesMinus[i].begin(); it != lanesMinus[i].end(); ++it) {
		  chroms.insert(it->first);
		}
	}

	for (set<string>::iterator pchr=chroms.begin(); pchr!=chroms.end(); pchr++)
	{
		string chr = *pchr;
		//cout << "start chr" << endl;

		//get all locations for current chrom
		set <int> starts;
		for(unsigned int i = 0; i < lanesPlus.size(); i++ ){
			for(map<int, cMeth>::iterator it = lanesPlus[i][chr].begin(); it != lanesPlus[i][chr].end(); ++it) {
				starts.insert(it->first);
			}
		}
		for(unsigned int i = 0; i < lanesMinus.size(); i++ ){
			for(map<int, cMeth>::iterator it = lanesMinus[i][chr].begin(); it != lanesMinus[i][chr].end(); ++it) {
				starts.insert(it->first);
			}
		}

		///////////////////////////////////////////
		////////////////////////copied from docomp function. So there are marks for nonR calculation and similar stuff 
		//mergeRatioFilesWorker();
		bool endOfChr = false;
		int pstart = 0;
		int batchSize = BATCHMAX;
		boost::thread_group mergeio;
		//this way of threading makes sure the output is still sorted
		while( ! endOfChr )
		{
			boost::thread_group g;
			vector < vector<string> > out;
			out.resize(option.threads);
			int thisStart = pstart;

			//start nonR calculation on threads
			//cout << "starting subthreads" << endl;
			int i = 0;
			for( ; i < option.threads;  )
			{
				if(starts.size() - pstart < BATCHMAX)
				{
					//cout << i << endl;
					batchSize = starts.size() - pstart;
					boost::thread *tp = new boost::thread( mergeRatioFilesWorker, boost::ref(lanesPlus), boost::ref(lanesMinus),chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]), option );
					g.add_thread(tp);
					//compBatchLoc(lane, chr, starts, pstart, batchSize, outStr[i]);
					pstart = pstart + batchSize;
					i++;
					//cout << i << endl;
					//cout << g.size() << endl;
					endOfChr = true;
					break;//no more location in current chrom. breaks out of inner most for loop
				} else
				{
					//cout << i << endl;
					boost::thread *tp = new boost::thread( mergeRatioFilesWorker, boost::ref(lanesPlus), boost::ref(lanesMinus),chr, boost::ref(starts), pstart, batchSize, boost::ref(out[i]), option );
					g.add_thread(tp);
					//compBatchLoc(lane, chr, starts, pstart, batchSize, outStr[i]);
					pstart = pstart + batchSize;
					i++;
					//cout << i << endl;
					//cout << g.size() << endl;
				}
			}

			//start R calculation on main thread while nonR calculation is on backgroup thread
			vector<string> rout;

			g.join_all();
			//cout << "sub trheads joined......waiting to join output thread" << endl;


			//put both nonR and R calculation to file in a backgroud thread
			mergeio.join_all();
			//cout << "joined output thread" << endl;

			boost::thread *tprout = new boost::thread( appendToFile, boost::ref(mergedLaneFile), out, rout);
			mergeio.add_thread(tprout);
			//cout << "added output thread" << endl;
		}
		///////////////////////////////////////////
		mergeio.join_all();
	}
	mergedLaneFile.close();
}

//x=(x1*x2***xn)^(1/n)
long double geometricMean(vector<long double> & v){
	long double result = 0.0;
	for(unsigned int i=0;i < v.size(); i++){
		if(v[i]>0){
			result += log10(v[i]);
		} else {
			cerr << "can not be negative " << v[i] << endl;
			exit(1);
		}
	}
	result = pow((long double)10.0, result / v.size() );
	return result;
}

//x=(x1*x2***xn)^(1/n). returned is log10(x)
long double log10GeometricMean(vector<long double> & v){
	long double result = 0.0;
	for(unsigned int i=0;i < v.size(); i++){
		if(v[i]>0){
			result += log10(v[i]);
		} else {
			cerr << "can not be negative " << v[i] << endl;
			exit(1);
		}
	}
	result = result / v.size() ;
	return result;
}
//returns x = x1*x2***xn
long double mutiP(vector<long double> & v){
	long double result = 0.0;
	for(unsigned int i=0;i < v.size(); i++){
		if(v[i]>0){
			result += log10(v[i]);
		} else {
			cerr << "can not be negative " << v[i] << endl;
			exit(1);
		}
	}
	result = pow((long double)10.0, result );
	return result;
}
//Strategy 1: search DMR over genome
//start with file F1 that has all the quantities calculated.
//dmc strategies to classify dmcs and specify the identifier for each record in a file. This file F2 could contain only dmcs.
//use file F1 and F2 to get all cpg and dmc info, search dmr
//Strategy 2: calculate quantities for predefined region
//File F1 is needed, DMC info is not needed.
//file F3 for the prefdefined region is needed.

void readCompToHashM1(string fileName, int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, Opts & option ){

	long double genomeSimP = 0.0;
	int cytosinesCountGenome = 0;

	stringstream fakeLaneId;
	fakeLaneId << "100" << i << "100" << j;
	int laneIdFake = string_to_int(fakeLaneId.str());
	int chrIndex = -1;
	int startIndex = -1;
	int endIndex = -1;
	int t1Index = -1;
	int r1Index = -1;
	int ci1Index = -1;
	int t2Index = -1;
	int r2Index = -1;
	int ci2Index = -1;
	int nDifIndex = -1;
	int cDifIndex = -1;
	int difCIIndex = -1;
	int simPIndex = -1;
	int fetPIndex = -1;

	string line = "";
	ifstream inputf(fileName.c_str(), ios::in);
	ofstream coveredCpgs;
	coveredCpgs.open((fileName + ".cvd.txt").c_str(), ios_base::out);
	while (inputf.good()) {
		getline(inputf, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));

		if( find(line.begin(), line.begin()+1, '#') == line.begin() )
		{
			chrIndex = find (fields.begin(), fields.end(), "#chrom") - fields.begin();
			startIndex = find (fields.begin(), fields.end(), "start") - fields.begin();
			endIndex = find (fields.begin(), fields.end(), "end") - fields.begin();
			t1Index = find (fields.begin(), fields.end(), "totalC_" + itos(i)) - fields.begin();
			r1Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(i)) - fields.begin();
			ci1Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(i)) - fields.begin();
			t2Index = find (fields.begin(), fields.end(), "totalC_" + itos(j)) - fields.begin();
			r2Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(j)) - fields.begin();
			ci2Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(j)) - fields.begin();
			nDifIndex = find (fields.begin(), fields.end(), "nominalDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			cDifIndex = find (fields.begin(), fields.end(), "credibleDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			difCIIndex = find (fields.begin(), fields.end(), "difCI_" + itos(j)+ "-" + itos(i)) - fields.begin();
			simPIndex = find (fields.begin(), fields.end(), "p_sim_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			fetPIndex = find (fields.begin(), fields.end(), "p_fet_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			//cout << i << "\t" << j <<"\t"<< chrIndex << "\t" << startIndex << "\t" <<t1Index<< "\t" <<t2Index <<endl;
			//cout << i << "\t" << j <<"\t"<< r1Index << "\t" << r2Index << "\t" <<ci1Index<< "\t" <<ci2Index <<"\t" <<fetPIndex<<endl;
		}
		else{
			if(fields[fetPIndex] == "NA"){continue;} //if fetP field is "NA" then next;
			//cout << i << ".." << j << ".." << endl;
			string chr = fields[chrIndex];
			int start = string_to_int(fields[startIndex]);
			int end = string_to_int(fields[endIndex]);
			//string loc = chr + "\t" + fields[startIndex] + "\t" + fields[endIndex];
			int t1 = string_to_int(fields[t1Index]);
			//cout << loc << "\t" << t1 ;
			double r1 = atof((fields[r1Index]).c_str());
			//cout << "\t" << r1 ;
			string ci1 = fields[ci1Index];
			//cout << "\t" << ci1 <<endl;
			int t2 = string_to_int(fields[t2Index]);
			double r2 = atof((fields[r2Index]).c_str());
			string ci2 = fields[ci2Index];
			double simP = atof((fields[simPIndex]).c_str());
			//not a bug but to be fixed
			if(simP == 0){simP = 0.00000000000000001;}
			//

			double fetP = atof((fields[fetPIndex]).c_str());
			double nDif = atof((fields[nDifIndex]).c_str());
			double cDif = atof((fields[cDifIndex]).c_str());
			string difCI = fields[difCIIndex];

			//vector <string> ab;
			//boost::split(ab, difCI, boost::is_any_of(","));
			//double cDif = abs(atof(ab[0].c_str())) < abs(atof(ab[1].c_str())) ? atof(ab[0].c_str()) : atof(ab[1].c_str());
			string cls = "insigOrLowDepth";
			if( (fetP < option.pFetDmc) && (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp) ){
				if(nDif > 0){
					cls = nDif > option.minNominalDif ? "strongHyper" : "hyper";
				} else {
					cls = abs(nDif) > option.minNominalDif ? "strongHypo" : "hypo";
				}
			}
			if( (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp) ){
				pairResult[laneIdFake][chr][start] = PairCompOut(chr, start, end, t1, r1, ci1, t2, r2, ci2, simP, fetP, nDif, difCI, cDif, cls);
				coveredCpgs << pairResult[laneIdFake][chr][start].ToString() << endl;
				genomeSimP += log10(simP);
				//cout << simP << "\t" << genomeSimP << endl;
				cytosinesCountGenome += 1;
			}
		}
	}
	cout << "done" << i << ".." <<j <<endl;
	inputf.close();
	coveredCpgs.close();
}

void readCompToHashM2(string fileName, int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, statSim & genomeSim, Opts & option ){

	long double genomeSimP = 0.0;
	int cytosinesCountGenome = 0;

	stringstream fakeLaneId;
	fakeLaneId << "100" << i << "100" << j;
	int laneIdFake = string_to_int(fakeLaneId.str());
	int chrIndex = -1;
	int startIndex = -1;
	int endIndex = -1;
	int t1Index = -1;
	int r1Index = -1;
	int ci1Index = -1;
	int t2Index = -1;
	int r2Index = -1;
	int ci2Index = -1;
	int nDifIndex = -1;
	int cDifIndex = -1;
	int difCIIndex = -1;
	int simPIndex = -1;
	int fetPIndex = -1;

	string line = "";
	ifstream inputf(fileName.c_str(), ios::in);
	ofstream coveredCpgs;
	coveredCpgs.open((fileName + ".cvd.txt").c_str(), ios_base::out);
	while (inputf.good()) {
		getline(inputf, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));

		if( find(line.begin(), line.begin()+1, '#') == line.begin() )
		{
			chrIndex = find (fields.begin(), fields.end(), "#chrom") - fields.begin();
			startIndex = find (fields.begin(), fields.end(), "start") - fields.begin();
			endIndex = find (fields.begin(), fields.end(), "end") - fields.begin();
			t1Index = find (fields.begin(), fields.end(), "totalC_" + itos(i)) - fields.begin();
			r1Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(i)) - fields.begin();
			ci1Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(i)) - fields.begin();
			t2Index = find (fields.begin(), fields.end(), "totalC_" + itos(j)) - fields.begin();
			r2Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(j)) - fields.begin();
			ci2Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(j)) - fields.begin();
			nDifIndex = find (fields.begin(), fields.end(), "nominalDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			cDifIndex = find (fields.begin(), fields.end(), "credibleDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			difCIIndex = find (fields.begin(), fields.end(), "difCI_" + itos(j)+ "-" + itos(i)) - fields.begin();
			simPIndex = find (fields.begin(), fields.end(), "p_sim_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			fetPIndex = find (fields.begin(), fields.end(), "p_fet_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			//cout << i << "\t" << j <<"\t"<< chrIndex << "\t" << startIndex << "\t" <<t1Index<< "\t" <<t2Index <<endl;
			//cout << i << "\t" << j <<"\t"<< r1Index << "\t" << r2Index << "\t" <<ci1Index<< "\t" <<ci2Index <<"\t" <<fetPIndex<<endl;
		}
		else{
			if(fields[fetPIndex] == "NA"){continue;} //if fetP field is "NA" then next;
			//cout << i << ".." << j << ".." << endl;
			string chr = fields[chrIndex];
			int start = string_to_int(fields[startIndex]);
			int end = string_to_int(fields[endIndex]);
			//string loc = chr + "\t" + fields[startIndex] + "\t" + fields[endIndex];
			int t1 = string_to_int(fields[t1Index]);
			//cout << loc << "\t" << t1 ;
			double r1 = atof((fields[r1Index]).c_str());
			//cout << "\t" << r1 ;
			string ci1 = fields[ci1Index];
			//cout << "\t" << ci1 <<endl;
			int t2 = string_to_int(fields[t2Index]);
			double r2 = atof((fields[r2Index]).c_str());
			string ci2 = fields[ci2Index];
			double simP = atof((fields[simPIndex]).c_str());
			//not a bug but to be fixed
			if(simP == 0){simP = 0.00000000000000001;}
			//

			double fetP = atof((fields[fetPIndex]).c_str());
			double nDif = atof((fields[nDifIndex]).c_str());
			double cDif = atof((fields[cDifIndex]).c_str());
			string difCI = fields[difCIIndex];

			//vector <string> ab;
			//boost::split(ab, difCI, boost::is_any_of(","));
			//double cDif = abs(atof(ab[0].c_str())) < abs(atof(ab[1].c_str())) ? atof(ab[0].c_str()) : atof(ab[1].c_str());
			string cls = "insigOrLowDepth";
			if( (simP < option.pSimDmc) && (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp)
			 //TODO && rank ? need put this cls part outside after the hash is created
				){
				if(cDif > 0){
					cls = cDif > option.minCredibleDif ? "strongHyper" : "hyper";
				} else {
					cls = abs(cDif) > option.minCredibleDif ? "strongHypo" : "hypo";
				}
			}
			if( (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp) ){
				pairResult[laneIdFake][chr][start] = PairCompOut(chr, start, end, t1, r1, ci1, t2, r2, ci2, simP, fetP, nDif, difCI, cDif, cls);
				if(option.doDmcScan) coveredCpgs << pairResult[laneIdFake][chr][start].ToString() << endl;
				genomeSimP += log10(simP);
				//cout << simP << "\t" << genomeSimP << endl;
				cytosinesCountGenome += 1;
			}
		}
	}
	genomeSim.cytosinesCountGenome = cytosinesCountGenome;
	genomeSim.genomeSimP = genomeSimP;
	cout << "done" << i << ".." <<j <<endl;
	inputf.close();
	coveredCpgs.close();
}


//call Dmc and Dmr according to method 1
//Dmc is defined as depth >= minDepth(suggest 10), p<p_fet_cut(suggest 0.05), nMethDif > nominalMethDifCutoff(0.3333)
//Dmr is defined as minC >=3, maxDist 300, p<p_fet_cut(suggest 0.05), nMethDif > nominalMethDifCutoff(0.3333)
void DmcDmrM1(int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, Opts & option){
	string label = option.labels[i] + "_vs_" + option.labels[j];
	ofstream dmcFile;
	string dmcName = ("dmc_M1_" + label );
	dmcFile.open((dmcName + ".txt").c_str(), ios_base::out);
	dmcFile << "#chrom\tstart\tend"
			<< "\t" << "totalC_" << i << "\t" << "nominalRatio_" << i << "\t" << "ratioCI_" << i
			<< "\t" << "totalC_" << j << "\t" << "nominalRatio_" << j << "\t" << "ratioCI_" << j
			<< "\t" << "nominalDif_" << j << "-" << i << "\t" << "credibleDif_" << j << "-" << i
			<< "\t" << "difCI_" << j << "-" << i
			<< "\t" << "p_sim_" << j << "_v_" << i << "\t" << "p_fet_" << j << "_v_" << i
			<< "\t" << "class" << endl;
	ofstream dmcTrack;
	dmcTrack.open((dmcName + ".bed").c_str(), ios_base::out);
	string trackHeadDmc = "track name=\"" + label + "_DMC\" description=\"" + label +"_DMC\" visibility=3 useScore=1 itemRgb=On";
	dmcTrack << trackHeadDmc << endl;

	map <int, map <string, map<int, int> > > dmcCombinedReg; //fakeLaneId->chr->start ==> end
																//note that the end is only the position of C on + strand
	map <int, map <string, map<int, vector<int> > > > allLocsInReg ; //the cytosines inside dmcCombinedReg;
	map <int, map <string, map<int, vector<int> > > > dmcLocsInReg ; //the DMCs inside dmcCombinedReg;
	for(map <int, map <string, map<int, PairCompOut> > >::iterator id = pairResult.begin();id != pairResult.end(); id++)
	{
		int laneId = id->first;
		map <string, map<int, PairCompOut> > thisLane = id->second;
		for(map <string, map<int, PairCompOut> >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr++)
		{
			string chr = pchr->first;
			map<int, PairCompOut> thisChr = pchr->second;

			int lastStart = 0;
			string lastClass = "";
			int consecDmc = 0;
			int startOfcombinedReg = 0;
			int endOfCombinedReg = 0;
			vector <int> allLocs;
			vector <int> dmcLocs;
			for(map<int, PairCompOut>::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++)
			{
				int start = pt->first;
				PairCompOut loc = pt->second;
				if(option.doDmcScan) dmcFile << loc.ToString() << endl;
				if( (loc.cls.find("strong") != string::npos) ) // isDmc.
				{
					if(option.doDmcScan) dmcTrack << loc.ToBedTrackString() << endl;
				}
				//cout << i << ".........." << j << "........." << start << endl;
				if( loc.ta >= option.minDepthForComp && loc.tb >= option.minDepthForComp)
				{
					if( (loc.cls.find("strong") != string::npos) ) // isDmc.
					{
						if(consecDmc == 0){		//initilize a combined region
							consecDmc = 1;
							startOfcombinedReg = start;
							endOfCombinedReg = start;
							allLocs.push_back(start);
							dmcLocs.push_back(start);
						}
						else	//check if this dmc belongs to current region or a new region
						{
							allLocs.push_back(start);
							dmcLocs.push_back(start);
							if( (loc.cls == lastClass) && (start - lastStart < option.maxDistConsDmcs)) //extend region//bug: what if the last element of current chr
							{
								endOfCombinedReg = start;
								consecDmc ++;
							} else {	//finish the combined region
								if(consecDmc >= option.minDmcsInDmr){ //then this is a candidate region
									dmcCombinedReg[laneId][chr][startOfcombinedReg] = endOfCombinedReg;
									//cout << startOfcombinedReg << "\t" << endOfCombinedReg <<endl;
									while (allLocs.back() != endOfCombinedReg) allLocs.pop_back();
									while (dmcLocs.back() != endOfCombinedReg) dmcLocs.pop_back();

									//copy(allLocs.begin(), allLocs.end(), ostream_iterator<int>(cout, ","));
									//cout << allLocs.size() << ".." << consecDmc << "" << endl;

									allLocsInReg[laneId][chr][startOfcombinedReg] = allLocs;
									dmcLocsInReg[laneId][chr][startOfcombinedReg] = dmcLocs;
									//vector <Dmr> dmrs;
									//getDmr(chr, allLocs, thisChr, dmrs );
								}

								consecDmc = 1; //new region starts
								allLocs.clear();
								dmcLocs.clear();
								startOfcombinedReg = start;
								endOfCombinedReg = start;
								allLocs.push_back(start);
								dmcLocs.push_back(start);
							}
						}
						lastStart = start;
						lastClass = loc.cls;
					} else {
						if(consecDmc >= 1)
						{
							allLocs.push_back(start);
						}
					}
					//finish this location
				}
			}
			//finish this chrom
			//may start call dmr now.
			//getDmr()
		}
	}
	dmcFile.close();
	dmcTrack.close();
	string chromsizes = option.reference.substr(0, option.reference.find_last_of(".")) + ".chrom.sizes"; 
	string dmcBigbed = "egrep -v '^track|^browser|^#' " + (dmcName + ".bed") + " > " + dmcName +".temp" + " && bedToBigBed " + dmcName + ".temp " + chromsizes + " " + dmcName + ".bb" + " && rm " + dmcName + ".temp";
	cout << "running: " << dmcBigbed << endl;
	system(dmcBigbed.c_str());
	cout << i << "...." << j <<endl;
	//done with DMC call

	int debug = 0;
	if(debug){
		//filter with coverage condition
		//this is actually testing . the coverage condition is already included in reading procedures
		map <int, map <string, map<int, vector<int> > > > goodLocsInReg ; //the cytosines inside dmcCombinedReg and satisfy the depth condition
		for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
			int laneId = id->first;
			map <string, map<int, vector<int> > > thisLane = id->second;
			for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
				string chr = pchr->first;
				map<int, vector<int> > thisChr = pchr->second;
				for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
					int start = pt->first;
					vector<int> allLocs = allLocsInReg[laneId][chr][start];
					for(vector<int>::iterator l = allLocs.begin(); l != allLocs.end(); l ++){
						int location = *l;
						if(pairResult[laneId][chr][location].ta < option.minDepthForComp || pairResult[laneId][chr][location].tb < option.minDepthForComp ){
							cerr << "ERRRRRRRRRRRR here" << endl;
							cout << chr << "\t" << start << "\t" << location << endl;
							//cout << pairResult[laneId][chr][location].ToString() << endl;
							exit(1);
						}
					}
					goodLocsInReg[laneId][chr][start] = allLocs;
				}
			}
		}
	}

	cout << "starting to call DMR" << endl;
	string dmrName = ("dmr_M1_" + label );
	ofstream dmrTable;
	dmrTable.open((dmrName + ".txt").c_str(), ios_base::out);
	dmrTable << "#chrom\tstart\tend"
			<< "\t" << "meanRatio_" << i << "\t" << "totalC_" << i << "\t" << "cSites_" << i
			<< "\t" << "meanRatio_" << j << "\t" << "totalC_" << j << "\t" << "cSites_" << j
			<< "\t" << "methDif_" << j << "-" << i
			<< "\t" << "p_" << j << "_v_" << i << "\t" << "class_" << j << "_v_" << i << endl;

	ofstream dmrTrack;
	dmrTrack.open((dmrName + ".bed").c_str(), ios_base::out);
	string trackHead = "track name=\"" + label + "_ScanDMR\" description=\"" + label +"_ScanDMR\" visibility=3 useScore=1 itemRgb=On";
	dmrTrack << trackHead << endl;
	//start to call DMR between lane i and j. //read only relevant info to save RAM.
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		//cout << laneId << endl;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			//cout << chr << endl;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> allLocs = allLocsInReg[laneId][chr][start];
				vector<int> dmcLocs = dmcLocsInReg[laneId][chr][start];
				//cout << "good=" << start << "\t" << dmcLocs.size() << endl;
				Dmr mostStrongDif;

				//start calling DMR for the region in laneId->chr->start
				int stop = 0;
				for(unsigned int dmcWindow = dmcLocs.size(); dmcWindow >= option.minDmcsInDmr && stop == 0; dmcWindow --){
					vector <Dmr> subReg;
					for(unsigned int dmcIndex = 0; dmcIndex + dmcWindow - 1 < dmcLocs.size(); dmcIndex ++){
						int windowLt = dmcLocs[dmcIndex];
						int windowRt = dmcLocs[dmcIndex + dmcWindow - 1];

						//cout << "this window" << endl;
						//get methdif and pvalue for this window from windowLt to windowRt
						double windowMethDif = 0;
						int ma = 0; //methylated cytosines in sample a
						int ua = 0; //unmethylated cytosines in sample b
						int mb = 0;
						int ub = 0;
						double mra = 0.0; //mean ratio in sample a
						double mrb = 0.0; //mean ratio in sample b
						int count = 0;
						for(unsigned int allIndex = 0; allIndex < allLocs.size(); allIndex ++){
							if(allLocs[allIndex] >= windowLt && allLocs[allIndex] <= windowRt){
								PairCompOut thisLoc = pairResult[laneId][chr][allLocs[allIndex]];
								windowMethDif += thisLoc.nDif;
								count += 1;
								mra += thisLoc.ra;
								mrb += thisLoc.rb;
								ma = ma + int (thisLoc.ta * thisLoc.ra + 0.5);
								ua = ua + int (thisLoc.ta * (1 - thisLoc.ra) + 0.5);
								mb = mb + int (thisLoc.tb * thisLoc.rb + 0.5);
								ub = ub + int (thisLoc.tb * (1 - thisLoc.rb) + 0.5);
							}
						}
						windowMethDif = windowMethDif / double( count );
						mra = mra / double( count );
						mrb = mrb / double( count );
						double p = 1.0;
						int varray[] = {ma, ua, mb, ub};
						//cout << "..." << ma << "\t" << ua << "\t" << mb << "\t" << ub << endl;
						fet_22( p, varray, 2, 2 );
						string cls = "insigOrLowDepth";
						if( (p < option.pFetDmr) ){
						//if(1){
							if(windowMethDif > 0){
								cls = windowMethDif > option.minNominalDif ? "strongHyper" : "hyper";
							} else {
								cls = abs(windowMethDif) > option.minNominalDif ? "strongHypo" : "hypo";
							}
						}
						//cout << "here " << p << endl;
						//Dmr(string Chr, int Start, int End, int AllSites, int DmcSites, int MC1, int UC1, int MC2, int UC2, double MethDif, double FetP)
						Dmr thisWindow(chr, windowLt, windowRt, count, dmcWindow, ma, ua, mb, ub, mra, mrb, windowMethDif, p, cls);
						subReg.push_back(thisWindow);
					}

					//determine the most strong region: #dmcWindow as 1st priority, then |#methDiff| >0.33
					sort(subReg.begin(), subReg.end(), sortDmrByMethDif);
					Dmr topReg = subReg[0];
					if(abs(topReg.methDif) >= option.minNominalDif && topReg.fetP < option.pFetDmr){
					//if(abs(topReg.methDif) >= option.minCredibleDif ){
						mostStrongDif = topReg;
						stop = 1;
					}
				}
				//dmrTable << "ss:" << start << endl;
				//done with call DMR for this candidate region
				if(stop == 1) {
					dmrTable << mostStrongDif.ToString() << endl;
					dmrTrack << mostStrongDif.ToBedTrackString() << endl;
				}
				//cout << "ss:" << start << endl;
			}

		}
	}
	dmrTable.close();
	dmrTrack.close();
	chromsizes = option.reference.substr(0, option.reference.find_last_of(".")) + ".chrom.sizes"; 
	string dmrBigbed = "egrep -v '^track|^browser|^#' " + (dmrName + ".bed") + " > " + dmrName +".temp" + " && bedToBigBed " + dmrName + ".temp " + chromsizes + " " + dmrName + ".bb";
	cout << "running: " << dmrBigbed << endl;
	system(dmrBigbed.c_str());
	//combine DMCs first


	debug = 0;
	if(debug){
	cout << "size 1 " << dmcCombinedReg.size() << endl;
	cout << "size 2 " << allLocsInReg.size() << endl;
	ofstream yFile;
	yFile.open("tempxxxyyy", ios_base::out);
	for(map <int, map <string, map<int, int> > >::iterator id = dmcCombinedReg.begin(); id != dmcCombinedReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, int> > thisLane = id->second;
		for(map <string, map<int, int> >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, int> thisChr = pchr->second;
			for(map<int,int>::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				int end = pt->second;

				yFile << chr << "\t" << start << "\t" << end << endl;
			}
		}
	}
	yFile.close();
	ofstream zFile;
	zFile.open("tempxxxzzz", ios_base::out);
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> ends = pt->second;
				for(unsigned int i = 0; i < ends.size(); i ++){
					PairCompOut tt = pairResult[laneId][chr][ends[i]];
					zFile << tt.ToString() << endl;
				}

			}
		}
	}
	zFile.close();
	ofstream aFile;
	aFile.open("tempxxxaaa", ios_base::out);
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = dmcLocsInReg.begin(); id != dmcLocsInReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> ends = pt->second;
				for(unsigned int i = 0; i < ends.size(); i ++){
					PairCompOut tt = pairResult[laneId][chr][ends[i]];
					aFile << tt.ToString() << endl;
				}

			}
		}
	}
	aFile.close();
	}

}

//call Dmc and Dmr according to method 2
//Dmc is defined as cMethDif > credibleMethDifCutoff(0.2)
//Dmr is defined as minC >=3, maxDist 300, cMethDif > credibleMethDifCutoff(0.2)
void DmcDmrM2(int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, statSim & genomeSim, Opts & option){
	long double genomeSimP = genomeSim.genomeSimP;
	int cytosinesCountGenome = genomeSim.cytosinesCountGenome;
	long double avgSim = genomeSimP / cytosinesCountGenome;
	
	string label = option.labels[i] + "_vs_" + option.labels[j];
	ofstream dmcFile;
	string dmcName = ("dmc_M2_" + label );
	dmcFile.open((dmcName + ".txt").c_str(), ios_base::out);
	dmcFile << "#chrom\tstart\tend"
			<< "\t" << "totalC_" << i << "\t" << "nominalRatio_" << i << "\t" << "ratioCI_" << i
			<< "\t" << "totalC_" << j << "\t" << "nominalRatio_" << j << "\t" << "ratioCI_" << j
			<< "\t" << "nominalDif_" << j << "-" << i << "\t" << "credibleDif_" << j << "-" << i
			<< "\t" << "difCI_" << j << "-" << i
			<< "\t" << "p_sim_" << j << "_v_" << i << "\t" << "p_fet_" << j << "_v_" << i
			<< "\t" << "class" << endl;
	ofstream dmcTrack;
	dmcTrack.open((dmcName + ".bed").c_str(), ios_base::out);
	string trackHeadDmc = "track name=\"" + label + "_DMC\" description=\"" + label +"_DMC\" visibility=3 useScore=1 itemRgb=On";
	dmcTrack << trackHeadDmc << endl;

	map <int, map <string, map<int, int> > > dmcCombinedReg; //fakeLaneId->chr->start ==> end
																//note that the end is only the position of C on + strand
	map <int, map <string, map<int, vector<int> > > > allLocsInReg ; //the cytosines inside dmcCombinedReg;
	map <int, map <string, map<int, vector<int> > > > dmcLocsInReg ; //the DMCs inside dmcCombinedReg;
	for(map <int, map <string, map<int, PairCompOut> > >::iterator id = pairResult.begin();id != pairResult.end(); id++)
	{
		int laneId = id->first;
		map <string, map<int, PairCompOut> > thisLane = id->second;
		for(map <string, map<int, PairCompOut> >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr++)
		{
			string chr = pchr->first;
			map<int, PairCompOut> thisChr = pchr->second;

			int lastStart = 0;
			string lastClass = "";
			int consecDmc = 0;
			int startOfcombinedReg = 0;
			int endOfCombinedReg = 0;
			vector <int> allLocs;
			vector <int> dmcLocs;
			for(map<int, PairCompOut>::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++)
			{
				int start = pt->first;
				PairCompOut loc = pt->second;
				if(option.doDmcScan) dmcFile << loc.ToString() << endl;
				if( (loc.cls.find("strong") != string::npos) ) // isDmc.
				{
					if(option.doDmcScan) dmcTrack << loc.ToBedTrackString() << endl;
				}
				//cout << i << ".........." << j << "........." << start << endl;
				//bug: 20121004: dsun: the use of cls(read from file) makes it a bug for recalculation of DMR using different criteria ?
				if( loc.ta >= option.minDepthForComp && loc.tb >= option.minDepthForComp)
				{
					if( (loc.cls.find("strong") != string::npos) ) // isDmc.
					{
						if(consecDmc == 0){		//initilize a combined region
							consecDmc = 1;
							startOfcombinedReg = start;
							endOfCombinedReg = start;
							allLocs.push_back(start);
							dmcLocs.push_back(start);
						}
						else	//check if this dmc belongs to current region or a new region
						{
							allLocs.push_back(start);
							dmcLocs.push_back(start);
							if( (loc.cls == lastClass) && (start - lastStart < option.maxDistConsDmcs)) //extend region//bug: what if the last element of current chr
							{
								endOfCombinedReg = start;
								consecDmc ++;
							} else {	//finish the combined region
								if(consecDmc >= option.minDmcsInDmr){ //then this is a candidate region
									dmcCombinedReg[laneId][chr][startOfcombinedReg] = endOfCombinedReg;
									//cout << startOfcombinedReg << "\t" << endOfCombinedReg <<endl;
									while (allLocs.back() != endOfCombinedReg) allLocs.pop_back();
									while (dmcLocs.back() != endOfCombinedReg) dmcLocs.pop_back();

									//copy(allLocs.begin(), allLocs.end(), ostream_iterator<int>(cout, ","));
									//cout << allLocs.size() << ".." << consecDmc << "" << endl;

									allLocsInReg[laneId][chr][startOfcombinedReg] = allLocs;
									dmcLocsInReg[laneId][chr][startOfcombinedReg] = dmcLocs;
									//vector <Dmr> dmrs;
									//getDmr(chr, allLocs, thisChr, dmrs );
								}

								consecDmc = 1; //new region starts
								allLocs.clear();
								dmcLocs.clear();
								startOfcombinedReg = start;
								endOfCombinedReg = start;
								allLocs.push_back(start);
								dmcLocs.push_back(start);
							}
						}
						lastStart = start;
						lastClass = loc.cls;
					} else {
						if(consecDmc >= 1)
						{
							allLocs.push_back(start);
						}
					}
					//finish this location
				}
			}
			//finish this chrom
			//may start call dmr now.
			//getDmr()
		}
	}
	dmcFile.close();
	dmcTrack.close();
	string chromsizes = option.reference.substr(0, option.reference.find_last_of(".")) + ".chrom.sizes"; 
	string dmcBigbed = "egrep -v '^track|^browser|^#' " + (dmcName + ".bed") + " > " + dmcName +".temp" + " && bedToBigBed " + dmcName + ".temp " + chromsizes + " " + dmcName + ".bb" + " && rm " + dmcName + ".temp";
	cout << "running: " << dmcBigbed << endl;
	system(dmcBigbed.c_str());
	cout << i << "...." << j <<endl;
	//done with DMC call

	int debug = 0;
	if(debug){
		//filter with coverage condition
		//this is actually testing . the coverage condition is already included in reading procedures
		map <int, map <string, map<int, vector<int> > > > goodLocsInReg ; //the cytosines inside dmcCombinedReg and satisfy the depth condition
		for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
			int laneId = id->first;
			map <string, map<int, vector<int> > > thisLane = id->second;
			for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
				string chr = pchr->first;
				map<int, vector<int> > thisChr = pchr->second;
				for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
					int start = pt->first;
					vector<int> allLocs = allLocsInReg[laneId][chr][start];
					for(vector<int>::iterator l = allLocs.begin(); l != allLocs.end(); l ++){
						int location = *l;
						if(pairResult[laneId][chr][location].ta < option.minDepthForComp || pairResult[laneId][chr][location].tb < option.minDepthForComp ){
							cerr << "ERRRRRRRRRRRR here" << endl;
							cout << chr << "\t" << start << "\t" << location << endl;
							//cout << pairResult[laneId][chr][location].ToString() << endl;
							exit(1);
						}
					}
					goodLocsInReg[laneId][chr][start] = allLocs;
				}
			}
		}
	}

	cout << "starting to call DMR" << endl;
	string dmrName = ("dmr_M2_" + label );
	ofstream dmrTable;
	dmrTable.open((dmrName + ".txt").c_str(), ios_base::out);
	dmrTable << "#chrom\tstart\tend"
			<< "\t" << "meanRatio_" << i << "\t" << "totalC_" << i << "\t" << "cSites_" << i
			<< "\t" << "meanRatio_" << j << "\t" << "totalC_" << j << "\t" << "cSites_" << j
			<< "\t" << "methDif_" << j << "-" << i
			<< "\t" << "p_" << j << "_v_" << i << "\t" << "class_" << j << "_v_" << i << endl;

	ofstream dmrTrack;
	dmrTrack.open((dmrName + ".bed").c_str(), ios_base::out);
	string trackHead = "track name=\"" + label + "_ScanDMR\" description=\"" + label +"_ScanDMR\" visibility=3 useScore=1 itemRgb=On";
	dmrTrack << trackHead << endl;
	//start to call DMR between lane i and j. //read only relevant info to save RAM.
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		//cout << laneId << endl;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			//cout << chr << endl;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> allLocs = allLocsInReg[laneId][chr][start];
				vector<int> dmcLocs = dmcLocsInReg[laneId][chr][start];
				//cout << "good" << endl;
				Dmr mostStrongDif;

				//start calling DMR for the region in laneId->chr->start
				int stop = 0;
				for(unsigned int dmcWindow = dmcLocs.size(); dmcWindow >= option.minDmcsInDmr && stop == 0; dmcWindow --){
					vector <Dmr> subReg;
					for(unsigned int dmcIndex = 0; dmcIndex + dmcWindow - 1 < dmcLocs.size(); dmcIndex ++){
						int windowLt = dmcLocs[dmcIndex];
						int windowRt = dmcLocs[dmcIndex + dmcWindow - 1];

						//cout << "this window" << endl;
						//get methdif and pvalue for this window from windowLt to windowRt
						double windowMethDif = 0;
						int ma = 0; //methylated cytosines in sample a
						int ua = 0; //unmethylated cytosines in sample b
						int mb = 0;
						int ub = 0;
						double mra = 0.0; //mean ratio in sample a
						double mrb = 0.0; //mean ratio in sample b
						int count = 0;
						for(unsigned int allIndex = 0; allIndex < allLocs.size(); allIndex ++){
							if(allLocs[allIndex] >= windowLt && allLocs[allIndex] <= windowRt){
								PairCompOut thisLoc = pairResult[laneId][chr][allLocs[allIndex]];
								windowMethDif += thisLoc.cDif;
								count += 1;
								mra += thisLoc.ra;
								mrb += thisLoc.rb;
								ma = ma + int (thisLoc.ta * thisLoc.ra + 0.5);
								ua = ua + int (thisLoc.ta * (1 - thisLoc.ra) + 0.5);
								mb = mb + int (thisLoc.tb * thisLoc.rb + 0.5);
								ub = ub + int (thisLoc.tb * (1 - thisLoc.rb) + 0.5);
							}
						}
						windowMethDif = windowMethDif / double( count );
						mra = mra / double( count );
						mrb = mrb / double( count );
						double p = 1.0;
						int varray[] = {ma, ua, mb, ub};
						//cout << "..." << ma << "\t" << ua << "\t" << mb << "\t" << ub << endl;
						fet_22( p, varray, 2, 2 );
						string cls = "insigOrLowDepth";
						//if( (p < option.pFetDmr) ){
						if(1){
							if(windowMethDif > 0){
								cls = windowMethDif > option.minCredibleDif ? "strongHyper" : "hyper";
							} else {
								cls = abs(windowMethDif) > option.minCredibleDif ? "strongHypo" : "hypo";
							}
						}
						//cout << "here " << p << endl;
						//Dmr(string Chr, int Start, int End, int AllSites, int DmcSites, int MC1, int UC1, int MC2, int UC2, double MethDif, double FetP)
						Dmr thisWindow(chr, windowLt, windowRt, count, dmcWindow, ma, ua, mb, ub, mra, mrb, windowMethDif, p, cls);
						subReg.push_back(thisWindow);
					}

					//determine the most strong region: #dmcWindow as 1st priority, then |#methDiff| >0.33
					sort(subReg.begin(), subReg.end(), sortDmrByMethDif);
					Dmr topReg = subReg[0];
					//if(abs(topReg.methDif) >= option.minNominalDif && topReg.fetP < option.pFetDmr){
					if(abs(topReg.methDif) >= option.minCredibleDif ){
						mostStrongDif = topReg;
						stop = 1;
					}
				}
				//dmrTable << "ss:" << start << endl;
				//done with call DMR for this candidate region
				if(stop == 1) {
//				//This extended_mostStrongDif is buggy by using tun0/comp.F1.MG.vs.F1.PG.txt
//
//					int start0 = mostStrongDif.start;
//					int end0 = mostStrongDif.end;
//					map<int, PairCompOut>::iterator ptA = pairResult[laneId][chr].find(start0);
//					map<int, PairCompOut>::iterator ptB = pairResult[laneId][chr].find(end0);
//					double methDif = mostStrongDif.methDif;
//					int allSites = mostStrongDif.allSites;
//
//					for(; ptA != pairResult[laneId][chr].begin(); )
//					{
//						ptA --;
//						int start = ptA->first;
//						PairCompOut loc = ptA->second;
//						if(chr == "chr17" && start > 12962600 && start< 12962799){cout << "::: " << start << " - " << log10(loc.simP) << " :: " << avgSim << endl;}
//						if( methDif*loc.cDif >=0 && abs( (methDif*allSites + loc.cDif) / (allSites+1) ) >= option.minCredibleDif && log10(loc.simP) < avgSim ){ //same direction of cDif, more sig than avgSim, minCDif condition satisfied
//							methDif = methDif*allSites + loc.cDif;
//							allSites = allSites + 1;
//							methDif = methDif / allSites;
//						}
//						else {
//							ptA ++;
//							break;
//						}
//					}
//
//					for(; ptB != pairResult[laneId][chr].end(); )
//					{
//						ptB ++;
//						int start = ptB->first;
//						PairCompOut loc = ptB->second;
//						if(chr == "chr17" && start > 12962600 && start< 12962799){cout << "::: " << start << " - " << log10(loc.simP) << " :: " << avgSim << endl;}
//						if( methDif*loc.cDif >=0 && abs( (methDif*allSites + loc.cDif) / (allSites+1) ) >= option.minCredibleDif && log10(loc.simP) < avgSim ){
//							methDif = methDif*allSites + loc.cDif;
//							allSites = allSites + 1;
//							methDif = methDif / allSites;
//						}
//						else {
//							//ptB --;
//							break;
//						}
//					}
//
//					int windowLt = ptA->first;
//					int windowRt = 0;
//
//					double windowMethDif = 0;
//					int ma = 0; //methylated cytosines in sample a
//					int ua = 0; //unmethylated cytosines in sample b
//					int mb = 0;
//					int ub = 0;
//					double mra = 0.0; //mean ratio in sample a
//					double mrb = 0.0; //mean ratio in sample b
//					int count = 0;
//
//					for(map<int, PairCompOut>::iterator pt = ptA; pt != ptB; pt ++){
//						PairCompOut thisLoc = pt->second;
//						windowRt = pt->first;
//						windowMethDif += thisLoc.cDif;
//						if(chr == "chr17" && start > 12962600 && start< 12962799){cout << ":::== " << pt->first << endl;}
//						count += 1;
//						mra += thisLoc.ra;
//						mrb += thisLoc.rb;
//						ma = ma + int (thisLoc.ta * thisLoc.ra + 0.5);
//						ua = ua + int (thisLoc.ta * (1 - thisLoc.ra) + 0.5);
//						mb = mb + int (thisLoc.tb * thisLoc.rb + 0.5);
//						ub = ub + int (thisLoc.tb * (1 - thisLoc.rb) + 0.5);
//					}
//
//					windowMethDif = windowMethDif / double( count );
//					mra = mra / double( count );
//					mrb = mrb / double( count );
//					double p = 1.0;
//					int varray[] = {ma, ua, mb, ub};
//					fet_1k( p, varray, 2, 2 );
//
//					string cls = "insigOrLowDepth";
//					if(1){
//						if(windowMethDif > 0){
//							cls = windowMethDif > option.minCredibleDif ? "strongHyper" : "hyper";
//						} else {
//							cls = abs(windowMethDif) > option.minCredibleDif ? "strongHypo" : "hypo";
//						}
//					}
//					Dmr extended_mostStrongDif(chr, windowLt, windowRt, count, mostStrongDif.dmcSites, ma, ua, mb, ub, mra, mrb, windowMethDif, p, cls);
//
//					dmrTable << extended_mostStrongDif.ToString() << endl;
//					dmrTrack << extended_mostStrongDif.ToBedTrackString() << endl;
//
					//comment out the following two lines if use the extended_mostStrongDif
					dmrTable << mostStrongDif.ToString() << endl;
					dmrTrack << mostStrongDif.ToBedTrackString() << endl;
				}
				//dmrTable << "ss:" << start << endl;
			}

		}
	}
	dmrTable.close();
	dmrTrack.close();
	chromsizes = option.reference.substr(0, option.reference.find_last_of(".")) + ".chrom.sizes"; 
	string dmrBigbed = "egrep -v '^track|^browser|^#' " + (dmrName + ".bed") + " > " + dmrName +".temp" + " && bedToBigBed " + dmrName + ".temp " + chromsizes + " " + dmrName + ".bb";
	cout << "running: " << dmrBigbed << endl;
	system(dmrBigbed.c_str());
	//combine DMCs first


	debug = 0;
	if(debug){
	cout << "size 1 " << dmcCombinedReg.size() << endl;
	cout << "size 2 " << allLocsInReg.size() << endl;
	ofstream yFile;
	yFile.open("tempxxxyyy", ios_base::out);
	for(map <int, map <string, map<int, int> > >::iterator id = dmcCombinedReg.begin(); id != dmcCombinedReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, int> > thisLane = id->second;
		for(map <string, map<int, int> >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, int> thisChr = pchr->second;
			for(map<int,int>::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				int end = pt->second;

				yFile << chr << "\t" << start << "\t" << end << endl;
			}
		}
	}
	yFile.close();
	ofstream zFile;
	zFile.open("tempxxxzzz", ios_base::out);
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> ends = pt->second;
				for(unsigned int i = 0; i < ends.size(); i ++){
					PairCompOut tt = pairResult[laneId][chr][ends[i]];
					zFile << tt.ToString() << endl;
				}

			}
		}
	}
	zFile.close();
	ofstream aFile;
	aFile.open("tempxxxaaa", ios_base::out);
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = dmcLocsInReg.begin(); id != dmcLocsInReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> ends = pt->second;
				for(unsigned int i = 0; i < ends.size(); i ++){
					PairCompOut tt = pairResult[laneId][chr][ends[i]];
					aFile << tt.ToString() << endl;
				}

			}
		}
	}
	aFile.close();
	}

}


void readCompToHashPred(string fileName, int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, statSim & genomeSim, Opts & option ){
	long double genomeSimP = 0.0;
	int cytosinesCountGenome = 0;

	stringstream fakeLaneId;
	fakeLaneId << "100" << i << "100" << j;
	int laneIdFake = string_to_int(fakeLaneId.str());
	int chrIndex = -1;
	int startIndex = -1;
	int endIndex = -1;
	int t1Index = -1;
	int r1Index = -1;
	int ci1Index = -1;
	int t2Index = -1;
	int r2Index = -1;
	int ci2Index = -1;
	int nDifIndex = -1;
	int cDifIndex = -1;
	int difCIIndex = -1;
	int simPIndex = -1;
	int fetPIndex = -1;

	string line = "";
	ifstream inputf(fileName.c_str(), ios::in);
	ofstream coveredCpgs;
	coveredCpgs.open((fileName + ".cvd.txt").c_str(), ios_base::out);
	while (inputf.good()) {
		getline(inputf, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));

		if( find(line.begin(), line.begin()+1, '#') == line.begin() )
		{
			chrIndex = find (fields.begin(), fields.end(), "#chrom") - fields.begin();
			startIndex = find (fields.begin(), fields.end(), "start") - fields.begin();
			endIndex = find (fields.begin(), fields.end(), "end") - fields.begin();
			t1Index = find (fields.begin(), fields.end(), "totalC_" + itos(i)) - fields.begin();
			r1Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(i)) - fields.begin();
			ci1Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(i)) - fields.begin();
			t2Index = find (fields.begin(), fields.end(), "totalC_" + itos(j)) - fields.begin();
			r2Index = find (fields.begin(), fields.end(), "nominalRatio_" + itos(j)) - fields.begin();
			ci2Index = find (fields.begin(), fields.end(), "ratioCI_" + itos(j)) - fields.begin();
			nDifIndex = find (fields.begin(), fields.end(), "nominalDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			cDifIndex = find (fields.begin(), fields.end(), "credibleDif_" + itos(j)+ "-" + itos(i)) - fields.begin();
			difCIIndex = find (fields.begin(), fields.end(), "difCI_" + itos(j)+ "-" + itos(i)) - fields.begin();
			simPIndex = find (fields.begin(), fields.end(), "p_sim_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			fetPIndex = find (fields.begin(), fields.end(), "p_fet_" + itos(j)+ "_v_" + itos(i)) - fields.begin();
			//cout << i << "\t" << j <<"\t"<< chrIndex << "\t" << startIndex << "\t" <<t1Index<< "\t" <<t2Index <<endl;
			//cout << i << "\t" << j <<"\t"<< r1Index << "\t" << r2Index << "\t" <<ci1Index<< "\t" <<ci2Index <<"\t" <<fetPIndex<<endl;
		}
		else{
			if(fields[fetPIndex] == "NA"){continue;} //if fetP field is "NA" then next;
			//cout << i << ".." << j << ".." << endl;
			string chr = fields[chrIndex];
			int start = string_to_int(fields[startIndex]);
			int end = string_to_int(fields[endIndex]);
			//string loc = chr + "\t" + fields[startIndex] + "\t" + fields[endIndex];
			int t1 = string_to_int(fields[t1Index]);
			//cout << loc << "\t" << t1 ;
			double r1 = atof((fields[r1Index]).c_str());
			//cout << "\t" << r1 ;
			string ci1 = fields[ci1Index];
			//cout << "\t" << ci1 <<endl;
			int t2 = string_to_int(fields[t2Index]);
			double r2 = atof((fields[r2Index]).c_str());
			string ci2 = fields[ci2Index];
			double simP = atof((fields[simPIndex]).c_str());
			//not a bug but to be fixed
			if(simP == 0){simP = 0.00000000000000001;}
			//

			double fetP = atof((fields[fetPIndex]).c_str());
			double nDif = atof((fields[nDifIndex]).c_str());
			double cDif = atof((fields[cDifIndex]).c_str());
			string difCI = fields[difCIIndex];

			//vector <string> ab;
			//boost::split(ab, difCI, boost::is_any_of(","));
			//double cDif = abs(atof(ab[0].c_str())) < abs(atof(ab[1].c_str())) ? atof(ab[0].c_str()) : atof(ab[1].c_str());
			string cls = "insigOrLowDepth";
			if( (fetP < option.pFetDmc) && (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp) ){
				if(nDif > 0){
					cls = nDif > option.minNominalDif ? "strongHyper" : "hyper";
				} else {
					cls = abs(nDif) > option.minNominalDif ? "strongHypo" : "hypo";
				}
			}
			if( (t1 >= option.minDepthForComp) && (t2 >= option.minDepthForComp) && ( abs(cDif) > option.filterCredibleDif ) ){
				pairResult[laneIdFake][chr][start] = PairCompOut(chr, start, end, t1, r1, ci1, t2, r2, ci2, simP, fetP, nDif, difCI, cDif, cls);
				//coveredCpgs << pairResult[laneIdFake][chr][start].ToString() << endl;
				genomeSimP += log10(simP);
				//cout << simP << "\t" << genomeSimP << endl;
				cytosinesCountGenome += 1;
			}
		}
	}
	genomeSim.cytosinesCountGenome = cytosinesCountGenome;
	genomeSim.genomeSimP = genomeSimP;
	cout << "done" << i << ".." <<j <<endl;
	inputf.close();
	coveredCpgs.close();
}

void isPredefinedRegionDmr(string regFile, string fileName, int i, int j, map <int, map <string, map<int, PairCompOut> > > & pairResult, statSim & genomeSim, Opts & option){
	//map <int, map <string, map<int, PairCompOut> > > pairResult ;
	//end of copy
	long double genomeSimP = genomeSim.genomeSimP;
	int cytosinesCountGenome = genomeSim.cytosinesCountGenome;

	stringstream fakeLaneId;
	fakeLaneId << "100" << i << "100" << j;
	int laneIdFake = string_to_int(fakeLaneId.str());

	string label = option.labels[i] + "_vs_" + option.labels[j];
	long double avgSim = genomeSimP / cytosinesCountGenome;
	cout << "total sim=" << genomeSimP << "\t" << pow((long double)10.0, genomeSimP) << endl;
	cout << "size=" << cytosinesCountGenome << " avg sim=" << avgSim << endl;

	map <int, map <string, map<int, int> > > predefinedReg; //fakeLaneId->chr->start ==> end
																//note that the end is only the position of C on + strand
	map <int, map <string, map<int, vector<int> > > > allLocsInReg ; //the cytosines inside predefinedReg; //fakeLaneId->chr->regStart ==> vector of sites

	map <int, map <string, map<int, string> > > predefinedRegAnno;

	//string interFile = regFile + ".Inter.AllCvdCpgs.bed";
	string interFile = regFile + ".Inter.AllCpgs.bed";
	//string regInterAllCpgsCmd = "intersectBed -a " + regFile + " -b " + fileName + ".cvd.txt" + " -wa -wb > " + interFile;
	string regInterAllCpgsCmd = "intersectBed -a " + regFile + " -b " + fileName + " -wa -wb > " + interFile;
	cout << "running: " << regInterAllCpgsCmd << endl;
	system(regInterAllCpgsCmd.c_str());
	//format: chr	start 	end	chr start 	end
	//format: chr	start 	end	name	score	strand chr start 	end

	string line = "";

	int numCols = 0;
	ifstream regFileStream(regFile.c_str(), ios_base::in);
	while(regFileStream.good()){
		getline(regFileStream, line);
		if( (line == "") || (line == "\n") || line.find("track")!=string::npos){continue;}
		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		numCols = fields.size();
		break;
	}
	cout << "numCols = " << numCols << endl;
	//exit(1);

	cout << "start reading "<< endl;
	int nc=0;
	ifstream regInterAllCpgs(interFile.c_str(), ios_base::in);
	while(regInterAllCpgs.good()){
		getline(regInterAllCpgs, line);
		if( (line == "") || (line == "\n") ){continue;}
		nc++;
		//cout << "processing line " << nc << endl;
		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		string chr = fields[0];
		int regStart = string_to_int(fields[1]);
		int regEnd = string_to_int(fields[2]);
		int cpgStart = string_to_int(fields[numCols+1]);

		predefinedReg[laneIdFake][chr][regStart] = regEnd;
		if(numCols > 3) {
			stringstream ss;
			copy(fields.begin()+3, fields.begin()+numCols, ostream_iterator<string>(ss, "\t"));
			predefinedRegAnno[laneIdFake][chr][regStart] = ss.str();
		}

		if( allLocsInReg.count(laneIdFake) != 0 && allLocsInReg[laneIdFake].count(chr) != 0 && allLocsInReg[laneIdFake][chr].count(regStart) != 0 ){
			vector<int> allLocs = allLocsInReg[laneIdFake][chr][regStart];
			allLocs.push_back(cpgStart);
			allLocsInReg[laneIdFake][chr][regStart] = allLocs;
		}else{
			vector <int> allLocs;
			allLocs.push_back(cpgStart);
			allLocsInReg[laneIdFake][chr][regStart] = allLocs;
		}

	}
	regInterAllCpgs.close();

	/*
	//cout << "start filtering"<< endl;
	//filter with coverage condition
	//this is actually testing . the coverage condition is already included in reading procedures
	map <int, map <string, map<int, vector<int> > > > goodLocsInReg ; //the cytosines inside dmcCombinedReg and satisfy the depth condition
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				vector<int> allLocs = allLocsInReg[laneId][chr][start];
				vector<int> goodLocs ;
				//cout << "processing "<< start << endl;

				for(vector<int>::iterator l = allLocs.begin(); l != allLocs.end(); l ++){
					int location = *l;
					if(pairResult[laneId][chr][location].ta < option.minDepthForComp || pairResult[laneId][chr][location].tb < option.minDepthForComp ){
						//cout << chr << "\t" << start << "\t" << location << "\t" << pairResult[laneId][chr][location].ta << "\t" <<  pairResult[laneId][chr][location].tb << endl;
						//allLocs.erase(l);
						cerr << "ERRRRRRRRRRRR here" << endl;
						//cout << chr << "\t" << start << "\t" << location << endl;
						//cout << pairResult[laneId][chr][location].ToString() << endl;
						exit(1);
					} else {
						goodLocs.push_back(location);
					}
				}
				allLocsInReg[laneId][chr][start] = goodLocs;
				goodLocsInReg[laneId][chr][start] = allLocs;
			}
		}
	}
	 */

	//stats related
	long long nCpgsInFeature = 0;
	double meanS1RatioOfCpgsInFeature = 0;
	double meanRatioOfRegionsInFeature1 = 0;
	double meanS2RatioOfCpgsInFeature = 0;
	double meanRatioOfRegionsInFeature2 = 0;
	double meanCDifOfCpgsInFeature = 0;
	int nRegionsInFeature = 0;
	int nHypoDmcsInFeatureM1 = 0;
	int nHyperDmcsInFeatureM1 = 0;
	int nHypoDmcsInFeatureM2 = 0;
	int nHyperDmcsInFeatureM2 = 0;

	cout << "starting to call DMR" << endl;
	string dmrName = ("dmr_" + regFile + "_" + label );
	ofstream dmrTable;
	dmrTable.open((dmrName + ".txt").c_str(), ios_base::out);

	dmrTable << "#chrom\tstart\tend";
	if(numCols > 3) {
		dmrTable << "\t" << "name\tscore\tstrand";
	}
	dmrTable
		<< "\t" << "meanRatio_" << i << "\t" << "totalC_" << i << "\t" << "cSites_" << i
		<< "\t" << "meanRatio_" << j << "\t" << "totalC_" << j << "\t" << "cSites_" << j
		<< "\t" << "meanOfNominalDif" << j << "-" << i
		<< "\t" << "p_fet_" << j << "_v_" << i << "\t" << "class_" << j << "_v_" << i
		<< "\t" << " sim, gm " << "\t"<< "sim" << "\t" << "gm"
		<< "\t" << " relativeSim, gm " << "\t" << "relativeSim" << "\t" << "relativeGm"
		<< "\t" << "meanOfCredibleDif" << "\t" << "meanOfAbsCredibleDif";
	if(option.minNominalDif >= 0){
		dmrTable << "\tminNominalDif\t" << "nHypoDmc" << "\t" << "nHyperDmc" ;
	}
	if(option.minCredibleDif >= 0){
		dmrTable << "\tminCredibleDif\t" << "nHypoDmc" << "\t" << "nHyperDmc" ;
	}
	dmrTable << endl;

	ofstream dmrTrack;
	dmrTrack.open((dmrName + ".bed").c_str(), ios_base::out);
	string trackHead = "track name=\"" + label + "_ScanDMR\" description=\"" + label +"_ScanDMR\" visibility=3 useScore=1 itemRgb=On";
	dmrTrack << trackHead << endl;
	//start to call DMR between lane i and j. //read only relevant info to save RAM.
	for(map <int, map <string, map<int, vector<int> > > >::iterator id = allLocsInReg.begin(); id != allLocsInReg.end(); id ++){
		int laneId = id->first;
		//cout << laneId << endl;
		map <string, map<int, vector<int> > > thisLane = id->second;
		for(map <string, map<int, vector<int> > >::iterator pchr = thisLane.begin(); pchr != thisLane.end(); pchr ++){
			string chr = pchr->first;
			//cout << chr << endl;
			map<int, vector<int> > thisChr = pchr->second;
			for(map<int,vector<int> >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
				int start = pt->first;
				int windowLt = start;								//regStart
				int windowRt = predefinedReg[laneId][chr][start];	//regEnd

				vector<int> allLocs = allLocsInReg[laneId][chr][start];
				//cout << "good" << endl;
				//cout << start <<"\t";
				//copy(allLocs.begin(), allLocs.end(), ostream_iterator<int>(cout, ","));cout << endl;
				//cout << allLocs.size()<< endl;

						double windowMethDif = 0;
						int ma = 0; //methylated cytosines in sample a
						int ua = 0; //unmethylated cytosines in sample b
						int mb = 0;
						int ub = 0;
						double mra = 0.0; //mean ratio in sample a
						double mrb = 0.0; //mean ratio in sample b
						int count = 0;
						int nHypoDmc = 0;
						int nHyperDmc = 0;
						int nHypoDmcM1 = 0;
						int nHyperDmcM1 = 0;
						
						nRegionsInFeature ++;
                                                
						vector <long double> simPVector;
						long double sumOfLogSimP = 0.0;
						long double sumOfCredibleMethDif = 0.0;
						long double sumOfAbsCredibleMethDif = 0.0;
						for(unsigned int allIndex = 0; allIndex < allLocs.size(); allIndex ++){
							//cout << allIndex << "\t" << windowLt << "\t" << windowRt << endl;
							if(allLocs[allIndex] >= windowLt - 1 && allLocs[allIndex] <= windowRt){ // this statment is not necessary for defined regions// -1 is to include the situation of a cpg is across the boundary of region.
								if(pairResult[laneId][chr].count(allLocs[allIndex]) == 0)continue; //this loc may not satisfy the filtering conditions.
								PairCompOut thisLoc = pairResult[laneId][chr][allLocs[allIndex]];
								windowMethDif += thisLoc.nDif;
								count += 1;
								nCpgsInFeature ++;
                                                                
								//cout << allIndex << "\t" << windowLt << "\t" << windowRt << endl;
								if(abs(thisLoc.cDif) > option.minCredibleDif){
									if(thisLoc.cDif>0) nHyperDmc ++;
									if(thisLoc.cDif<0) nHypoDmc ++;
								}
								if(abs(thisLoc.nDif) > option.minNominalDif){
									if(thisLoc.nDif>0) nHyperDmcM1 ++;
									if(thisLoc.nDif<0) nHypoDmcM1 ++;
								}								
								mra += thisLoc.ra;
								mrb += thisLoc.rb;

								meanS1RatioOfCpgsInFeature += thisLoc.ra;
								meanS2RatioOfCpgsInFeature += thisLoc.rb;

								ma = ma + int (thisLoc.ta * thisLoc.ra + 0.5);
								ua = ua + int (thisLoc.ta * (1 - thisLoc.ra) + 0.5);
								mb = mb + int (thisLoc.tb * thisLoc.rb + 0.5);
								ub = ub + int (thisLoc.tb * (1 - thisLoc.rb) + 0.5);
								simPVector.push_back(thisLoc.simP);
								sumOfLogSimP += log10(thisLoc.simP);
								sumOfCredibleMethDif += thisLoc.cDif;
								sumOfAbsCredibleMethDif += abs(thisLoc.cDif);
							}
						}
						if(count ==0)continue;
						windowMethDif = windowMethDif / double( count );
						mra = mra / double( count );
						mrb = mrb / double( count );
                        meanRatioOfRegionsInFeature1 += mra;
                        meanRatioOfRegionsInFeature2 += mrb;
                        meanCDifOfCpgsInFeature += sumOfCredibleMethDif;
						double p = 1.0;
						int varray[] = {ma, ua, mb, ub};
						//cout << "..." << ma << "\t" << ua << "\t" << mb << "\t" << ub << endl;
						fet_22( p, varray, 2, 2 );

						//long double geometricMeanSimP = geometricMean(simPVector);
						//cout << sumOfLogSimP << endl;
						long double relativeSim = sumOfLogSimP - count * avgSim;
						string cls = "insigOrLowDepth";
						if( (p < option.pFetDmr) ){
							if(windowMethDif > 0){
								cls = windowMethDif > option.minNominalDif ? "strongHyper" : "hyper";
							} else {
								cls = abs(windowMethDif) > option.minNominalDif ? "strongHypo" : "hypo";
							}
						}
						//cout << "here " << p << endl;
						//Dmr(string Chr, int Start, int End, int AllSites, int DmcSites, int MC1, int UC1, int MC2, int UC2, double MethDif, double FetP)
						int dmcWindow = -1;
						Dmr thisWindow(chr, windowLt, windowRt, count, dmcWindow, ma, ua, mb, ub, mra, mrb, windowMethDif,  p, cls);
						//dmrTable << thisWindow.ToString() << endl;
						stringstream ss;
						ss 	<< setprecision(3) << chr << "\t" << windowLt << "\t" << windowRt << "\t";
						if(numCols > 3){
							ss << predefinedRegAnno[laneIdFake][chr][windowLt] ;
						}
						ss	<< mra << "\t" << ma+ua << "\t" << count << "\t"
							<< mrb << "\t" << mb+ub << "\t" << count << "\t"
							<< windowMethDif << "\t" << p << "\t" << cls << "\t"
							<< " sim, gm " << "\t"<< sumOfLogSimP << "\t" << sumOfLogSimP / count << "\t"
							<< "relativeSim, gm\t" << relativeSim << "\t" << relativeSim / count << "\t"
							<< sumOfCredibleMethDif/count << "\t" << sumOfAbsCredibleMethDif/count;
						if(option.minNominalDif >= 0){
							ss << "\t" << "minNominalDif\t" << nHypoDmcM1 << "\t" << nHyperDmcM1;
						}
						if(option.minCredibleDif >= 0){
							ss << "\t" << "minCredibleDif\t" << nHypoDmc << "\t" << nHyperDmc;
						}
						dmrTable << ss.str() << endl;
						dmrTrack << thisWindow.ToBedTrackString() << endl;

						nHyperDmcsInFeatureM1 += nHyperDmcM1;
						nHyperDmcsInFeatureM2 += nHyperDmc;
						nHypoDmcsInFeatureM1 += nHypoDmcM1;
						nHypoDmcsInFeatureM2 += nHypoDmc;
			}
		}

	}
	dmrTable.close();
	dmrTrack.close();

	long long nDmcsInFeatureM1 = (nHypoDmcsInFeatureM1 + nHyperDmcsInFeatureM1);
    long long nDmcsInFeatureM2 = (nHypoDmcsInFeatureM2 + nHyperDmcsInFeatureM2);
    meanS1RatioOfCpgsInFeature = meanS1RatioOfCpgsInFeature / nCpgsInFeature;
    meanS2RatioOfCpgsInFeature = meanS2RatioOfCpgsInFeature / nCpgsInFeature;
    meanRatioOfRegionsInFeature1 = meanRatioOfRegionsInFeature1 / nRegionsInFeature;
    meanRatioOfRegionsInFeature2 = meanRatioOfRegionsInFeature2 / nRegionsInFeature;


    ofstream dmrStat;
    dmrStat.open((dmrName + ".stat").c_str(), ios_base::out);

    dmrStat    << "nCpgsInFeature" <<"\t"<< "meanS1RatioOfCpgsInFeature" <<"\t"<< "meanS2RatioOfCpgsInFeature" <<"\t"<< "meanRatioOfCpgsInS2-meanRatioOfCpgsInS1"
	<<"\t"	<< "meanCDifOfCpgsInFeature"
	<<"\t"  << "nDmcsInFeatureM1" <<"\t"<< "nHypoDmcsInFeatureM1" <<"\t"<< "nHyperDmcsInFeatureM1"
	<<"\t"  << "nDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nHyperDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM1/nDmcsInFeatureM1"
	<<"\t"  << "nDmcsInFeatureM2" <<"\t"<< "nHypoDmcsInFeatureM2" <<"\t"<< "nHyperDmcsInFeatureM2"
	<<"\t"  << "nDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nHyperDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM2/nDmcsInFeatureM2"
	<<"\t"  << "nRegionsInFeature" <<"\t"<< "meanS1RatioOfRegionsInFeature" <<"\t"<< "meanS2RatioOfRegionsInFeature" <<"\t"<< "meanRatioOfRegionsInS2 - meanRatioOfRegionsInS1" <<endl;

    dmrStat    << nCpgsInFeature <<"\t"<< meanS1RatioOfCpgsInFeature <<"\t"<< meanS2RatioOfCpgsInFeature <<"\t"<< meanS2RatioOfCpgsInFeature - meanS1RatioOfCpgsInFeature
    <<"\t"	<< meanCDifOfCpgsInFeature / nCpgsInFeature
	<<"\t"  << nDmcsInFeatureM1 <<"\t"<< nHypoDmcsInFeatureM1 <<"\t"<< nHyperDmcsInFeatureM1
    <<"\t"  << nDmcsInFeatureM1/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM1)/double(nCpgsInFeature) <<"\t"<< (nHyperDmcsInFeatureM1)/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM1)/double (nDmcsInFeatureM1)
    <<"\t"  << nDmcsInFeatureM2 <<"\t"<< nHypoDmcsInFeatureM2 <<"\t"<< nHyperDmcsInFeatureM2
    <<"\t"  << nDmcsInFeatureM2/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM2)/double(nCpgsInFeature) <<"\t"<< (nHyperDmcsInFeatureM2)/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM2)/double (nDmcsInFeatureM2)
    <<"\t"  << nRegionsInFeature <<"\t"<< meanRatioOfRegionsInFeature1 <<"\t"<< meanRatioOfRegionsInFeature2 <<"\t"<< meanRatioOfRegionsInFeature2 - meanRatioOfRegionsInFeature1 <<endl;
    dmrStat.close();

    string chromsizes = option.reference.substr(0, option.reference.find_last_of(".")) + ".chrom.sizes"; 
	string dmrBigbed = "egrep -v '^track|^browser|^#' " + (dmrName + ".bed") + " > " + dmrName +".temp" + " && bedToBigBed " + dmrName + ".temp " + chromsizes + " " + dmrName + ".bb";
	cout << "running: " << dmrBigbed << endl;
	system(dmrBigbed.c_str());
	//combine DMCs first


	//string interFile = regFile + ".Inter.AllCvdCpgs.bed";
	string regAllCpgsFile = regFile + ".AllCpgs.bed";
	//string regInterAllCpgsCmd = "intersectBed -a " + regFile + " -b " + fileName + ".cvd.txt" + " -wa -wb > " + interFile;
	string regAllCpgsCmd = "intersectBed -a " + fileName + " -b " + regFile + " -wa | uniq > " + regAllCpgsFile;
	cout << "running: " << regAllCpgsCmd << endl;
	system(regAllCpgsCmd.c_str());

	//stats related
	long long nCpgsInAll = 0;
	double meanS1RatioOfCpgsInAll = 0;
	double meanS2RatioOfCpgsInAll = 0;
	double meanCDifOfCpgsInAll = 0;

	int nHypoDmcsInAllM1 = 0;
	int nHyperDmcsInAllM1 = 0;
	int nHypoDmcsInAllM2 = 0;
	int nHyperDmcsInAllM2 = 0;

	double meanS1RatioOfHyperDmcsM1InAll = 0;
	double meanS2RatioOfHyperDmcsM1InAll = 0;
	double meanS1RatioOfHypoDmcsM1InAll = 0;
	double meanS2RatioOfHypoDmcsM1InAll = 0;

	double meanS1RatioOfHyperDmcsM2InAll = 0;
	double meanS2RatioOfHyperDmcsM2InAll = 0;
	double meanS1RatioOfHypoDmcsM2InAll = 0;
	double meanS2RatioOfHypoDmcsM2InAll = 0;

	nCpgsInFeature = 0;
	meanS1RatioOfCpgsInFeature = 0;
	meanS2RatioOfCpgsInFeature = 0;
	meanCDifOfCpgsInFeature = 0;

	nHypoDmcsInFeatureM1 = 0;
	nHyperDmcsInFeatureM1 = 0;
	nHypoDmcsInFeatureM2 = 0;
	nHyperDmcsInFeatureM2 = 0;

	double meanS1RatioOfHyperDmcsM1InFeature = 0;
	double meanS2RatioOfHyperDmcsM1InFeature = 0;
	double meanS1RatioOfHypoDmcsM1InFeature = 0;
	double meanS2RatioOfHypoDmcsM1InFeature = 0;

	double meanS1RatioOfHyperDmcsM2InFeature = 0;
	double meanS2RatioOfHyperDmcsM2InFeature = 0;
	double meanS1RatioOfHypoDmcsM2InFeature = 0;
	double meanS2RatioOfHypoDmcsM2InFeature = 0;

	map <int, map <string, map<int, int> > > allLocsInFeature; //fakeLaneId->chr->start ==> 1
	line = "";
	cout << "start reading "<< endl;
	ifstream regAllCpgs(regAllCpgsFile.c_str(), ios_base::in);
	while(regAllCpgs.good()){
		getline(regAllCpgs, line);
		if( (line == "") || (line == "\n") ){continue;}

		vector <string> fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		string chr = fields[0];
		int start = string_to_int(fields[1]);

		allLocsInFeature[laneIdFake][chr][start] = 1;
	}
	regAllCpgs.close();


	for(map<string, map<int, PairCompOut> >::iterator pchr = pairResult[laneIdFake].begin(); pchr != pairResult[laneIdFake].end(); pchr ++){
		string chr = pchr->first;
		map<int, PairCompOut> thisChr = pchr->second;
		for(map<int, PairCompOut>::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
			int start = pt->first;
			if(pairResult[laneIdFake][chr].count(start) == 0)continue; //this loc may not satisfy the filtering conditions.
			nCpgsInAll ++;
			PairCompOut thisLoc = pt->second;//pairResult[laneIdFake][chr][start];

			if(abs(thisLoc.cDif) > option.minCredibleDif){
				if(thisLoc.cDif>0) {
					nHyperDmcsInAllM2 ++;
					meanS1RatioOfHyperDmcsM2InAll += thisLoc.ra;
					meanS2RatioOfHyperDmcsM2InAll += thisLoc.rb;
				}
				if(thisLoc.cDif<0) {
					nHypoDmcsInAllM2 ++;
					meanS1RatioOfHypoDmcsM2InAll += thisLoc.ra;
					meanS2RatioOfHypoDmcsM2InAll += thisLoc.rb;
				}
			}
			if(abs(thisLoc.nDif) > option.minNominalDif){
				if(thisLoc.nDif>0) {
					nHyperDmcsInAllM1 ++;
					meanS1RatioOfHyperDmcsM1InAll += thisLoc.ra;
					meanS2RatioOfHyperDmcsM1InAll += thisLoc.rb;
				}
				if(thisLoc.nDif<0) {
					nHypoDmcsInAllM1 ++;
					meanS1RatioOfHypoDmcsM1InAll += thisLoc.ra;
					meanS2RatioOfHypoDmcsM1InAll += thisLoc.rb;
				}
			}

			meanS1RatioOfCpgsInAll += thisLoc.ra;
			meanS2RatioOfCpgsInAll += thisLoc.rb;
			meanCDifOfCpgsInAll += thisLoc.cDif;
		}
	}


	for(map <string, map<int, int > >::iterator pchr = allLocsInFeature[laneIdFake].begin(); pchr != allLocsInFeature[laneIdFake].end(); pchr ++){
		string chr = pchr->first;
		map<int, int > thisChr = pchr->second;
		for(map<int,int >::iterator pt = thisChr.begin(); pt != thisChr.end(); pt ++){
			int start = pt->first;
			if(pairResult[laneIdFake][chr].count(start) == 0)continue; //this loc may not satisfy the filtering conditions.
			nCpgsInFeature ++;
			PairCompOut thisLoc = pairResult[laneIdFake][chr][start];

			if(abs(thisLoc.cDif) > option.minCredibleDif){
				if(thisLoc.cDif>0) {
					nHyperDmcsInFeatureM2 ++;
					meanS1RatioOfHyperDmcsM2InFeature += thisLoc.ra;
					meanS2RatioOfHyperDmcsM2InFeature += thisLoc.rb;
				}
				if(thisLoc.cDif<0) {
					nHypoDmcsInFeatureM2 ++;
					meanS1RatioOfHypoDmcsM2InFeature += thisLoc.ra;
					meanS2RatioOfHypoDmcsM2InFeature += thisLoc.rb;
				}
			}
			if(abs(thisLoc.nDif) > option.minNominalDif){
				if(thisLoc.nDif>0) {
					nHyperDmcsInFeatureM1 ++;
					meanS1RatioOfHyperDmcsM1InFeature += thisLoc.ra;
					meanS2RatioOfHyperDmcsM1InFeature += thisLoc.rb;
				}
				if(thisLoc.nDif<0) {
					nHypoDmcsInFeatureM1 ++;
					meanS1RatioOfHypoDmcsM1InFeature += thisLoc.ra;
					meanS2RatioOfHypoDmcsM1InFeature += thisLoc.rb;
				}
			}

			meanS1RatioOfCpgsInFeature += thisLoc.ra;
			meanS2RatioOfCpgsInFeature += thisLoc.rb;
			meanCDifOfCpgsInFeature += thisLoc.cDif;
		}
	}
	long long nDmcsInAllM1 = (nHypoDmcsInAllM1 + nHyperDmcsInAllM1);
    long long nDmcsInAllM2 = (nHypoDmcsInAllM2 + nHyperDmcsInAllM2);
    meanS1RatioOfCpgsInAll = meanS1RatioOfCpgsInAll / nCpgsInAll;
    meanS2RatioOfCpgsInAll = meanS2RatioOfCpgsInAll / nCpgsInAll;
    meanCDifOfCpgsInAll = meanCDifOfCpgsInAll / nCpgsInAll;

	meanS1RatioOfHyperDmcsM1InAll /= nHyperDmcsInAllM1;
	meanS2RatioOfHyperDmcsM1InAll /= nHyperDmcsInAllM1;
	meanS1RatioOfHypoDmcsM1InAll /= nHypoDmcsInAllM1;
	meanS2RatioOfHypoDmcsM1InAll /= nHypoDmcsInAllM1;

	meanS1RatioOfHyperDmcsM2InAll /= nHyperDmcsInAllM2;
	meanS2RatioOfHyperDmcsM2InAll /= nHyperDmcsInAllM2;
	meanS1RatioOfHypoDmcsM2InAll /= nHyperDmcsInAllM2;
	meanS2RatioOfHypoDmcsM2InAll /= nHyperDmcsInAllM2;

	nDmcsInFeatureM1 = (nHypoDmcsInFeatureM1 + nHyperDmcsInFeatureM1);
    nDmcsInFeatureM2 = (nHypoDmcsInFeatureM2 + nHyperDmcsInFeatureM2);

    meanS1RatioOfCpgsInFeature = meanS1RatioOfCpgsInFeature / nCpgsInFeature;
    meanS2RatioOfCpgsInFeature = meanS2RatioOfCpgsInFeature / nCpgsInFeature;
    meanCDifOfCpgsInFeature = meanCDifOfCpgsInFeature / nCpgsInFeature;

	meanS1RatioOfHyperDmcsM1InFeature /= nHyperDmcsInFeatureM1;
	meanS2RatioOfHyperDmcsM1InFeature /= nHyperDmcsInFeatureM1;
	meanS1RatioOfHypoDmcsM1InFeature /= nHypoDmcsInFeatureM1;
	meanS2RatioOfHypoDmcsM1InFeature /= nHypoDmcsInFeatureM1;

	meanS1RatioOfHyperDmcsM2InFeature /= nHyperDmcsInFeatureM2;
	meanS2RatioOfHyperDmcsM2InFeature /= nHyperDmcsInFeatureM2;
	meanS1RatioOfHypoDmcsM2InFeature /= nHyperDmcsInFeatureM2;
	meanS2RatioOfHypoDmcsM2InFeature /= nHyperDmcsInFeatureM2;

	string dmcName = ("dmc_" + regFile + "_" + label );
	ofstream dmcStat;
	dmcStat.open((dmcName + ".stat").c_str(), ios_base::out);

	dmcStat    << "Feature" <<"\t"<< "nCpgsInFeature" <<"\t"<< "nCpgsInFeature/nCpgsInAll" <<"\t"<< "meanS1RatioOfCpgsInFeature" <<"\t"<< "meanS2RatioOfCpgsInFeature"
    <<"\t"	<< "meanNDifOfCpgsInFeature" <<"\t"<< "meanCDifOfCpgsInFeature"

    <<"\t"  << "nDmcsInFeatureM1" <<"\t"<< "nDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nDmcsInFeatureM1/nDmcsInAllM1"
    <<"\t"	<< "nHypoDmcsInFeatureM1" <<"\t"<< "nHyperDmcsInFeatureM1"
	<<"\t"	<< "nHypoDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nHyperDmcsInFeatureM1/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM1/nDmcsInFeatureM1"
	<<"\t"	<< "meanS1RatioOfHypoDmcsM1InFeature" <<"\t"<< "meanS2RatioOfHypoDmcsM1InFeature" <<"\t"<< "meanS1RatioOfHyperDmcsM1InFeature" <<"\t"<< "meanS2RatioOfHyperDmcsM1InFeature"

	<<"\t"  << "nDmcsInFeatureM2" <<"\t"<< "nDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nDmcsInFeatureM2/nDmcsInAllM2"
	<<"\t"	<< "nHypoDmcsInFeatureM2" <<"\t"<< "nHyperDmcsInFeatureM2"
	<<"\t"  << "nHypoDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nHyperDmcsInFeatureM2/nCpgsInFeature" <<"\t"<< "nHypoDmcsInFeatureM2/nDmcsInFeatureM2"
	<<"\t"	<< "meanS1RatioOfHypoDmcsM2InFeature" <<"\t"<< "meanS2RatioOfHypoDmcsM2InFeature" <<"\t"<< "meanS1RatioOfHyperDmcsM2InFeature" <<"\t"<< "meanS2RatioOfHyperDmcsM2InFeature"
	<<"\t"  << "nRegionsInFeature" <<"\t"<< "meanS1RatioOfRegionsInFeature" <<"\t"<< "meanS2RatioOfRegionsInFeature" <<"\t"<< "meanS2RatioOfRegionsInFeature-meanS1RatioOfRegionsInFeature"
	<< endl;

	dmcStat    << "All" <<"\t"<< nCpgsInAll <<"\t"<< nCpgsInAll/double(nCpgsInAll) <<"\t"<< meanS1RatioOfCpgsInAll <<"\t"<< meanS2RatioOfCpgsInAll
	<<"\t"	<< meanS2RatioOfCpgsInAll - meanS1RatioOfCpgsInAll <<"\t"<< meanCDifOfCpgsInAll

	<<"\t"  << nDmcsInAllM1 <<"\t"<< nDmcsInAllM1/double(nCpgsInAll) <<"\t"<< nDmcsInAllM1/double(nDmcsInAllM1)
	<<"\t"	<< nHypoDmcsInAllM1 <<"\t"<< nHyperDmcsInAllM1
    <<"\t"	<< (nHypoDmcsInAllM1)/double(nCpgsInAll) <<"\t"<< (nHyperDmcsInAllM1)/double(nCpgsInAll) <<"\t"<< (nHypoDmcsInAllM1)/double (nDmcsInAllM1)
    <<"\t"	<< meanS1RatioOfHypoDmcsM1InAll <<"\t"<< meanS2RatioOfHypoDmcsM1InAll <<"\t"<< meanS1RatioOfHyperDmcsM1InAll <<"\t"<< meanS2RatioOfHyperDmcsM1InAll

    <<"\t"  << nDmcsInAllM2 <<"\t"<< nDmcsInAllM2/double(nCpgsInAll) <<"\t"<< nDmcsInAllM2/double(nDmcsInAllM2)
    <<"\t"	<< nHypoDmcsInAllM2 <<"\t"<< nHyperDmcsInAllM2
    <<"\t"	<< (nHypoDmcsInAllM2)/double(nCpgsInAll) <<"\t"<< (nHyperDmcsInAllM2)/double(nCpgsInAll) <<"\t"<< (nHypoDmcsInAllM2)/double (nDmcsInAllM2)
    <<"\t"	<< meanS1RatioOfHypoDmcsM2InAll <<"\t"<< meanS2RatioOfHypoDmcsM2InAll <<"\t"<< meanS1RatioOfHyperDmcsM2InAll <<"\t"<< meanS2RatioOfHyperDmcsM2InAll
	<<"\t"  << "nan" <<"\t"<< "nan" <<"\t"<< "nan" <<"\t"<< "nan"
    << endl;

	dmcStat    << regFile <<"\t"<< nCpgsInFeature <<"\t"<< nCpgsInFeature/double(nCpgsInAll) <<"\t"<< meanS1RatioOfCpgsInFeature <<"\t"<< meanS2RatioOfCpgsInFeature
	<<"\t"	<< meanS2RatioOfCpgsInFeature - meanS1RatioOfCpgsInFeature <<"\t"<< meanCDifOfCpgsInFeature

	<<"\t"  << nDmcsInFeatureM1 <<"\t"<< nDmcsInFeatureM1/double(nCpgsInFeature) <<"\t"<< nDmcsInFeatureM1/double(nDmcsInAllM1)
	<<"\t"	<< nHypoDmcsInFeatureM1 <<"\t"<< nHyperDmcsInFeatureM1
    <<"\t"	<< (nHypoDmcsInFeatureM1)/double(nCpgsInFeature) <<"\t"<< (nHyperDmcsInFeatureM1)/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM1)/double (nDmcsInFeatureM1)
    <<"\t"	<< meanS1RatioOfHypoDmcsM1InFeature <<"\t"<< meanS2RatioOfHypoDmcsM1InFeature <<"\t"<< meanS1RatioOfHyperDmcsM1InFeature <<"\t"<< meanS2RatioOfHyperDmcsM1InFeature

    <<"\t"  << nDmcsInFeatureM2 <<"\t"<< nDmcsInFeatureM2/double(nCpgsInFeature) <<"\t"<< nDmcsInFeatureM2/double(nDmcsInAllM2)
    <<"\t"	<< nHypoDmcsInFeatureM2 <<"\t"<< nHyperDmcsInFeatureM2
    <<"\t"	<< (nHypoDmcsInFeatureM2)/double(nCpgsInFeature) <<"\t"<< (nHyperDmcsInFeatureM2)/double(nCpgsInFeature) <<"\t"<< (nHypoDmcsInFeatureM2)/double (nDmcsInFeatureM2)
    <<"\t"	<< meanS1RatioOfHypoDmcsM2InFeature <<"\t"<< meanS2RatioOfHypoDmcsM2InFeature <<"\t"<< meanS1RatioOfHyperDmcsM2InFeature <<"\t"<< meanS2RatioOfHyperDmcsM2InFeature
	<<"\t"  << nRegionsInFeature <<"\t"<< meanRatioOfRegionsInFeature1 <<"\t"<< meanRatioOfRegionsInFeature2 <<"\t"<< meanRatioOfRegionsInFeature2 - meanRatioOfRegionsInFeature1

    << endl;

	dmcStat.close();
}



int load_lut(int num_threads, string exep)
{

	//build table for pdiffInRegion
	if( boost::filesystem::exists( exep + "/lut_pdiffInRegion.dat" ) ){
		//read
		cout << "Reading " << exep << "/lut_pdiffInRegion.dat" << endl;

		int tableMax;
		FILE *FpLut=fopen((exep + "/lut_pdiffInRegion.dat").c_str(), "rb");

		double value;
		fread( &value, 1, sizeof(value), FpLut);
		tableMax = (int) value;

		for(int n1 = 1; n1 <= tableMax; n1 ++){
			for(int k1 = 0; k1 <= n1; k1 ++){
				for(int n2 = 1; n2 <= tableMax; n2 ++){
					for(int k2 = 0; k2 <= n2; k2 ++){
						MultiKey combi(n1,k1,n2,k2);
						fread( &value, 1, sizeof(value), FpLut);
						lut_pdiffInRegion[combi] = value;
						//cout << std::setprecision(20);
						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << lut_pdiffInRegion[combi] << endl;
					}
				}
			}
		}
		fclose(FpLut);

	} else {
		cout << "Building " << exep << "/lut_pdiffInRegion.dat" << endl;
		cout << "You should not see this message unless you are trying to build the database" << endl;
		cout << "Check if lut_pdiffInRegion.dat is at the same location with mcomp" << endl;

		int tableMax = 30;
		//build
		build_lut_pdiffInRegion(tableMax, num_threads);

		//write
		FILE *fpLut=fopen((exep + "/lut_pdiffInRegion.dat").c_str(), "wb");
		double value=(double)tableMax;
		fwrite(&value, 1, sizeof(value), fpLut); // table starts with tableMax;
		for(int n1 = 1; n1 <= tableMax; n1 ++){
			for(int k1 = 0; k1 <= n1; k1 ++){
				for(int n2 = 1; n2 <= tableMax; n2 ++){
					for(int k2 = 0; k2 <= n2; k2 ++){
						//value = pdiffInRegion(k1,n1,k2,n2,DIFFDELTA);
						MultiKey combi(n1,k1,n2,k2);
						value = lut_pdiffInRegion[combi];
						fwrite(&value, 1, sizeof(value), fpLut);
						//cout << std::setprecision(20);
						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << value << endl;
					}
				}
			}
		}
		fclose(fpLut);
	}

	//build table for pdiffCI
	if( boost::filesystem::exists( exep + "/lut_pdiffCI.dat" ) ){
		//read
		cout << "Reading " << exep << "/lut_pdiffCI.dat" << endl;
		int tableMax;
		FILE *FpLut=fopen((exep + "/lut_pdiffCI.dat").c_str(), "rb");

		double value;
		fread( &value, 1, sizeof(value), FpLut);
		tableMax = (int) value;

		CI ci; ci.a = -1; ci.b = -1;

		for(int n1 = 1; n1 <= tableMax; n1 ++){
			for(int k1 = 0; k1 <= n1; k1 ++){
				for(int n2 = 1; n2 <= tableMax; n2 ++){
					for(int k2 = 0; k2 <= n2; k2 ++){
						MultiKey combi(n1,k1,n2,k2);
						fread( &ci, 1, sizeof(ci), FpLut);
						lut_pdiffCI[combi] = ci;
						//cout << std::setprecision(20);
						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << lut_pdiffCI[combi].a << "\t" << lut_pdiffCI[combi].b << endl;
					}
				}
			}
		}
		fclose(FpLut);

	} else {
		cout << "Building " << exep << "/lut_pdiffCI.dat" << endl;
		int tableMax = 30;
		//build
		build_lut_pdiffCI(tableMax, num_threads);

		//write
		FILE *fpLut=fopen((exep + "/lut_pdiffCI.dat").c_str(), "wb");
		double value=(double)tableMax;
		fwrite(&value, 1, sizeof(value), fpLut); // table starts with tableMax;

		CI ci; ci.a = -1; ci.b = -1;

		for(int n1 = 1; n1 <= tableMax; n1 ++){
			for(int k1 = 0; k1 <= n1; k1 ++){
				for(int n2 = 1; n2 <= tableMax; n2 ++){
					for(int k2 = 0; k2 <= n2; k2 ++){
						//ci = pdiffCI(k1, n1, k2, n2, ALPHA, 1);
						MultiKey combi(n1,k1,n2,k2);
						ci = lut_pdiffCI[combi];
						fwrite(&ci, 1, sizeof(ci), fpLut);
						//cout << std::setprecision(20);
						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << ci.a << "\t" << ci.b << endl;
					}
				}
			}
		}
		fclose(fpLut);
	}

// Do not load `lut_fet` for now, since it is not used.

// 	//build table for fisher exact test
// 	if( boost::filesystem::exists( exep + "/lut_fet.dat" ) ){
// 		cout << "Reading " << exep << "/lut_fet.dat" << endl;
// 		//read
// 		int tableMax;
// 		FILE *FpLut=fopen((exep + "/lut_fet.dat").c_str(), "rb");
// 
// 		double value;
// 		fread( &value, 1, sizeof(value), FpLut);
// 		tableMax = (int) value;
// 
// 		for(int n1 = 1; n1 <= tableMax; n1 ++){
// 			for(int k1 = 0; k1 <= n1; k1 ++){
// 				for(int n2 = 1; n2 <= tableMax; n2 ++){
// 					for(int k2 = 0; k2 <= n2; k2 ++){
// 						MultiKey combi(n1,k1,n2,k2);
// 						fread( &value, 1, sizeof(value), FpLut);
// 						lut_fet[combi] = value;
// 						//cout << std::setprecision(20);
// 						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << lut_fet[combi] << endl;
// 					}
// 				}
// 			}
// 		}
// 		fclose(FpLut);
// 
// 	} else {
// 		cout << "Building " << exep << "/lut_fet.dat" << endl;
// 		//build_lut_fet_singleThread(50);//takes 102 seconds on a 3.5Ghz CPU
// 		//build_lut_fet_singleThread(60);//takes 203 seconds on a 3.5Ghz CPU
// 		//build_lut_fet(100, 24); //using the CC version 'fisher.c', it takes 23 seconds on a 3.5Ghz CPU
// 		int tableMax = 30;
// 		//build
// 		build_lut_fet(tableMax, num_threads); //The fexat.c is not thread safe.
// 		//build_lut_fet_singleThread(tableMax);
// 		//cout << "finish here" << endl;
// 		//write
// 		FILE *fpLut=fopen((exep + "/lut_fet.dat").c_str(), "wb");
// 		double value=(double)tableMax;
// 		fwrite(&value, 1, sizeof(value), fpLut); // table starts with tableMax;
// 		for(int n1 = 1; n1 <= tableMax; n1 ++){
// 			for(int k1 = 0; k1 <= n1; k1 ++){
// 				for(int n2 = 1; n2 <= tableMax; n2 ++){
// 					for(int k2 = 0; k2 <= n2; k2 ++){
// 						//int varray[] = {k1, n1-k1, k2, n2-k2};
// 						//fet_1k( value, varray, 2, 2 );
// 						MultiKey combi(n1,k1,n2,k2);
// 						value = lut_fet[combi];
// 						fwrite(&value, 1, sizeof(value), fpLut);
// 						//cout << std::setprecision(20);
// 						//cout << n1 << "\t" << k1 << "\t" << n2 << "\t" << k2 << "\t" << value << endl;
// 					}
// 				}
// 			}
// 		}
// 		fclose(fpLut);
// 	}




	return 0;
}




int main(int argc, char *argv[])
{
	//cout << boost::filesystem::system_complete(argv[0]) << endl; return 0;
	//double p_acc = pdiffInRegion(10,23,5,38,0.01);
	//cout << p_acc<< endl;
//	double value = -1;
//	int varray1[] = {0, 1, 1, 3};
//	fet_1k( value, varray1, 2, 2 );
//	cout << std::setprecision(20) << "0,1,1,3 = " << value << endl;
//	int varray2[] = {0, 1, 3, 0};
//	fet_1k( value, varray2, 2, 2 );
//	cout << std::setprecision(20) << "0,1,3,0 = " << value << endl;
//	cout << std::setprecision(20) << "0,1,1,3 = " << fet2x2(0,1,1,3) << endl;
//	cout << std::setprecision(20) << "0,1,3,0 = " << fet2x2(0,1,3,0) << endl;
//	0,1,1,3 = 1.0000000000000004441
//	0,1,3,0 = 0.24999999999999994449
//	0,1,1,3 = 1
//	0,1,3,0 = 0.25

//	build_lut_singleCI(100, 14);// this takes 1 second on 3.5Ghz cpu and makes it not import if lookup table is created.
//	build_lut_singleCI_singleThread(200);// this takes 1.2 second on 3.5Ghz cpu and makes it not import if lookup table is created.
//	build_lut_pdiffInRegion(30, 14);//this takes 2 minutes on 3.5Ghz cpu.
//	build_lut_pdiffInRegion(30, 24);//this takes 96 seconds on 3.5Ghz cpu.
//	build_lut_pdiffInRegion(40, 24);//this takes 331 seconds on 3.5Ghz cpu.
//	build_lut_pdiffCI(30, 14);//this takes 10 minutes on 3.5Ghz Cpu.
//	build_lut_pdiffCI(30, 24);//this takes 430 seconds on 3.5Ghz Cpu.
//	build_lut_pdiffInRegion_singleThread(20);
//	exit(0);
//	int v[3] = { 8,12, 16 };
//	std::vector<int> k(&v[0], &v[0]+3);
//	int w[3] = { 25,25, 40 };
//	std::vector<int> n(&w[0], &w[0]+3);
//	BiKey fits(-1, -1);
//	BetaBinomialFit(n, k, fits); //fits = <int> n1, <int> k1
//	cout << "fits=" << fits.n1 << ", " << fits.k1 << endl;
//	exit(0);
try{
	parse_options(argc, argv);
	Opts GO(opts);//global options


	//GO.out();

	if(GO.doMergeRatioFiles){
		cout << "#######################################" << endl;
		cout << "Start    merging of ratio files for withVariance=" << GO.withVariance << endl;
		for(unsigned int i = 0; i < GO.ratiosFiles.size(); i++ ){
			vector <string> toMergeFiles;
			boost::split(toMergeFiles, GO.ratiosFiles[i], boost::is_any_of(","));
			// mergeRatioFiles(toMergeFiles, GO.mergedRatiosFiles[i], GO);
			mergeRatioFilesByChrom(toMergeFiles, GO.mergedRatiosFiles[i], GO);
			cout << " done merging of ratio files " << GO.ratiosFiles[i] << " into " << GO.mergedRatiosFiles[i] << endl;
		}
		cout << "Finished merging of ratio files for withVariance=" << GO.withVariance << endl;
		cout << "#######################################" << endl << endl;
	}


	if(GO.doStrandSpecifiMeth){
		cout << "#######################################" << endl;
		cout << "Start strand specific methylation" << endl;
		for(unsigned int i = 0; i < GO.mergedRatiosFiles.size(); i++ ){
			strandSpecificMethForLane(GO.mergedRatiosFiles[i]);
		}
		cout << "Finished strand specific methylation" << endl;
		cout << "#######################################" << endl;
	}



	if(GO.doComp){

		if(GO.mergedRatiosFiles.size() < 2){
			cout << "Please input 2 or more sample files for comparison between samples" << endl;
			return 0;
		}

		cout << "Starting to do comparisons." << endl;

// Refactored to `get_exepath()` in `types.h` and `types.cpp`
// 		//string exep(  (char *)getauxval(AT_EXECFN) );
// 		string exep = do_readlink("/proc/self/exe");
// 		vector<string> splits;
// 		boost::split(splits, exep, boost::is_any_of("/"));
// 		splits.pop_back();
// 		exep = boost::algorithm::join(splits, "/");
// 		std::cout << exep << std::endl;
//		boost::filesystem::path exepath(argv[0]);
//		string exep = exepath.parent_path().string();
//		if(exep == ""){ //that means argv[0] is just a command without path infomation
//			string cmd ("which ");
//			cmd = cmd + argv[0];
//			exep = std::system(cmd.c_str());
//			vector<string> splits;
//			boost::split(splits, exep, boost::is_any_of("/"));
//			splits.pop_back();
//			exep = boost::algorithm::join(splits, "/");
//		}
		//this works fine if the command is a complete path like /path/bin/mcomp or ./mcomp
		//it breaks if it is just directly called like mcomp -r a.bed -r b.bed ...
		//in this situation, exepath is "mcomp" then "" (ie null)


		//boost::filesystem::path full_path( boost::filesystem::initial_path <boost::filesystem::path> () );
	    //full_path = boost::filesystem::system_complete( boost::filesystem::path( argv[0] ) );
	    //string exep = full_path.parent_path().string();
	    //this doesnot improve the situation like mcomp -r a.bed ...
	    //if your current working dir is /pathA and mcomp is /pathB/mcomp, it returns "/pathA/mcomp" then "/pathA"

		string exep = get_exepath();

		load_lut(GO.threads, exep);
		cout << "Finished loading lookup tables" << endl;
		//return 0;
//		build_lut_pdiffInRegion(30, GO.threads);
//		build_lut_pdiffCI(30, GO.threads);

		string fileName = GO.compFile;

		//read each mapped file into a hash which does not record which strand has whah info.
		map <int, map <string, map<int, cMeth> > > lane ; //lane->chrom->loc->(totalC, methC, strand, nextN);
		for(unsigned int i = 0; i < GO.mergedRatiosFiles.size(); i++ ){
			map <string, map<int, cMeth> > meth ; //chrom->loc->(totalC, methC, strand, nextN);
			readLaneToHash( GO.mergedRatiosFiles[i], GO, meth );
			lane[i] = meth;
			cout << " finish reading file " << i <<endl;
		}
		cout << " finish reading" << endl;

		doComp(lane,fileName, GO);
		cout << "Finished comparison." << endl;

	}


	if(1){

		cout << "Start to call DMCs" << endl;
		//
		//start to call DMC and DMR between lane i and j. //read only relevant info to save RAM.
		//this part is detached from the previous comparison, for the posibility that it may be separted later on.

		double minCredibleDif = 0;

		for( unsigned int i = 0; i < GO.labels.size(); i++ )
		{
			for( unsigned int j = i + 1; j < GO.labels.size(); j++ )	//compare i with j;
			{
				if(GO.doDmcScan || GO.doDmrScan){
					//for(int tempMethod = 1; tempMethod <=2; tempMethod ++){
						if(GO.dmrMethods & 0x1 ){
							//map <int, map <string, map<int, cMeth> > > lane ; //lane->chrom->loc->(totalC, methC, strand, nextN);
							map <int, map <string, map<int, PairCompOut> > > pairResult ; //lane->chrom->loc->(totalC, methC, strand, nextN);
							//hash function: lane i vs lane j ==> lane 100i100j // nobody would run with 100+ samples
							readCompToHashM1(GO.compFile, i, j, pairResult, GO);
							if(GO.doDmrScan){
								DmcDmrM1(i, j, pairResult, GO);
							}
						}
						if( GO.dmrMethods & 0x2 ){
							statSim genomeSim;
							//map <int, map <string, map<int, cMeth> > > lane ; //lane->chrom->loc->(totalC, methC, strand, nextN);
							map <int, map <string, map<int, PairCompOut> > > pairResult ; //lane->chrom->loc->(totalC, methC, strand, nextN);
							//hash function: lane i vs lane j ==> lane 100i100j // nobody would run with 100+ samples
							readCompToHashM2(GO.compFile, i, j, pairResult, genomeSim, GO);
							if(GO.doDmrScan){
								DmcDmrM2(i, j, pairResult, genomeSim, GO);
							}
						}
						if( GO.dmrMethods & 0x4 ){
							//string inputFile = "dmc_M2_" + itos(i) + "_vs_" + itos(j) + ".txt";
							string inputFile  = "dmc_M2_" + GO.labels[i] + "_vs_" + GO.labels[j] + ".txt";
							string outputFile = "dmr_M3_" + GO.labels[i] + "_vs_" + GO.labels[j] + ".txt";
							//string outputFile = "dmr_M3_" + itos(i) + "_vs_" + itos(j);
							HMMCONFIG conf;
							hmm(inputFile, outputFile, conf);
						}

				}

				cout << "start to check predefined regions" << endl;

				//TODO:the output file name need change
				if(GO.predefinedFeature.size() > 0){
					map <int, map <string, map<int, PairCompOut> > > pairResult ; //lane->chrom->loc->(totalC, methC, strand, nextN);
					statSim genomeSim;
					readCompToHashPred(GO.compFile, i, j, pairResult, genomeSim, GO);
					for(unsigned int m = 0; m < GO.predefinedFeature.size(); m++){
						vector<string> features;
						boost::split(features, GO.predefinedFeature[m], boost::is_any_of(","));
						for(unsigned int n = 0; n < features.size(); n++){
							cout << "start to process predefined region " << features[n] << endl;
							isPredefinedRegionDmr(features[n], GO.compFile, i, j, pairResult, genomeSim, GO);
						}
					}

				}
			}
		}


	}
	// */


	cout << "Program successfully finished"<<endl;
}
catch (exception& e)
{
  cout << e.what() << endl;
}
catch(...)
{
	std::cerr <<"unknown error" <<std::endl;
}

return 0;

/*TODO
 * 1.
 */




















}
