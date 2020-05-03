/* todo: need make opts quality score local because different files may have different quality score base */
// change next to nextBase, totalC to pTotalC, methC to mMethC, ...
// There might be a bug for bisulfite conversion ratio. See the 0 in strand_specific_meth.xlsx->call.
// todo: add an option to set cut size automatically by model
// done: merge stats from two+ files
// done: bug when dealing with overlapped PE reads
// another project: overhang size and remove duplicate reads and use of (random?) adpaters to avoid removal of duplicate reads
// another project: methylated regions that deviate the average: that's predefined region. Like methylated CGI/PROMEOTER. UNMETHYLATED GENE BODY.
// todo: call mcomp at end of program
// todo: test dense_map to see if there's improvement on speed
// todo: read whole dna into ram and drop XR dependence; drop ZS dependence by flag operations
// for my wgbs project: suggest to set --minMMFragSize 103
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
//#include <cmath>
//#include <math.h>
#include <string>
//#include <unistd.h>
#include <vector>
#include <algorithm>
#include <map>
//#include <valarray>
//#include <numeric>
#include <iterator>
////#include <typeinfo>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
////#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
//#include <boost/date_time.hpp>

//depth for saturation on statistics
#define DEPTHSAT 40

namespace po = boost::program_options;

using namespace std;

#include <sys/types.h>
#include <unistd.h>

#include <sstream>
#include <zlib.h>
#include "sam.h"
//#include "bam.h"

//this line maybe putin to function parse_options
po::variables_map options;


string itos(int i)
{
    stringstream s;
    s << i;
    return s.str();
}

string fastaheader2id(const string & header) {
	string id;
	if (header.empty() || header[0]!='>') return id;
	vector< string > fields;
	boost::split(fields, header, boost::is_any_of(" \t")); // space before possible descriptions
	if (!fields.empty()) {
		id=fields[0].substr(1);
	}
	return id;
}

//read DNA file once and create map chrom -> sequence takes about 2.5GB for mm9.
//to save some RAM is why I do it by chrom.
void dnaByChr(string FILE_NAME, string chrom, string & fatable, vector <int> & numc)
{
	int found = 0;
	string line ;
	ifstream inputf;
	inputf.open(FILE_NAME.c_str(), ifstream::in);
	int nc 		= 0;
	int ncg 	= 0;
	int nch 	= 0;
	int nchg 	= 0;
	int nchh 	= 0;
	while (inputf.good()) {
		getline(inputf, line); //line already chomps the trailing "\n"
		if( (line == "") || (line == "\n") ){continue;} //line == "" means it's blank line or EOF
		if(found){
			if(line[0] == '>'){ //start of next chrom;finish reading
				for(unsigned int i = 0; i < fatable.size(); i++){
					if(fatable[i] == 'C'){
						nc ++;
						if(fatable[i+1] == 'G'){
							ncg ++;
						} else {
							nch ++;
							if(fatable[i+2] == 'G'){
								nchg ++;
							} else {
								nchh ++;
							}
						}
					}
				}
				numc.push_back(nc);
				numc.push_back(ncg);
				numc.push_back(nch);
				numc.push_back(nchg);
				numc.push_back(nchh);
				break;
			}
			boost::to_upper(line);

			fatable += line;
		}
		if( found == 0 && line[0]=='>' && fastaheader2id(line)==chrom ){ //start reading
			found = 1;
		}
	}
	inputf.close();

	//##process the last chrom;
	for(unsigned int i = 0; i < fatable.size(); i++){
		if(fatable[i] == 'C'){
			nc ++;
			if(fatable[i+1] == 'G'){
				ncg ++;
			} else {
				nch ++;
				if(fatable[i+2] == 'G'){
					nchg ++;
				} else {
					nchh ++;
				}
			}
		}
	}
	numc.push_back(nc);
	numc.push_back(ncg);
	numc.push_back(nch);
	numc.push_back(nchg);
	numc.push_back(nchh);
}


class Opts //Options
{
public:
	int threads;
	vector <string> mappedFiles;
	string outputDir;
	string webOutputDir;
	string sampleName;
	string genome;
	string reference;
	int cytosineMinScore;
	int nextBaseMinScore;

	//this option is absorbed by fullMode.
	//int reportSkippedBase;
	int qualityScoreBase;
	int trimWGBSEndRepairPE2Seq;
	int trimWGBSEndRepairPE1Seq;
	int processPEOverlapSeq;
	int trimRRBSEndRepairSeq;
	int skipRandomChrom; 	//skip random chrom?
	int requiredFlag;
	int excludedFlag;
	int minFragSize;
	int minMMFragSize;
	int reportCpX; 			//'G': output CG methy file
	int reportCHX; 			//'G': output CHG methy file
	int fullMode;			//if off, only *.G.bed and *_stat.txt are generated; if on, *.bed, *_skip.bed, and *_strand.bed are generated.
	int statsOnly;
	int keepTemp;
	int reportStrand;		//output *_strand.bed? //this file has been incorporated in *bam.bed file and is only for debuging purpose now.
	int reportSkip; 		//output *_skip.bed?
	int reportCombined; 	//output *.bed?
	Opts(){
		threads = 1;
		outputDir = "./";
		webOutputDir = "./";
		sampleName = "mSuite";
		genome = "";
		reference = "";
		cytosineMinScore = 20;
		nextBaseMinScore = 3;
		//reportSkippedBase = 1;
		qualityScoreBase = 0; //0 means autodetection; Sanger=>33;Solexa=>59;Illumina=>64; See wiki "FASTQ_format" for details
		trimWGBSEndRepairPE2Seq = 3;
		trimWGBSEndRepairPE1Seq = 3;
		processPEOverlapSeq = 1;
		trimRRBSEndRepairSeq = 2;
		skipRandomChrom = 1;
		requiredFlag = 0;
		excludedFlag = 0;
		reportCpX = 'G';
		reportCHX = 'X';
		fullMode = 0;
		statsOnly = 0;
		keepTemp = 0;
		reportStrand = 0;
		reportSkip = 0;
		reportCombined = 0;
	};
} opts;


int parse_options(int ac, char * av[]){
	po::options_description desc("Allowed options;");
	desc.add_options()
	("help,h", 								"Produce help message. Common options are provided with single letter format. Parameter defaults are in brackts. Example command: mCall -m Ko.bam; mCall -m wt_r1.bam -m wt_r2.bam -sampleName Wt; See doc for more details.")
	("mappedFiles,m",						po::value< vector<string> >()->multitoken(), "Specify the names of RRBS/WGBS alignment files for methylation calling. Multiple files can be provided to combine them(eg. lanes or replicates) into a single track;")
	("sampleName", 							po::value<string>(), "If two or more mappedFiles are specifed, this option generates a merged result; Ignored for one input file;")
	("outputDir",							po::value<string>(), "The name of the output directory;")
	("webOutputDir",						po::value<string>(), "The name of the web-accessible output directory for UCSC Genome Browser tracks;")
	("genome,g", 							po::value<string>(), "The UCSC Genome Browser identifier of source genome assembly; mm9 for example;")
	("reference,r", 						po::value<string>(), "Reference DNA fasta file; It's required if CHG methylation is wanted;")
	("cytosineMinScore", 					po::value<int>()->default_value(20), "Threshold for cytosine quality score (default: 20). Discard the base if threshold is not reached;")
	("nextBaseMinScore", 					po::value<int>()->default_value(3), "Threshold for the next base quality score(default: 3,ie, better than B or #); Possible values: -1 makes the program not to check if next base matches reference; any positive integer or zero makes the program to check if next base matches reference and reaches this score threshold;")
	//("reportSkippedBase", 				po::value<int>()->default_value(1), "Specify if bases that are not accepted for methylation analysis should be written to an extra output file;")
	("qualityScoreBase", 					po::value<int>()->default_value(0), "Specify quality score system: 0 means autodetection; Sanger=>33;Solexa=>59;Illumina=>64; See wiki FASTQ_format for details;")
	("trimWGBSEndRepairPE2Seq", 			po::value<int>()->default_value(3), "How to trim end-repair sequence from begin of +-/-- reads from Pair End WGBS Sequencing; 0: no trim; n(positive integer): trim n bases from begin of +-/-- reads; -2: model determined n; -1: trim from beginning to before 1st methylated C; Suggest 3; n>readLen is equivalent to use PE1 reads;")
	("trimWGBSEndRepairPE1Seq", 			po::value<int>()->default_value(3), "How to trim end-repair sequence from end   of ++/-+ reads from Pair End WGBS Sequencing; 0: no trim; n(positive integer): trim n + NM bases from end of ++/-+ reads if fragSize <= maxReadLen; -2: model determined n; Suggest 3;")
	("processPEOverlapSeq", 				po::value<int>()->default_value(1), "1/0 makes the program count once/twice the overlap seq of two pairs;")
	("trimRRBSEndRepairSeq",				po::value<int>()->default_value(2), "How to trim end-repair sequence for RRBS reads; RRBS or WGBS protocol can be automatically detected; 0: no trim; 2: trim the last CG at exactly end of ++/-+ reads and trim the first CG at exactly begin of +-/-- reads like the WGBS situation;")
	("skipRandomChrom", 					po::value<int>()->default_value(1), "Specify whether to skip random and hadrop chrom;")
	("requiredFlag,f", 						po::value<int>()->default_value(0), "Requiring samtools flag; 0x2(properly paried), 0x40(PE1), 0x80(PE2), 0x100(not unique), r=0x10(reverse); Examples: -f 0x10 <=> +-/-+ (Right) reads; -f 0x40 <=> ++/-+ (PE1) reads; -f 0x50 <=> -+ read; -f 0x90 <=> +- read;")
	("excludedFlag,F", 						po::value<int>()->default_value(0), "Excluding samtools flag; Examples: -f 0x2 -F 0x100 <=> uniquely mapped pairs; -F 0x10 <=> ++/-- (Left) reads; -F 0x40 <=> -f 0x80 +-/-- (PE2) reads; -f 0x40 -F 0x10 <=> ++ read; -f 0x80 -F 0x10 <=> -- read; ")
	("minFragSize", 						po::value<int>()->default_value(0), "Requiring min fragment size, the 9th field in sam file; Since non-properly-paired read has 0 at 9th field, setting this option is requiring properly paired and large enough fragment size;")
	("minMMFragSize", 						po::value<int>()->default_value(0), "Requiring min fragment size for multiply matched read; Same as option above but only this option is only applicable to reads with flag 0x100 set as 1;")
	("reportCpX", 							po::value<char>()->default_value('G'), "X=G generates a file for CpG methylation; A/C/T generates file for CpA/CpC/CpT meth;")
	("reportCHX", 							po::value<char>()->default_value('X'), "X=G generates a file for CHG methylation; A/C/T generates file for CHA/CHC/CHT meth; This file is large;")
	("fullMode,a", 							po::value<int>()->default_value(0), "Specify whether to turn on full mode. Off(0): only *.G.bed, *.HG.bed and *_stat.txt are allowed to be generated. On(1): file *.HG.bed, *.bed, *_skip.bed, and *_strand.bed are forced to be generated. Extremely large files will be generated at fullMode.")
	("statsOnly", 							po::value<int>()->default_value(0), "Off(0): no effect. On(1): only *_stat.txt is generated.")
	("keepTemp", 							po::value<int>()->default_value(0), "Specify whether to keep temp files;")
	("threads,p",							po::value<int>()->default_value(1),"Number of threads on all mapped file. Suggest 1~8 on EACH input file depending RAM size and disk speed.")
	;
	//toadd
	//minDepthForSingleSampleProfiling 10//use this for various sample profiling and etc.

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
	//int cur_pid = pid_t ;
	stringstream configFileName;
	configFileName << "run.config." << getpid();
	configFile.open(configFileName.str().c_str(), ios_base::app);

	for(std::map<string,po::variable_value>::iterator iter = options.begin(); iter != options.end(); ++iter)
	{
		string k =  (*iter).first;

		cout		<<k<<"=";
		configFile 	<<k<<"=";


		if( k == "threads"){
			opts.threads 	= 	options[k].as<int>();
			configFile 		<<	options[k].as<int>();
			cout 			<<	options[k].as<int>();
		}
		else if( k == "mappedFiles"){
			opts.mappedFiles = options[k].as< vector<string> >();
			for(unsigned int i = 0 ; i < opts.mappedFiles.size(); i++)
			{
				cout 		<< opts.mappedFiles[i]<<" ";//copy(myvector.begin(), myvector.end(), ostream_iterator<int>(cout, "\n"));
				configFile 	<< opts.mappedFiles[i]<<" ";

				if ( !boost::filesystem::exists( opts.mappedFiles[i] ) )
				{
				  cout << endl << "Can't find file " << opts.mappedFiles[i] << endl;
				  exit(1);
				}

			}
		}
		else if( k == "webOutputDir"){
			opts.webOutputDir 	= 	options[k].as<string>();
			configFile 			<<	options[k].as<string>();
			cout 				<<	options[k].as<string>();
		}
		else if( k == "outputDir"){
			opts.outputDir 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k== "sampleName"){
			opts.sampleName 	= 	options[k].as<string>();
			configFile 			<<	options[k].as<string>();
			cout 				<<	options[k].as<string>();
//			if(opts.mappedFiles.size() != 1){
//				cerr 			<< "sampleName ignormed because multiple mappedFiles are provided!" <<endl;
//				break;
//			}
		}
		else if( k == "genome"){
			opts.genome 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k == "reference"){
			opts.reference 	= 	options[k].as<string>();
			configFile 		<<	options[k].as<string>();
			cout 			<<	options[k].as<string>();
		}
		else if( k == "cytosineMinScore"){
			opts.cytosineMinScore 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "nextBaseMinScore"){
			opts.nextBaseMinScore 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
//		else if( k == "reportSkippedBase"){
//			opts.reportSkippedBase 	= 	options[k].as<int>();
//			configFile 				<<	options[k].as<int>();
//			cout 					<<	options[k].as<int>();
//		}
		else if( k == "qualityScoreBase"){
			opts.qualityScoreBase 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "trimWGBSEndRepairPE2Seq"){
			opts.trimWGBSEndRepairPE2Seq 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "trimWGBSEndRepairPE1Seq"){
			opts.trimWGBSEndRepairPE1Seq 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "processPEOverlapSeq"){
			opts.processPEOverlapSeq 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "trimRRBSEndRepairSeq"){
			opts.trimRRBSEndRepairSeq 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "skipRandomChrom"){
			opts.skipRandomChrom 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "requiredFlag"){
			opts.requiredFlag 		= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "excludedFlag"){
			opts.excludedFlag 		= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "minFragSize"){
			opts.minFragSize 		= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "minMMFragSize"){
			opts.minMMFragSize 	= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "reportCpX"){
			opts.reportCpX 	= 	options[k].as<char>();
			configFile 				<<	options[k].as<char>();
			cout 					<<	options[k].as<char>();
		}
		else if( k == "reportCHX"){
			opts.reportCHX 	= 	options[k].as<char>();
			configFile 				<<	options[k].as<char>();
			cout 					<<	options[k].as<char>();
		}
		else if( k == "fullMode"){
			opts.fullMode 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "statsOnly"){
			opts.statsOnly 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else if( k == "keepTemp"){
			opts.keepTemp 			= 	options[k].as<int>();
			configFile 				<<	options[k].as<int>();
			cout 					<<	options[k].as<int>();
		}
		else{
			cerr << "Please modify code for this new type."<<endl;
		}

		cout			<<endl;
		configFile 		<<endl;
	}
	configFile.close();
	configFile.clear();

	if( options["fullMode"].as<int>() == 0 )
	{
		opts.reportCombined = 0;
		opts.reportSkip = 0;
		opts.reportStrand = 0;
	}

	if( options["mappedFiles"].empty() )
	{
		cout << "Mandatory parameters missing. Program will terminate now."<<endl;
		exit(1);
	}
	return 0;
};

//not used;
void workerFunc(const char* msg, unsigned delaySecs)
{
    boost::posix_time::seconds workTime(delaySecs);

    cout << msg << endl;

    // Pretend to do something useful...
    boost::this_thread::sleep(workTime);

    cout << "Worker: finished after delaying " << delaySecs <<  endl;
}

//not used;
map<string,int> create_chrIndex()
{
  map<string,int> chri;

  chri["chr0"] = 0;
  chri["chr1"] = 1;
  chri["chr2"] = 2;
  chri["chr3"] = 3;
  chri["chr4"] = 4;
  chri["chr5"] = 5;
  chri["chr6"] = 6;
  chri["chr7"] = 7;
  chri["chr8"] = 8;
  chri["chr9"] = 9;
  chri["chr10"] = 10;
  chri["chr11"] = 11;
  chri["chr12"] = 12;
  chri["chr13"] = 13;
  chri["chr14"] = 14;
  chri["chr15"] = 15;
  chri["chr16"] = 16;
  chri["chr17"] = 17;
  chri["chr18"] = 18;
  chri["chr19"] = 19;
  chri["chr20"] = 20;
  chri["chr21"] = 21;
  chri["chr22"] = 22;
  chri["chr23"] = 23;
  chri["chrL"] = 76;
  chri["chrM"] = 77;
  chri["chrX"] = 88;
  chri["chrY"] = 89;
  chri["chrZ"] = 90;
  chri["chrUn"]= 85;

  chri["chr0_random"] = 100;
  chri["chr1_random"] = 101;
  chri["chr2_random"] = 102;
  chri["chr3_random"] = 103;
  chri["chr4_random"] = 104;
  chri["chr5_random"] = 105;
  chri["chr6_random"] = 106;
  chri["chr7_random"] = 107;
  chri["chr8_random"] = 108;
  chri["chr9_random"] = 109;
  chri["chr10_random"] = 110;
  chri["chr11_random"] = 111;
  chri["chr12_random"] = 112;
  chri["chr13_random"] = 113;
  chri["chr14_random"] = 114;
  chri["chr15_random"] = 115;
  chri["chr16_random"] = 116;
  chri["chr17_random"] = 117;
  chri["chr18_random"] = 118;
  chri["chr19_random"] = 119;
  chri["chr20_random"] = 120;
  chri["chr21_random"] = 121;
  chri["chr22_random"] = 122;
  chri["chr23_random"] = 123;
  chri["chrM_random"] = 177;
  chri["chrX_random"] = 188;
  chri["chrY_random"] = 189;
  chri["chrZ_random"] = 190;
  chri["chrUn_random"] = 185;

  return chri;
}

//not used; //global variable chri
map <string, int> chri = create_chrIndex();

//not used;
string chrIndexToName(int x){
	stringstream chrom;
	chrom << "chr";
	if(x % 100 >=77){ //M/X/Y/Z/ random
		chrom << char( x % 100);
	} else {
		chrom << (x % 100);
		if( (x%100)== 85 ){chrom << 'n';}
	}
	if(x >= 100){
		chrom << "_random";
	}
	return chrom.str();
}


class fileProperty
{
public:
	string format; //bam|sam|bsp
	int    qualScoreBase;
	long long   nAllReads;
	long long  nMappedReads;
	string protocol; //RRBS or WGBS
	int    readLen;

	fileProperty(){
		format = "";
		qualScoreBase = 0;
		nAllReads = 0;
		nMappedReads = 0;
		protocol = "WGBS";
		readLen = 100;
	};
	fileProperty(string f, int q, long long na, long long nm, string p, int rl) : format(f), qualScoreBase(q), nAllReads(na), nMappedReads(nm), protocol(p), readLen(rl) {};
};


//global variable chroms;
vector <string> chroms; //index may not be tid
map <string, map<int, string> > chromsByLane; //chromsByLane[laneName][tid] = chromName
map <string, fileProperty> propertiesByLane; //propertiesByLane[laneName] = fileProperty
vector <string> chromNames; //merged chromNames from all lanes


void getCharSetFromString(set <char> & chars, string & s)
{
	for(unsigned int i = 0; i < s.size();i++){
		chars.insert(s[i]);
	}
}

void set_insert(set<char> & set1, set <char> & set2)
{
	for(set<char>::iterator it = set2.begin(); it != set2.end(); it++){
		set1.insert(*it);
	}
};


int string_to_int( string s){return atoi(s.c_str());}

char rctable(char m)//reverse complimentary of bases
{
	switch(m){
	case 'A':
		return 'T';
	case 'C':
		return 'G';
	case 'G':
		return 'C';
	case 'T':
		return 'A';
	case 'X':
		return 'Y';
	case 'N':
		return 'M';
	default:
		cerr << "unexpected sequence " << m <<endl;
		exit(1);
	}
};

char bsPosTable(char m) //bisulfite conversion table of bases on Postive strand
{
	if(m == 'C'){
		return 'T';
	}else{
		return m;
	}
}

char bsNegTable(char m) //bisulfite conversion table of bases on Negative strand
{
	if(m == 'G'){
		return 'A';
	}else{
		return m;
	}
}

class alignedRead{
public:
	string name;
	string chrom;
	int pos;
	char strand; 		//1st char in strand
	char strandCopy; 	//2nd char in strand
	string seq;
	string ref;
	string qual;
	alignedRead(){
		name = "*";
		chrom = "*";
		pos = 0;
		strand = '*';
		strandCopy = '*';
		seq = "*";
		ref = "*";
		qual = "*";
	};
};

//not currently used

alignedRead getAlignment(string line, string format)
{
	vector <string> fields;
	boost::split(fields, line, boost::is_any_of("\t"));
	alignedRead read;
	if(format == "SAM"){
		if(fields[1].find('u') != string::npos )return read;
		read.chrom = fields[2];
		read.pos = string_to_int(fields[3]) - 1;
		read.seq = fields[9];
		read.qual = fields[10];

		size_t strandPos = line.find("ZS:Z:");
		if(strandPos != string::npos){
			read.strand = line[strandPos+5];
			read.strandCopy = line[strandPos+6];
		}
		else{
			cerr << "ZS:Z:<str> field not found; Is it BAM/SAM format!\n";
			exit(1);
		}

		size_t refPos = line.find("R:Z:");
		if(refPos != string::npos){
			string ref = "";
			for(unsigned int i = 4; i < line.size(); i++ ){
				if(line[refPos+i] == 'A' || line[refPos+i] == 'C' ||line[refPos+i] == 'G' ||line[refPos+i] == 'T' ){
					ref += line[refPos+i] ;
				}
				else{		//end of ref field;
					break;
				}
			}
		}
		else{
			cerr << "R:Z:<str> field not found; Please enable -R option in BSMAP!\n";
			exit(1);
		}

		return read;

	}
	else if( format == "BSP" ){
		if(fields[3] == "NM" || fields[3] == "QC") return read;
		read.chrom = fields[4];
		read.pos = string_to_int(fields[5]);
		read.strand = fields[6][0];
		read.strandCopy = fields[6][1];
		read.seq = fields[1];
		read.ref = fields[8];
		return read;
		//there's no quality output in bsp format
	}
	else{
		cerr << "only SAM/BAM/BSP format are supported\n";
		exit(1);
	}
}

//not used
void performCallThreadMode(string file)//not finished; not currently used..
{
	cout << file <<endl;

	ofstream skipBpFile;
	skipBpFile.open((file + "_skip.bed").c_str(), ios_base::out);
	skipBpFile << "#chrom\tstart\tend\tratio\ttotalC\tmeanC\tstrand\tnext" <<endl;


	string format = file.substr(file.find_last_of(".") + 1);
	boost::to_upper(format);
	if( format == "SAM") {
		cout << "SAM..." << endl;
	}
	else if( format == "BAM" ) {
		cout << "BAM..." << endl;
	}
	else if( format == "BSP" ) {
		cout << "BSP..." << endl;
	}

	int BatchNum = 100;
	vector<alignedRead> reads;
	reads.resize(BatchNum);
	vector<alignedRead>::iterator p=reads.begin();

	samfile_t * BAM_fp;
	bam1_t *BAM_b;
	BAM_fp = samopen(file.c_str(), "rb", 0);
	BAM_b = bam_init1();

	while(samread(BAM_fp,BAM_b) > 0 && p != reads.end() )
	{
		p->name = string((char*)bam1_qname(BAM_b));
		p->chrom = BAM_fp->header->target_name[BAM_b->core.tid];
		p->ref = string( (char *) bam_aux2Z(bam_aux_get(BAM_b, "XR") ) );
		string zs = string( (char *) bam_aux2Z(bam_aux_get(BAM_b, "ZS") ) );
		p->strand = zs[0];
		p->strandCopy = zs[1];

		//following are just retrieving seq and qual.
		size_t l_seq=BAM_b->core.l_qseq;
   		p->seq.assign(l_seq,0);
   		p->qual.assign(l_seq,0);
		char * s = (char*) bam1_seq(BAM_b);
		char * t = (char*) bam1_qual(BAM_b);
		//string a = string( (char*) bam1_aux(BAM_b) ); //why doesn't this work?
		//int f2 = bam_aux2i(bam_aux_get(BAM_b, "NM") ) ;
		for(size_t i=0;i<l_seq;i++){
			p->seq[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
			p->qual[i]=t[i]+33;
		}
	}
}

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

class stat
{
public:
	long long sites; //gcc defines long and long long  as 8 bytes but visual studio defines long and long long as 4 and 8 bytes.
	double mean;
	long long methC;
	long long totalC;
	double depth;
	stat(){
		sites = -1;
		mean = -1;
		methC = -1;
		totalC = -1;
		depth = -1;
	};
	stat(long long s, double mn, long long mC, long long t, double d) : sites(s), mean(mn), methC(mC), totalC(t), depth(d) {};
};


class statByDepth
{
public:
	long long ncg;
	long long nch;
	long long nchg;
	double rcg;
	double rch;
	double rchg;

	//these can be derived from above;
	long long nc; 		//nc = ncg + nch;
	long long nchh;		//nchh = nch - nchg
	double rc;		//rc = ( rcg*ncg + rch*nch ) / nc;
	double rchh;	//rchh = (rch*nch - rchg*nchg) / nchh;

	statByDepth(){
		ncg = 0;
		nch = 0;
		nchg = 0;
		rcg = 0.0;
		rch = 0.0;
		rchg = 0.0;

		nc = 0;
		nchh =0;
		rc = 0.0;
		rchh = 0.0;
	}
};

//global variables for stats file
//vector < map<char, stat> >  statbVec; //each element is for each chrom; //statbVec[chromIndex][next]=statObj
//vector < map<char, map<char, stat> > >  statsVec; //each element is for each chrom; //statsVec[chromIndex][strand][next]=statObj
//vector < map<int, statByDepth> > statbByDepthVec; //each element is for each chrom;

//todo: bug: very very tiny chance of thread conflict;
map <string, vector <int> > numcByChrom; //each key is for each chrom; //numcByChrom[chromName] = v[];

map <string, map < string, map<char, stat> > > statbByLane; //statbByLane[laneName][chromName][next] = statObj
map <string, map < string, map<char, map<char, stat> > > > statsByLane; //statsByLane[laneName][chromName][strand][next]=statObj
map <string, map < string, map<int, statByDepth> > > statbByDepthByLane; //statbByLane[laneName][chromName][depth] = statByDepthObj



void qualityScoreFormatDetection(string file, string format, int & qscorebase) //set opts.qualityScoreBase here during auto detection
{
	string ASCII033to058 = "!\"#$%&'()*+,-./0123456789:"; //should be !"#$%&'()*+,-./0123456789: there's double quotation.//SangerOnly
	string ASCII059to063 = "<=>?"; //SolexaOnly illumina v1.0.
	string ASCII064to097 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`a"; //should be @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh there's backward slash;
	string ASCII067to073 = "CDEFGHI";
	string ASCII098to104 = "bcdefgh";
	string ASCII105to126 = "ijklmnopqrstuvwxyz{|}~"; //most probably IlluminaOnly.
	//033-058 => SangerFormat;
	//No 033-058, but 059to063 => SolexasFormat;
	//No 033-063, but 064to096 => IlluminaFormat;
	//097-104 ==> IlluminaFormat;
	//105-126 ==> IlluminaFormat;
	set <char> S033to058;
	set <char> S059to063;
	set <char> S067to073;
	set <char> S098to104;

	getCharSetFromString(S033to058, ASCII033to058);
	getCharSetFromString(S059to063, ASCII059to063);
	getCharSetFromString(S067to073, ASCII067to073);
	getCharSetFromString(S098to104, ASCII098to104);

    int count = 0;
    string line = "";
    set <char> scores;

    if( format == "SAM" )	//open SAM format and read at most 999 lines to push quality scores into set scores;
    {
		ifstream inputf(file.c_str(), ios::in);
		while (inputf.good()) {
			getline(inputf, line);
	        if( line == "" || find(line.begin(), line.begin()+1, '@') == line.begin() ){continue;}

			if(count > 999){
				break;
			}
			count += 1 ;

			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			set <char> score;
			getCharSetFromString(score, fields[10]);
			set_insert(scores, score);
		}
		inputf.close();
	}

    if(format == "BAM" )
    {

		samfile_t * BAM_fp;
		bam1_t *b;
		BAM_fp = samopen(file.c_str(), "rb", 0);
		b = bam_init1();

		alignedRead read;

		while(samread(BAM_fp,b) > 0 )
		{
			if(b->core.tid < 0)continue; //if chr?? is *, then tid = -1;
			if(count > 999){
				break;
			}
			count += 1 ;

			size_t l_seq = b->core.l_qseq;
			read.qual.assign(l_seq,0);
			char * t = (char*) bam1_qual(b);
			for(size_t i = 0; i < l_seq; i++){
				read.qual[i]=t[i]+33;
			}

			set <char> score;
			getCharSetFromString(score, read.qual);
			set_insert(scores, score);
		}
		samclose(BAM_fp);
    }


    //process set scores to guess the quality score system
	vector <char> inter(126);
	//copy(inter.begin(), inter.end(), ostream_iterator<char>(cout, "5"));
	set_intersection (scores.begin(), scores.end(), S033to058.begin(), S033to058.end(), inter.begin());

	cout << "For file " <<file << ", the ";
	if(inter.size() > 0){
		cout << "quality score format is Sanger format based at 33!" <<endl;
		//opts.qualityScoreBase = 33;
		qscorebase = 33;
	}
	else
	{
		set_intersection (scores.begin(), scores.end(), S059to063.begin(), S059to063.end(), inter.begin());
		if(inter.size() > 0){
			cout << "quality score format is Solexa format based at 59!" <<endl;
			//opts.qualityScoreBase = 59;
			qscorebase = 59;
		}
		else
		{
			set_intersection (scores.begin(), scores.end(), S098to104.begin(), S098to104.end(), inter.begin());
			int goodScoreIllumina 	= int(inter.size());
			set_intersection (scores.begin(), scores.end(), S067to073.begin(), S067to073.end(), inter.begin());
			int goodScoreSanger 	= int(inter.size());
			if(goodScoreSanger / goodScoreIllumina >1000){ // is it possible that there are very few scores in "abcdefgh" region for sanger format?
				cout << "quality score format is Sanger format based at 33!" <<endl;
				//opts.qualityScoreBase = 33;
				qscorebase = 33;
			}
			else{
				cout << "quality score format is Sanger format based at 64!" <<endl;
				//opts.qualityScoreBase = 64;
				qscorebase = 64;
			}
		}
	}

}


//todo: need sam version
//it also detects that what is the max read length and if the 1st 2 bases are CG or TG if more than half then it is RRBS protocol. //should be 100% for -- read
int detectChromsReadLenProtocol(string file, string format, int & readLen, string & protocol)
{
	int countCGorTG = 0;
	int maxReadLen = 0;
	protocol = "WGBS";
	int nline = 0;

	if(format == "BAM")
	{
		samfile_t * BAM_fp;
		bam1_t *b;
		BAM_fp = samopen(file.c_str(), "rb", 0);
		b = bam_init1();

		for(int x = 0; x < BAM_fp->header->n_targets; x++){
			chromsByLane[file][x] = BAM_fp->header->target_name[x];
		}

		while(samread(BAM_fp,b) > 0 && nline < 999)
		{
			if(b->core.tid < 0)continue;
			if( b->core.l_qseq > maxReadLen ){
				maxReadLen = b->core.l_qseq;
			}

			if( b->core.flag & 0x80 ) continue; //if this read is PE2 read
			if( b->core.flag & 0x10 ) continue; //if this read is reverse strand
			//left is read mapped to ++ strand

			nline ++;

			//following are just retrieving seq and qual.
			string seqs;
			size_t l_seq=b->core.l_qseq;
			seqs.assign(l_seq,0);
			char * s = (char*) bam1_seq(b);
			for(size_t i=0;i<2;i++){
				seqs[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
			}
			if( ((seqs[0] == 'C') || (seqs[0] == 'T')) && seqs[1] == 'G' ){
				countCGorTG ++;
			}
		}
		samclose(BAM_fp);
	} // end of BAM format

	if(format == "SAM")
	{
		//get chroms
		string line = "";
		ifstream input(file.c_str(), ios::in);
		while( input.good() ) {
			getline(input, line);
			if( find(line.begin(), line.begin()+1, '@') == line.begin() ){
				vector <string> fields;
				boost::split(fields, line, boost::is_any_of("\t"));
				vector <string> fs;
				boost::split(fs, fields[1], boost::is_any_of(":"));

				if(fs[0] == "SN" ){
					chromsByLane[file][nline] = fs[1];
					nline ++;
				}
			} else {
				break;
			}
		}
		input.close();

		//get readLen and protocol
		nline = 0;
		line = "";
		ifstream inputf(file.c_str(), ios::in);
		while (inputf.good() && nline < 999) {
			getline(inputf, line);
			if( line == "" || find(line.begin(), line.begin()+1, '@') == line.begin() ){continue;}

			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));

			if(fields[2] == "*")continue;
			string seqs = fields[9];
			if( seqs.size() > maxReadLen ){
				maxReadLen = seqs.size();
			}

			unsigned int flag = string_to_int(fields[1]);
			if( flag & 0x80 ) continue; //if this read is PE2 read
			if( flag & 0x10 ) continue; //if this read is reverse strand

			nline ++;

			if( ((seqs[0] == 'C') || (seqs[0] == 'T')) && seqs[1] == 'G' ){
				countCGorTG ++;
			}
		}
		inputf.close();
	}	//end of SAM foramt


	if( double(countCGorTG)/nline > 5.0/8.0){
		protocol = "RRBS";
	}
	readLen = maxReadLen;

	cout << "Protocol and read length are detected as " << protocol <<" and " << readLen << " bases for file " << file <<endl;
	return 0;
}


void getReadCounts(string file, string format, long long & nAllReads, long long & nMappedReads)
{
	nAllReads = 0;
	nMappedReads = 0;

	if(format == "BAM"){
		samfile_t * BAM_fp;
		bam1_t *b;
		BAM_fp = samopen(file.c_str(), "rb", 0);
		if (BAM_fp == 0) {
			fprintf(stderr, "Fail to open BAM file %s\n", file.c_str());
			exit(1);
		}


		b = bam_init1();
		while(samread(BAM_fp,b) > 0 )
		{
			nAllReads++;
			if(b->core.tid < 0)continue; //if chr?? is *, then tid = -1; //it's bug that -2 %5 =-2. -2 %1 = 0; -1%1 =0;
			if(opts.skipRandomChrom && (chromsByLane[file][b->core.tid].find("random") != string::npos || chromsByLane[file][b->core.tid].find("had") != string::npos ) ){
				continue;
			}
			nMappedReads ++;
		}
		samclose(BAM_fp);
	}

	if(format == "SAM"){
		string line = "";
		ifstream inputf(file.c_str(), ios::in);
		while (inputf.good() ) {
			getline(inputf, line);
			if( line == "" || find(line.begin(), line.begin()+1, '@') == line.begin() ){continue;}

			nAllReads ++;

			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			if(fields[2] == "*")continue;
			if(opts.skipRandomChrom && (fields[2].find("random") != string::npos || fields[2].find("had") != string::npos ) ){
				continue;
			}
			nMappedReads ++;
		}
		inputf.close();
	}

	cout << "For file " << file << " the number of all reads is " << nAllReads << " and the number of mapped reads is " << nMappedReads << endl;

}


bool found_XR_and_ZS(string file, string format)
{
	bool found = true;

	int nline = 0;

	if(format == "BAM")
	{
		samfile_t * BAM_fp;
		bam1_t *b;
		BAM_fp = samopen(file.c_str(), "rb", 0);
		b = bam_init1();

		while(samread(BAM_fp,b) > 0 && nline < 5)
		{
			if(b->core.tid < 0)continue;
			nline ++;
			uint8_t *ss = bam1_aux(b);
			char XR[2]= {'X', 'R'};
			char ZS[2]= {'Z', 'S'};
			bool foundXR = false;
			bool foundZS = false;

			while (ss < b->data + b->data_len) {
					uint8_t type, key[2];
					key[0] = ss[0];
					key[1] = ss[1];
					ss += 2;
					type = *ss;
					++ss;
					//printf( "\t%c%c:", key[0], key[1]);
					//cout <<endl;
					if(XR[0] == key[0] && XR[1] == key[1]){foundXR = true;}
					if(ZS[0] == key[0] && ZS[1] == key[1]){foundZS = true;}

					if (type == 'A') {
						//printf( "A:%c", *ss);
						++ss;
					}
					else if (type == 'C') {
						//printf( "i:%u", *ss);
						++ss;
					}
					else if (type == 'c') {
						//printf( "i:%d", *ss);
						++ss;
					}
					else if (type == 'S') {
						//printf( "i:%u", *(uint16_t*)ss);
						ss += 2;
					}
					else if (type == 's') {
						//printf( "i:%d", *(int16_t*)ss);
						ss += 2;
					}
					else if (type == 'I') {
						//printf( "i:%u", *(uint32_t*)ss);
						ss += 4;
					}
					else if (type == 'i') {
						//printf( "i:%d", *(int32_t*)ss);
						ss += 4;
					}
					else if (type == 'f') {
						//printf( "f:%g", *(float*)ss);
						ss += 4;
					}
					else if (type == 'd') {
						//printf( "d:%lg", *(double*)ss);
						ss += 8;
					}
					else if (type == 'Z' || type == 'H') {
						//printf( "%c:", type);
						//while (*ss) printf( "%c", *ss++);
						while (*ss) ss++;
						++ss;
					}
			}


			if(foundXR && foundZS){ //for the nline(5) mapped alignments, as long as one alginment does not contain XR or ZS, found will be set to false.
				//continue;
			}else{
				found = false;
				return found;
			}
		}

		samclose(BAM_fp);
	}
	if(format == "SAM")
	{
		string line = "";
		ifstream input(file.c_str(), ios::in);
		while( input.good() && nline < 5) {
			getline(input, line);
			if( line == "" || find(line.begin(), line.begin()+1, '@') == line.begin() ){continue;}

			vector <string> fields;
			boost::split(fields, line, boost::is_any_of("\t"));

			if(fields[2] == "*")continue;

			nline ++;
			bool foundXR = (line.find("XR") != string::npos);
			bool foundZS = (line.find("ZS") != string::npos);
			if(foundXR && foundZS){ //for the nline(5) mapped alignments, as long as one alginment does not contain XR or ZS, found will be set to false.
				//continue;
			}else{
				found = false;
				return found;
			}
		}
		input.close();

	}

	cout << "XR and ZS fields are found in first " <<nline <<" mapped alignment" <<endl;
	return found;
}



int strandCombinedBed( string laneName, string ofile, string chromName, map<int, cMeth> & meth, ofstream & statsFile )
{
	cout << "Writing strand combined statistics for file " << laneName << " and chrom " << chromName << endl;
	if(meth.size() == 0) return 0;

	ofstream bBedFile;
	bBedFile.open((ofile + ".bed").c_str(), ios_base::out);
	bBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
	//if(opts.reference != "" && opts.reportCHX != 'X') bBedFile << "\tlocalSeq";
	if(opts.reference != "" ) bBedFile << "\tlocalSeq";
	bBedFile << endl;
	bBedFile << setprecision(3);

	ofstream bGBedFile;
	bGBedFile.open((ofile + ".G.bed").c_str(), ios_base::out);
	bGBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
	//if(opts.reference != "" && opts.reportCHX != 'X') bGBedFile << "\tlocalSeq";
	if(opts.reference != "" ) bGBedFile << "\tlocalSeq";
	bGBedFile << endl;
	bGBedFile << setprecision(3);

	ofstream bHGBedFile;
	bHGBedFile.open((ofile + ".HG.bed").c_str(), ios_base::out);
	bHGBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
	//if(opts.reference != "" && opts.reportCHX != 'X') bHGBedFile << "\tlocalSeq";
	if(opts.reference != "" ) bHGBedFile << "\tlocalSeq";
	bHGBedFile << endl;
	bHGBedFile << setprecision(3);

	vector<int> numc;
	string dna = "";
	if(opts.reference != ""){
		dnaByChr(opts.reference, chromName, dna, numc);
	}
	numcByChrom[chromName] = numc;

	map<char, stat>  statb;

	map<int, statByDepth> statbByDepth;

	char nexts[] ={'C', 'G', 'A', 'T'};
	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
		statb[ nexts[n] ] = stat(0, 0, 0, 0, 0);
	}

	for(map<int,cMeth>::iterator l = meth.begin(); l != meth.end(); l ++){
		cMeth loci = (*l).second;
		int pos = (*l).first;
		int tot = loci.totalC;
		int me = loci.methC;
		char st = loci.strand;
		char next = loci.nextN;

		char nn = 'Z';
		string localSeq;
		if(opts.reference != ""){
			nn = dna[pos+2];
			if(st == '-' && pos>1){nn = rctable(dna[pos-2]);}
			stringstream localSeqss;

			//localSeqss << dna[pos-2] << dna[pos-1] << 'C' << next << nn;
			//this is error because when it is not CG methylation, that happened asymmetric.
			char cpos1=(pos>1)?dna[pos-2]:'N';
			char cpos2=(pos>0)?dna[pos-1]:'N';
			localSeqss << cpos1 << cpos2 << dna[pos] << dna[pos+1] << dna[pos+2];

			localSeq = localSeqss.str();
		}

		double ratioBoth = 0.0;
		int totBoth = 0;
		int meBoth = 0;

		stringstream bBedFileToPrint;
		bBedFileToPrint << setprecision(3);

		if( st == '+' and next == 'G'){ 			//need combine this CG on + strand with the CG on - strand
			if(meth.count(pos + 1) > 0){			// there is coverage on C on both strands
				cMeth Minus = meth[pos + 1];
				int totMinus = Minus.totalC;
				int meMinus = Minus.methC;

				totBoth = tot + totMinus;
				meBoth = me + meMinus;
				ratioBoth = double (meBoth) / totBoth;
				statb[next].sites += 1;
				statb[next].totalC += totBoth;
				statb[next].methC += meBoth;
				statb[next].mean += ratioBoth;
				bBedFileToPrint << chromName <<"\t"<< pos <<"\t"<< pos + 2 <<"\t"<< ratioBoth <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< 'B' <<"\t"<< next <<"\t"<< '+' <<"\t"<< tot <<"\t"<< me <<"\t"<< '-' <<"\t"<< totMinus <<"\t"<< meMinus;
			} else {								//##coverag only on + strand, no coverage on the C on the - strand
				totBoth = tot;
				meBoth = me;
				ratioBoth = double(meBoth) / totBoth;
				statb[next].sites += 1;
				statb[next].totalC += totBoth;
				statb[next].methC += meBoth;
				statb[next].mean += ratioBoth;
				bBedFileToPrint << chromName <<"\t"<< pos <<"\t"<< pos + 2 <<"\t"<< ratioBoth <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< st <<"\t"<< next <<"\t"<< '+' <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< '-' <<"\t"<< '0' <<"\t"<< '0';
			}
		}
		else if( st == '-' and next == 'G'){//		##this CG is already combined with the CG on + strand if it exists
			if( meth.count(pos -1 ) > 0){
				continue;					//##if such key exists then this key has already been combined in previous plus strand result.
			}
			else{//except KeyError: 				//##only coverage on - strand, no coverage on the C on the + strand
				totBoth = tot;
				meBoth = me;
				ratioBoth = double(meBoth) / totBoth;
				statb[next].sites += 1;
				statb[next].totalC += totBoth;
				statb[next].methC += meBoth;
				statb[next].mean += ratioBoth;
				bBedFileToPrint << chromName <<"\t"<< pos - 1 <<"\t"<< pos + 1 <<"\t"<< ratioBoth <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< st <<"\t"<< next <<"\t"<< '+' <<"\t"<< '0' <<"\t"<< '0' <<"\t"<< '-' <<"\t"<< totBoth <<"\t"<< meBoth;
			}
		}
		else{//							##next != G : no need and not possible to expand to 2 bases for CCCCC like sequence
			totBoth = tot;
			meBoth = me;
			ratioBoth = double(meBoth) / totBoth;
			statb[next].sites += 1;
			statb[next].totalC += totBoth;
			statb[next].methC += meBoth;
			statb[next].mean += ratioBoth;

			if( st == '+' ){
				//fprintf(bBedFile, "%s\t%d\t%d\t%.3f\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n", chrom, pos, pos + 1, ratioBoth, totBoth, meBoth,'+',next, '+', totBoth, meBoth, '-', 0, 0);
				bBedFileToPrint << chromName <<"\t"<< pos <<"\t"<< pos + 1 <<"\t"<< ratioBoth <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< '+' <<"\t"<< next <<"\t"<< '+' <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< '-' <<"\t"<< '0'  <<"\t"<< '0';
			}
			else if( st == '-' ){
				//fprintf(bBedFile, "%s\t%d\t%d\t%.3f\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\n", chrom, pos, pos + 1, ratioBoth, totBoth, meBoth,'-',next, '+', 0, 0, '-', totBoth, meBoth);
				bBedFileToPrint << chromName <<"\t"<< pos <<"\t"<< pos + 1 <<"\t"<< ratioBoth <<"\t"<< totBoth <<"\t"<< meBoth <<"\t"<< '-' <<"\t"<< next <<"\t"<< '+' <<"\t"<< '0' <<"\t"<< '0' <<"\t"<< '-' <<"\t"<< totBoth <<"\t"<< meBoth;
			}
		}

		if(opts.reference != ""){
		//if(opts.reference != "" && opts.reportCHX != 'X'){
			bBedFileToPrint << "\t" << localSeq;
		}

		bBedFileToPrint << endl;

		if(opts.statsOnly){
			//nothing
		} else {
			if( opts.fullMode == 1 ) bBedFile << bBedFileToPrint.str();

			if( opts.reportCpX == next ){ 							//CpG meth
				if( opts.skipRandomChrom == 1 && chromName.find("random") != string::npos ){
					//continue;
				} else {
					bGBedFile << bBedFileToPrint.str();
				}
			}
			if( opts.reportCHX == nn && opts.reportCpX != next ){ 	//CHG meth
				if( opts.skipRandomChrom == 1 && chromName.find("random") != string::npos ){
					//continue;
				} else {
					bHGBedFile << bBedFileToPrint.str();
				}
			}
		}

		for(int d = 1; d <= DEPTHSAT; d++){
			if(totBoth >= d){
				statbByDepth[d].nc ++ ;
				statbByDepth[d].rc += ratioBoth;
				if(next == 'G'){
					statbByDepth[d].ncg ++ ;
					statbByDepth[d].rcg += ratioBoth;
				} else {
					statbByDepth[d].nch ++ ;
					statbByDepth[d].rch += ratioBoth;
					if(nn == 'G'){
						statbByDepth[d].nchg ++ ;
						statbByDepth[d].rchg += ratioBoth;
					} else {
						statbByDepth[d].nchh ++ ;
						statbByDepth[d].rchh += ratioBoth;
					}
				}
			}
		}
	}

	bBedFile.close();
	bGBedFile.close();
	bHGBedFile.close();
	for(int d = 1; d <= DEPTHSAT; d++){
		statbByDepth[d].rc /= statbByDepth[d].nc ;
		statbByDepth[d].rcg /= statbByDepth[d].ncg ;
		statbByDepth[d].rch /= statbByDepth[d].nch ;
		statbByDepth[d].rchg /= statbByDepth[d].nchg ;
		statbByDepth[d].rchh /= statbByDepth[d].nchh ;
	}
	statbByDepthByLane[laneName][chromName]=statbByDepth;
	//statbByDepthVec.push_back(statbByDepth);
	////START:This part is different than output in strandSpecificBed() for the CG dimers;
	statsFile << "next\tstrand\tsites\tmean\ttotalC\tmethC\tglobal\tdepth\n";
	//start to output statistics on strand specific bed file

	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
		char nt = nexts[n];
		double globalRatio = 0;
	    if( statb[ nt ].sites == 0 ){
			statb[nt].mean = 0;
			statb[nt].depth = 0;
	    }else{
			statb[nt].mean = statb[nt].mean / double( statb[nt].sites );
			statb[nt].depth = statb[nt].totalC / double( statb[nt].sites );
	    	globalRatio = double(statb[nt].methC)/statb[nt].totalC;
	    }
	    //statbVec.push_back(stat(statb[nt].sites, statb[nt].mean, statb[nt].methC, statb[nt].totalC, statb[nt].depth));
	    statsFile << nt <<"\t"<< 'B' <<"\t"<< statb[nt].sites <<"\t"<< statb[nt].mean <<"\t"<< statb[nt].totalC <<"\t"<< statb[nt].methC <<"\t"<< globalRatio <<"\t"<< statb[nt].depth << endl;
	}
	statbByLane[laneName][chromName]=statb;
	//statbVec.push_back(statb);

	double bisulfiteConversion = 0;
	if( (statb['C'].totalC + statb['A'].totalC + statb['T'].totalC) > 0 ){
		bisulfiteConversion = double(statb['C'].methC + statb['A'].methC + statb['T'].methC) / (statb['C'].totalC + statb['A'].totalC + statb['T'].totalC);
	}
	statsFile << "bisulfiteConversionFail: " << bisulfiteConversion << ", or bisulfite conversion ratio = " << 1 - bisulfiteConversion << endl;
	////END:This part is different than output in strandSpecificBed() for the CG dimers;

	statsFile << "depth\tNumC\tNumCG\tNumCH\tNumCHG\tNumCHH\tMeanRatioC\tMeanRatioCG\tMeanRatioCH\tMeanRatioCHG\tMeanRatioCHH" << endl;
	if(opts.reference != ""){
		statsFile << "0(genome)" <<"\t"<< numc[0] <<"\t"<< numc[1] <<"\t"<< numc[2] <<"\t"<< numc[3] <<"\t"<< numc[4] <<"\t"<< "NA\tNA\tNA\tNA\tNA" << endl;
	}
	for(int d = 1; d <= DEPTHSAT; d++){
		statsFile << d <<"\t"<< statbByDepth[ d ].nc
				<<"\t"<< statbByDepth[ d ].ncg <<"\t"<< statbByDepth[ d ].nch
				<<"\t"<< statbByDepth[ d ].nchg <<"\t"<< statbByDepth[ d ].nchh
				<<"\t"<< statbByDepth[ d ].rc
				<<"\t"<< statbByDepth[ d ].rcg <<"\t"<< statbByDepth[ d ].rch
				<<"\t"<< statbByDepth[ d ].rchg <<"\t"<< statbByDepth[ d ].rchh
				<< endl;
	}

	return 0;
};

void mergeStatbVec( string laneName, ofstream & statsFile )
{
	statsFile << endl;
	statsFile << "Strand combined stats where CGs on two strands are regarded as one Cytosine." << endl;
	statsFile << "next\tstrand\tsites\tmean\ttotalC\tmethC\tglobal\tdepth\n";

	map<char, stat>  x;

	//int chrs = statbByLane[laneName].size();
	char nexts[] ={'C', 'G', 'A', 'T'};
	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++)
	{
		long long xsites = 0;
		double xmean = 0.0;
		long long xtotalC = 0;
		long long xmethC = 0;
		double globalRatio = 0.0;
		double xdepth = 0.0;
		for(map < string, map<char, stat> >::iterator it = statbByLane[laneName].begin(); it != statbByLane[laneName].end(); it++)
		{
			string chromName = it->first;
			xsites 	+= statbByLane[laneName][chromName][ nexts[n] ].sites;
			xmean 	+= statbByLane[laneName][chromName][ nexts[n] ].sites * statbByLane[laneName][chromName][ nexts[n] ].mean;
			xtotalC += statbByLane[laneName][chromName][ nexts[n] ].totalC;
			xmethC 	+= statbByLane[laneName][chromName][ nexts[n] ].methC;
		}
		xmean = xmean / xsites;
		globalRatio = double(xmethC) / xtotalC;
		xdepth = double(xtotalC) / xsites;

		statsFile << nexts[n] <<"\t"<< 'B' <<"\t"<< xsites <<"\t"<< xmean <<"\t"<< xtotalC <<"\t"<< xmethC <<"\t"<< globalRatio <<"\t"<< xdepth << endl;

		x[nexts[n]] = stat(xsites, xmean, xmethC, xtotalC, xdepth);
	}

	double bisulfiteConversion = 0;
	if( (x['C'].totalC + x['A'].totalC + x['T'].totalC) > 0 ){
		bisulfiteConversion = double(x['C'].methC + x['A'].methC + x['T'].methC) / (x['C'].totalC + x['A'].totalC + x['T'].totalC);
	}
	statsFile << "bisulfiteConversionFail: " << bisulfiteConversion << ", or bisulfite conversion ratio = " << 1 - bisulfiteConversion << endl;

	//statsFile << endl << "All reads = " << propertiesByLane[file].nAllReads << "; Mapped reads = " << propertiesByLane[file].nMappedReads << endl << endl;

	statsFile << "depth\tNumC\tNumCG\tNumCH\tNumCHG\tNumCHH\tMeanRatioC\tMeanRatioCG\tMeanRatioCH\tMeanRatioCHG\tMeanRatioCHH" << endl;
	long long nc = 0; long long ncg = 0; long long nch = 0; long long nchg = 0; long long nchh = 0;
	if(opts.reference != ""){
		for(map<string, vector<int> >::iterator it = numcByChrom.begin(); it != numcByChrom.end(); it++){
			vector <int> numc = it->second;
			nc += numc[0];
			ncg += numc[1];
			nch += numc[2];
			nchg += numc[3];
			nchh += numc[4];
		}
		statsFile << "0(genome)" <<"\t"<< nc <<"\t"<< ncg <<"\t"<< nch <<"\t"<< nchg <<"\t"<< nchh <<"\t"<< "NA\tNA\tNA\tNA\tNA" << endl;
	}
	for(int d = 1; d <= DEPTHSAT; d++){
		nc = 0; ncg = 0; nch = 0; nchg = 0; nchh = 0;
		double rc = 0; double rcg = 0; double rch = 0; double rchg = 0; double rchh = 0;
		//map <string, map < string, map<int, statByDepth> > > statbByDepthByLane; //statbByLane[laneName][chromName][depth] = statByDepthObj
		for(map < string, map<int, statByDepth> >::iterator it = statbByDepthByLane[laneName].begin(); it != statbByDepthByLane[laneName].end(); it++)
		{
			string chromName = it->first;

			nc 		+= statbByDepthByLane[laneName][chromName][ d ].nc;
			ncg 	+= statbByDepthByLane[laneName][chromName][ d ].ncg;
			nch 	+= statbByDepthByLane[laneName][chromName][ d ].nch;
			nchg 	+= statbByDepthByLane[laneName][chromName][ d ].nchg;
			nchh 	+= statbByDepthByLane[laneName][chromName][ d ].nchh;

			if( statbByDepthByLane[laneName][chromName][ d ].nc   > 0 ) rc 		+= statbByDepthByLane[laneName][chromName][ d ].nc 	 * statbByDepthByLane[laneName][chromName][ d ].rc;
			if( statbByDepthByLane[laneName][chromName][ d ].ncg  > 0 ) rcg 	+= statbByDepthByLane[laneName][chromName][ d ].ncg	 * statbByDepthByLane[laneName][chromName][ d ].rcg;
			if( statbByDepthByLane[laneName][chromName][ d ].nch  > 0 ) rch 	+= statbByDepthByLane[laneName][chromName][ d ].nch  * statbByDepthByLane[laneName][chromName][ d ].rch;
			if( statbByDepthByLane[laneName][chromName][ d ].nchg > 0 ) rchg	+= statbByDepthByLane[laneName][chromName][ d ].nchg * statbByDepthByLane[laneName][chromName][ d ].rchg;
			if( statbByDepthByLane[laneName][chromName][ d ].nchh > 0 ) rchh 	+= statbByDepthByLane[laneName][chromName][ d ].nchh * statbByDepthByLane[laneName][chromName][ d ].rchh;
		}

		rc 	 = rc 	/ nc;
		rcg  = rcg 	/ ncg;
		rch  = rch 	/ nch;
		rchg = rchg / nchg;
		rchh = rchh / nchh;
		statsFile << d <<"\t"<< nc <<"\t"<< ncg <<"\t"<< nch <<"\t"<< nchg <<"\t"<< nchh <<"\t"<< rc <<"\t"<< rcg <<"\t"<< rch <<"\t"<< rchg <<"\t"<< rchh << endl;
	}

}

int strandSpecificBed(string laneName, string ofile, string chromName, map<int, cMeth> & meth, ofstream & statsFile )
{
	cout << "Writing strand specific statistics for file " << laneName << " and chrom " << chromName << endl;
	if(meth.size() == 0) return 0;

	map<char, map<char, stat> > stats;

	char strands[] = { '+','-','B' };
	char nexts[] ={'C', 'G', 'A', 'T'};
	for(unsigned int s = 0; s < sizeof(strands)/sizeof(char); s++){
		for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
			stats[ strands[s] ] [ nexts[n] ] = stat(0, 0, 0, 0, 0);
		}
	}

	ofstream strandBedFile;
	strandBedFile.open((ofile + "_strand.bed").c_str(), ios_base::out);
	strandBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext" <<endl;
	strandBedFile << setprecision(3);


	for(map<int,cMeth>::iterator l = meth.begin(); l != meth.end(); l ++){
		int start = l->first;
		cMeth loci = (*l).second;

		double ratio = double(loci.methC)/loci.totalC;
		stats[loci.strand][loci.nextN].sites += 1;
		stats[loci.strand][loci.nextN].totalC += loci.totalC;
		stats[loci.strand][loci.nextN].methC += loci.methC;
		stats[loci.strand][loci.nextN].mean += ratio;
		//fprintf(strandBedFile, "%s\t%d\t%d\t%.3f\t%d\t%d\t%s\t%s\n",chrom, (*l).first, (*l).first + 1, ratio, loci.totalC, loci.methC,loci.strand, loci.nextN);
		if((opts.fullMode == 1) && (opts.statsOnly == 0))
			strandBedFile << chromName <<"\t"<< start <<"\t"<< start + 1 <<"\t"<< ratio <<"\t"<< loci.totalC <<"\t"<< loci.methC <<"\t"<< loci.strand <<"\t"<< loci.nextN << endl;

	}
	strandBedFile.close();
	//finished output of _strand.bed file

	//start to output statistics on strand specific bed file
	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
		stats['B'][ nexts[n] ].sites = stats['+'][ nexts[n] ].sites + stats['-'][ nexts[n] ].sites ;
		stats['B'][ nexts[n] ].totalC = stats['+'][ nexts[n] ].totalC + stats['-'][ nexts[n] ].totalC ;
		stats['B'][ nexts[n] ].methC = stats['+'][ nexts[n] ].methC + stats['-'][ nexts[n] ].methC ;
		stats['B'][ nexts[n] ].mean = stats['+'][ nexts[n] ].mean + stats['-'][ nexts[n] ].mean ;
	}
	statsFile << "next\tstrand\tsites\tmean\ttotalC\tmethC\tglobal\tdepth\n";
	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
		for(unsigned int s = 0; s < sizeof(strands)/sizeof(char); s++){
			char next = nexts[n];
			char strand = strands[s];
			double globalRatio = 0;
			if(stats[strand][next].sites == 0 ){
				stats[strand][next].mean = 0;
				stats[strand][next].depth = 0;
			} else {
				stats[strand][next].mean = stats[strand][next].mean / double( stats[strand][next].sites );
				stats[strand][next].depth = stats[strand][next].totalC / double( stats[strand][next].sites );
				globalRatio = stats[strand][next].methC / double( stats[strand][next].totalC );
			}
			//statsVec.push_back(stat(stats[strand][next].sites, stats[strand][next].mean, stats[strand][next].methC, stats[strand][next].totalC, stats[strand][next].depth));
			statsFile << next <<"\t"<< strand <<"\t"<< stats[strand][next].sites <<"\t"<< stats[strand][next].mean <<"\t"<< stats[strand][next].totalC <<"\t"<< stats[strand][next].methC <<"\t"<< globalRatio <<"\t"<< stats[strand][next].depth << endl;
		}
	}
	statsByLane[laneName][chromName] = stats;
	//statsVec.push_back(stats);

	double bisulfiteConversion = 0;
	if( stats['B']['C'].totalC + stats['B']['A'].totalC + stats['B']['T'].totalC > 0 ){
		bisulfiteConversion = double(stats['B']['C'].methC + stats['B']['A'].methC + stats['B']['T'].methC) / (stats['B']['C'].totalC + stats['B']['A'].totalC + stats['B']['T'].totalC);
	}
	//printf( "%s\t%.3f\n", "bisulfiteConversionFail: ", bisulfiteConversion);
	statsFile << "bisulfiteConversionFail: " << bisulfiteConversion << ", or bisulfite conversion ratio = " << 1 - bisulfiteConversion << endl;

	//finished generating the strand specific bed stats
	return 0;
};

void mergeStatsVec( string laneName, ofstream & statsFile )
{
	statsFile << "Strand specific stats where CGs on two strands are regarded as two Cytosines." << endl;
	statsFile << "next\tstrand\tsites\tmean\ttotalC\tmethC\tglobal\tdepth\n";

	char nexts[] ={'C', 'G', 'A', 'T'};
	char strands[] = { '+','-','B' };
	//int chrs = statsVec.size();
	map<char, map<char, stat> >  x; //merged stats table

	for(unsigned int n = 0; n < sizeof(nexts)/sizeof(char); n++){
		for(unsigned int s = 0; s < sizeof(strands)/sizeof(char); s++){
			char next = nexts[n];
			char strand = strands[s];
			long long xsites = 0;
			double xmean = 0.0;
			long long xtotalC = 0;
			long long xmethC = 0;
			double globalRatio = 0.0;
			double xdepth = 0.0;

			//map <string, map < string, map<char, map<char, stat> > > > statsByLane; //statsByLane[laneName][chromName][strand][next]=statObj
			for(map < string, map< char, map<char, stat> > >::iterator it = statsByLane[laneName].begin(); it != statsByLane[laneName].end(); it++)
			{
				string chromName = it->first;
				xsites 	+= statsByLane[laneName][chromName][strand][ next ].sites;
				xmean 	+= statsByLane[laneName][chromName][strand][ next ].sites * statsByLane[laneName][chromName][strand][ next ].mean;
				xtotalC += statsByLane[laneName][chromName][strand][ next ].totalC;
				xmethC 	+= statsByLane[laneName][chromName][strand][ next ].methC;
			}
			xmean = xmean / xsites;
			globalRatio = double(xmethC) / xtotalC;
			xdepth = double(xtotalC) / xsites;

			statsFile << next <<"\t"<< strand <<"\t"<< xsites <<"\t"<< xmean <<"\t"<< xtotalC <<"\t"<< xmethC <<"\t"<< globalRatio <<"\t"<< xdepth << endl;

			x[strand][next] = stat(xsites, xmean, xmethC, xtotalC, xdepth);

		}
	}
	double bisulfiteConversion = 0;
	if( x['B']['C'].totalC + x['B']['A'].totalC + x['B']['T'].totalC > 0 ){
		bisulfiteConversion = double(x['B']['C'].methC + x['B']['A'].methC + x['B']['T'].methC) / (x['B']['C'].totalC + x['B']['A'].totalC + x['B']['T'].totalC);
	}
	statsFile << "bisulfiteConversionFail: " << bisulfiteConversion << ", or bisulfite conversion ratio = " << 1 - bisulfiteConversion << endl;
}


//The ZS field is equivalent to the following way
//if( ((b->core.flag & 0x10) == 0) && ((b->core.flag & 0x40) != 0) ){
//	cout << "++" << "..." << zs << endl;
//}
//if( ((b->core.flag & 0x10) != 0) && ((b->core.flag & 0x40) == 0) ){
//	cout << "+-" << "..." << zs << endl;
//}
//if( ((b->core.flag & 0x10) != 0) && ((b->core.flag & 0x40) != 0) ){
//	cout << "-+" << "..." << zs << endl;
//}
//if( ((b->core.flag & 0x10) == 0) && ((b->core.flag & 0x40) == 0) ){
//	cout << "--" << "..." << zs << endl;
//}

void createHashFromRead(alignedRead & read,  map<int, cMeth> & meth, ofstream & skipBpFile )
{
	//start creating hash of cpgs found from current line of read
	char match = 'C';
	char convert = 'T';
	if(read.strand == '-'){
		match = 'G';
		convert = 'A';
	}
	
	//TODO:make this an option later
	//do removal of end repair sequence. 
	//0: no trim; -1: trim from beginning to 1st methylated C in the read in process; -2: identify the position of methylated C from ++/-+ reads; n>=0: trim n bases from beginning;
	unsigned int minIndex = 0;
	unsigned int maxIndex = read.ref.size();
	if(opts.trimWGBSEndRepairPE2Seq == -1)
	{
		if( read.strandCopy == '-' )
		{
			if( read.strand == '+' ){	//trim from last C to the end in the sequence
				maxIndex = 0;		//if 'C' is not found, then this read is disregarded
				char methy[] = {'C'};
				string::iterator it = find_end(read.seq.begin(),read.seq.end(), methy, methy+1);
				if (it!=read.seq.end()) maxIndex = int( it - read.seq.begin() );
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< maxIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim+-FromMaxIndex" <<"\t"<< "withTrimSize" <<"\t"<< read.ref.size() - maxIndex << endl;
			}
			if( read.strand == '-' ){	//trim from beginning to the first G in the sequence
				minIndex = read.ref.size(); 	//if 'G' is not found, then this read is disregarded
				string::iterator it = find(read.seq.begin(), read.seq.end(), 'G');
				if (it!=read.seq.end()) minIndex = int( it - read.seq.begin() );
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< minIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim--BeforeMinIndex" <<"\t"<< "withTrimSize" <<"\t"<< minIndex << endl;
			}
		}
	}
	if(opts.trimWGBSEndRepairPE2Seq > 0){
		if( read.strandCopy == '-' )
		{
			if( read.strand == '+' ){
				maxIndex = read.seq.size() - opts.trimWGBSEndRepairPE2Seq;
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< maxIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim+-FromMaxIndex" <<"\t"<< "withTrimSize" <<"\t"<< read.ref.size() - maxIndex <<"\t"<< read.name << endl;
			}
			if( read.strand == '-' ){
				minIndex = opts.trimWGBSEndRepairPE2Seq;
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< minIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim--BeforeMinIndex" <<"\t"<< "withTrimSize" <<"\t"<< minIndex <<"\t"<< read.name << endl;
			}
		}
	}
	//size_t index = read.ref.find(match); //note that string::find is problematic, it always gives an index out of range.
	string::iterator p = find(read.ref.begin(), read.ref.end(), match);
	while( p != read.ref.end() ){

		unsigned int index = p - read.ref.begin();
		if( index < minIndex || index > maxIndex ) {	//do removal of end repair sequence
			p = find( p + 1, read.ref.end(),match);
			continue;
		}

		int loc = read.pos + index;

		char nextN = 'Z';
		if(read.strand == '+'){						//find nextN; //code could be simpler
			if(index + 1 < read.ref.size()){
				nextN = read.ref[index + 1];
			}else{
				nextN = 'X';
			}
		}
		else if(read.strand == '-'){
			if(index > 0){
				nextN = read.ref[index - 1];
			}else{
				nextN = 'X';
			}
			//cout <<"HHH" << nextN << ".."<< read.pos << "\t"<< index <<"\t"<<nextN<<"\t"<< read.seq<<"\t"<< read.ref <<"\t"<< read.strand<<"\tHHH\n";
			nextN = rctable(nextN);
		}

		//if C falls on the last base of a read. This looks like a new option but it's part of opts.nextBaseMinScore.
		//if one does not care whether nextbase matches, then opts.nextBaseMinScore == -1. Then we do not need check what's nextN.
		//if one does     care whether nextbase matches, then opts.nextBaseMinScore  > -1. Then we do     need check what's nextN.
		//This part is separated with the situation when C is not the last base of a read for easier reading of code.
		if( opts.nextBaseMinScore > -1 ){
			if( nextN == 'X' || nextN == 'Y' ){			//skip when C falls on the last base.
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< index <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\tLast\n";
				//fprintf (skipBpFile, "%s\t%d\t%d\t%c\t%s\t%s\t%s\t%s\n",read.chrom.c_str(),read.pos, int(index), read.strand, read.seq.c_str(), read.ref.c_str(), read.qual.c_str(), "Last");
				p = find( p + 1, read.ref.end(),match);
				//index = read.ref.find(match,index+1);
				continue;
			}
		}

		if( opts.cytosineMinScore ){					//skip when C does not satisfies the score threshold
			if(read.qual[index] - opts.qualityScoreBase < opts.cytosineMinScore ){
				//fprintf(skipBpFile, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\n",read.chrom.c_str(),read.pos, int(index), read.strand, read.seq.c_str(), read.ref.c_str(), read.qual.c_str(), "CMinScore", read.qual[index]-opts.qualityScoreBase);
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< index <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "CMinScore" <<"\t"<< read.qual[index] - opts.qualityScoreBase << endl;
				p = find( p + 1, read.ref.end(),match);
				//index =read.ref.find(match, index + 1);
				continue;
			}
		}

		if( opts.nextBaseMinScore > -1 ){ 			// -1 -> does not check whether next base matches reference or whether the score exceeds minimum.
			if(read.strand == '+'){
				if( (read.seq[index+1] != bsPosTable(read.ref[index+1]) && read.seq[index+1] != read.ref[index+1] ) || read.qual[index+1]-opts.qualityScoreBase <opts.nextBaseMinScore){
					//fprintf(skipBpFile, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n",read.chrom,read.pos, index, read.strand, read.seq, read.ref, read.qual, "NextBaseMinScore", read.seq[index+1], read.ref[index+1], read.qual[index+1]-opts.qualityScoreBase);
					if(opts.fullMode) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< index <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "NextBaseMinScore" <<"\t"<< read.seq[index+1] <<"\t"<< read.ref[index+1] <<"\t"<< read.qual[index+1]-opts.qualityScoreBase << endl;
					p = find( p + 1, read.ref.end(),match);
					//index = read.ref.find(match, index + 1);
					continue;
				}
			}
			else if (read.strand == '-'){
				if( (read.seq[index  - 1] != bsNegTable(read.ref[index - 1]) && read.seq[index-1] != read.ref[index-1] ) || read.qual[index-1]-opts.qualityScoreBase <opts.nextBaseMinScore){
					//fprintf(skipBpFile, "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\n",read.chrom,read.pos, index, read.strand, read.seq, read.ref, read.qual, "NextBaseMinScore", read.seq[index-1], read.ref[index-1], read.qual[index-1]-opts.qualityScoreBase);
					if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< index <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "NextBaseMinScore" <<"\t"<< read.seq[index-1] <<"\t"<< read.ref[index-1] <<"\t"<< read.qual[index-1]-opts.qualityScoreBase << endl;
					p = find( p + 1, read.ref.end(),match);
					//index = read.ref.find(match, index + 1);
					continue;
				}
			}
		}

		if(1){ 										//process this cytosine on this read
			if(read.seq[index] == convert){ //unmethylated
				//if(chri[read.chrom] == 0){cout << chri[read.chrom] << endl;cout << read.chrom << endl; exit(1);}
				if( meth.count(loc) == 0){
					meth[loc] = cMeth(1, 0, read.strand, nextN);
				}else {
					meth[loc] = cMeth(meth[loc].totalC + 1, meth[loc].methC, read.strand, nextN);
				}
			}
			else if( read.seq[index] == match){ //methylated
				if( meth.count(loc) == 0){
					meth[loc] = cMeth(1, 1, read.strand, nextN);
				}else {
					meth[loc] = cMeth(meth[loc].totalC + 1, meth[loc].methC + 1, read.strand, nextN);
				}
			} else { 							//SNP
				if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< index <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "SNP" <<"\t"<< read.qual[index] - opts.qualityScoreBase << endl;
			}
		}

		//if(index + 1 >= read.ref.size() ){break;} //string::find function does not check this:(
		//index = read.ref.find(match, index+1);
		p = find( p + 1, read.ref.end(),match);
	}
	//finished creating hash of cpgs found from one line

};


int createHashFromFile(string file, string chromName, map<int, cMeth> & meth,  string & dna)
{
	cout <<"Start processing file "<< file << " on chrom " << chromName << endl;

	string format = file.substr(file.find_last_of(".") + 1);
	boost::to_upper(format);

	//read file by line and assign line to read; then process read object to create hash;
	//if use multithreads on one file, then skipBpFile need be locked and unlocked frequently and meth hash need be merged from multithreads.
	//I did not test but I suspect it's not going to be much faster.
	long nline = 0;  // number of all reads

	alignedRead read;
	//map<int, cMeth> meth ; //loc->(totalC, methC, strand, nextN);

	string ofile = file + "."+ chromName;
	ofstream skipBpFile;
	skipBpFile.open((ofile + "_skip.bed").c_str(), ios_base::out);
	skipBpFile << "#chrom\tstart\tend\tratio\ttotalC\tmeanC\tstrand\tnext" <<endl;


	if(format == "BAM")
	{
		samfile_t * BAM_fp;
		bam1_t *b;
		BAM_fp = samopen(file.c_str(), "rb", 0);
		if (BAM_fp == 0) {
			fprintf(stderr, "Fail to open BAM file %s\n", file.c_str());
			exit(1);
		}

		b = bam_init1();

		int found = 0;
		while(samread(BAM_fp,b) > 0 )
		{

			if( chromsByLane[file][b->core.tid] != chromName ) {
				if(found) break; //for sorted bam, otherwise it's a bug;
				continue;
			} else {
				found = 1;
			}

			read.name = string((char*)bam1_qname(b));

			if( abs(b->core.isize) < opts.minFragSize ){
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< b->core.isize <<"\t"<<  "minFragSize" << endl;
				continue;
			}
			if( abs(b->core.isize) < opts.minMMFragSize && (b->core.flag & 0x100) ){
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< b->core.isize <<"\t"<< "minMMFragSize" << endl;
				continue;
			}
			if( (opts.requiredFlag == 0) || ( (b->core.flag & opts.requiredFlag) == opts.requiredFlag ) ){
				//no requirement on flag or this read has this requiredFlag
			} else {
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< (b->core.flag & opts.requiredFlag) <<"\t"<< "requiredFlag" << endl;
				continue;
			}
			if( (opts.excludedFlag > 0) && ( (b->core.flag & opts.excludedFlag) == opts.excludedFlag ) ){
				//this read has this excludedFlag
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< (b->core.flag & opts.excludedFlag) <<"\t"<< "excludedFlag" << endl;
				continue;
			}



			read.chrom = BAM_fp->header->target_name[b->core.tid];

			read.pos = b->core.pos;

			//retrieving seq and qual.
			size_t l_seq=b->core.l_qseq;
			read.seq.assign(l_seq,0);
			read.qual.assign(l_seq,0);
			char * s = (char*) bam1_seq(b);
			char * t = (char*) bam1_qual(b);
			for(size_t i=0;i<l_seq;i++){
				read.seq[i]=bam_nt16_rev_table[bam1_seqi(s,i)];
				read.qual[i]=t[i]+33;
			}


			string zs = string( (char *) bam_aux2Z(bam_aux_get(b, "ZS") ) );
			//read.ref = string( (char *) bam_aux2Z(bam_aux_get(b, "XR") ) );

			if(b->core.n_cigar >3) {continue;} //currently I do not allow multiple I/Ds. All allowed patterns are like: 50M, 10M2I38M, 10M2D40M.
			//seq and qual may be cut or inserted to according to cigar operation
			string newSeq = "";
			string newQual = "";
			int refLen = 0;
			int spos = 0;
			int okToProceed = 1;
			uint32_t* cigar= bam1_cigar(b);
			for( int k=0;k< b->core.n_cigar;++k)
			{
				int cop =cigar[k] & BAM_CIGAR_MASK; // operation
				int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
				switch(cop)
				{
					case BAM_CMATCH:
						//printf("M");
						newSeq  += read.seq.substr(spos, cl );
						newQual += read.qual.substr(spos, cl );
						spos += cl;
						refLen += cl;
						break;
					case BAM_CINS:
						//printf("I");
						//for(int itemp = 0; itemp < cl; itemp++){refStr += 'X';}
						spos += cl;
						//refLen += 0;
						break;
					case BAM_CDEL:
						//printf("D");
						for(int itemp = 0; itemp < cl; itemp++){newSeq += 'N'; newQual += 'B';}
						//spos += 0;
						refLen += cl;
						break;
					case BAM_CREF_SKIP:
						okToProceed = 0;
						//printf("N");
						break;
					case BAM_CSOFT_CLIP:
						okToProceed = 0;
						//printf("S");
						break;
					case BAM_CHARD_CLIP:
						okToProceed = 0;
						//printf("R");
						break;
					case BAM_CPAD:
						okToProceed = 0;
						//printf("P");
						break;
					default:
						okToProceed = 0;
						//printf("?");
						break;
				}
			}

			if(! okToProceed ){continue;} //unexpected CIGAR
			read.seq  = newSeq;
			read.qual = newQual;
			read.ref  = dna.substr(read.pos, refLen);
			//cout << refStr << endl;
			//exit(0);
			read.strand = zs[0];
			read.strandCopy = zs[1];


			//will drop off the dependence of ZS and XR;
//			if(read.ref[0] > 96) { // that means there are two lower case letters // in new format
//				//cout << read.ref << endl;
//				read.ref = read.ref.substr(2, read.seq.size() );
//				//cout << read.ref << endl << endl;
//			}
			boost::to_upper(read.ref);
			//cout << read.ref << endl;
			//exit(0);
			//check if PE reads have overlap
			if(opts.processPEOverlapSeq){
				// isize >0 means it's properly paired and it's the leftmost end and this number means the fragsize. 0 means not properly paired, though sometimes it's mapped. negative number means the right end of the frag.
				if( b->core.isize > 0 && b->core.mpos - read.pos < l_seq) //current read is left end && there's overlap between the two ends.
				{
				//this does not handle the ==0 situation
				//if( b->core.mpos - read.pos > 0 && b->core.mpos - read.pos < l_seq ) //current read is left end && there's overlap between the two ends.
					read.seq  = read.seq.substr (0, b->core.mpos - read.pos);
					read.ref  = read.ref.substr (0, b->core.mpos - read.pos);
					read.qual = read.qual.substr(0, b->core.mpos - read.pos);
					l_seq = b->core.mpos - read.pos;
					if(opts.fullMode) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< l_seq <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "OverlapPEReads" << endl;
				}
			}

			//trim the end of PE1 read if maxReadLen is larger than fragment // example fragment = 95bp for 100mer sequencing, then trim read from end;
			if(opts.trimWGBSEndRepairPE1Seq && abs(b->core.isize) < propertiesByLane[file].readLen)
			{
				int bpsToTrim = 0;
				if( (b->core.flag & 0x2) ){ //if this read is properly paired, just trim opts.trimWGBSEndRepairPE1Seq bases from end of read
					bpsToTrim = opts.trimWGBSEndRepairPE1Seq;
				}
				if( (b->core.flag & 0x2) == 0 ){ // if this read is not properly paired, trim opts.trimWGBSEndRepairPE1Seq + NM bases from end of read
					size_t nm = bam_aux2i(bam_aux_get(b,"NM"));
					bpsToTrim = opts.trimWGBSEndRepairPE1Seq + nm;
				}
				if(opts.fullMode) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "trimWGBSEndRepairPE1Seq" << endl;

				if(l_seq < bpsToTrim) continue;
				if( read.strandCopy == '+' )
				{
					if( read.strand == '+' ){
						read.seq  = read.seq.substr (0, l_seq - bpsToTrim);
						read.ref  = read.ref.substr (0, l_seq - bpsToTrim);
						read.qual = read.qual.substr(0, l_seq - bpsToTrim);
						//if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< maxIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim+-FromMaxIndex" <<"\t"<< "withTrimSize" <<"\t"<< read.ref.size() - maxIndex << endl;
					}
					if( read.strand == '-' ){
						read.seq  = read.seq.substr (bpsToTrim, l_seq - bpsToTrim);
						read.ref  = read.ref.substr (bpsToTrim, l_seq - bpsToTrim);
						read.qual = read.qual.substr(bpsToTrim, l_seq - bpsToTrim);
						read.pos  = read.pos + bpsToTrim; //todo : put this trimWGBSEndRepairPE1Seq stuff in createHashFromRead()
						//if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< minIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim--BeforeMinIndex" <<"\t"<< "withTrimSize" <<"\t"<< minIndex << endl;
					}
				}
			}
			//cout << read.ref << endl;
			//exit(0);
			//finished retrieving fields from BAM file for current line
			createHashFromRead(read, meth, skipBpFile);
		}

		//finished creating hash of cpgs found from all lines
		samclose(BAM_fp);

	}	//end of BAM format



	if(format == "SAM")
	{
	    string line = "";
	    ifstream inputf(file.c_str(), ios::in);
	    while (inputf.good()) {
	    	getline(inputf, line);
	        if( line == "" || find(line.begin(), line.begin()+1, '@') == line.begin() ){continue;}

			//nline ++;
			//if(nline % 1000000 == 0){cout << "reading " << nline << " lines in file "<<file<<endl;}

	    	vector <string> fields;
	    	boost::split(fields, line, boost::is_any_of("\t"));

	    	read.chrom = fields[2];
			//if(read.chrom == "*")continue;
			if( read.chrom != chromName ) continue;


			int isize = string_to_int(fields[8]);
			unsigned int flag = string_to_int(fields[1]);
			if( abs(isize) < opts.minFragSize ){
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< isize <<"\t"<<  "minFragSize" << endl;
				continue;
			}
			if( abs(isize) < opts.minMMFragSize && (flag & 0x100) ){
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< isize <<"\t"<< "minMMFragSize" << endl;
				continue;
			}
			if( (opts.requiredFlag == 0) || ( (flag & opts.requiredFlag) == opts.requiredFlag ) ){
				//no requirement on flag or this read has this requiredFlag
			} else {
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< (flag & opts.requiredFlag) <<"\t"<< "requiredFlag" << endl;
				continue;
			}
			if( (opts.excludedFlag > 0) && ( (flag & opts.excludedFlag) == opts.excludedFlag ) ){
				//this read has this excludedFlag
				if(opts.fullMode) skipBpFile << read.name <<"\t"<< (flag & opts.excludedFlag) <<"\t"<< "excludedFlag" << endl;
				continue;
			}




			read.name = fields[0];

			read.pos = string_to_int(fields[3]) - 1;

			read.seq = fields[9];
			read.qual = fields[10];


			size_t strandPos = line.find("ZS:Z:");
			if(strandPos != string::npos){
				read.strand = line[strandPos+5];
				read.strandCopy = line[strandPos+6];
			}
			else{
				cerr << line << endl;
				cerr << strandPos << endl;
				cerr << "ZS:Z:<str> field not found; Is it BAM/SAM format?\n";
				exit(1);
			}

			size_t refPos = line.find("R:Z:");
			if(refPos != string::npos){
				unsigned int startAt = 4;
				if(line[refPos + startAt] > 96){
					//cout << line[refPos + startAt] << "==" << line[refPos + startAt] -96 << endl;
					startAt = 6;
				}
				string ref = "";
				for(unsigned int i = startAt; i < line.size(); i++ ){
					//cout << i << ",,,,";
					if(line[refPos+i] == 'A' || line[refPos+i] == 'C' ||line[refPos+i] == 'G' ||line[refPos+i] == 'T' ){
						ref += line[refPos+i] ;
					}
					else{		//end of ref field;
						break;
					}
				}
				read.ref = ref;
			}
			else{
				cerr << "R:Z:<str> field not found; Please enable -R option in BSMAP!\n";
				exit(1);
			}
			//finished retrieving fields from BAM file for current line
			//cout << line << endl;
			//cout << read.ref<< endl;

			int l_seq = read.seq.size();
			//check if PE reads have overlap
			if(opts.processPEOverlapSeq){
				int mpos = string_to_int(fields[7]) - 1;
				// isize >0 means it's properly paired and it's the leftmost end and this number means the fragsize. 0 means not properly paired, though sometimes it's mapped. negative number means the right end of the frag.
				if( isize > 0 && mpos - read.pos < l_seq) //current read is left end && there's overlap between the two ends.
				{
				//this does not handle the ==0 situation
				//if( b->core.mpos - read.pos > 0 && b->core.mpos - read.pos < l_seq ) //current read is left end && there's overlap between the two ends.
					read.seq  = read.seq.substr (0, mpos - read.pos);
					read.ref  = read.ref.substr (0, mpos - read.pos);
					read.qual = read.qual.substr(0, mpos - read.pos);
					l_seq = mpos - read.pos;
					if(opts.fullMode) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< l_seq <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "OverlapPEReads" << endl;
				}
			}

			//trim the end of PE1 read if maxReadLen is larger than fragment // example fragment = 95bp for 100mer sequencing, then trim read from end;
			if(opts.trimWGBSEndRepairPE1Seq && abs(isize) < propertiesByLane[file].readLen)
			{
				int bpsToTrim = 0;
				if( (flag & 0x2) ){ //if this read is properly paired, just trim opts.trimWGBSEndRepairPE1Seq bases from end of read
					bpsToTrim = opts.trimWGBSEndRepairPE1Seq;
				}
				if( (flag & 0x2) == 0 ){ // if this read is not properly paired, trim opts.trimWGBSEndRepairPE1Seq + NM bases from end of read
					string nmstr = "";
					size_t refPos = line.find("NM:i:");
					if(refPos != string::npos){
						for(unsigned int i = 5; i < line.size(); i++){
							if(line[refPos+i] >= 48 && line[refPos+i] <= 57 ){
								nmstr += line[refPos+i] ;
							} else {
								break;
							}
						}
					}
					int nm = string_to_int(nmstr);
					//cout << nm << endl;
					bpsToTrim = opts.trimWGBSEndRepairPE1Seq + nm;
				}
				if(opts.fullMode) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< read.strand <<"\t"<< read.name <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< read.qual <<"\t"<< "trimWGBSEndRepairPE1Seq" << endl;

				if(l_seq < bpsToTrim) continue;
				if( read.strandCopy == '+' )
				{
					if( read.strand == '+' ){
						read.seq  = read.seq.substr (0, l_seq - bpsToTrim);
						read.ref  = read.ref.substr (0, l_seq - bpsToTrim);
						read.qual = read.qual.substr(0, l_seq - bpsToTrim);
						//if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< maxIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim+-FromMaxIndex" <<"\t"<< "withTrimSize" <<"\t"<< read.ref.size() - maxIndex << endl;
					}
					if( read.strand == '-' ){
						read.seq  = read.seq.substr (bpsToTrim, l_seq - bpsToTrim);
						read.ref  = read.ref.substr (bpsToTrim, l_seq - bpsToTrim);
						read.qual = read.qual.substr(bpsToTrim, l_seq - bpsToTrim);
						read.pos  = read.pos + bpsToTrim; //todo : put this trimWGBSEndRepairPE1Seq stuff in createHashFromRead()
						//if(opts.fullMode == 1) skipBpFile << read.chrom <<"\t"<< read.pos <<"\t"<< minIndex <<"\t"<< read.strand <<"\t"<< read.seq <<"\t"<< read.ref <<"\t"<< "Trim--BeforeMinIndex" <<"\t"<< "withTrimSize" <<"\t"<< minIndex << endl;
					}
				}
			}
//			cout << read.seq << endl;
//			cout << read.ref << endl;
//			cout << read.qual << endl;

			createHashFromRead(read, meth, skipBpFile);

	    }
		//finished creating hash of cpgs found from all lines
	    if(nline == 0){
	    	cerr << "No read is found in file " << file <<endl;
	    }
		inputf.close();

	}	//end of SAM foramt


	skipBpFile.close();

	ofstream statsFile;
	statsFile.open((ofile + "_stat.txt").c_str(), ios_base::out);
	strandSpecificBed(file, ofile, chromName, meth, statsFile);
	strandCombinedBed(file, ofile, chromName, meth, statsFile);
	statsFile.close();
	// need remove the blank files
	if( opts.fullMode == 0 ){
		system( ("rm -rf " + ofile + "_strand.bed").c_str() );
		system( ("rm -rf " + ofile + "_skip.bed").c_str() );
		system( ("rm -rf " + ofile + ".bed").c_str() );
	}



	cout <<"Finished processing file "<< file << " on chrom " << chromName << endl;
	return 0;
};



void performCallChrom(vector <string> files, string chromName)
{

	map<int, map<int, cMeth> > meth ; //fileIndex->loc->(totalC, methC, strand, nextN);

	vector<int> numc;
	string dna = "";
	if(opts.reference != ""){
		dnaByChr(opts.reference, chromName, dna, numc);
	}

	boost::thread_group g;
	for(unsigned int i = 0; i< files.size(); i++)
	{
		boost::thread *tp = new boost::thread( createHashFromFile, files[i], chromName, boost::ref(meth[i]), boost::ref(dna));
		g.add_thread(tp);

	}
	g.join_all();


	//merge lanes into one file for current chromName
	if(opts.mappedFiles.size()>1)
	{
		set <int> starts;
		for(unsigned int i = 0; i < meth.size(); i++ ){
			//cout << "For file " << files[i] << ", number of cytosines processed: " << meth[i].size() <<endl;
			for(map<int, cMeth>::iterator it = meth[i].begin(); it != meth[i].end(); ++it) {
				//cout << it->first << endl;
				starts.insert(it->first);
			}
		}
		cout << "For file " << opts.sampleName << ", number of cytosines processed: " << starts.size() << " for chrom " << chromName << endl;
		map<int, cMeth> mergedMeth;
		for(set<int>::iterator pstart = starts.begin(); pstart != starts.end(); pstart++ )
		{
			int start = *pstart;
			cMeth mergedLoci;
			for(unsigned int i = 0; i< files.size(); i++)
			{
				if(meth[i].count(start) > 0){
					mergedLoci.totalC += meth[i][start].totalC;
					mergedLoci.methC  += meth[i][start].methC;
					mergedLoci.strand  = meth[i][start].strand;
					mergedLoci.nextN   = meth[i][start].nextN;
				}
			}
			//if(mergedLoci.totalC>5)cout << chromName << "\t" << start << "\t" << mergedLoci.totalC << endl;
			mergedMeth[start] = mergedLoci;
		}

		string ofile = opts.sampleName + "."+ chromName;
		ofstream statsFile;
		statsFile.open((ofile + "_stat.txt").c_str(), ios_base::out);
		strandSpecificBed(opts.sampleName, ofile, chromName, mergedMeth, statsFile);
		strandCombinedBed(opts.sampleName, ofile, chromName, mergedMeth, statsFile);
		statsFile.close();
//		// need remove the blank files
//		if( opts.fullMode == 0 ){
//			system( ("rm " + ofile + "_strand.bed").c_str() );
//			system( ("rm " + ofile + ".bed").c_str() );
//		}

	}

	meth.clear();
	cout << endl << "All files done for chrom " << chromName << endl;
};


void performCallChromStep(vector <string> files, int chromOffset, int chromStep)
{
	for(unsigned int i = chromOffset; i < chromNames.size(); i += chromStep){
		if(opts.skipRandomChrom && (chromNames[i].find("random") != string::npos || chromNames[i].find("had") != string::npos ) ){
			//skip this chrom
		} else {
			performCallChrom( files, chromNames[i]);
		}
	}
};

void mergeChroms(string file){
	string sysCmd = "";

	ofstream bGBedFile;
	bGBedFile.open((file + ".G.bed").c_str(), ios_base::out);
	bGBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
	if(opts.reference != "") bGBedFile << "\tlocalSeq";
	bGBedFile << endl;
	bGBedFile.close();
	sysCmd = "cat " + file + ".*.G.bed | grep -v chrom >> " + file + ".G.bed";
	system(sysCmd.c_str());


	ofstream bHGBedFile;
	bHGBedFile.open((file + ".HG.bed").c_str(), ios_base::out);
	bHGBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
	if(opts.reference != "") bHGBedFile << "\tlocalSeq";
	bHGBedFile << endl;
	bHGBedFile.close();
	sysCmd = "cat " + file + ".*.HG.bed | grep -v chrom >> " + file + ".HG.bed";
	system(sysCmd.c_str());

	if(opts.keepTemp){
	} else {
		sysCmd = "rm " + file + ".*.G.bed ";
		system(sysCmd.c_str());

		sysCmd = "rm " + file + ".*.HG.bed ";
		system(sysCmd.c_str());

		sysCmd = "rm " + file + ".*_stat.txt ";
		system(sysCmd.c_str());
	}

	ofstream statsFile;
	statsFile.open((file + "_stat.txt").c_str(), ios_base::out);
	statsFile << endl << "All reads = " << propertiesByLane[file].nAllReads << "; Mapped reads = " << propertiesByLane[file].nMappedReads << endl << endl;
	mergeStatsVec(file, statsFile);
	mergeStatbVec(file, statsFile);
	statsFile.close();

	if(opts.fullMode == 1){
			ofstream bBedFile;
			bBedFile.open((file + ".bed").c_str(), ios_base::out);
			bBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext\tPlus\ttotalC\tmethC\tMinus\ttotalC\tmethC";
			if(opts.reference != "") bBedFile << "\tlocalSeq";
			bBedFile << endl;
			bBedFile.close();
			sysCmd = "cat " + file + ".chr??.bed " + file + ".chr?.bed | grep -v chrom >> " + file + ".bed";
			system(sysCmd.c_str());

			ofstream skipBpFile;
			skipBpFile.open((file + "_skip.bed").c_str(), ios_base::out);
			skipBpFile << "#chrom\tstart\tend\tratio\ttotalC\tmeanC\tstrand\tnext" <<endl;
			skipBpFile.close();
			sysCmd = "cat " + file + ".*_skip.bed | grep -v chrom >> " + file + "_skip.bed";
			system(sysCmd.c_str());

			ofstream strandBedFile;
			strandBedFile.open((file + "_strand.bed").c_str(), ios_base::out);
			strandBedFile << "#chrom\tstart\tend\tratio\ttotalC\tmethC\tstrand\tnext" <<endl;
			strandBedFile.close();
			sysCmd = "cat " + file + ".*_strand.bed | grep -v chrom >> " + file + "_strand.bed";
			system(sysCmd.c_str());
		}
}

void initPrep(string file)
{
	//string file = opts.mappedFiles[i];

	////START:format detection; qualityScore detection; reference detection;
	string format = file.substr(file.find_last_of(".") + 1);
	boost::to_upper(format);
	if( format == "SAM") {
		cout << "From the extension of file " << file << ", program is parsing file according to SAM foramt" << endl;
	}
	else if( format == "BAM" ) {
		cout << "From the extension of file " << file << ", program is parsing file according to BAM foramt" << endl;
	}
	else if( format == "BSP" ) {
		cerr << "BSP format does not contain quality string and is not supported by this program" << endl;
		exit(1);
	}
	else{
		cerr << "Please use file extension .sam or .bam; other formats are not supported." <<endl;
		exit(1);
	}

	if(! found_XR_and_ZS(file, format)){
		cerr << "XR:Z or ZR:Z:, or ZS:Z: field not found, that is totally fine." <<endl;
		//exit(1);
		//todo:make this an internal decision
	}

	int qscorebase = opts.qualityScoreBase;
	qualityScoreFormatDetection(file, format, qscorebase);
	if(opts.qualityScoreBase == 0) {
		opts.qualityScoreBase = qscorebase;
	}
	//cout << opts.qualityScoreBase << endl;

	int readLen = 0;
	string protocol = "";
	detectChromsReadLenProtocol(file, format, readLen, protocol);
	////END:format detection; qualityScore detection; reference detection;

	////START: initialize nAllReads, nMappedReads
	long long nAllReads = 0; long long nMappedReads = 0;
	getReadCounts(file, format, nAllReads, nMappedReads);

	////END: initialize nAllReads, nMappedReads

	fileProperty property(format, qscorebase, nAllReads, nMappedReads, protocol, readLen);
	propertiesByLane[file] = property;

}

int main( int ac, char * av[])
{
	try
	{
		parse_options(ac, av);
		cout << "Program started" <<endl;

		set <string> chromsFromLanes;

		//initial preparation. reading through each file to get nuber of reads and mapped reads
		boost::thread_group gi;
		for(unsigned int i = 0; i< opts.mappedFiles.size(); i++){
			boost::thread *tp = new boost::thread( initPrep, opts.mappedFiles[i]);
			gi.add_thread(tp);
		}
		gi.join_all();

		if(opts.mappedFiles.size()>1){
			long long na = 0; long long nm = 0;
			for(unsigned int i = 0; i< opts.mappedFiles.size(); i++){
				string file = opts.mappedFiles[i];
				na += propertiesByLane[file].nAllReads;
				nm += propertiesByLane[file].nMappedReads;
			}
			propertiesByLane[opts.sampleName] = fileProperty(opts.sampleName, 0, na, nm, "WGBS", 100);//only na and nm are used
		}

		for(map<string, map<int, string> >::iterator it=chromsByLane.begin(); it!=chromsByLane.end(); it++){
			for(map<int, string>::iterator fit=it->second.begin(); fit!=it->second.end(); fit++){
				chromsFromLanes.insert(fit->second);
			}
		}
		chromNames.resize(chromsFromLanes.size());
		copy(chromsFromLanes.begin(), chromsFromLanes.end(), chromNames.begin());

		//perform methylation call chrom by chrom
		boost::thread_group g;
		for(int i = 0; i< opts.threads; i++){
			boost::thread *tp = new boost::thread( performCallChromStep, opts.mappedFiles, i, opts.threads);
			g.add_thread(tp);
		}
		g.join_all();

		//merge results from all chroms for each lane file
		boost::thread_group g1;
		for(unsigned int i = 0; i< opts.mappedFiles.size(); i++)
		{
			//string file = (i < opts.mappedFiles.size() ) ? opts.mappedFiles[i] : opts.sampleName;
			string file = opts.mappedFiles[i];
			boost::thread *tp = new boost::thread( mergeChroms, file);
			g1.add_thread(tp);
		}

		//merge results from all chroms for sample file
		if(opts.mappedFiles.size()>1){
			boost::thread *tp = new boost::thread( mergeChroms, opts.sampleName);
			g1.add_thread(tp);
		} else {
			if(opts.mappedFiles[0] == opts.sampleName){
				//
			} else {
				string sysCmd = "ln -s " + opts.mappedFiles[0] + ".G.bed " + opts.sampleName + ".G.bed";
				system(sysCmd.c_str());
			}
		}

		g1.join_all();

		cout << "Program successfully finished"<<endl;
	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
	}

	return 0;
};
