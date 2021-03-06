#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <thread>
#include <cstdio>
#include "sam.h"

using namespace boost::program_options;
using namespace std;

class Opts {
	public:
		string infile;
		string outfile;
		string reference;
		int length;
		int numreads;
		string rscript;
	public:
		Opts():
			infile("")
			, outfile("")
			, reference("")
			, rscript("")
			, length(150)
			, numreads(2000000) { }
	public:
		void out() {
			cout << "infile: " << infile << endl;
			cout << "outfile: " << outfile << endl;
			cout << "reference: " << reference << endl;
			cout << "length: " << length << endl;
			cout << "numreads: " << numreads << endl;
			cout << "rscript: " << rscript << endl;
		}
} opts;

int parse_options(int ac, const char ** av) {
	try
	{
		options_description desc{"Allowed options"};
		desc.add_options()
			("help,h", "Produce help message.")
			("infile,i", value<string>()->default_value(""), "Input BAM file.")
			("outfile,o", value<string>()->default_value(""), "Output statistics.")
			("reference,r", value<string>()->default_value(""), "Reference FASTA file. This option is required.")
			("length,l", value<int>()->default_value(150), "Read length. Length of the query sequence in the BAM file may be shorter than the read length, but the read length should ensure the longest mapping. Default: 150.")
			("numreads,n", value<int>()->default_value(20000000), "Number of reads. First `n` reads will be examined. Be aware of extremely low CpG methylation levels when chrM is the first chromosome. Default: 20000000.")
			("rscript,s", value<string>()->default_value(""), "Rscript for mbias plot. Default: `$bindir/bseqc2mbiasplot.R`.")
			;

		variables_map vm;
		store(parse_command_line(ac, av, desc), vm);
		notify(vm);

		if (vm.count("help")) {
			cout << desc << endl;
			cout << "Examples:" <<endl;
			cout << "  " << av[0] << " -i in.bam -o result.txt -r hg38.fa" << endl;
			cout << endl;
			cout << "Date: 2020/05/17" << endl;
			cout << "Authors: Jin Li <lijin.abc@gmail.com>" << endl;
			exit(1);
		}

		for(map<string, variable_value>::iterator it=vm.begin(); it!=vm.end(); ++it) {
			string k=it->first;
			if(k=="infile"){
				opts.infile=vm[k].as<string>();
			} else if(k=="length"){
				opts.length=vm[k].as<int>();
			} else if(k=="numreads"){
				opts.numreads=vm[k].as<int>();
			} else if(k=="reference"){
				opts.reference=vm[k].as<string>();
			} else if(k=="outfile"){
				opts.outfile=vm[k].as<string>();
			} else if(k=="rscript"){
				opts.rscript=vm[k].as<string>();
			} else {
				cerr << "Error: invalid option " << k << endl;
				exit(1);
			}
		}
		if (opts.infile.empty()) {
			cerr << "Error: -i|--infile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.outfile.empty()) {
			cerr << "Error: -o|--outfile must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.reference.empty()) {
			cerr << "Error: -r|--reference must be specified." << endl;
			cout << desc << endl;
			exit(1);
		}
		if (opts.rscript.empty()) {
			opts.rscript=boost::filesystem::canonical("/proc/self/exe").parent_path().string()+"/bseqc2mbiasplot.R";
		}
		opts.out();
	} catch (const error &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}

// PE: true, SE: false
bool estimatelayout(string & infile) {
	samfile_t *in=0;
	if ((in=samopen(infile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << infile << endl;
		exit(-1);
	}
	int r=0;
	int count=0;
	bam1_t *b=bam_init1();
	bool layoutpe=false;
	while (count<1 && (r=samread(in, b))>=0) {
		uint32_t flag=b->core.flag;
		if (flag & 0x100) continue;
		if (flag & 0x1) layoutpe=true;
		count++;
	}
	samclose(in);
	return layoutpe;
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

int refbychr(string infile, string chr, string & ref)
{
	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty()) continue;
		if (line[0]=='>' && fastaheader2id(line)==chr) {
			while ((fin.good() && !fin.eof()) && getline(fin, line)) {
				if (line.empty()) continue;
				if (line[0]=='>') break;
				ref+=line;
			}
		}
	}
	fin.close();
	return 0;
}

int addtag(const bam1_t *b, string & ref, vector< vector< int > > & readcounts) {
	if(b->core.n_cigar > 3) return 1;

	size_t qlen=b->core.l_qseq;
	if (opts.length<qlen) {
		cerr << "Error: the input length (-l|--length) is less than the length of query sequence. "
			<< opts.length << "<" << qlen
			<< ". Please increase the input length and run the program again." << endl;
		exit(1);
	}

	char * s = (char*) bam1_seq(b);
	string seq(qlen, '\0');
	for (int i = 0; i<qlen; ++i) {
		seq[i] = bam_nt16_rev_table[bam1_seqi(s, i)];
	}

	string newSeq;
	int refLen=0;
	int spos=0;
	int okToProceed=1;
	uint32_t* cigar=bam1_cigar(b);
	for( int k=0; k< b->core.n_cigar; ++k)
	{
		int cop =cigar[k] & BAM_CIGAR_MASK; // operation
		int cl = cigar[k] >> BAM_CIGAR_SHIFT; // length
		switch(cop)
		{
			case BAM_CMATCH:
				newSeq  += seq.substr(spos, cl );
				spos += cl;
				refLen += cl;
				break;
			case BAM_CINS:
				spos += cl;
				break;
			case BAM_CDEL:
				for(int itemp = 0; itemp < cl; itemp++) { newSeq += 'N'; }
				refLen += cl;
				break;
			case BAM_CREF_SKIP:
				okToProceed = 0;
				break;
			case BAM_CSOFT_CLIP:
				okToProceed = 0;
				break;
			case BAM_CHARD_CLIP:
				okToProceed = 0;
				break;
			case BAM_CPAD:
				okToProceed = 0;
				break;
			default:
				okToProceed = 0;
				break;
		}
	}

	if(! okToProceed) return 1;
	seq=newSeq;
	string tag=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
	string refseq;
	bool forwardstrand=(tag[0]=='+');
	bool forwardsave=(tag=="++" || tag=="--");
	if (forwardstrand) {
		refseq=ref.substr(b->core.pos, refLen+1);
	} else {
		if (b->core.pos==0) {
			refseq="N"+ref.substr(0, refLen);
		} else {
			refseq=ref.substr(b->core.pos-1, refLen+1);
		}
	}
	boost::to_upper(refseq);

	if (forwardstrand) {
		if (forwardsave) {
			for (int i=0; i<qlen; ++i) {
				if (refseq[i]=='C') {
					if (refseq[i+1]=='G') {
						readcounts[0][i]+=(seq[i]=='C');
						readcounts[1][i]+=(seq[i]=='T');
					} else {
						readcounts[2][i]+=(seq[i]=='C');
						readcounts[3][i]+=(seq[i]=='T');
					}
				}
			}
		} else {
			for (int i=0; i<qlen; ++i) {
				if (refseq[i]=='C') {
					if (refseq[i+1]=='G') {
						readcounts[0][qlen-1-i]+=(seq[i]=='C');
						readcounts[1][qlen-1-i]+=(seq[i]=='T');
					} else {
						readcounts[2][qlen-1-i]+=(seq[i]=='C');
						readcounts[3][qlen-1-i]+=(seq[i]=='T');
					}
				}
			}
		}
	} else {
		if (forwardsave) {
			for (int i=0; i<qlen; ++i) {
				if (refseq[i+1]=='G') {
					if (refseq[i]=='C') {
						readcounts[0][i]+=(seq[i]=='G');
						readcounts[1][i]+=(seq[i]=='A');
					} else {
						readcounts[2][i]+=(seq[i]=='G');
						readcounts[3][i]+=(seq[i]=='A');
					}
				}
			}
		} else {
			for (int i=0; i<qlen; ++i) {
				if (refseq[i+1]=='G') {
					if (refseq[i]=='C') {
						readcounts[0][qlen-1-i]+=(seq[i]=='G');
						readcounts[1][qlen-1-i]+=(seq[i]=='A');
					} else {
						readcounts[2][qlen-1-i]+=(seq[i]=='G');
						readcounts[3][qlen-1-i]+=(seq[i]=='A');
					}
				}
			}
		}
	}
	return 0;
}

int file2readtags(string & infile, map< string, vector< string > > & readtags) {
	samfile_t *in=0;
	if ((in=samopen(infile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << infile << endl;
		return 1;
	}
	int r=0;
	int count=0;
	bam1_t *b=bam_init1();
	while (count<opts.numreads && (r=samread(in, b))>=0) {
		uint32_t flag=b->core.flag;
		if (flag & 0x100) continue;
		string qname=string((char*)bam1_qname(b));
		string zs=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
		map< string, vector< string > > :: iterator it=readtags.find(qname);
		if (readtags.end()!=it) {
			if (flag & 0x40) {
				it->second[0]=zs;
			} else if (flag & 0x80) {
				it->second[1]=zs;
			}
		} else {
			vector< string > tag(2, "N");
			if (flag & 0x40) {
				tag[0]=zs;
			} else if (flag & 0x80) {
				tag[1]=zs;
			}
			readtags[qname]=tag;
		}
		count++;
	}
	samclose(in);
	return 0;
}

int readtags2tagstats(map< string, vector< string > > & readtags, map< string, int > & tagstats) {
	for (map< string, vector< string > > :: iterator it=readtags.begin(); readtags.end()!=it; ++it) {
		string tag=it->second[0] + "," + it->second[1];
		map< string, int > :: iterator tit=tagstats.find(tag);
		if (tagstats.end()!=tit) {
			tit->second++;
		} else {
			tagstats[tag]=1;
		}
	}
	return 0;
}

int file2tagreadcounts(string & infile, map< string, vector< string > > & readtags, map< string, vector< vector< vector< int > > > > & tagreadcounts) {
	samfile_t *in=0;
	if ((in=samopen(infile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << infile << endl;
		return 1;
	}
	int r=0;
	int count=0;
	bam1_t *b=bam_init1();
	map< string, string > chrref; // chr->ref
	while (count<opts.numreads && (r=samread(in, b))>=0) {
		uint32_t flag=b->core.flag;
		if (flag & 0x100) continue;
		count++;

		string qname=string((char*)bam1_qname(b));
		string tag=readtags[qname][0] + "," + readtags[qname][1];
		string chr=in->header->target_name[b->core.tid];
		map< string, string > :: iterator itref=chrref.find(chr);
		if (chrref.end()==itref) {
			string ref;
			refbychr(opts.reference, chr, ref);
			chrref[chr]=ref;
		}
		if (flag & 0x40) {
			vector< vector< int > > &readcounts=tagreadcounts[tag][0];
			addtag(b, chrref[chr], readcounts);
		} else if (flag & 0x80) {
			vector< vector< int > > &readcounts=tagreadcounts[tag][1];
			addtag(b, chrref[chr], readcounts);
		}
	}
	samclose(in);
	return 0;
}

int tagreadcounts2tagrcstrand(map< string, vector< vector< vector< int > > > > &tagreadcounts, map< string, vector< vector< int > > > &tagrcstrand) {
	for (map< string, vector< vector< vector< int > > > > :: iterator it=tagreadcounts.begin(); it!=tagreadcounts.end(); ++it) {
		vector< string > tags;
		boost::split(tags, it->first, boost::is_any_of(","));
		for (int e=0; e<tags.size(); e++) {
			string tag=tags[e];
			if (tag!="N") {
				map< string, vector< vector< int > > > :: iterator ittag=tagrcstrand.find(tag);
				if (ittag!=tagrcstrand.end()) {
					for (int i=0; i<ittag->second.size(); i++) {
						for (int j=0; j<ittag->second[i].size(); j++) {
							ittag->second[i][j]+=it->second[e][i][j];
						}
					}
				} else {
					vector< vector< int > > rc(4, vector< int > (opts.length, 0));
					for (int i=0; i<rc.size(); i++) {
						for (int j=0; j<rc[i].size(); j++) {
							rc[i][j]=it->second[e][i][j];
						}
					}
					tagrcstrand[tag]=rc;
				}
			}
		}
	}
	return 0;
}

int tagrcs2methbsr(map< string, vector< vector< vector< int > > > > & tagreadcounts, map< string, vector< double > > & methbsr)
{
	for (map< string, vector< vector< vector< int > > > > :: iterator it=tagreadcounts.begin(); tagreadcounts.end()!=it; ++it) {
		string tag=it->first;
		vector< double > bsr(4, -1);
		vector< int > end1sum(4, 0);
		for (int i=0; i<it->second[0].size(); ++i) {
			for (int j=0; j<it->second[0][i].size(); ++j) {
				end1sum[i]+=it->second[0][i][j];
			}
		}
		if (end1sum[0]+end1sum[1]>0) {
			bsr[0]=1.0*end1sum[0]/(end1sum[0]+end1sum[1]);
		}
		if (end1sum[2]+end1sum[3]>0) {
			bsr[1]=1.0*end1sum[2]/(end1sum[2]+end1sum[3]);
		}

		vector< int > end2sum(4, 0);
		for (int i=0; i<it->second[1].size(); ++i) {
			for (int j=0; j<it->second[1][i].size(); ++j) {
				end2sum[i]+=it->second[1][i][j];
			}
		}
		if (end2sum[0]+end2sum[1]>0) {
			bsr[2]=1.0*end2sum[0]/(end2sum[0]+end2sum[1]);
		}
		if (end2sum[2]+end2sum[3]>0) {
			bsr[3]=1.0*end2sum[2]/(end2sum[2]+end2sum[3]);
		}
		methbsr[tag]=bsr;
	}
	return 0;
}

int tagrcs2methbsrstrand(map< string, vector< vector< int > > > &tagrcstrand, map< string, vector< double > > &methbsrstrand)
{
	for (map< string, vector< vector< int > > > :: iterator it=tagrcstrand.begin(); tagrcstrand.end()!=it; ++it) {
		string tag=it->first;
		vector< double > bsr(2, -1);
		vector< int > sum(4, 0);
		for (int i=0; i<it->second.size(); ++i) {
			for (int j=0; j<it->second[i].size(); ++j) {
				sum[i]+=it->second[i][j];
			}
		}
		if (sum[0]+sum[1]>0) {
			bsr[0]=1.0*sum[0]/(sum[0]+sum[1]);
		}
		if (sum[2]+sum[3]>0) {
			bsr[1]=1.0*sum[2]/(sum[2]+sum[3]);
		}
		methbsrstrand[tag]=bsr;
	}
	return 0;
}

int estimateprotocol(map< string, int > & tagstats, int & protocol, double & errorrate) {
	// 0: traditional library based on shotgun approach
	// 1: traditional library based on Nextera transposase approach
	// 2: pico
	protocol=0;
	if ((tagstats.end()!=tagstats.find("++,+-")
				&& tagstats.end()!=tagstats.find("+-,++")
				&& tagstats["++,+-"]<10*tagstats["+-,++"]
				&& tagstats["+-,++"]<10*tagstats["++,+-"])
			|| (tagstats.end()!=tagstats.find("-+,--")
				&& tagstats.end()!=tagstats.find("--,-+")
				&& tagstats["-+,--"]<10*tagstats["--,-+"]
				&& tagstats["--,-+"]<10*tagstats["-+,--"]
				))
	{
		protocol=2;
	} else if ((tagstats.end()!=tagstats.find("++,+-")
				&& tagstats.end()!=tagstats.find("+-,++")
				&& tagstats["+-,++"]>tagstats["++,+-"])
			|| (tagstats.end()!=tagstats.find("-+,--")
				&& tagstats.end()!=tagstats.find("--,-+")
				&& tagstats["--,-+"]>tagstats["-+,--"])
			|| (tagstats.end()!=tagstats.find("+-,++")
				&& tagstats.end()==tagstats.find("++,+-"))
			|| (tagstats.end()!=tagstats.find("--,-+")
				&& tagstats.end()==tagstats.find("-+,--"))
			)
	{
		protocol=1;
	}

	int total=0;
	for(map< string, int > :: iterator it=tagstats.begin(); it!=tagstats.end(); ++it) {
		total+=it->second;
	}
	if (total==0) return 1;
	int postivenumber=0;
	if (protocol==2) {
		vector < string > picotags {
			"++,+-", "+-,++", "-+,--", "--,-+"
				, "++,N", "N,++", "+-,N", "N,+-"
				, "-+,N", "N,-+", "--,N", "N,--"
		};
		for (string &tag : picotags) {
			map< string, int > :: iterator it=tagstats.find(tag);
			if (tagstats.end()!=it) {
				postivenumber+=it->second;
			}
		}
	} else if (protocol==0) {
		vector < string > tradtags {
			"++,+-", "-+,--"
				, "++,N", "N,+-"
				, "-+,N", "N,--"
		};
		for (string &tag : tradtags) {
			map< string, int > :: iterator it=tagstats.find(tag);
			if (tagstats.end()!=it) {
				postivenumber+=it->second;
			}
		}
	} else if (protocol==1) {
		vector < string > tradtags {
			"+-,++", "--,-+"
				, "+-,N", "N,++"
				, "--,N", "N,-+"
		};
		for (string &tag : tradtags) {
			map< string, int > :: iterator it=tagstats.find(tag);
			if (tagstats.end()!=it) {
				postivenumber+=it->second;
			}
		}
	}
	errorrate=1.0-1.0*postivenumber/total;
	return 0;
}

int tagrc2file(string & outfile, map< string, vector< vector< vector< int > > > > & tagrc)
{
	ofstream fout(outfile);
	fout << "tag" << '\t' << "end" << '\t' << "position" << '\t' << "numC_CG" << '\t' << "numT_CG" << '\t' << "numC_CH" << '\t' << "numT_CH" << endl;
	for (map< string, vector< vector< vector< int > > > > :: iterator it=tagrc.begin(); it!=tagrc.end(); ++it) {
		for (int e=0; e<it->second.size(); e++) {
			for (int i=0; i<opts.length; i++) {
				fout << it->first << '\t' << e << '\t' << i+1 << '\t' << it->second[e][0][i] << '\t' << it->second[e][1][i] << '\t' << it->second[e][2][i] << '\t' << it->second[e][3][i] << endl;
			}
		}
	}
	return 0;
}

int tagrcstrand2file(string & outfile, map< string, vector< vector< int > > > & tagrcstrand)
{
	ofstream fout(outfile);
	fout << "tag" << '\t' << "position" << '\t' << "numC_CG" << '\t' << "numT_CG" << '\t' << "numC_CH" << '\t' << "numT_CH" << endl;
	for (map< string, vector< vector< int > > > :: iterator it=tagrcstrand.begin(); it!=tagrcstrand.end(); ++it) {
		for (int i=0; i<opts.length; i++) {
			fout << it->first << '\t' << i+1 << '\t' << it->second[0][i] << '\t' << it->second[1][i] << '\t' << it->second[2][i] << '\t' << it->second[3][i] << endl;
		}
	}
	return 0;
}

string errorstatus(double meth1ch, double meth2ch) {
	int status1=-1;
	if (meth1ch>0) {
		status1=0;
		if (meth1ch>0.05) {
			status1=1;
		} else if (meth1ch>0.02) {
			status1=2;
		}
	}
	int status2=-1;
	if (meth2ch>0) {
		status2=0;
		if (meth2ch>0.05) {
			status2=1;
		} else if (meth2ch>0.02) {
			status2=2;
		}
	}

	string status;
	if (status1==1 || status2==1) {
		status="Error";
	} else if (status1==0 && status2==0) {
		status="OK";
	} else {
		if (status1>=0) {
			status=(status1==0)?"OK":"Alert";
			if (status2>=0) {
				status+= ",";
				status+=(status2==0)?"OK":"Alert";
			}
		} else {
			if (status2>=0) {
				status=(status2==0)?"OK":"Alert";
			}
		}
	}
	return status;
}

int methbsr2file(ostream & out, int protocol, map< string, int > & tagstats, map< string, vector< double > > & methbsr, map< string, vector< double > > & methbsrstrand)
{
	set< string > validtags;
	if (protocol==2) {
		validtags={
			"++,+-", "+-,++", "-+,--", "--,-+"
				, "++,N", "N,++", "+-,N", "N,+-"
				, "-+,N", "N,-+", "--,N", "N,--"
		};
	} else if (protocol==0) {
		validtags={
			"++,+-", "-+,--"
				, "++,N", "N,+-"
				, "-+,N", "N,--"
		};
	} else if (protocol==1) {
		validtags={
			"+-,++", "--,-+"
				, "N,++", "+-,N"
				, "N,-+", "--,N"
		};
	}

	// 24 PE pairs
	out << "tag" << '\t' << "numreads" << '\t' << "end1CG" << '\t' << "end1CH" << '\t' << "end2CG" << '\t' << "end2CH" << '\t' << "bcrstatus" << endl;
	for (map< string, int > :: iterator it=tagstats.begin(); it!=tagstats.end(); ++it) {
		string tag=it->first;
		int numreads=it->second;
		double meth1cg=methbsr[tag][0];
		double meth1ch=methbsr[tag][1];
		double meth2cg=methbsr[tag][2];
		double meth2ch=methbsr[tag][3];
		string status;
		if (validtags.end()!=validtags.find(tag)) {
			status=errorstatus(meth1ch, meth2ch);
		}
		out << tag << '\t' << numreads << '\t' << meth1cg << '\t' << meth1ch << '\t' << meth2cg << '\t' << meth2ch << '\t' << status << endl;
	}

	// 4 single strands
	out << "tag" << '\t' << "methCG" << '\t' << "methCH" << endl;
	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		out << it->first;
		for (int i=0; i<it->second.size(); i++) {
			out << '\t' << it->second[i];
		}
		out << endl;
	}
	return 0;
}

int mbiasplot(string & infile, bool pico, string & outfile) {
	string cmd = "R --slave --no-save --no-restore --no-init-file";
	cmd+=" -e width=6";
	cmd+=" -e height=3";
	if (pico) {
		cmd+=" -e pico=T";
	} else {
		cmd+=" -e pico=F";
	}
	cmd+=" -e \"xlab='Position in read (bp)'\"";
	cmd+=" -e \"ylab='Methylation level'\"";
	cmd+=" -e \"infile='" + infile + "'\"";
	cmd+=" -e \"outfile='" + outfile + "'\"";
	cmd+=" -e \"source('" + opts.rscript + "')\"";
	cout << cmd << endl;
	FILE *fp;
	char info[10240];
	fp = popen(cmd.c_str(), "r");
	if (fp==NULL) {
		fprintf(stderr, "popen error.\n");
		return EXIT_FAILURE;
	}
	while (fgets(info, 10240, fp) != NULL) {
		fprintf(stderr, "%s", info);
	}
	pclose(fp);
	return 0;
}

int failureqc(ostream & out, double errorrate, int protocol, map< string, vector< double > > &methbsrstrand, map< string, vector< double > > &methbsr) {
	int failure=0; // 0: pass, 1: failure, 2: warning
	if (errorrate>0.05) {
		out << "We detected library error (failure): too high of the mapping error rate " << errorrate << endl;
		failure=1;
	} else if (errorrate>0.02) {
		out << "We detected library error (warning): a little high of the mapping error rate " << errorrate << endl;
		failure=2;
	}

	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		double bcr=1.0-it->second[1];
		if (bcr<0.95) {
			out << "We detected library error (failure): too low of the bisulfite conversion rate " << bcr << " for strand " << it->first << endl;
			failure=1;
		} else if (bcr<0.98) {
			out << "We detected library error (warning): a little low of the bisulfite conversion rate " << bcr << " for strand " << it->first << endl;
			if (failure==0) {
				failure=2;
			}
		}
	}

	double minmeth=1.0;
	double maxmeth=0;
	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		if (it->second[0]>maxmeth) {
			maxmeth=it->second[0];
		}
		if (it->second[0]<minmeth) {
			minmeth=it->second[0];
		}
	}
	if (maxmeth-minmeth>0.05) {
		out << "We detected library error (failure): inconsistent average methylation level in four strands" << endl;
		failure=1;
	} else if (maxmeth-minmeth>0.02) {
		out << "We detected library error (warning): relatively inconsistent average methylation level in four strands" << endl;
		if (failure==0) {
			failure=2;
		}
	}

	if (protocol==2) {
		vector< string > tags{"++,+-", "+-,++", "-+,--", "--,-+"};
		for (string &tag: tags) {
			map< string, vector< double > > :: iterator it=methbsr.find(tag);
			if (methbsr.end()!=it) {
				if ((it->second[0]-it->second[2]>0.05) || (it->second[2]-it->second[0]>0.05)) {
					out << "We detected library error (failure): inconsistent average methylation levels of two ends in " << tag << endl;
					failure=1;
				}
			}
		}
	} else if (protocol==0) {
		vector< string > tags{"++,+-", "-+,--"};
		for (string &tag: tags) {
			map< string, vector< double > > :: iterator it=methbsr.find(tag);
			if (methbsr.end()!=it) {
				if ((it->second[0]-it->second[2]>0.05) || (it->second[2]-it->second[0]>0.05)) {
					out << "We detected library error (failure): inconsistent average methylation levels of two ends in " << tag << endl;
					failure=1;
				}
			}
		}
	} else if (protocol==1) {
		vector< string > tags{"+-,++", "--,-+"};
		for (string &tag: tags) {
			map< string, vector< double > > :: iterator it=methbsr.find(tag);
			if (methbsr.end()!=it) {
				if ((it->second[0]-it->second[2]>0.05) || (it->second[2]-it->second[0]>0.05)) {
					out << "We detected library error (failure): inconsistent average methylation levels of two ends in " << tag << endl;
					failure=1;
				}
			}
		}
	}
	if (failure==0) {
		out << "Info: successful library." << endl;
	}
	return failure;
}

int bseqc_pe() {
	map< string, vector< string > > readtags; // qname->[tag1,tag2]
	file2readtags(opts.infile, readtags);

	map< string, int > tagstats; // tag->number
	readtags2tagstats(readtags, tagstats);
	int protocol=0;
	double errorrate=-1;
	estimateprotocol(tagstats, protocol, errorrate);

	boost::filesystem::create_directories(boost::filesystem::absolute(opts.outfile).parent_path());
	ofstream fout(opts.outfile);
	if (protocol==2) {
		fout << "Pico library construction detected. 12 positive PE mapping pairs:\n(++,+-), (+-,++), (-+,--), (--,-+), (++,N), (N,++), (+-,N), (N,+-), (-+,N), (N,-+), (--,N), (N,--)" << endl;
	} else if (protocol==0) {
		fout << "Traditional library construction based on shotgun approach detected. 6 positive PE mapping pairs:\n(++,+-), (-+,--), (++,N), (N,+-), (-+,N), (N,--)" << endl;
	} else if (protocol==1) {
		fout << "Traditional library construction based on Nextera transposase approach detected. 6 positive PE mapping pairs:\n(+-,++), (--,-+), (+-,N), (N,++), (--,N), (N,-+)" << endl;
	}
	if (errorrate>=0) {
		fout << "Estimated error rate: " << errorrate << ", positive rate: " << 1-errorrate << endl;
	}

	map< string, vector< vector< vector< int > > > > tagreadcounts; // tag->[end1,end2]; end[12]=[CG, CH]=[[c1,...,clength],[t1,...,tlength]]
	for (map< string, int > :: iterator it=tagstats.begin(); tagstats.end()!=it; ++it) {
		vector< vector< vector< int > > > twoendsreadcounts(2, vector< vector< int > > (4, vector< int > (opts.length, 0)));
		tagreadcounts[it->first]=twoendsreadcounts;
	}
	file2tagreadcounts(opts.infile, readtags, tagreadcounts);

	string obname=boost::filesystem::path(opts.outfile).replace_extension().string();
	string outtagrcfile=obname+"_mbias_pe.txt";
	tagrc2file(outtagrcfile, tagreadcounts);
	string outtagrcplot=obname+"_mbias_pe.pdf";
	mbiasplot(outtagrcfile, protocol==2, outtagrcplot);

	map< string, vector< double > > methbsr; // tag->[end1,end2]; end[12]=[CG, CH]
	tagrcs2methbsr(tagreadcounts, methbsr);

	map< string, vector< vector< int > > > tagrcstrand;
	tagreadcounts2tagrcstrand(tagreadcounts, tagrcstrand);

	string outtagrcstrandfile=obname+"_mbias_strand.txt";
	tagrcstrand2file(outtagrcstrandfile, tagrcstrand);
	string outtagrcstrandplot=obname+"_mbias_strand.pdf";
	mbiasplot(outtagrcstrandfile, protocol==2, outtagrcstrandplot);

	map< string, vector< double > > methbsrstrand; // tag->[CG, CH]
	tagrcs2methbsrstrand(tagrcstrand, methbsrstrand);

	methbsr2file(fout, protocol, tagstats, methbsr, methbsrstrand);
	int exit_code=failureqc(fout, errorrate, protocol, methbsrstrand, methbsr);
	return exit_code;
}

int file2tagstats(string & infile, map< string, int > & tagstats) {
	samfile_t *in=0;
	if ((in=samopen(infile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << infile << endl;
		return 1;
	}
	int r=0;
	int count=0;
	bam1_t *b=bam_init1();
	while (count<opts.numreads && (r=samread(in, b))>=0) {
		uint32_t flag=b->core.flag;
		if (flag & 0x100) continue;
		string zs=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
		tagstats[zs]++;
		count++;
	}
	samclose(in);
	return 0;
}

int file2tagreadcounts(string & infile, map< string, vector< vector< int > > > & tagreadcounts) {
	samfile_t *in=0;
	if ((in=samopen(infile.c_str(), "rb", 0))==0) {
		cerr << "Error: not found " << infile << endl;
		return 1;
	}
	int r=0;
	int count=0;
	bam1_t *b=bam_init1();
	map< string, string > chrref; // chr->ref
	while (count<opts.numreads && (r=samread(in, b))>=0) {
		uint32_t flag=b->core.flag;
		if (flag & 0x100) continue;
		string chr=in->header->target_name[b->core.tid];
		map< string, string > :: iterator itref=chrref.find(chr);
		if (chrref.end()==itref) {
			string ref;
			refbychr(opts.reference, chr, ref);
			chrref[chr]=ref;
		}
		string tag=string((char *)bam_aux2Z(bam_aux_get(b, "ZS")));
		vector< vector< int > > &readcounts=tagreadcounts[tag];
		addtag(b, chrref[chr], readcounts);
		count++;
	}
	samclose(in);
	return 0;
}

int methbsr2file(ostream & out, map< string, int > & tagstats, map< string, vector< double > > & methbsrstrand)
{
	out << "tag" << '\t' << "numreads" << '\t' << "methCG" << '\t' << "methCH" << endl;
	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		string tag=it->first;
		out << tag << '\t' << tagstats[tag];
		for (int i=0; i<it->second.size(); i++) {
			out << '\t' << it->second[i];
		}
		out << endl;
	}
	return 0;
}

int failureqc_se(ostream & out, double errorrate, int protocol, map< string, vector< double > > &methbsrstrand) {
	int failure=0; // 0: pass, 1: failure, 2: warning

	if (errorrate>0.05) {
		out << "We detected library error (failure): too high of the mapping error rate " << errorrate << endl;
		failure=1;
	} else if (errorrate>0.02) {
		out << "We detected library error (warning): a little high of the mapping error rate " << errorrate << endl;
		failure=2;
	}

	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		if (protocol==0 && (it->first=="+-" || it->first=="--")) continue; // traditional library by shotgun should not check +- and --
		if (protocol==1 && (it->first=="++" || it->first=="-+")) continue; // traditional library by Nextera transposase should not check ++ and -+
		double bcr=1.0-it->second[1];
		if (bcr<0.95) {
			out << "We detected library error (failure): too low of the bisulfite conversion rate " << bcr << " for strand " << it->first << endl;
			failure=1;
		} else if (bcr<0.98) {
			out << "We detected library error (warning): a little low of the bisulfite conversion rate " << bcr << " for strand " << it->first << endl;
			if (failure==0) {
				failure=2;
			}
		}
	}

	double minmeth=1.0;
	double maxmeth=0;
	for (map< string, vector< double > > :: iterator it=methbsrstrand.begin(); it!=methbsrstrand.end(); ++it) {
		if (protocol==0 && (it->first=="+-" || it->first=="--")) continue; // traditional library by shotgun should not check +- and --
		if (protocol==1 && (it->first=="++" || it->first=="-+")) continue; // traditional library by Nextera transposase should not check ++ and -+
		if (it->second[0]>maxmeth) {
			maxmeth=it->second[0];
		}
		if (it->second[0]<minmeth) {
			minmeth=it->second[0];
		}
	}
	if (maxmeth-minmeth>0.05) {
		if (protocol==2) {
			out << "We detected library error (failure): inconsistent average methylation level in four strands" << endl;
		} else {
			out << "We detected library error (failure): inconsistent average methylation level in two strands" << endl;
		}
		failure=1;
	} else if (maxmeth-minmeth>0.02) {
		if (protocol==2) {
			out << "We detected library error (warning): relatively inconsistent average methylation level in four strands" << endl;
		} else {
			out << "We detected library error (warning): relatively inconsistent average methylation level in two strands" << endl;
		}
		if (failure==0) {
			failure=2;
		}
	}

	if (failure==0) {
		out << "Info: successful library." << endl;
	}
	return failure;
}

int estimateprotocol_se(map< string, int > & tagstats, int & protocol, double & errorrate) {
	// 0: traditional library based on shotgun approach
	// 1: traditional library based on Nextera transposase approach
	// 2: pico
	protocol=0;
	if ((tagstats.end()!=tagstats.find("++")
				&& tagstats.end()!=tagstats.find("+-")
				&& tagstats["++"]<10*tagstats["+-"]
				&& tagstats["+-"]<10*tagstats["++"])
			|| (tagstats.end()!=tagstats.find("-+")
				&& tagstats.end()!=tagstats.find("--")
				&& tagstats["-+"]<10*tagstats["--"]
				&& tagstats["--"]<10*tagstats["-+"]
				))
	{
		protocol=2;
	} else if ((tagstats.end()!=tagstats.find("++")
				&& tagstats.end()!=tagstats.find("+-")
				&& tagstats["+-"]>tagstats["++"])
			|| (tagstats.end()!=tagstats.find("-+")
				&& tagstats.end()!=tagstats.find("--")
				&& tagstats["--"]>tagstats["-+"])
			|| (tagstats.end()!=tagstats.find("+-")
				&& tagstats.end()==tagstats.find("++"))
			|| (tagstats.end()!=tagstats.find("--")
				&& tagstats.end()==tagstats.find("-+"))
			)
	{
		protocol=1;
	}

	int total=0;
	for(map< string, int > :: iterator it=tagstats.begin(); it!=tagstats.end(); ++it) {
		total+=it->second;
	}
	if (total==0) return 1;
	// error rate:
	// 1. pico library
	//   four valid strands: ++, +-, -+, --
	//   no more evidence to exclude any of 4 strands
	// 2. traditional library
	//   Two valid strands for shotgun approach: ++ and -+
	//   Two valid strands for Nextera transposase: +- and --
	if (protocol==0) { // shotgun
		int postivenumber=0;
		vector < string > tradtags {"++", "-+"};
		for (string &tag : tradtags) {
			map< string, int > :: iterator it=tagstats.find(tag);
			if (tagstats.end()!=it) {
				postivenumber+=it->second;
			}
		}
		errorrate=1.0-1.0*postivenumber/total;
	} else if (protocol==1) { // Nextera transposase
		int postivenumber=0;
		vector < string > tradtags {"+-", "--"};
		for (string &tag : tradtags) {
			map< string, int > :: iterator it=tagstats.find(tag);
			if (tagstats.end()!=it) {
				postivenumber+=it->second;
			}
		}
		errorrate=1.0-1.0*postivenumber/total;
	}
	return 0;
}

int bseqc_se() {
	boost::filesystem::create_directories(boost::filesystem::absolute(opts.outfile).parent_path());
	ofstream fout(opts.outfile);
	fout << "Single-end mapping detected. Strand QC will be examined only." << endl;

	map< string, int > tagstats; // tag->number
	file2tagstats(opts.infile, tagstats);

	int protocol=0;
	double errorrate=-1;
	estimateprotocol_se(tagstats, protocol, errorrate);

	if (protocol==2) {
		fout << "Pico library construction detected. 4 strands mapping are positive: ++, +-, -+, --" << endl;
	} else if (protocol==0) {
		fout << "Traditional library construction based on shotgun approach detected. 2 strands are positive: ++ and -+" << endl;
	} else if (protocol==1) {
		fout << "Traditional library construction based on Nextera transposase approach detected. 2 strands are positive: +- and --" << endl;
	}
	if (errorrate>=0) {
		fout << "Estimated error rate: " << errorrate << ", positive rate: " << 1-errorrate << endl;
	}

	map< string, vector< vector< int > > > tagreadcounts; // tag->[CG_C, CG_T, CH_C, CH_T]=[[c1,...,clength],[t1,...,tlength]]
	for (map< string, int > :: iterator it=tagstats.begin(); tagstats.end()!=it; ++it) {
		vector< vector< int > > readcounts(4, vector< int > (opts.length, 0));
		tagreadcounts[it->first]=readcounts;
	}
	file2tagreadcounts(opts.infile, tagreadcounts);

	string obname=boost::filesystem::path(opts.outfile).replace_extension().string();
	string outtagrcstrandfile=obname+"_mbias_strand.txt";
	tagrcstrand2file(outtagrcstrandfile, tagreadcounts);
	string outtagrcstrandplot=obname+"_mbias_strand.pdf";
	mbiasplot(outtagrcstrandfile, false, outtagrcstrandplot);

	map< string, vector< double > > methbsrstrand; // tag->[CG, CH]
	tagrcs2methbsrstrand(tagreadcounts, methbsrstrand);

	methbsr2file(fout, tagstats, methbsrstrand);
	int exit_code=failureqc_se(fout, errorrate, protocol, methbsrstrand);
	return exit_code;
}

int bseqc() {
	bool layoutpe=estimatelayout(opts.infile);
	if (layoutpe) {
		return bseqc_pe();
	}
	return bseqc_se();
}

int main(int argc, const char ** argv)
{
	parse_options(argc, argv);
	int exit_code=bseqc();
	return exit_code;
}
