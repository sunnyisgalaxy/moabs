#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <time.h>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>

using namespace std;

int string_to_int( string s){return atoi(s.c_str());}
string itos(int i)
{
    stringstream s;
    s << i;
    return s.str();
}



std::string do_readlink(std::string const& path) {
    char buff[1024];
    ssize_t len = ::readlink(path.c_str(), buff, sizeof(buff)-1);
    if (len != -1) {
      buff[len] = '\0';
      return std::string(buff);
    } else {
     /* handle error condition */
    }
}


std::string get_exepath(){
		//string exep(  (char *)getauxval(AT_EXECFN) );
		std::string exep = do_readlink("/proc/self/exe");
		std::vector<std::string> splits;
		boost::split(splits, exep, boost::is_any_of("/"));
		splits.pop_back();
		exep = boost::algorithm::join(splits, "/");
		
		//std::cout << exep << std::endl;
		
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
	
		return exep;
}

bool isconsensus(std::vector< int > & tcs, std::vector< int > & mcs, int depth=10) {
	double maxratio=0.0;
	double minratio=1.0;
	for (int locus=0; locus<tcs.size(); locus++) {
		if (tcs[locus]<depth) return false;
		double locusratio = 1.0*mcs[locus]/tcs[locus];
		maxratio=std::max(maxratio, locusratio);
		minratio=std::min(minratio, locusratio);
	}
	return maxratio<minratio+0.1;
}
