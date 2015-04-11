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


