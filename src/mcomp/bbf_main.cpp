#include "bbf.h"

int main(int arc, const char* argv[]){

// 	std::cout << arc << std::endl;
	std::vector<int> n;
	std::vector<int> k;

	for(int i = 1; i < arc; i++){
		if(i<= arc /2){
			n.push_back( atoi (argv[i]));
		} else {
			k.push_back( atoi (argv[i]));
		}
	}
// 	std::copy(n.begin(), n.end(), std::ostream_iterator<int>(std::cout, "~"));
// 	std::cout << "====================" << std::endl;
// 	std::copy(k.begin(), k.end(), std::ostream_iterator<int>(std::cout, "~"));
// 	std::cout << "====================" << std::endl;
	BiKey bk(-1,-1);
	BetaBinomialFit(n, k, bk);
	std::cout << bk.n1 << "," << bk.k1 << std::endl;
	return 0;
}
