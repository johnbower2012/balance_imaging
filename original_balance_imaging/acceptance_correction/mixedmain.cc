#include "bal.h"
using namespace std;

int main(int argc, char * const argv[]){
	if(argc!=2){
		printf("USAGE: balance parsdirname\n");
		exit(1);
	}
	string parsdirname=argv[1];
	CBalance *balance=new CBalance(parsdirname);

	balance->InitArrays(); // Call only once, after parameters set
												 // and after acceptance is instantiated (needs to know ETA range to set bins)

	balance->CalcWeights();
	//balance->WriteWeightTables();
	//balance->WriteGQGPTables(); // output in tex format in results/directory
	balance->CalcMixed();
	return 0;
}
