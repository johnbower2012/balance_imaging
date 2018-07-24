#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <armadillo>
#include <random>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

int main(int argc,char *argv[]){
	double Itest=0.0,dx=0.001,x;
	for(x=0.5*dx;x<40;x+=dx){
		Itest+=2.0*dx/(cosh(x)*cosh(2*x));
	}
	printf("Itest= %g =? %g\n",Itest,(sqrt(2.0)-1.0)*PI);
	
	return 0;
}

