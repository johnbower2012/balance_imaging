#include "coral.h"
using namespace std;

int main(){
	complex<double> cx,ceta,dcx,ci(0.0,1.0);
	const int ellmax=2;
	int ell;
	double FL,GL,FLprime,GLprime,FLa,FLb,GLa,GLb;
	double q=2.0,cxr=0.0,cxi,dx=0.00001;
	
	ceta=70.0/(137.036*ci*q);
	printf("ceta=(%10.3e,%10.3e)\n",real(ceta),imag(ceta));
	dcx=dx*ci;
	for(cxi=0.1;cxi<5.0;cxi+=1.0){
		cx=ci*cxi;
		printf("_____________________________________ cxi=%g ________________________________________________\n",cxi);
		for(ell=0;ell<=ellmax;ell++){
			CoulWave::GetFGprime_ComplexQ(ell,cx,ceta,&FL,&GL,&FLprime,&GLprime);
			CoulWave::GetFGprime_ComplexQ(ell,cx-dcx,ceta,&FLa,&GLa,&FLprime,&GLprime);
			CoulWave::GetFGprime_ComplexQ(ell,cx+dcx,ceta,&FLb,&GLb,&FLprime,&GLprime);
			printf("______ ell=%d _______\n",ell);
			printf("ell=%d, cx=(%g,%g), Fl=%g, Gl=%g, Flprime=%g, Glprime=%g\n", ell,real(cx),imag(cx),FL,GL,FLprime,GLprime);
			printf("ell=%d, cx=(%g,%g), Fl=%g, Gl=%g, Flprime=%g, Glprime=%g\n",ell,real(cx),imag(cx),FL,GL,0.5*(FLb-FLa)/dx,0.5*(GLb-GLa)/dx);
			printf("FL test, ZERO =? %g\n",(FLb+FLa-2.0*FL)/(dx*dx) +(ell+1.0)*ell*FL/real(cx*cx) +2.0*FL*real(ceta/cx) -FL);
			printf("GL test, ZERO =? %g\n",(GLb+GLa-2.0*GL)/(dx*dx) +(ell+1.0)*ell*GL/real(cx*cx) +2.0*GL*real(ceta/cx) -GL);
		}
	}
	exit(1);
	return 0;
}
