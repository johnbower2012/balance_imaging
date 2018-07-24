#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
	double Rx=4,Ry=3,Rz=6;
	int Nsample=1000000;
	double delq,q,r,x,y,z,qz,qperp=0.0,ctheta,g1,g2,g3,g4;
	int isample,Nq,iq;
 	CWaveFunction *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/wfparameters.dat");
	CRandom *randy=new CRandom(-1234);
	
	Nq=wf->GetNQMAX();
	delq=wf->GetDELQ();
	printf("check Nq=%d, delq=%g\n",Nq,delq);
	double *C=new double[Nq];
	for(iq=0;iq<Nq;iq++)
		C[iq]=0.0;
	
	for(isample=0;isample<Nsample;isample++){
		randy->gauss2(&g1,&g2); randy->gauss2(&g3,&g4);
		x=Rx*g1; y=Ry*g2; z=Rz*g3; r=sqrt(x*x+y*y+z*z);
		for(iq=0;iq<Nq;iq++){
			q=(0.5+iq)*delq;
			qz=q*(1.0-2.0*randy->ran());
			qperp=sqrt(q*q-qz*qz);
			ctheta=(qz*z+qperp*x)/(q*r);
			//printf("q=%g, r=%g, ctheta=%g\n",q,r,ctheta);
			C[iq]+=wf->CalcPsiSquared(q,r,ctheta);
		}
	}
	for(iq=0;iq<Nq;iq++){
		q=(0.5+iq)*delq;
		C[iq]=C[iq]/double(Nsample);
		printf("%6.1f %g\n",q,C[iq]);
	}
	
}

