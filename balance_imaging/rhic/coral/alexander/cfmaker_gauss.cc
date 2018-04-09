#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	char sdirname[120],dirname[120],dataroot[120];
	double Rx,Ry,Rz,lambda,roff[3],r2[3][3];
	double q,hbarc=197.326;
	int i;

	CSourceCalc_Gaussian *scalc;
	scalc=new CSourceCalc_Gaussian();
	Rx=14.639; Ry=5.418; Rz=6.787; lambda=0.5190;
	scalc->SetSPars(lambda,Rx,Ry,Rz);
	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
	
	sprintf(dataroot,"/Users/pratt/data");
	
	// Calc Using CCHArray Objects
	CCHArray *source=new CCHArray("parameters/apars_sf.dat");
	scalc->CalcS(source);
	source->FillRemainderX();
	printf("SF projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	source->PrintProjections();
	source->WriteProjections("sf_gauss.dat");
	CKernel *kernel=new CKernel("parameters/kparameters.dat");
	kernel->Calc(wf);
	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");	
	S2CF::s2c(source,kernel,cf);
	cf->FillRemainderX();
	printf("CF projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	cf->PrintProjections();
	cf->WriteProjections("cf_gauss.dat");

}

