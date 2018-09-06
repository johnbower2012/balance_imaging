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
	Rx=3.0; Ry=4.0; Rz=5.0; lambda=1.0;
	scalc->SetSPars(lambda,Rx,Ry,Rz);
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");
	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	//CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
	
	sprintf(dataroot,"/Users/pratt/data");
	
	// Calc Using 3D Objects
	C3DArray *s3d=new C3DArray("parameters/apars3d_sf.dat");
	CKernelWF *kernelwf=new CKernelWF("parameters/kparameters.dat");
	kernelwf->Calc(wf);
	scalc->CalcS(s3d);
	s3d->PrintMoments();
	S2CF::s2c(s3d,kernelwf,c3d);
	//S2CF::s2c(s3d,wf,c3d);
	printf("_____________ c3d _________________\n");
	c3d->PrintProjections();
	//c3d->PrintMoments();
	printf("exponential guess\n");
	for(i=0;i<30;i++){
		q=(i+0.5)*2.0;
		printf("%3d %10.3e %10.3e %10.3e\n",i,exp(-4.0*q*q*9.0/(hbarc*hbarc)),exp(-4.0*q*q*25.0/(hbarc*hbarc)),exp(-4.0*q*q*49.0/(hbarc*hbarc)));
	}
	sprintf(dirname,"%s/cfdata/gauss_Rx%g_Ry%g_Rz%g_lambda%g/3d",dataroot,Rx,Ry,Rz,lambda);
	c3d->WriteArray(dirname);	

	// Calc Using CCHArray Objects
	CCHArray *source=new CCHArray("parameters/apars_sf.dat");
	scalc->CalcS(source);
	source->FillRemainderX();
	//printf("source projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	//source->PrintProjections();
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	CKernel *kernel=new CKernel("parameters/kparameters.dat");
	kernel->Calc(wf);
	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");	
	S2CF::s2c(source,kernel,cf);
	cf->FillRemainderX();
	//cf->PrintProjections();
	ArrayCalc::Calc3DArrayFromAExpArray(cf,c3d);
	//c3d->PrintProjections();
	//c3d->PrintMoments();
	sprintf(dirname,"%s/cfdata/gaussCH_Rx%g_Ry%g_Rz%g_lambda%g/3d",dataroot,Rx,Ry,Rz,lambda);
	c3d->WriteArray(dirname);
	
}

