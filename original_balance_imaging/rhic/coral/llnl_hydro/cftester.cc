#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	char sdirname[120],dirname[120],dataroot[120];
	double Rx,Ry,Rz,lambda;
	sprintf(dataroot,"/Users/pratt/data");
	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");
	Rx=3.0; Ry=5.0; Rz=7.0; lambda=1.0;
		CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");

// Calc Using CCHArray Objects
	sprintf(dirname,"%s/cfdata/gaussCH_Rx%g_Ry%g_Rz%g_lambda%g",dataroot,Rx,Ry,Rz,lambda);
/*
	CSourceCalc_Gaussian *scalc;
	scalc=new CSourceCalc_Gaussian();
	scalc->SetSPars(lambda,Rx,Ry,Rz);
	CCHArray *source=new CCHArray("parameters/apars_sf.dat");
	scalc->CalcS(source);
	source->FillRemainderX();
	printf("source projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	source->PrintProjections();
	printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	CKernel *kernel=new CKernel("parameters/kparameters.dat");
	kernel->Calc(wf);
	S2CF::s2c(source,kernel,cf);
	cf->WriteAX(dirname);
*/
	cf->ReadAX(dirname);

	cf->FillRemainderX();
	cf->PrintProjections();
	
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");
	sprintf(dirname,"%s/cfdata/gauss_Rx%g_Ry%g_Rz%g_lambda%g/3d",dataroot,Rx,Ry,Rz,lambda);
	c3d->ReadArray(dirname);
	c3d->PrintProjections();

	C3DArray *c3dCH=new C3DArray("parameters/apars3d_cf.dat");
	sprintf(dirname,"%s/cfdata/gaussCH_Rx%g_Ry%g_Rz%g_lambda%g/3d",dataroot,Rx,Ry,Rz,lambda);
	c3dCH->ReadArray(dirname);
	c3dCH->PrintProjections();
	
	C3DArray *c3dprime=new C3DArray("parameters/apars3d_cf.dat");
	ArrayCalc::Calc3DArrayFromAExpArray(cf,c3dprime);
	c3dprime->PrintProjections();
	
	int i,n=c3d->GetNXMAX();
	double *q=new double[n];
	double *c1=new double[n];
	double *c2=new double[n];
	double *c3=new double[n];
	double *c4=new double[n];
	double ctheta,phi,stheta;
	double qx,qy,qz,dq;
	dq=c3d->GetDELX();
	int icontinue,imc,nmc=100000;
	CRandom *randy=new CRandom(-12345);
	double x,y,z,cbar,root2=sqrt(2.0),ctheta_wf,r;
	do{
		printf("Enter ctheta & phi: ");
		scanf("%lf %lf",&ctheta,&phi);
		stheta=sqrt(1.0-ctheta*ctheta);
		for(i=0;i<n;i++){
			q[i]=(i+0.5)*dq;
			qx=q[i]*stheta*cos(phi);
			qy=q[i]*stheta*sin(phi);
			qz=q[i]*ctheta;
			c1[i]=cf->AExpand(qx,qy,qz);
			c2[i]=c3dCH->GetElement(qx,qy,qz);
			c3[i]=c3d->GetElement(qx,qy,qz);
			cbar=0.0;
			for(imc=0;imc<nmc;imc++){
				x=root2*Rx*randy->gauss();
				y=root2*Ry*randy->gauss();
				z=root2*Rz*randy->gauss();
				r=sqrt(x*x+y*y+z*z);
				ctheta_wf=(qx*x+qy*y+qz*z)/(q[i]*r);
				cbar+=(wf->GetPsiSquared(q[i],r,ctheta_wf)-1.0);
			}
			c4[i]=cbar/double(nmc);
			printf("%6.2f  %7.3f  %7.3f  %7.3f %7.3f\n",q[i],c1[i],c2[i],c3[i],c4[i]);
		}
		printf("Enter 0 to quit, other to continue: ");
		scanf("%d",&icontinue);
	}while(icontinue!=0);
}

