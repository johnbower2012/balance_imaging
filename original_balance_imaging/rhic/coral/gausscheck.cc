#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  int ir;
  double x,y,z,ctheta,phi,guess,r;
  //const double PI=4.0*atan(1.0);
	double Rx=15,Ry=4,Rz=3;
	double Xoff=0.0,Yoff=0.0,Zoff=0.0,lambda=1.0;
  char sdirname[160];
  char ArrayParsFilename[120];
  CSourceCalc_Gaussian *scalc;

  // Create Array to store Cart. Harmonic source info
  sprintf(ArrayParsFilename,"parameters/apars_sf.dat");
	CCHArray *Asource=new CCHArray(ArrayParsFilename);
  Asource->PrintPars();

  // Create Source Calc Object
  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(lambda,Rx,Ry,Rz,Xoff,Yoff,Zoff);
  parameter::PrintPars(scalc->spars);

  // Have Source Calc Object Fill Array
  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);
	Asource->PrintProjections();
	
  CWaveFunction *wf;
  CKernel *kernel;
  wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
  kernel=new CKernel("parameters/kparameters.dat");
  kernel->Calc(wf);

	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");
	S2CF::s2c(Asource,kernel,cf);
	cf->FillRemainderX();
	cf->Print(0,0,0);
	cf->WriteAX("/Users/scottepratt/data/cfdata/fakedata_kt250");
	cf->PrintProjections();

/*
  // Make 3D Cartesian array from Cart. Harmonic Data
  sprintf(ArrayParsFilename,"parameters/apars3d.dat\0");
  C3DArray *threed=new C3DArray(ArrayParsFilename);
  ArrayCalc::Calc3DArrayFromAExpArray(Asource,threed);
  printf("Created 3darray\n");

  // Calc S for some random angles and compare to correct answer
  for(ir=1;ir<Asource->GetNRADIAL();ir+=2){
    ctheta=1.0-2.0*randy.ran();
    phi=2.0*PI*randy.ran();
    r=(0.5+ir)*Asource->GetRADSTEP();
    z=r*ctheta;
    x=r*sqrt(1.0-ctheta*ctheta);
    y=x*sin(phi);
    x=x*cos(phi);
    guess=exp(-0.25*( ((x-Xoff)*(x-Xoff)/(Rx*Rx))
		      +((y-Yoff)*(y-Yoff)/(Ry*Ry))
		      +((z-Zoff)*(z-Zoff)/(Rz*Rz))));
    guess=lambda*guess/(Rx*Ry*Rz*pow(4.0*PI,1.5));
    printf("r=%g, S(%g,%g,%g)=%g =? %g =? %g\n",r,x,y,z,
	   Asource->AExpand(x,y,z),guess,
	   threed->GetElement(x,y,z));
  }
  
  // Write array info to file
  sprintf(sdirname,"sdata/gauss_R%g,%g,%g_off%g,%g,%g",
	  Rx,Ry,Rz,Xoff,Yoff,Zoff);
  Asource->WriteAllA(sdirname); 
*/
}

