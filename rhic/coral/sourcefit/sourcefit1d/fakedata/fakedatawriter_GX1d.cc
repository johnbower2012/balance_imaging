#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <ctime>

using namespace std;

#include "coral.h"

int main(){
  CRandom *randy;
  CCHArray *Cfake;
  CCHArray *Ctheory;
  CCHArray *Efake;
  CCHArray *Sgx;
  CSourceCalc_GX1D *scalc;
  CWaveFunction *wf;
  CKernel *kernel;
  double R=5.0,X=10.0,lambda=0.6,Xfrac=0.5,a=5.0;
  double error=0.05,chisquare,emag;
  int lx=0,ly=0,lz=0,lxcalc=0,lycalc=0,lzcalc=0,iq,dlx=1,dly=1,dlz=1;
  string pardirname="../parameters/pp";
  randy=new CRandom(-time(NULL));
	
  // Create Arrays
  string apars_filename;
	apars_filename=pardirname+"/aparsCH_source.dat";
  Sgx=new CCHArray(apars_filename);
  apars_filename=pardirname+"/aparsCH_CF.dat";
  Cfake=new CCHArray(apars_filename);
  Ctheory=new CCHArray(apars_filename);
  Efake=new CCHArray(apars_filename);
	
  string fakedatadirname="GX1d/pp";
  string fakeerrordirname="GX1d/pp_errors";
	
  string wfparsfilename=pardirname+"/kparameters.dat";
  string kdatadirname="../../kdata/pp";
  kernel=new CKernel(wfparsfilename);
  //kernel->ReadData(kdatadirname);
  wf=new CWaveFunction_pp_schrod(wfparsfilename);
  kernel->Calc(wf);
	kernel->WriteData(kdatadirname);
	
  scalc=new CSourceCalc_GX1D();
  scalc->SetSPars(lambda,Xfrac,R,X,a);
  scalc->CalcS(Sgx);
  scalc->NormCheck(Sgx);
  S2CF::s2c(Sgx,kernel,Ctheory);
	
	
  // Make error matrix
  if(Efake->GetXSYM()) dlx=2;
  if(Efake->GetYSYM()) dly=2;
  if(Efake->GetZSYM()) dlz=2;
  int Lmax=Efake->GetLMAX();
  for(iq=0;iq<Efake->GetNRADIAL();iq++){
    for(lx=0;lx<2;lx+=dlx){
      for(ly=0;ly<=Lmax-lx;ly+=dly){
				for(lz=0;lz<=Lmax-lx-ly;lz+=dlz){
					emag=error/(iq+1.0);
					Cfake->SetElement(lx,ly,lz,iq,emag*randy->gauss()
						+Ctheory->GetElement(lx,ly,lz,iq));
					Efake->SetElement(lx,ly,lz,iq,emag);
				}
      }
    }
  }
	
  chisquare=CFCalc::GetChiSquared(lxcalc,lycalc,lzcalc,Cfake,Efake,Ctheory);
  printf("chisquare=%g\n",chisquare);
  
  // Write Arrays
  Cfake->WriteAX(fakedatadirname);
  Efake->WriteAX(fakeerrordirname);
	
  return 0;
  
}

