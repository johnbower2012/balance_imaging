#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CCHArray *Cfake;
  C3DArray *Cfake3D;
  C3DArray *Ctheory3D;
  C3DArray *Efake3D;
  CCHArray *S;
  CSourceCalc_Blast *scalc;
  CWaveFunction *wf;
  CKernel *kernel;
  double lambda=0.6,R=13.0,Tau=12.0,DelTau=5.0;
  double error=0.05,chisquare;
  string pardirname="../parameters/pp";

  // Create Arrays
  string apars_filename=pardirname+"/aparsCH_source.dat";
  S=new CCHArray(apars_filename);
  apars_filename=pardirname+"/aparsCH_CF.dat";
  Cfake=new CCHArray(apars_filename);
  apars_filename=pardirname+"/apars3D_CF.dat";
  Cfake3D=new C3DArray(apars_filename);
  Ctheory3D=new C3DArray(apars_filename);
  Efake3D=new C3DArray(apars_filename);

  string fakedatadirname="blast/pp";
  string fakeerrordirname="blast/pp_errors";

  string wfparsfilename=pardirname+"/kparameters.dat";
  string kdatadirname="../../kdata/pp";

  kernel=new CKernel(wfparsfilename);
  //kernel->ReadData(kdatadirname);
  wf=new CWaveFunction_pp_schrod(wfparsfilename);
  kernel->Calc(wf);
  kernel->WriteData(kdatadirname);

  scalc=new CSourceCalc_Blast();
  scalc->SetSPars(lambda,R,Tau,DelTau,0.7,110.0,800.0,1.5,938.28,938.28);
  scalc->CalcS(S);
  scalc->NormCheck(S);
  S2CF::s2c(S,kernel,Cfake);
  ArrayCalc::Calc3DArrayFromAExpArray(Cfake,Ctheory3D);

  //Cfake3D->ReadArray(fakedatadirname);
  //Efake3D->ReadArray(fakeerrordirname);

  // Add random errors to Cfake3D
  Efake3D->RandomizeGaussian(error);
  ArrayCalc::AddArrays(Ctheory3D,Efake3D,Cfake3D);
  // Make error matrix and calculate chi square
  Efake3D->MakeConstant(error);
  chisquare=CFCalc::GetChiSquared(Cfake3D,Efake3D,Ctheory3D);
  printf("chisquare=%g\n",chisquare);
  
  // Write Arrays
  Cfake3D->WriteArray(fakedatadirname);
  Efake3D->WriteArray(fakeerrordirname);

  return 0;
  
}

