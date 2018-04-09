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
  CCHArray *Sgauss;
  CSourceCalc_Gaussian *scalc;
  CWaveFunction *wf;
  CKernel *kernel;
  double lambda=0.5,Rx=7.0,Ry=5.0,Rz=3.0,Xoff=0.0;
  double error=0.005,chisquare;
  string pardirname="../parameters/pipi";

  // Create Arrays
  string apars_filename=pardirname+"/aparsCH_source.dat";
  Sgauss=new CCHArray(apars_filename);
  apars_filename=pardirname+"/aparsCH_CF.dat";
  Cfake=new CCHArray(apars_filename);
  apars_filename=pardirname+"/apars3D_CF.dat";
  Cfake3D=new C3DArray(apars_filename);
  Ctheory3D=new C3DArray(apars_filename);
  Efake3D=new C3DArray(apars_filename);

  string fakedatadirname="3dgaussian/pipi";
  string fakeerrordirname="3dgaussian/pipi_errors";

  string wfparsfilename;
  wfparsfilename=pardirname+"/kparameters.dat";
  string kdatadirname="../../kdata/pipi";

  kernel=new CKernel(wfparsfilename);
  //kernel->ReadData(kdatadirname);
  wf=new CWaveFunction_pipluspiplus_sqwell(wfparsfilename);
  kernel->Calc(wf);
  kernel->WriteData(kdatadirname);

  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(lambda,Rx,Ry,Rz,Xoff,0.0,0.0);
  scalc->CalcS(Sgauss);
  scalc->NormCheck(Sgauss);
  S2CF::s2c(Sgauss,kernel,Cfake);
  ArrayCalc::Calc3DArrayFromAExpArray(Cfake,Ctheory3D);

  //Cfake3D->ReadArray(fakedatadirname);
  //Efake3D->ReadArray(fakeerrordirname);

  // Add random errors to Cfake3D
  Efake3D->RandomizeGaussian(error);
  //Efake3D->ScaleArray(error);
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

