#include "coral.h"

int main (int argc, char *argv[]){
  double lambda=0.6,R=13.0,Tau=12.0,DelTau=5.0;
  double *xcontour,*ycontour;
  int nfound_contour=0,npts_contour=20;
  xcontour=new double[npts_contour];
  ycontour=new double[npts_contour];

  CCF2SFit *fitter;
  CCHArray *Ctheory;
  CCHArray *S;
  C3DArray *Ctheory3D;
  C3DArray *Cdata3D;
  C3DArray *Cerror3D;
  CSourceCalc_Blast *scalc;
  CKernel *kernel;
  char filename[200];
  char dirname[200];
  double chisquare;
  char pardirname[100];
  sprintf(pardirname,"parameters/pp\0");

  sprintf(filename,"%s/aparsCH_source.dat\0",pardirname);
  S=new CCHArray(filename);
  sprintf(filename,"%s/aparsCH_CF.dat\0",pardirname);
  Ctheory=new CCHArray(filename);
  sprintf(filename,"%s/apars3D_CF.dat\0",pardirname);
  Ctheory3D=new C3DArray(filename);
  
  Cdata3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/blast/pp\0");
  Cdata3D->ReadArray(dirname);
  Cerror3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/blast/pp_errors\0");
  Cerror3D->ReadArray(dirname);

  sprintf(filename,"%s/kparameters.dat\0",pardirname);
  kernel=new CKernel(filename);
  sprintf(dirname,"../kdata/pp\0");
  kernel->Read(dirname);
    
  scalc=new CSourceCalc_Blast();
  scalc->SetSPars(lambda,R,Tau,DelTau);
  parameter::set(scalc->spars,"Nsample",100);
  scalc->CalcS(S);
  scalc->NormCheck(S);

  fitter=new CCF2SFit_Blast(scalc,Cdata3D,Cerror3D,Ctheory3D,
				      Ctheory,S,kernel);
  //fitter->FixPar(4);

  // Find minimum
  fitter->PrintPars();

  fitter->Metropolis(20);
  fitter->UpdateStepMatrix();
  fitter->PrintPars();
  fitter->PrintStepMatrix();

  fitter->Metropolis(40);
  fitter->UpdateStepMatrix();
  fitter->Metropolis(60);
  fitter->UpdateStepMatrix();

  fitter->PrintPars();

  for(int iq=0;iq<20;iq++)
    printf("C(%d)=%g\n",iq,Ctheory->GetElement(0,0,0,iq));

  return 0;
}
