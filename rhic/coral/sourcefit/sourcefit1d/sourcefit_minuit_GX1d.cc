#include "coral.h"
// I have compiler problems with g++ if I try to make minuit.o, then link to it,
// thus I simply include minuit.cc here to be compiled alongside the main prog.
#include "sfit_minuit.cc"

int main (int argc, char *argv[]){
  int lx,ly,lz;
  double *xcontour,*ycontour;
  int nfound_contour=0,npts_contour=20;
  xcontour=new double[npts_contour];
  ycontour=new double[npts_contour];

  CCF2S_Minuit_GX1D *mninfo;
  CCHArray *Ctheory;
  CCHArray *S;
  CCHArray *Cdata;
  CCHArray *Cerror;
  CSourceCalc_GX1D *scalc;
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
  Cdata=new CCHArray(filename);
  sprintf(dirname,"fakedata/GX1d/pp\0");
  Cdata->ReadAX(dirname);
  Cerror=new CCHArray(filename);
  sprintf(dirname,"fakedata/GX1d/pp_errors\0");
  Cerror->ReadAX(dirname);

  sprintf(filename,"%s/kparameters.dat\0",pardirname);
  kernel=new CKernel(filename);
  sprintf(dirname,"../kdata/pp\0");
  kernel->Read(dirname);
  
  lx=ly=lz=0;
  scalc=new CSourceCalc_GX1D();
  scalc->SetSPars(0.6,0.5,5.0,9.0,5.0);
  scalc->CalcS(S);
  scalc->NormCheck(S);

  mninfo=new CCF2S_Minuit_GX1D(scalc,Cdata,Cerror,Ctheory,S,kernel);

  mninfo->lx=lx; mninfo->ly=ly; mninfo->lz=lz;
  mninfo->FixPar(3);
  mninfo->FixPar(1);

  // Find minimum
  printf("Beginning Parameters : \n");
  mninfo->ViewPars();
  mninfo->StratLevel(1);
  mninfo->Minimize();

  mninfo->ErrorMatrix();
  // Plot contour for 1st and 2nd parameter where fcn = fcnmin + 0.25
  mninfo->Contour(2,3,npts_contour,xcontour,ycontour);
  //mninfo->Minos();
  //mninfo->Simplex();

  return 0;
}
