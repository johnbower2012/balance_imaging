#include "coral.cc"
// I have compiler problems with g++ if I try to make minuit.o, then link to it,
// thus I simply include minuit.cc here to be compiled alongside the main prog.
#include "sfit_minuit.cc"

int main (int argc, char *argv[]){
  double *xcontour,*ycontour;
  int nfound_contour=0,npts_contour=20;
  xcontour=new double[npts_contour];
  ycontour=new double[npts_contour];

  CCF2S_Minuit *mninfo;
  CCHArray *Ctheory;
  CCHArray *Sgauss;
  C3DArray *Ctheory3D;
  C3DArray *Cdata3D;
  C3DArray *Cerror3D;
  CSourceCalc_Gaussian *scalc;
  CKernel *kernel;
  char filename[200];
  char dirname[200];
  double chisquare;
  char pardirname[100];
  sprintf(pardirname,"parameters/pp\0");

  sprintf(filename,"%s/aparsCH_source.dat\0",pardirname);
  Sgauss=new CCHArray(filename);
  sprintf(filename,"%s/aparsCH_CF.dat\0",pardirname);
  Ctheory=new CCHArray(filename);
  sprintf(filename,"%s/apars3D_CF.dat\0",pardirname);
  Ctheory3D=new C3DArray(filename);
  
  Cdata3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/3dgaussian/pp\0");
  Cdata3D->ReadArray(dirname);
  Cerror3D=new C3DArray(filename);
  sprintf(dirname,"fakedata/3dgaussian/pp_errors\0");
  Cerror3D->ReadArray(dirname);

  sprintf(filename,"%s/kparameters.dat\0",pardirname);
  kernel=new CKernel(filename);
  sprintf(dirname,"../kdata/pp\0");
  kernel->Read(dirname);
    
  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(0.5,7,5,3,0,0,0);
  scalc->CalcS(Sgauss);
  scalc->NormCheck(Sgauss);

  mninfo=new CCF2S_Minuit_3DGaussian(scalc,Cdata3D,Cerror3D,Ctheory3D,
				Ctheory,Sgauss,kernel);
  mninfo->FixPar(4);

  // Find minimum
  mninfo->ViewPars();
  mninfo->StratLevel(0);
  mninfo->Minimize();

  for(int iq=0;iq<20;iq++)
    printf("C(%d)=%g\n",iq,Ctheory->GetElement(0,0,0,iq));

  mninfo->ErrorMatrix();
  // Plot contour for 1st and 2nd parameter where fcn = fcnmin + 0.25
  //mninfo->Contour(1,2,npts_contour,xcontour,ycontour);
  //mninfo->Minos();

  return 0;
}
