#include "coral.h"

int main (int argc, char *argv[]){
  int lx,ly,lz;

  CCF2SFit fitter;
  string filename;
  string dirname;
  string pardirname;
  double chisquare;

  pardirname="parameters/pp";

	filename=pardirname+"/aparsCH_source.dat";
  fitter.sourceCH=new CCHArray(filename);

	filename=pardirname+"/aparsCH_CF.dat";
  fitter.ctheoryCH=new CCHArray(filename);
  fitter.cexpCH=new CCHArray(filename);
	printf("test a\n");

  dirname="fakedata/GX1d/pp";
  fitter.cexpCH->ReadAX(dirname);
  fitter.cerrorCH=new CCHArray(filename);
  dirname="fakedata/GX1d/pp_errors";
  fitter.cerrorCH->ReadAX(dirname);

	filename=pardirname+"/kparameters.dat";
  fitter.kernel=new CKernel(filename);
  dirname="../kdata/pp";
  fitter.kernel->ReadData(dirname);
	printf("test d\n");
  
  fitter.sourcecalc=new CSourceCalc_GX1D();

  lx=ly=lz=0;
  fitter.SetL(lx,ly,lz);
  fitter.AddPar("lambda",0.5,0.02,0.0,1.0);
  fitter.AddPar("Xfrac",0.4,0.02,0.0,1.0);
  fitter.AddPar("R",4,0.1,1.0,25.0);
  fitter.AddPar("X",10.0,0.1,0.0,30.0);
  fitter.AddPar("a",5.0,0.02,0.0,2.0);
  fitter.FixPar("a");
  fitter.FixPar("X");
  fitter.SetCalcFlag(1);
  fitter.SetMCSourceFlag(0);

  // Find minimum
  fitter.PrintPars();
  fitter.SteepestDescent(5);
  fitter.Newton(5);
  fitter.PrintPars();
  fitter.PrintErrorMatrix();

  fitter.Metropolis(400);
  fitter.UpdateStepMatrix();
  fitter.Metropolis(400);
  fitter.PrintPars();
  fitter.PrintErrorMatrix();

  return 0;
}
