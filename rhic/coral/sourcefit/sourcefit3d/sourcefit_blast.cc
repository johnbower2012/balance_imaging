#include "coral.h"

int main (int argc, char *argv[]){

  CCF2SFit fitter;
  string filename;
  string dirname;
  string pardirname="parameters/pp";

	filename=pardirname+"/apars3D_CF.dat\0";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);

  dirname="fakedata/blast/pp";
  fitter.cexp3D->ReadArray(dirname);
  dirname="fakedata/blast/pp_errors";
  fitter.cerror3D->ReadArray(dirname);

  filename=pardirname+"/kparameters.dat";
  fitter.kernel=new CKernel(filename);
  dirname="../kdata/pp";
  fitter.kernel->ReadData(dirname);

  fitter.lista=new CMCList(0);
  fitter.listb=fitter.lista;

  fitter.sourcecalc=new CSourceCalc_Blast();
  parameter::set(fitter.sourcecalc->spars,"Nsample",40);

  fitter.AddPar("lambda",0.6,0.02,0.1,2.0);
  fitter.AddPar("R",13.0,0.1,2.0,20.0);
  fitter.AddPar("Tau",12.0,0.1,2.0,25.0);
  fitter.AddPar("DelTau",5.0,0.1,0.0,15.0);
  fitter.SetMCSourceFlag(1);
  fitter.SetCalcFlag(3);

  // Find minimum
  fitter.PrintPars();

  //fitter.SteepestDescent(10);
  fitter.Metropolis(50);
  fitter.UpdateStepMatrix();
  fitter.Metropolis(50);

  fitter.PrintPars();

  return 0;
}
