#include "coral.h"

int main (int argc, char *argv[]){

  CCF2SFit fitter;
  string filename;
  string dirname;
  string pardirname="parameters/pipi";

  filename=pardirname+"/aparsCH_source.dat";
  fitter.sourceCH=new CCHArray(filename);

  filename=pardirname+"/aparsCH_CF.dat";
  fitter.ctheoryCH=new CCHArray(filename);

  filename=pardirname+"/apars3D_CF.dat";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);

  dirname="fakedata/3dgaussian/pipi";
  fitter.cexp3D->ReadArray(dirname);
  dirname="fakedata/3dgaussian/pipi_errors";
  fitter.cerror3D->ReadArray(dirname);

  filename=pardirname+"/kparameters.dat";
  fitter.kernel=new CKernel(filename);
  dirname="../kdata/pipi";
  fitter.kernel->ReadData(dirname);

  fitter.sourcecalc=new CSourceCalc_Gaussian();
	
  fitter.AddPar("lambda",0.6,0.2,0.1,2.0);
  fitter.AddPar("Rx",5.0,2.0,0.0,15.0);
  fitter.AddPar("Ry",5.0,2.0,0.0,15.0);
  fitter.AddPar("Rz",5.0,2.0,0.0,15.0);
  fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
  fitter.FixPar("Xoff");
	//fitter.FixPar("lambda");

  fitter.sourcecalc->CalcS(fitter.sourceCH);
  fitter.sourcecalc->NormCheck(fitter.sourceCH);

  fitter.SetCalcFlag(2);

  // Find minimum
  fitter.PrintPars();

	printf("_____________ BEGIN STEEPEST DESCENT ______________\n");
  fitter.SteepestDescent(50);

	//printf("_____________ BEGIN CONJUGATE GRADIENT ______________\n");
	//fitter.ConjugateGradient(5);

	//printf("_____________ BEGIN NEWTON ______________\n");
	//fitter.Newton(4);
	//fitter.UpdateStepMatrix();
	//fitter.Newton(4);
	//fitter.UpdateStepMatrix();
	//fitter.Newton(4);

	/*printf("______________ BEGIN METROPOLIS __________________\n");
	fitter.UpdateStepMatrix();
  fitter.Metropolis(100);
  fitter.UpdateStepMatrix();
  fitter.Metropolis(100); */

  fitter.PrintPars();

  for(int iq=0;iq<fitter.ctheoryCH->GetNRADIAL();iq++)
    printf("C(%d)=%g =? %g\n",iq,fitter.ctheoryCH->GetElement(0,0,0,iq));

  return 0;
}
