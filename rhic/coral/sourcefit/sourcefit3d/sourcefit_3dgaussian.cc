#include "coral.h"

int main (int argc, char *argv[]){
  CCF2SFit fitter;
	int iq;

  fitter.sourceCH=new CCHArray("parameters/pipi/aparsCH_source.dat");
  fitter.ctheoryCH=new CCHArray("parameters/pipi/aparsCH_CF.dat");
  string filename="parameters/pipi/apars3D_CF.dat";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);
	fitter.cexpCH=new CCHArray("parameters/pipi/aparsCH_CF.dat");

  fitter.cexpCH->ReadAX("/Users/scottepratt/data/cfdata/fakedata_kt250");
	fitter.cexpCH->FillRemainderX();
	ArrayCalc::Calc3DArrayFromAExpArray(fitter.cexpCH,fitter.cexp3D);
  fitter.cerror3D->MakeConstant(0.001);
	int ix,iy,iz;
	for(ix=0;ix<3;ix++){
		for(iy=0;iy<3;iy++){
			for(iz=0;iz<3;iz++){
				fitter.cerror3D->SetElement(ix,iy,iz,10.0);
			}
		}
	}
	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/pipi/kparameters.dat");
	fitter.kernel=new CKernel("parameters/pipi/kparameters.dat");
  //fitter.kernel->ReadData("/Users/scottepratt/data/kdata/pipluspiplus");
	fitter.kernel->Calc(wf);

  fitter.sourcecalc=new CSourceCalc_Gaussian();
	
  fitter.AddPar("lambda",0.31,0.01,0.6,1.4);
  fitter.AddPar("Rx",10.42,0.05,4.0,20.0);
  fitter.AddPar("Ry",4.00,0.05,2.0,8.0);
  fitter.AddPar("Rz",3.15,0.05,2.0,15.0);
  fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
  fitter.FixPar("Xoff");
	//fitter.FixPar("lambda");

  fitter.sourcecalc->CalcS(fitter.sourceCH);
  fitter.sourcecalc->NormCheck(fitter.sourceCH);

  fitter.SetCalcFlag(2);

  // Find minimum
  fitter.PrintPars();

	fitter.SteepestDescent(6);

	//fitter.Newton(1);
	/*fitter.UpdateStepMatrix();
	fitter.Newton(4);
	fitter.UpdateStepMatrix();
	fitter.Newton(4);
	*/
	
	printf("___________________________________________\n");
	//fitter.UpdateStepMatrix();
	//fitter.SteepestDescent(3);
	//printf("FINISHED STEEPEST DESCENT : NOW TRY CG\n");
	//fitter.ConjugateGradient(5);

	/* fitter.Metropolis(100);
  fitter.UpdateStepMatrix();
  fitter.Metropolis(100);
	fitter.UpdateStepMatrix();
  fitter.Metropolis(100);
	*/

	fitter.PrintPars();
	fitter.ctheoryCH->FillRemainderX();
	printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
	fitter.cexpCH->PrintProjections();
	printf("+++++++++++   GAUSSIAN FIT  +++++++++++++++\n");
	fitter.ctheoryCH->PrintProjections();
	

  return 0;
}
