#include "coral.h"

int main (int argc, char *argv[]){
	string parsdirname="parameters/";
	char qchar[100];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	string dataroot="/Users/pratt/data/cfdata/";
	string dirname[4];
	double kt[4]={200,300,400,500};
	dirname[0]=dataroot+"oscar_"+qualifier+"/kt200";
	dirname[1]=dataroot+"oscar_"+qualifier+"/kt300";
	dirname[2]=dataroot+"oscar_"+qualifier+"/kt400";
	dirname[3]=dataroot+"oscar_"+qualifier+"/kt500";
	FILE *fptr=fopen("results.dat","w");
	fprintf(fptr,"! kt  lambda  Rout  Routlab Rside   Rlong\n");
	CCF2SFit fitter;
	int iq;
  fitter.sourceCH=new CCHArray(parsdirname+"apars_cf.dat");
  fitter.ctheoryCH=new CCHArray(parsdirname+"apars_cf.dat");
  fitter.cexpCH=new CCHArray(parsdirname+"apars_cf.dat");
  fitter.cerrorCH=new CCHArray(fitter.cexpCH->GetLMAX(),fitter.cexpCH->GetNRADIAL(),fitter.cexpCH->GetRADSTEP(),fitter.cexpCH->GetXSYM(),fitter.cexpCH->GetYSYM(),fitter.cexpCH->GetZSYM());

	for(int iq=0;iq<fitter.cerrorCH->GetNRADIAL();iq++) fitter.cerrorCH->SetElement(0,0,0,iq,0.1);
	fitter.cerrorCH->SetElement(0,0,0,0,1000.0);
	fitter.cerrorCH->SetElement(0,0,0,1,1000.0);
	fitter.cerrorCH->SetElement(0,0,0,2,1000.0);
	
	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);

	fitter.kernel=new CKernel("parameters/kparameters.dat");
	fitter.kernel->Calc(wf);

	fitter.sourcecalc=new CSourceCalc_Gaussian();

	fitter.AddPar("lambda",0.65,0.1,0.1,1.4);
	fitter.AddPar("Rx",5.0,0.3,2.0,30.0);
	fitter.AddPar("Ry",5.0,0.3,1.0,20.0);
	fitter.AddPar("Rz",5.0,0.3,1.0,20.0);
  fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
  fitter.FixPar("Xoff");

  fitter.sourcecalc->CalcS(fitter.sourceCH);
  fitter.sourcecalc->NormCheck(fitter.sourceCH);

  fitter.SetCalcFlag(7);

	for(int ikt=0;ikt<4;ikt++){
		fitter.ResetChiSquared();
		fitter.cexpCH->ReadAX(dirname[ikt]);
		fitter.cexpCH->FillRemainderX();
		//fitter.SteepestDescent(2);
		fitter.Newton(8);
		//fitter.Metropolis(500);
		fitter.PrintPars();
		printf("%6.1f %6.4f %6.3f %6.3f %6.3f %6.3f\n", kt[ikt],fitter.GetPar("lambda"),fitter.GetPar("Rx"),	fitter.GetPar("Rx")*139.57/sqrt(kt[ikt]*kt[ikt]+139.57*139.57),fitter.GetPar("Ry"),fitter.GetPar("Rz"));
		fprintf(fptr,"%6.1f %6.4f %6.3f %6.3f %6.3f %6.3f\n", kt[ikt],fitter.GetPar("lambda"),fitter.GetPar("Rx"),fitter.GetPar("Rx")*139.57/sqrt(kt[ikt]*kt[ikt]+139.57*139.57),fitter.GetPar("Ry"),fitter.GetPar("Rz"));
		printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
		fitter.cexpCH->PrintProjections();
		fitter.cexpCH->PrintMoments();
		printf("+++++++++++   GAUSSIAN FITS  +++++++++++++++\n");
		fitter.ctheoryCH->PrintProjections();
	}
	fclose(fptr);
  return 0;
}
