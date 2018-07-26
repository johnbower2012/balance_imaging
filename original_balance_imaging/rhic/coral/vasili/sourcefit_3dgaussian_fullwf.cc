#include "coral.h"

int main (int argc, char *argv[]){
	CCF2SFit fitter;
	char qchar[100],ktchars[10];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	int iq;
	double kt[10]={50,100,150,200,250,300,350,400,450,500};
	double Rout,Rside,Rlong;
	string dirname,dataroot="/Users/vasili/data/";
	string filename="results/results3d_"+qualifier+".dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"! kt  lambda  Rout  Routlab Rside   Rlong\n");

	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");

	filename="parameters/apars3D_cf.dat";
	fitter.ctheory3D=new C3DArray(filename);
	fitter.cexp3D=new C3DArray(filename);
	fitter.cerror3D=new C3DArray(filename);
	fitter.kernelwf=new CKernelWF("parameters/kparameters.dat");
	fitter.kernelwf->Calc(wf);
	fitter.source3D=new C3DArray("parameters/apars3d_sf.dat");

	fitter.cerror3D->MakeConstant(1.0);
	int ix,iy,iz;

	fitter.sourcecalc=new CSourceCalc_Gaussian();
	parameter::set(fitter.sourcecalc->spars,"NMC",1000);

	fitter.AddPar("lambda",0.22844,0.05,0.1,1.4);
	fitter.AddPar("Rx",6.9312,0.5,2.0,30.0);
	fitter.AddPar("Ry",6.2819,0.5,2.0,20.0);
	fitter.AddPar("Rz",13.083,0.5,2.0,20.0);
	fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
	fitter.FixPar("Xoff");
	fitter.SetCalcFlag(9); // Set flag to 9 for full WF fitting(slow)
	// or to 10 for Bowler-Sinyukov(fast). The BS is more like what experiments do

	for(int ikt=0;ikt<10;ikt++){
		fitter.ResetChiSquared();
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/oscar_"+qualifier+"_fullwf/kt"+ktchars+"_3d";
		printf("dirname=%s\n",dirname.c_str());
		fitter.cexp3D->ReadArray(dirname);
		printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
		fitter.cexp3D->PrintProjections();
		//fitter.SteepestDescent(2);
		//fitter.PrintPars();
		fitter.Newton(2);
		fitter.PrintPars();
		Rout=fitter.GetPar("Rx")*139.57/sqrt(kt[ikt]*kt[ikt]+139.57*139.57);
		Rside=fitter.GetPar("Ry");
		Rlong=fitter.GetPar("Rz");
		fprintf(fptr,"%6.1f %6.4f %6.3f %6.3f %6.3f %6.3f %6.3f\n", kt[ikt],fitter.GetPar("lambda"),fitter.GetPar("Rx"),Rout,Rside,Rlong,Rout/Rside);
		printf("+++++++++++   GAUSSIAN FITS  +++++++++++++++\n");
		fitter.ctheory3D->PrintProjections();
	}
	fclose(fptr);

	return 0;
}