#include "coral.h"

int main (int argc, char *argv[]){
  CCF2SFit fitter;
	char qchar[100],ktchars[10];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	double mass=493.677;
	//double mass=139.58;
	int iq;
	double kt[9]={100,200,300,400,500,600,700,800,900};
	double Rout,Rside,Rlong;
	string dirname,dataroot="/Users/pratt/data/";
	string filename="results/results3d_bowlersinyukov_"+qualifier+".dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"! kt  lambda  Rout  Routlab Rside   Rlong\n");


  filename="parameters/apars3D_cf.dat";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters_bowlersinyukov.dat",1.0,mass,mass,0.5);
	CKernel *kernel=new CKernel("parameters/kparameters_bowlersinyukov.dat");
	kernel->Calc(wf);

	fitter.cerror3D->MakeConstant(0.1);
	int ix,iy,iz;

	/*
	for(ix=0;ix<2;ix++){
		for(iy=0;iy<2;iy++){
			for(iz=0;iz<2;iz++){
				fitter.cerror3D->SetElement(ix,iy,iz,1000.0);
			}
		}                                                                                            
	}
	*/

	fitter.sourcecalc=new CSourceCalc_Gaussian();
	fitter.kernel=kernel;

	fitter.AddPar("lambda",1.0,0.05,0.1,1.4);
	fitter.AddPar("Rx",3.0,0.25,2.0,30.0);
	fitter.AddPar("Ry",5.0,0.25,2.0,15.0);
	fitter.AddPar("Rz",7.0,0.25,2.0,25.0);
	fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
	fitter.FixPar("Xoff");
	fitter.SetCalcFlag(10); // for BS

	for(int ikt=0;ikt<1;ikt++){
		fitter.ResetChiSquared();
		sprintf(ktchars,"%g",kt[ikt]);
		//dirname=dataroot+"cfdata/oscar_"+qualifier+"/kt"+ktchars+"_3d";
		dirname=dataroot+"cfdata/"+qualifier+"/3d";
		printf("dirname=%s\n",dirname.c_str());
		fitter.cexp3D->ReadArray(dirname);
		printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
		fitter.cexp3D->PrintProjections();
		fitter.SteepestDescent(2);
		fitter.Newton(4);
		fitter.PrintPars();
		Rout=fitter.GetPar("Rx")*mass/sqrt(kt[ikt]*kt[ikt]+mass*mass);
		Rside=fitter.GetPar("Ry");
		Rlong=fitter.GetPar("Rz");
		fprintf(fptr,"%6.1f %6.4f %6.3f %6.3f %6.3f %6.3f %6.3f\n", kt[ikt],fitter.GetPar("lambda"),fitter.GetPar("Rx"),Rout,Rside,Rlong,Rout/Rside);
		printf("+++++++++++   GAUSSIAN FITS  +++++++++++++++\n");
		fitter.ctheory3D->PrintProjections();
	}
	fclose(fptr);

  return 0;
}
