#include "coral.h"

int main (int argc, char *argv[]){
  CCF2SFit fitter;
	char qchar[100],ktchars[10];
	int iq;
	double kt[10]={50,100,150,200,250,300,350,400,450,500};
	double Rout,Rside,Rlong;
	string dirname,dataroot="/Users/pratt/data/";
	string filename="results/results3d_steffen.dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"! kt  lambda  Rout  Routlab Rside   Rlong\n");


  filename="parameters/apars3D_cf.dat";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	//CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);

	fitter.cerror3D->MakeConstant(0.001);
	int ix,iy,iz;
	/*for(ix=0;ix<2;ix++){
		for(iy=0;iy<2;iy++){
			for(iz=0;iz<2;iz++){
				fitter.cerror3D->SetElement(ix,iy,iz,1000.0);
			}
		}                                                                                            
	}*/

	fitter.sourcecalc=new CSourceCalc_Gaussian();

	fitter.AddPar("lambda",0.5,0.02,0.1,1.4);
	fitter.AddPar("Rx",6.0,0.1,2.0,30.0);
	fitter.AddPar("Ry",6.0,0.1,2.0,20.0);
	fitter.AddPar("Rz",6.0,0.1,2.0,20.0);
  fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
  fitter.FixPar("Xoff");
  fitter.SetCalcFlag(8); // for pure boson wf
  
	for(int ikt=0;ikt<10;ikt++){
		fitter.ResetChiSquared();
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/steffen/kt"+ktchars+"_3d";
		printf("dirname=%s\n",dirname.c_str());
		fitter.cexp3D->ReadArray(dirname);
		printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
		fitter.cexp3D->PrintProjections();
		fitter.Newton(6);
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
