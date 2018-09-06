#include "coral.h"

int main (int argc, char *argv[]){
	CCF2SFit fitter;
	char qchar[100],ktchars[10];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	int iq;
	double kt[9]={100,200,300,400,500,600,700,800,900};
	double Rout,Rside,Rlong,mass=493.677;
	string dirname,dataroot="/Users/pratt/data/";
	string filename="results/results3d_"+qualifier+".dat";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"! kt  lambda  Rout  Routlab Rside   Rlong\n");

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",1.0,mass,mass,1.0);

	filename="parameters/apars3D_cf.dat";
	fitter.ctheory3D=new C3DArray(filename);
	fitter.cexp3D=new C3DArray(filename);
	fitter.cerror3D=new C3DArray(filename);
	fitter.kernelwf=new CKernelWF("parameters/kparameters.dat");
	fitter.kernelwf->Calc(wf);
	fitter.source3D=new C3DArray("parameters/apars3d_sf.dat");

	fitter.cerror3D->MakeConstant(1.0);
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
	parameter::set(fitter.sourcecalc->spars,"NMC",1000);

	/*fitter.AddPar("lambda",0.8767,0.05,0.1,1.4);
	fitter.AddPar("Rx",4.779,0.25,2.0,30.0);
	fitter.AddPar("Ry",4.324,0.25,2.0,20.0);
	fitter.AddPar("Rz",5.545,0.25,2.0,20.0);*/
	fitter.AddPar("lambda",1,0.05,0.1,1.4);
	fitter.AddPar("Rx",3,0.25,2.0,30.0);
	fitter.AddPar("Ry",5,0.25,2.0,20.0);
	fitter.AddPar("Rz",7,0.25,2.0,20.0);
	fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
	fitter.FixPar("Xoff");
	fitter.SetCalcFlag(9); 

	for(int ikt=0;ikt<1;ikt++){
		fitter.ResetChiSquared();
		sprintf(ktchars,"%g",kt[ikt]);
		//dirname=dataroot+"cfdata/oscar_"+qualifier+"/kt"+ktchars+"_3d";
		dirname=dataroot+"cfdata/"+qualifier+"/3d";
		printf("dirname=%s\n",dirname.c_str());
		fitter.cexp3D->ReadArray(dirname);
		printf("++++++++++++ C(q) from S(x) +++++++++++++++\n");
		fitter.cexp3D->PrintProjections();
		//fitter.SteepestDescent(2);
		//fitter.PrintPars();
		//fitter.SteepestDescent(4);
		//fitter.PrintPars();
		fitter.Newton(2);
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
