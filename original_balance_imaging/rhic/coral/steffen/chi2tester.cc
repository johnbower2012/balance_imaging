#include "coral.h"

int main (int argc, char *argv[]){
  CCF2SFit fitter;
	int iq;
	char dataroot[120],dirname[120];
	sprintf(dataroot,"/Users/pratt/data");
	sprintf(dirname,"%s/cfdata/gaussCH_Rx3_Ry4_Rz5_lambda1/3d",dataroot);

  fitter.sourceCH=new CCHArray("parameters/apars_cf.dat");
  fitter.ctheoryCH=new CCHArray("parameters/apars_cf.dat");
  string filename="parameters/apars3D_cf.dat";
  fitter.ctheory3D=new C3DArray(filename);
  fitter.cexp3D=new C3DArray(filename);
  fitter.cerror3D=new C3DArray(filename);
	fitter.cexpCH=new CCHArray("parameters/apars_cf.dat");
	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	//CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);

//fitter.cexpCH->ReadAX("/Users/pratt/data/cfdata/fakedata_kt250");
//fitter.cexpCH->FillRemainderX();
//ArrayCalc::Calc3DArrayFromAExpArray(fitter.cexpCH,fitter.cexp3D);

	fitter.cexp3D->ReadArray(dirname);
	
	fitter.cexp3D->PrintProjections();
	fitter.cerror3D->MakeConstant(0.001);
	int ix,iy,iz,nmax;
	nmax=fitter.cerror3D->GetNXMAX();            
	for(ix=0;ix<nmax;ix++){           
		for(iy=0;iy<nmax;iy++){
			for(iz=0;iz<nmax;iz++){  
				fitter.cerror3D->SetElement(ix,iy,iz,exp(double(ix*ix+iy*iy+iz*iz)/1000.0));  
			}
		}                                                                                            
	}
	fitter.kernel=new CKernel("parameters/kparameters.dat");
	fitter.kernel->Calc(wf);

	fitter.sourcecalc=new CSourceCalc_Gaussian();
	
	fitter.AddPar("lambda",1.0,0.05,0.6,1.4);
	fitter.AddPar("Rx",3,0.25,2.0,20.0);
	fitter.AddPar("Ry",4,0.25,2.0,2.0);
	fitter.AddPar("Rz",5,0.25,2.0,20.0);
  fitter.AddPar("Xoff",0.0,0.0,0.0,15.0);
  fitter.FixPar("Xoff");
	//fitter.FixPar("lambda");
	
	parameter::PrintPars(fitter.sourcecalc->spars);
  fitter.sourcecalc->CalcS(fitter.sourceCH);
	fitter.sourceCH->FillRemainderX();
  fitter.sourcecalc->NormCheck(fitter.sourceCH);
	S2CF::s2c(fitter.sourceCH,fitter.kernel,fitter.ctheoryCH);
	fitter.ctheoryCH->FillRemainderX();
	ArrayCalc::Calc3DArrayFromAExpArray(fitter.ctheoryCH,fitter.ctheory3D);

	CFCalc::GetChiSquared(fitter.cexp3D,fitter.cerror3D,fitter.ctheory3D);
	
	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
	
	fitter.ctheoryCH->PrintProjections();

	printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");

	fitter.cexp3D->PrintProjections();
	fitter.ctheory3D->PrintProjections();

  return 0;
}
