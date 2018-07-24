#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: cfmaker3d_oscar qualifier\n");
		exit(-1);
	}
	string dataroot="/Users/pratt/data/llnl_hydro/";
	string qualifier=argv[1];
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	string dirname,filename;
	double kt[10]={50,100,150,200,250,300,350,400,450,500};
	//double kt[10]={100,200,300,400,500};
	char ktchars[5];
	int ikt;

	CMCList *list=NULL;
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	//CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
	//CKernelWF *kernel=new CKernelWF("parameters/kparameters.dat");
	//kernel->Calc(wf);

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",dataroot+qualifier+"/bb_oscar.dat");
	//parameter::set(scalc->spars,"OSCARfilename",dataroot+qualifier+"/hydro_oscar.dat");
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(ikt=0;ikt<10;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		//S2CF::s2c(list,list,kernel,c3d,1000000);
		S2CF::s2c_bosons(list,c3d);
		sprintf(ktchars,"%g",kt[ikt]);
		//dirname=dataroot+qualifier+"/cfdata_nocascade_bosons/kt"+ktchars+"_3d";
		//dirname=dataroot+qualifier+"/cfdata_bosons/kt"+ktchars+"_3d";
		filename=dataroot+qualifier+"/cfdata_bosons/kt"+ktchars+"_3d.dat";
		//c3d->WriteArray(dirname);
		c3d->WriteArray_PHENIX(filename);
		c3d->PrintProjections();
	}
}

