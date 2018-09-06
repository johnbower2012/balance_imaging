#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	// this assumes OSCARfilename is of form hydro_qualifier.dat and is in directory 
	// dataroot/hydrodata/
	// You should set dataroot below to correct value. "qualifier" is read in.
	//It could be something like "b4" for impact parameter, b=4.
	// 3d CFs are stored in dataroot/cfdata/oscar_qualifier_fullwf/kt%g_3d/
	// The files store the 3-d info for a given ix and iy into separate files
	// It will also write a file 3Darraypars.dat into this directory for future reading
	string dataroot="/Users/vasili/data/";
	char qchar[100];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	string dirname,filename;
	double kt[6]={100,200,300,400,500,600};
	char ktchars[9];
	int ikt;

	CMCList *list=NULL;
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CKernelWF *kernel=new CKernelWF("parameters/kparameters.dat");
	kernel->Calc(wf);

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	filename=dataroot+"hydrodata/hydro_"+qualifier+".dat";
	parameter::set(scalc->spars,"OSCARfilename",filename);
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(ikt=0;ikt<6;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		S2CF::s2c(list,list,kernel,c3d,100000);
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/oscar_"+qualifier+"_fullwf/kt"+ktchars+"_3d";
		c3d->WriteArray(dirname);
		c3d->PrintProjections();
	}
}

