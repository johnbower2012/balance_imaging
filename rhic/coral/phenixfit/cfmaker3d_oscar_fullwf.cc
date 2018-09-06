#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	string dataroot="/Users/pratt/data/";
	char qchar[100];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	//int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int idlista[2]={321,-321},idlistb[2]={321,-321};
	int nida=2,nidb=2;
	string dirname;
	double kt[9]={100,200,300,400,500,600,700,800,900};
	char ktchars[9];
	int ikt;

	CMCList *list=NULL;
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",1,493.677,493.677,1.0);
	//CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
	CKernelWF *kernel=new CKernelWF("parameters/kparameters.dat");
	kernel->Calc(wf);

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",dataroot+"bb/"+qualifier+"/oscar.dat");
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(ikt=0;ikt<9;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		S2CF::s2c(list,list,kernel,c3d,100000);
		//S2CF::s2c_bosons(list,c3d);
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/oscar_kk_"+qualifier+"_fullwf/kt"+ktchars+"_3d";
		c3d->WriteArray(dirname);
		c3d->PrintProjections();
	}
}

