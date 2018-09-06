#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	string dirname,filename,ktstring;
	double kt[6]={100,200,300,400,500,600};
	char ktchars[9];
	int ikt;

	CMCList *list=NULL;
	CCHArray *Asource=new CCHArray("parameters/apars_sf.dat");
	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic("parameters/kparameters.dat",0,139.58,139.58,1.0);
	CKernel *kernel=new CKernel("parameters/kparameters.dat");
	kernel->Calc(wf);

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",string("oscardata/bb_oscar.dat"));
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(int ikt=0;ikt<6;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		scalc->CombineMCLists(list,list,Asource);
		scalc->NormCheck(Asource);
		Asource->FillRemainderX();
		//Asource->PrintProjections();
		//scalc->CalcEffGaussPars(Asource);

		S2CF::s2c(Asource,kernel,cf);
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dirname="results/kt"+ktstring;
		cf->WriteAX(dirname);
		cf->FillRemainderX();
		cf->PrintProjections();

	}
}

