#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	char qchar[100];
	printf("Enter qualifer : ");
	scanf("%s",qchar);
	string qualifier=qchar;
	double kt[4]={200,300,400,500};
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	string dataroot="/Users/pratt/data/bb/";
	string parsdirname="parameters/";
	string dirname;
	char ktchars[10];

	CMCList *list=NULL;
	CCHArray *Asource=new CCHArray(parsdirname+"apars_sf.dat");
	CCHArray *cf=new CCHArray(parsdirname+"apars_cf.dat");
	//C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

	//CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/kparameters.dat");
	CWaveFunction_generic *wf=new CWaveFunction_generic(parsdirname+"kparameters.dat",0,139.58,139.58,1.0);
	CKernel *kernel=new CKernel(parsdirname+"kparameters.dat");
	kernel->Calc(wf);

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR(parsdirname+"spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",dataroot+qualifier+"/oscar.dat");
	scalc->SetIDs(idlista,nida,idlistb,nidb);


	dataroot="/Users/pratt/data/";
	for(int ikt=1;ikt<2;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		//list->PrintMoments(60.0);

		scalc->CombineMCLists(list,list,Asource,10000);
		scalc->NormCheck(Asource);
		Asource->FillRemainderX();
		Asource->PrintProjections();
		Asource->WriteProjections("sf_kt300.dat");
		scalc->CalcEffGaussPars(Asource);

		S2CF::s2c(Asource,kernel,cf);
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/oscar_"+qualifier+"/kt"+ktchars;
		cf->WriteAX(dirname);
		cf->FillRemainderX();
		cf->PrintProjections();
		cf->WriteProjections("cf_kt300.dat");


	}
}

