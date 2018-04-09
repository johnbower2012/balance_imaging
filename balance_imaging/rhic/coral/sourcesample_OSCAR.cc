#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	char sdirname[120],cfdirname[120];
	double kt;
	CMCList *list=NULL;
	// Create Array for storing source info
	CCHArray *Asource=new CCHArray("parameters/apars_oscar_sf.dat");

	CCHArray *cf=new CCHArray("parameters/apars_oscar_cf.dat");

	CWaveFunction_pipluspiplus_sqwell *wf=new CWaveFunction_pipluspiplus_sqwell("parameters/wfparameters.dat");
	CKernel *kernel=new CKernel("parameters/wfparameters.dat");
  kernel->Calc(wf);
	//kernel->ReadData("/Users/scottepratt/data/kdata/pipluspiplus");

	// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(double kt=250.0; kt<251; kt+=100){
		printf("____________ Calculating S for kt = %g __________________\n",kt);
		parameter::set(scalc->spars,"PT",2.0*kt);
	// Calculate source array
		scalc->CalcS(list,list);
		scalc->CombineMCLists(list,list,Asource);
		scalc->NormCheck(Asource);
		Asource->FillRemainderX();
		scalc->CalcEffGaussPars(Asource);  // fits l=0,1,2 moments of source function

  // Write source info to file
		sprintf(sdirname,"/Users/scottepratt/data/sdata/fakedata_kt%g",kt);
		Asource->WriteAX(sdirname);
		Asource->PrintProjections();

		S2CF::s2c(Asource,kernel,cf);
		cf->FillRemainderX();
		//cf->Print(0,0,0);
		//sprintf(cfdirname,"/Users/scottepratt/data/cfdata/OSCAR_pionsonly_kt%g",kt);
		//cf->WriteAllA(cfdirname);
		cf->PrintProjections();
	}
}

