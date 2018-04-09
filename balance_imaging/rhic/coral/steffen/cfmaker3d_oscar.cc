#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	string dataroot="/Users/pratt/data/";
	int idlista[3]={211,-211,111},idlistb[3]={211,-211,111};
	int nida=3,nidb=3;
	string dirname;
	double kt[10]={50,100,150,200,250,300,350,400,450,500};
	char ktchars[10];
	int ikt;

	CMCList *list=NULL;
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",dataroot+"steffen/hydro_input_a.f19");
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(ikt=0;ikt<10;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		//S2CF::s2c(list,list,kernel,c3d,1000000);
		S2CF::s2c_bosons(list,c3d);
		sprintf(ktchars,"%g",kt[ikt]);
		dirname=dataroot+"cfdata/steffen/kt"+ktchars+"_3d";
		c3d->WriteArray(dirname);
		c3d->PrintProjections();
	}
}

