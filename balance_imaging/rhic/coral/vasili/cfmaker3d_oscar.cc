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
	C3DArray *c3d=new C3DArray("parameters/apars3d_cf.dat");

// Initialize Source Calc Object
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");
	filename="oscardata/bb_oscar.dat";
	parameter::set(scalc->spars,"OSCARfilename",filename);
	parameter::PrintPars(scalc->spars);
	scalc->SetIDs(idlista,nida,idlistb,nidb);

	for(ikt=0;ikt<6;ikt++){
		printf("____________ Calculating S for kt = %g __________________\n",kt[ikt]);
		parameter::set(scalc->spars,"PT",2.0*kt[ikt]);
		scalc->CalcS(list,list);
		S2CF::s2c_bosons(list,c3d);
		sprintf(ktchars,"%g",kt[ikt]);
		ktstring=ktchars;
		dirname="results/kt"+ktstring+"_3d";
		c3d->WriteArray(dirname);
		c3d->PrintProjections();
	}
}


