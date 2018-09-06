#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
	int idlist[3]={211,-211,111};
	int nid=3;
	string dataroot="/Users/pratt/data/bb/";
	string parsdirname="parameters/";

	CMCList *list=NULL;
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR(parsdirname+"spars_OSCAR.dat");
	parameter::set(scalc->spars,"OSCARfilename",dataroot+"default_short/oscar.dat");
	scalc->SetIDs(idlist,nid,idlist,nid);

	parameter::set(scalc->spars,"PT",300.0);
	scalc->CalcS(list,list);
}

