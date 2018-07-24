#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf=new CWaveFunction_pp_schrod("parameters/wfparameters.dat");;
  CKernel *kernel=new CKernel("parameters/wfparameters.dat");
  kernel->Calc(wf);

	CCHArray *source=new CCHArray("parameters/apars_sf.dat");
	string dirname="/usr/users/micha/mysourcedata";
	// you need file lx0_ly0_lz0.tmp
	source.ReadArray(dirname.c_str());
	
	CCHArray *cf=new CCHArray("parameters/apars_cf.dat");
	
	S2CF::s2c(source,kernel,cf);
	cf->Print(0,0,0);
	cf->WriteAX("/usr/users/micha/mycfdata");

}

