#include "coral.h"
using namespace std;

// This will read a source from file and calculate correlation function

int main(){
	CSourceCalc_OSCAR *scalc;
	scalc=new CSourceCalc_OSCAR("parameters/spars_OSCAR.dat");

	CCHArray *source=new CCHArray("parameters/apars_oscar_sf.dat");
	source->ReadAX("/Users/scottepratt/data/sdata/fakedata_kt250");
	source->FillRemainderX();
	source->PrintProjections();
	scalc->CalcEffGaussPars(source);
	
  CKernel *kernel=new CKernel("parameters/wfparameters.dat");
  kernel->ReadData("/Users/scottepratt/data/kdata/pipluspiplus");

	CCHArray *cf=new CCHArray("parameters/apars_oscar_cf.dat");
	S2CF::s2c(source,kernel,cf);
	cf->FillRemainderX();
	cf->Print(0,0,0);
	cf->WriteAX("/Users/scottepratt/data/cfdata/fakedata_kt250");
	cf->PrintProjections();
	
	C3DArray *c3d=new C3DArray("parameters/apars3d.dat");
	ArrayCalc::Calc3DArrayFromAExpArray(cf,c3d);
	c3d->WriteArray("/Users/scottepratt/data/cfdata/fakedata_kt250/3d");

}

