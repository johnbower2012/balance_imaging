#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf;
  CKernel *kernel;

  string wfparsfilename="parameters/pipi/kparameters.dat";
	string kdatadirname="/Users/scottepratt/data/kdata/pipluspiplus";

  kernel=new CKernel(wfparsfilename);

  // Either read or calc
  //kernel->Read(kdatadirname);
  wf=new CWaveFunction_pipluspiplus_sqwell(wfparsfilename);
  kernel->Calc(wf);
  //kernel->Calc_ClassCoul(938.28,493.677,1);

  kernel->WriteData(kdatadirname);

}

