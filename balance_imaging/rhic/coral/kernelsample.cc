#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wf;
  CKernel *kernel;

  wf=new CWaveFunction_pipluspiplus_sqwell("parameters/wfparameters.dat");
  kernel=new CKernel("parameters/wfparameters.dat");
  kernel->Calc(wf);
  kernel->WriteData("/Users/scottepratt/data/kdata/pipluspiplus");

}

