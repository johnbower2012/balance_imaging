#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>

using namespace std;

#include "coral.h"

int main(){
  CWaveFunction *wfpp;
  CKernel *kernelpp;

  wfpp=new CWaveFunction_kpluspiplus_sqwell("parameters/wfparameters.dat");
  kernelpp=new CKernel("parameters/wfparameters.dat");
  kernelpp->Calc(wfpp);
  kernelpp->WriteData("kernel/pipluskplus");
	
  CWaveFunction *wfpm;
  CKernel *kernelpm;

  wfpm=new CWaveFunction_kpluspiminus_sqwell("parameters/wfparameters.dat");
  kernelpm=new CKernel("parameters/wfparameters.dat");
  kernelpm->Calc(wfpm);
  kernelpm->WriteData("kernel/pipluskminus");
	
	exit(1);

}

