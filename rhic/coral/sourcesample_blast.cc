#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"

int main(){
  int nsample;
  CCHArray *Asource,*Xsource;
  Asource=new CCHArray("parameters/apars_blast.dat");

  CSourceCalc_Blast *scalc=new CSourceCalc_Blast();
  //printf("Enter nsample (NMC=nsample^2): ");
  //scanf("%d",&nsample);
  nsample=1000;
  parameter::set(scalc->spars,"Nsample",nsample);
  parameter::set(scalc->spars,"etaG",1.0);
  parameter::set(scalc->spars,"Ma",139.59);
  parameter::set(scalc->spars,"Mb",139.59);
  parameter::set(scalc->spars,"DelTau",5.0);
  parameter::set(scalc->spars,"Pt",100.0);
  parameter::PrintPars(scalc->spars);

  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

}

