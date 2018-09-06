//
//  acceptance.cc
//  

#include "bal.h"
using namespace std;

CAcceptance::CAcceptance(){
	ETAMIN=-1.0;
	ETAMAX=1.0;
}

void CAcceptance::CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca){
  accept=true;
	efficiency=1.0;
}