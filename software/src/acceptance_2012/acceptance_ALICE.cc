	//
	//  acceptance_CHEAP.cc
	//

#include "bal.h"
#include "TALICEAcceptanceCint.h"

using namespace std;

void SetupParameterization(const char* filename);
void SetCentrality(Int_t iCentrality);
Double_t IsAccepted(Int_t gPdgCode, Double_t gPt, Double_t gRapidity, Double_t gPhi, double *dca);

CAcceptance_ALICE::CAcceptance_ALICE() : CAcceptance(){
	ETAMIN=-0.8; // Don't bother calling Acceptance Routine if outside these boundaries.
	ETAMAX=0.8;
	ptmin=200.0; 
	ptmax=2000.0;
	alice=new TALICEAcceptance();
	filename="../alice/efficiencyALICE.root";
	alice->SetupParameterization(filename.c_str());
	alice->SetCentrality(0);
}

void CAcceptance_ALICE::CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca){
	double eta,pt,pmag,gammav,phi,*p=part->p,y=part->y;
	int pid=part->resinfo->code;
	
	/*
	accept=false;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	if(eta>ETAMIN && eta<ETAMAX){
		if(pt>ptmin && pt<ptmax){
			accept=true;
			efficiency=0.8;
		}
	}*/
	
		
	accept=false;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	if(eta>ETAMIN && eta<ETAMAX){
		if(pt>ptmin && pt<ptmax){
			accept=true;
			phi=atan2(part->p[2],part->p[1]);
			efficiency=alice->IsAccepted(pid,0.001*pt,y,phi,dca);
		}
	}
}