	//
	//  acceptance_CHEAP.cc
	//

#include "bal.h"
#include "TALICEAcceptanceCint.h"

using namespace std;

CAcceptance_STAR::CAcceptance_STAR() : CAcceptance(){
	ETAMIN=-1.0; // Don't bother calling Acceptance Routine if outside these boundaries.
	ETAMAX=1.0;
	ptmin=150.0; 
	ptmax=3000.0;
	centrality=4; // 4=central (0-5%), 3=(10-20%), 2=(30-40%), 1=(50-60%), 0=(70-80%)
}

void CAcceptance_STAR::CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca){
	double eta,pt,pmag,gammav,phi,*p=part->p,y=part->y;
	int pid=part->resinfo->code,starpid;
	if(abs(pid)==211) starpid=1;
	else if(abs(pid)==321) starpid=2;
	else if(pid==-2212) starpid=3;
	else if(pid==2212) starpid=4;
	else{
		printf("CAcceptance_STAR::CalcAcceptance, pid=%d isn't in STAR list\n",pid);
		exit(1);
	}
	
	accept=false;
	efficiency=0.0;
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	pmag=sqrt(pt*pt+p[3]*p[3]);
	eta=atanh(p[3]/pmag);
	if(eta>ETAMIN && eta<ETAMAX){
		if(pt>ptmin && pt<ptmax){
			accept=true;
			phi=atan2(part->p[2],part->p[1]);
			star_acc_eff(starpid,0.001*pt,eta,phi,centrality,accept,efficiency);
		}
	}
}


void CAcceptance_STAR::star_acc_eff(int pid,double pt,double eta,double phi,int cen,bool &accept,double &eff){
		// Routine to return acceptance and efficiency
		// Gary D. Westfall
		// October, 2011
		// Based on results published by STAR PRC 79 034909 2009
		// Input
		//  pid=1 pion (either sign)
		//  pid=2 kaon (with sign)
		//  pid=3 anti-proton
		//  pid=4 proton
		//  cen=4 0-5% (Central)
		//  cen=3 10-20%
		//  cen=2 30-40%
		//  cen=1 50-60%
		//  cen=0 70-80% (Peripheral)
	
		//  basic acceptance cuts are
		//   |eta| < 1.0
		//   full phi acceptance
		//   calculate the efficiency based on particle ID, pt, and centrality
	
		// This routine is set up for full STAR TPC+TOF
	
	int i;
	double f1,f2,f3,f4,f5;
	double P0,P1,P2,P3,P4,P5,P6;
		// Pion parameters
	double P0_pion[5]={0.840,0.818,0.809,0.781,0.759};
	double P1_pion[5]={0.129,0.111,0.109,0.097,0.070};
	double P2_pion[5]={4.661,3.631,3.224,2.310,1.373};
		// Kaon parameters
	double P0_kaon[5]={0.608,0.585,0.503,0.494,0.450};
	double P1_kaon[5]={0.238,0.234,0.231,0.229,0.229};
	double P2_kaon[5]={2.425,3.034,3.968,3.492,3.925};
	double P3_kaon[5]={0.099,0.085,0.152,0.139,0.149};
		// pbar parameters
	double P0_pbar[5]={0.317,0.233,0.227,0.245,0.246};
	double P1_pbar[5]={0.303,0.303,0.300,0.305,0.301};
	double P2_pbar[5]={19.473,13.480,11.183,15.567,13.054};
	double P3_pbar[5]={-0.004,0.028,0.041,0.026,0.039};
	double P4_pbar[5]={10.156,4.160,2.031,1.895,0.889};
	double P5_pbar[5]={0.006,0.104,0.153,0.107,0.160};
		// p parameters
	double P0_p[5]={0.189,0.201,0.193,0.187,0.221};
	double P1_p[5]={0.306,0.310,0.303,0.308,0.307};
	double P2_p[5]={10.643,16.825,9.565,12.826,14.488};
	double P3_p[5]={0.023,0.026,0.041,0.042,0.033};
	double P4_p[5]={15.212,9.331,3.585,3.380,1.944};
	double P5_p[5]={0.131,0.127,0.194,0.194,0.181};
	
	accept=false;
	eff=0.0;
	
		// Do some data quality checks
	if(fabs(eta) > 1.0){
		return;
	}
	if((cen < 0) || (cen > 4)){
		return;
	}
	if(pt < 0.15){
		return;
	}
	
		// The basic parameters are OK, check the details
	
		// Pions
	if(pid == 1){
		if(pt < 1.6){
			accept=true;
			P0=P0_pion[cen];
			P1=P1_pion[cen];
			P2=P2_pion[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2));
			if(pt > 0.60){
				eff=0.64*eff;
				if(fabs(eta)>0.8) eff=0.0;
			}
		}
			//printf("eff=%g, pt=%g, eta=%g, cen=%d\n",eff,1000.0*pt,eta,cen);
	}
	
		// Kaons
	if(pid == 2){
		if((pt > 0.20) && (pt < 1.6)){
			accept=true;
			P0=P0_kaon[cen];
			P1=P1_kaon[cen];
			P2=P2_kaon[cen];
			P3=P3_kaon[cen];
			f1=P1/pt;
			eff=P0*exp(-pow(f1,P2))+P3*pt;
			if(pt > 0.60){
				eff=0.64*eff;
				if(fabs(eta)>0.8) eff=0.0;
			}
		}
	}
	
		// pbar
	if(pid == 3){
		if((pt > 0.4) && (pt < 3.0)){ // Note that STAR uses a lower cut of 0.4 GeV to suppress background protons
			accept=true;
			P0=P0_pbar[cen];
			P1=P1_pbar[cen];
			P2=P2_pbar[cen];
			P3=P3_pbar[cen];
			P4=P4_pbar[cen];
			P5=P5_pbar[cen];
			f1=pow(P1/pt,P2);
			f2=P0*exp(-f1);
			f3=P3*pt;
			f4=pow(P4/pt,P5);
			f5=exp(f4);
			eff=(f2+f3)*f5;
			if(pt > 1.0){
				eff=0.64*eff;
				if(fabs(eta)<0.8) eff=0.0;
			}
		}
	}
	
		// p
	if(pid == 4){
		if((pt > 0.4) && (pt < 3.0)){ // Note that STAR uses a lower cut of 0.4 GeV to suppress background protons
			accept=true;
			P0=P0_p[cen];
			P1=P1_p[cen];
			P2=P2_p[cen];
			P3=P3_p[cen];
			P4=P4_p[cen];
			P5=P5_p[cen];
			f1=pow(P1/pt,P2);
			f2=P0*exp(-f1);
			f3=P3*pt;
			f4=pow(P4/pt,P5);
			f5=exp(f4);
			eff=(f2+f3)*f5;
			if(pt > 1.0){
				eff=0.64*eff;
				if(fabs(eta)<0.8) eff=0.0;
			}
		}
	}

}