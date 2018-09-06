#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "coral.h"
class CBlast{
public:
	double T,tau,Rx,Ry,etaG,uxmax,uymax;
	void GetUR(double *u,double *r);
	void GetP(double m,double *u,double *p);
	CRandom *randy;
	CBlast(double Tset,double tauset,double Rxset,double Ryset,double etaGset,double uxmaxset,double uymaxset){
		T=Tset; tau=tauset; Rx=Rxset; Ry=Ryset; etaG=etaGset; uxmax=uxmaxset; uymax=uymaxset;
		randy=new CRandom(-1234);
	}
};

int main(){
	long long int nsample=5000,ipair;
	int iphia,iphib,iphi;
	const long long int NPHI=36,NPAIRS=100000000;
	long double C[NPHI]={0.0},D[NPHI]={0.0};
	double phia,phib,phic,dphia,dphib;
	double psi2_ac_pp,psi2_bc_pp,psi2_ac_pm,psi2_bc_pm,weight,sum;
	double T=100.0,Rx=13.0,Ry=13.0,uxmax=1.0,uymax=1.0,R=13.0,etaG=1.7,Pt,tau=12.0,ymax=1.0;
	double Ma=139.58,Mb=139.58,Mc=139.58;
	double rab[4],rc[4],uab[4],uc[4];
	double pa[4],pb[4],pc[4],ya,yb,yc;
	FILE *fptr=fopen("figs/Bhbt.dat","w");
	CBlast *blast=new CBlast(T,tau,Rx,Ry,etaG,uxmax,uymax);
	
	CWaveFunction_pipluspiplus_sqwell *cw_pp=new CWaveFunction_pipluspiplus_sqwell(string("parameters/wfparameters.dat"));
	CWaveFunction_pipluspiminus_sqwell *cw_pm=new CWaveFunction_pipluspiminus_sqwell(string("parameters/wfparameters.dat"));

//________________________________________________________________________________//
	for(ipair=0;ipair<NPAIRS;ipair++){
		do{
			blast->GetUR(uab,rab);
			blast->GetP(Ma,uab,pa);
			ya=atanh(pa[3]/pa[0]);
			blast->GetP(Mb,uab,pb);
			yb=atanh(pb[3]/pb[0]);
			blast->GetUR(uc,rc);
			blast->GetP(Mc,uc,pc);
			yc=atanh(pb[3]/pb[0]);
		}while(fabs(yc)>ymax || (fabs(ya)>ymax && fabs(yb)>ymax));
		//printf("rab=(%g,%g,%g,%g) rc=(%g,%g,%g,%g)\n",rab[0],rab[1],rab[2],rab[3],rc[0],rc[1],rc[2],rc[3]);
		phia=atan2(pa[2],pa[1]);
		phib=atan2(pb[2],pb[1]);
		phic=atan2(pc[2],pc[1]);
		dphia=fabs(phia-phic);
		if(dphia>PI) dphia=2.0*PI-dphia;
		iphia=(dphia/PI)*double(NPHI);
		dphib=fabs(phib-phic);
		if(dphib>PI) dphib=2.0*PI-dphib;
		iphib=(dphib/PI)*double(NPHI);
		if(iphia>NPHI || iphib>NPHI){
			printf("iphi out of bounds, iphia=%d, iphib=%d\n",iphia,iphib);
			exit(1);
		}
		psi2_ac_pm=cw_pm->GetPsiSquared(pa,rab,pc,rc);
		psi2_ac_pp=cw_pp->GetPsiSquared(pa,rab,pc,rc);
		psi2_bc_pm=cw_pm->GetPsiSquared(pb,rab,pc,rc);
		psi2_bc_pp=cw_pp->GetPsiSquared(pb,rab,pc,rc);
		weight=psi2_ac_pm-psi2_bc_pm-psi2_ac_pp+psi2_bc_pp;
		if(fabs(ya)<ymax){
			C[iphia]+=weight;
			D[iphia]+=2.0;
		}
		if(fabs(yb)<ymax){
			C[iphib]-=weight;
			D[iphib]+=2.0;
		}
		if((1+ipair)%(NPAIRS/10)==0) printf("finished %g percent\n",100*(1+ipair)/double(NPAIRS));
	}
	for(iphi=0;iphi<NPHI;iphi++){
		printf("B(iphi=%d)/MULT=%g, D[iphi]=%g\n",iphi,double(C[iphi]/D[iphi]),double(D[iphi]));
		fprintf(fptr,"%6.2f %10.3e %10.3e\n",(0.5+iphi)*180.0/double(NPHI),double(C[iphi]/D[iphi]),double(D[iphi]));
	}
	fclose(fptr);	
}

void CBlast::GetUR(double *u,double *r){
	double eta,uperp2;
	do{
		u[1]=uxmax*(1.0-2.0*randy->ran());
		u[2]=uymax*(1.0-2.0*randy->ran());
	}while(pow(u[1]/uxmax,2)+pow(u[2]/uymax,2)<1.0);
	uperp2=u[1]*u[1]+u[2]*u[2];
	r[1]=Rx*u[1]/uxmax;
	r[2]=Ry*u[2]/uymax;
	eta=etaG*randy->gauss();
	u[0]=sqrt(1.0+uperp2)*cosh(eta);
	u[3]=u[0]*tanh(eta);
	r[3]=tau*u[3];
	r[0]=tau*u[1];
}

void CBlast::GetP(double m,double *u,double *p){
	double p0[4];
	randy->generate_boltzmann(m,T,p0);
	Misc::Boost(u,p0,p);
}
	
	

