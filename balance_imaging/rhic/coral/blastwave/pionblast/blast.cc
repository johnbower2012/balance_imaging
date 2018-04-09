#include "blast.h"
using namespace std;

CBlast::CBlast(double tauset,double RHBTset,double etaGset){
	tau=tauset; RHBT=RHBTset; etaG=etaGset;
	randy=new CRandom(-time(NULL));
	printf("tau=%g, RHBT=%g, etaG=%g\n",tau,RHBT,etaG);
}

void CBlast::GetUR(double *u,double *r){
	double eta,uperp2,rtilde,rho,phis,phib,yt,uperp,gammamax,u0,g,gv;
	gammamax=cosh(ytmax);
	do{
		do{
			r[1]=(1.0-2.0*randy->ran());
			r[2]=(1.0-2.0*randy->ran());
			rtilde=sqrt(r[1]*r[1]+r[2]*r[2]);
		}while(rtilde>1.0);
		r[1]*=Rx;
		r[2]*=Ry;
		phis=atan2(r[2],r[1]);
		phib=atan2(r[2]/(Ry*Ry),r[1]/(Rx*Rx));
		yt=(p0+p2*cos(2.0*phib)+p4*cos(4.0*phib))*rtilde;
		uperp=sinh(yt);
		u[1]=uperp*cos(phib);
		u[2]=uperp*sin(phib);	
		u[0]=sqrt(1.0+uperp*uperp);
		u[3]=0.0;
		r[3]=0.0;
		r[0]=tau;
	}while(u[0]/gammamax<randy->ran());

	//printf("check, Rx=%g, Ry=%g\n",Rx,Ry);
	//printf("r=(%g,%g,%g,%g) u=(%g,%g,%g,%g), yt=%g, phib=%g, phis=%g\n",r[0],r[1],r[2],r[3],u[0],u[1],u[2],u[3],yt,phib*180.0/PI,phis*180.0/PI);
}

void CBlast::GetP(double m,double *u,double *p){
	double p0[4],weight,yt,udotp0,wmax=1.0+tanh(ytmax);
	do{
		randy->generate_boltzmann(m,T,p0);
		udotp0=p0[0]*u[0]+p0[1]*u[1]+p0[2]*u[2]+p0[3]*u[3];
		Misc::Boost(u,p0,p);
		weight=udotp0/(u[0]*p0[0]);
		//printf("weight=%g\n",weight);
		if(weight>wmax){
			printf("weight>1.0!!!, =%g\n",weight);
			exit(1);
		}
	} while(randy->ran()>weight/wmax);
}

// -----------------------------------------------------------------------------------------------------//

bool acceptance(double *p){
	double y=atanh(p[3]/p[0]);
	double pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	double ymax=1.0,ptmin=150.0,ptmax=2000.0;
	if(pt>ptmin && pt<ptmax && fabs(y)<ymax) return true;
	else return false;
}

double GetPhiEff2_bosons(double *pa,double *ra,double *pb,double *rb){
		double qdotr=(pa[0]-pb[0])*(ra[0]-rb[0])
			-(pa[1]-pb[1])*(ra[1]-rb[1])-(pa[2]-pb[2])*(ra[2]-rb[2])-(pa[3]-pb[3])*(ra[3]-rb[3]);
		return 1.0+cos(qdotr/HBARC);
}

