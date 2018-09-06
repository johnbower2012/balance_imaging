#include "blast.h"
using namespace std;

CBlast::CBlast(double tauset,double RHBTset,double etaGset){
	tau=tauset; RHBT=RHBTset; etaG=etaGset;
	randy=new CRandom(-time(NULL));
	printf("tau=%g, RHBT=%g, etaG=%g\n",tau,RHBT,etaG);
}

void CBlast::GetUR(double *u,double *r){
	double eta,uperp2,rtilde,rho,phis,phib,yt,uperp,gammamax;
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
		eta=0.0;
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

bool acceptance(double pt,double y){
	double ymax=1.0,ptmin=150.0,ptmax=2000.0;
	if(pt>ptmin && pt<ptmax && fabs(y)<ymax) return true;
	else return false;
}

double GetM(CRandom *randy){
	double m,xran;
	xran=randy->ran();
	/*
	if(xran<0.6) m=139.58;
	else{
	if(xran<0.8) m=493.677;
	else m=938.28;
	}
	*/
	m=139.58;
	return m;
}

void decompose(double *pa,double *pb,double *ra,double *rb,double &qlong,double &qout,double &qside,double &qinv,double &rinv,double &kt){
	double pt,s=0.0,Mlong,roots,r2=0.0,pdotr=0.0,pdotq=0.0;
	double ptot[4],q[4];
	const int g[4]={1,-1,-1,-1};
	int alpha;
	qinv=0.0;
	for(alpha=0;alpha<4;alpha++){
		ptot[alpha]=pa[alpha]+pb[alpha];
		s+=g[alpha]*ptot[alpha]*ptot[alpha];
		q[alpha]=pa[alpha]-pb[alpha];
		pdotq+=g[alpha]*ptot[alpha]*q[alpha];
		qinv-=g[alpha]*q[alpha]*q[alpha];
		r2+=g[alpha]*pow(ra[alpha]-rb[alpha],2);
		pdotr+=g[alpha]*ptot[alpha]*(ra[alpha]-rb[alpha]);
	}
	r2=r2-pdotr*pdotr/s;
	if(r2>0.0){
		printf("r2=%g, but should be <0",r2);
		exit(1);
	}
	rinv=sqrt(-r2);
	pt=sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
	Mlong=sqrt(s+pt*pt);
	roots=sqrt(s);

	qside=(ptot[1]*q[2]-ptot[2]*q[1])/pt;
	qlong=(ptot[0]*q[3]-ptot[3]*q[0])/Mlong;
	qout=(roots/Mlong)*(ptot[1]*q[1]+ptot[2]*q[2])/pt;
	qinv=sqrt(qinv+pdotq*pdotq/s);
	if(qinv<0.0){
		printf("qinv^2<0, =%g\n",qinv);
		exit(1);
	}
	qinv=0.5*qinv;
	qlong=0.5*qlong;
	qside=0.5*qside;
	qout=0.5*qout;
	kt=0.5*sqrt(ptot[1]*ptot[1]+ptot[2]*ptot[2]);
}

double GetPhiEff2(double *pa,double ma,double *ra,double *pb,double mb,double *rb,int sign,CBlast *blast){
	double qlong,qout,qside,qinv,rinv,mt,kt;
	decompose(pa,pb,ra,rb,qlong,qout,qside,qinv,rinv,kt);
	double lambdaHBT=0.5;
	double answer,V,E,alpha=197.323/137.036,arg,mu,Ea,Eb,eta,gamow;
	int iq;
	Ea=sqrt(qinv*qinv+ma*ma);
	Eb=sqrt(qinv*qinv+mb*mb);
	E=Ea+Eb;
	mu=Ea*Eb/E;
	eta=sign*137.036*mu/qinv;
	gamow=2.0*PI*eta/(exp(2.0*PI*eta)-1.0);
	if(eta>20) gamow=0.0;
	if(eta<-20) gamow=-2.0*PI*eta;
	if(fabs(eta<1.0E-6)) gamow=1.0;
	V=sign*alpha/rinv;
	//V=0.0;
	answer=1.0-V*( (mu/(qinv*qinv))+(0.5/mu)+1.5*(Ea-Eb)*(Ea-Eb)/(Ea*Eb*E) );
	if(sign==1){
		if(answer<gamow) answer=gamow;
	}
	else{
		if(answer>gamow) answer=gamow;
	}
	//answer=1.0;
	
	if(sign==1 && fabs(ma-mb)<1.0E-3 && fabs(ma)<500.0){
		double mt=sqrt(ma*ma+kt*kt);
		double rout,rside,rlong;
		rout=0.5*blast->Rx*sqrt(250.0/mt);
		rside=blast->Ry*rout/blast->Rx;
		rlong=0.6*blast->tau*sqrt(250.0/mt);
		rout=(mt/ma)*rout;
		arg=4.0*(rout*rout*qout*qout+rside*rside*qside*qside+rlong*rlong*qlong*qlong)/(HBARC*HBARC);
		if(arg<80) answer*=(1.0+exp(-arg));
		//printf("mt=%g, rout=%g, rside=%g, rlong=%g\n",mt,rout,rside,rlong);
	}
	answer=lambdaHBT*answer+(1.0-lambdaHBT);
	return answer;
}

double GetPhiEff2_bosons(double *pa,double ma,double *ra,double *pb,double mb,double *rb){
	if(fabs(ma-mb)<1.0 && ma<600){
		double qdotr=(pa[0]-pb[0])*(ra[0]-rb[0])
			-(pa[1]-pb[1])*(ra[1]-rb[1])-(pa[2]-pb[2])*(ra[2]-rb[2])-(pa[3]-pb[3])*(ra[3]-rb[3]);
		return cos(qdotr/HBARC);
	}
	else return 0;
}

