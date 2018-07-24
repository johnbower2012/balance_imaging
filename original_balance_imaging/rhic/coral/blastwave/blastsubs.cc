#include "blast.h"
using namespace std;

CBlast::CBlast(double Tset,double tauset,double Rxset,double Ryset,double etaGset,double uxmaxset,double uymaxset){
	T=Tset; tau=tauset; Rx=Rxset; Ry=Ryset; etaG=etaGset; uxmax=uxmaxset; uymax=uymaxset;
	randy=new CRandom(-time(NULL));
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
	eta=0.0;//etaG*randy->gauss();
	u[0]=sqrt(1.0+uperp2);//*cosh(eta);
	u[3]=0.0;//u[0]*tanh(eta);
	r[3]=0.0;//tau*sinh(eta);
	r[0]=tau;//*cosh(eta);
	//printf("r=(%g,%g,%g,%g), tau=%g=?%g\n",r[0],r[1],r[2],r[3],tau,sqrt(r[0]*r[0]-r[3]*r[3]));
}

void CBlast::GetP(double m,double *u,double *p){
	double p0[4];
	randy->generate_boltzmann(m,T,p0);
	Misc::Boost(u,p0,p);
}

// -----------------------------------------------------------------------------------------------------//

bool acceptance(double pt,double y){
	double ymax=1.0,ptmax=2000.0,ptmin=150.0;
	if(pt>ptmin && pt< ptmax && fabs(y)<ymax) return true;
	else return false;
}

double GetM(CRandom *randy){
	double m,xran;
	xran=randy->ran();
	if(xran<0.6) m=139.58;
	else{
		if(xran<0.8) m=493.677;
		else m=938.28;
	}
	//	m=139.58;
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

void CalcDelr(double *ra,double *rb,double *pa,double *pb){
	double ptot[4],u[4],rboost[4],r[4],s=0.0,pdotr=0.0,r2=0.0,g[4]={1.0,-1.0,-1.0,-1.0};
	int alpha;
	for(alpha=0;alpha<4;alpha++){
		ptot[alpha]=pa[alpha]+pb[alpha];
		r[alpha]=ra[alpha]-rb[alpha];
		pdotr+=ptot[alpha]*r[alpha]*g[alpha];
		s+=ptot[alpha]*ptot[alpha]*g[alpha];
		r2+=r[alpha]*r[alpha]*g[alpha];
	}
	for(alpha=0;alpha<4;alpha++){
		r[alpha]=r[alpha]-ptot[alpha]*pdotr/s;
		u[alpha]=ptot[alpha]/sqrt(s);
	}
	Misc::BoostToCM(u,r,rboost);
	printf("rboost=(%g,%g,%g,%g)\n",rboost[0],rboost[1],rboost[2],rboost[3]);
	r2=-r2+pdotr*pdotr/s;
	printf("r=%g\n",sqrt(r2));
}

double GetPhiEff2(double *pa,double ma,double *ra,double *pb,double mb,double *rb,int sign,parameterMap *parmap){
	double qlong,qout,qside,qinv,rinv,mt,kt;
	decompose(pa,pb,ra,rb,qlong,qout,qside,qinv,rinv,kt);
	double lambdaHBT=0.5;
	double answer,V,E,alpha=197.323/137.036,arg,mu,Ea,Eb,eta,gamow,scale=parameter::getD(*parmap,"HBTSCALE");
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

	if(sign==1 && fabs(ma-mb)<1.0E-3 && fabs(ma)<500.0){
		double mt=sqrt(ma*ma+kt*kt);
		double rout,rside,rlong;
		rout=scale*6.0*sqrt(250.0/mt);
		rside=scale*rout;
		rlong=scale*7.0*sqrt(250.0/mt);
		rout=(mt/ma)*rout;
		arg=4.0*(rout*rout*qout*qout+rside*rside*qside*qside+rlong*rlong*qlong*qlong)/(HBARC*HBARC);
		if(arg<80) answer*=(1.0+exp(-arg));
	}
	answer=lambdaHBT*answer+(1.0-lambdaHBT);
	return answer;
}




