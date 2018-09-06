#include <cstdlib>
#include <cmath>
#include <cstdio>
#include "coral.h"
using namespace std;

class CBlast{
public:
	double T,tau,RHBT,etaG,Rx,Ry,Mult,Npart,Nch;
	int pcent;
	void GetUR(double *u,double *r);
	void GetP(double m,double *u,double *p);
	CRandom *randy;
	CBlast(double tauset,double RHBTset,double etaGset);
	double eps,p0,p2,p4,ytmax;
	void Set0to5();
	void Set5to10();
	void Set10to20();
	void Set20to30();
	void Set30to40();
	void Set40to50();
	void Set50to60();
	void Set60to70();
	void Set70to80();
	void SetPars(int pcent);
};

bool acceptance(double pt,double y);

double GetM(CRandom *randy);

void decompose(double *pa,double *pb,double *ra,double *rb,double &qlong,double &qout,double &qside,double &qinv,double &rinv,double &kt);

void CalcDelr(double *ra,double *rb,double *pa,double *pb);

double GetPhiEff2(double *pa,double ma,double *ra,double *pb,double mb,double *rb,int sign,CBlast *blast);

double GetPhiEff2_bosons(double *pa,double ma,double *ra,double *pb,double mb,double *rb);
