#include "blast.h"
using namespace std;

void CBlast::SetPars(int pcentSet){
	pcent=pcentSet;
	if(pcent==0) Set0to5();
	else if(pcent==5) Set5to10();
	else if(pcent==10) Set10to20();
	else if(pcent==20) Set20to30();
	else if(pcent==30) Set30to40();
	else if(pcent==40) Set40to50();
	else if(pcent==50) Set50to60();
	else if(pcent==60) Set60to70();
	else if(pcent==70) Set70to80();
	else{
		printf("percent centrality =%d, not on recognized list\n",pcent);
		exit(1);
	}
	p4=0.0;
	//eps*=1.25;
	double scale=pow(Npart/394.0,1.0/3.0);
	//printf("SCALE=%g\n",scale);
	Mult=1000*pow(scale,3);
	RHBT*=scale;
	tau*=scale;
	Ry=sqrt((1.0+eps)*RHBT*RHBT);
	Rx=sqrt((1.0-eps)*RHBT*RHBT);
	ytmax=p0+fabs(p2)+fabs(p4);
	printf("IN CBLAST::SetPARS, RHBT=%g, tau=%g, Rx=%g, Ry=%g\n",RHBT,tau,Rx,Ry);
}

//Pion v_2{2}
void CBlast::Set0to5(){
	p2=0.031;
//p4=0.0022;
	p4=0;
	eps=0.0325;
	p0=1.03;
	T=96.5;
	Npart=347.3;
	Nch=406.1;
}

//Pion v_2
void CBlast::Set5to10(){
	p2=0.0406;
	p4=0.00407;
	eps=0.0482;
	p0=0.999;
	T=99.6;
	Npart=293.3;
	Nch=343;
}

//Pion v_2
void CBlast::Set10to20(){
	p2=0.0532;
	p4=0.00514;
	eps=0.078;
	p0=0.982;
	T=99.5;
	Npart=229.0;
	Nch=267.5;
}

//Pion v_2{2}
void CBlast::Set20to30(){
	p2=0.0691;
//p4=0.00514;
	p4=0;
	eps=0.127;
	p0=0.938;
	T=101.0;
	Npart=162.0;
	Nch=188.0;
}

//Pion v_2
void CBlast::Set30to40(){
	p2=0.0699;
	p4=0.00844;
	eps=0.151;
	p0=0.894;
	T=104.0;
	Npart=112.0;
	Nch=128.0;
}


//Pion v_2{2}
void CBlast::Set40to50(){
	p2=0.0882;
//p4=0.00514;
	p4=0;
	eps=0.181;
	p0=0.839;
	T=107.0;
	Npart=74.2;
	Nch=83.5;
}

//Pion v_2
void CBlast::Set50to60(){
	p2=0.0623;
	p4=0.0113;
	eps=0.221;
	p0=0.788;
	T=110.0;
	Npart=45.8;
	Nch=51.5;
}

//Pion v_2{2}
void CBlast::Set60to70(){
	p2=0.08;
//p4=0.0224;
	p4=0;
	eps=0.284;
	p0=0.688;
	T=120.0;
	Npart=25.9;
	Nch=29;
}

//Pion v_2
void CBlast::Set70to80(){
	p2=0.0477;
	p4=0.0396;
	eps=0.378;
	p0=0.608;
	T=125.0;
	Npart=13.0;
	Nch=14.5;
}

