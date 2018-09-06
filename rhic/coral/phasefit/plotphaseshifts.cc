#include "coral.h"
#include "xgraph.h"
#include "phaseshiftdata/read_phase.cc"
using namespace std;

int main(){
  CWaveFunction *wf;
  CReadPhaseShift *phasereader;
  CXGraph xgraph(400,400,50,50);
  xgraph.setaxes(0.0,0.0,400.0,180.0);
  xgraph.drawaxes();
  double q,delta;
  int iq,nqmax,ichannel,nchannels;
  char parsfilename[100];
  sprintf(parsfilename,"parameters/wfparameters.dat\0");
  wf=new CWaveFunction_kpluspiminus_sqwell(parsfilename);
  const double PI=4.0*atan(1.0);

  nqmax=wf->GetNQMAX();
  nchannels=wf->GetNCHANNELS();
  for(ichannel=0;ichannel<nchannels;ichannel++){
    if(ichannel==0) xgraph.setcolor("red");
    if(ichannel==1) xgraph.setcolor("green");
    if(ichannel==2) xgraph.setcolor("blue");
    if(ichannel==3) xgraph.setcolor("cyan");
    if(ichannel==4) xgraph.setcolor("yellow");
    for(iq=0;iq<nqmax;iq++){
      q=wf->GetQ(iq);
      delta=(180.0/PI)*wf->GetDELTA(ichannel,iq);
      if(delta<0.0) delta+=180.0;
      if(delta>180.0) delta-=180.0;
      xgraph.drawdiamond(q,delta,0.01);
    }
  }

  double *qarray,*darray;
  qarray=new double[nqmax];
  darray=new double[nqmax];
  xgraph.setcolor("black");

	/*
  // FOR PPi
  for(ichannel=0;ichannel<4;ichannel++){
    if(ichannel==0) phasereader=new CReadPhaseShift("PPi","s31");
    if(ichannel==1) phasereader=new CReadPhaseShift("PPi","p31");
    if(ichannel==2) phasereader=new CReadPhaseShift("PPi","p33");
    if(ichannel==3) phasereader=new CReadPhaseShift("PPi","s11");
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=wf->GetQ(iq);
      darray[iq]=(180.0/PI)*phasereader->ReadPhaseShift(qarray[iq]);
      if(darray[iq]<0.0) darray[iq]+=180.0;
      if(darray[iq]>180.0) darray[iq]-=180.0;
    }
    xgraph.plotline(qarray,darray,nqmax);
    delete (phasereader);
  }
	*/
  
	/*
  // FOR PK
  for(ichannel=0;ichannel<3;ichannel++){
    if(ichannel==0) phasereader=new CReadPhaseShift("PK","s11");
    if(ichannel==1) phasereader=new CReadPhaseShift("PK","p11");
    if(ichannel==2) phasereader=new CReadPhaseShift("PK","p13");
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=wf->GetQ(iq);
      darray[iq]=(180.0/PI)*phasereader->ReadPhaseShift(qarray[iq]);
      if(darray[iq]<0.0) darray[iq]+=180.0;
      if(darray[iq]>180.0) darray[iq]-=180.0;
    }
    xgraph.plotline(qarray,darray,nqmax);
    delete (phasereader);
  }
	*/
  
	/*
  //For pipi
  double ddeltadq;
  int ell,I;
  for(ichannel=0;ichannel<5;ichannel++){
    if(ichannel==0){I=0; ell=0;}
    if(ichannel==1){I=0; ell=2;}
    if(ichannel==2){I=1; ell=1;}
    if(ichannel==3){I=2; ell=0;}
    if(ichannel==4){I=2; ell=2;}
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=wf->GetQ(iq);
      WaveFunctionRoutines::getphaseshift_pipi(I,ell,qarray[iq],&delta,&ddeltadq);
      darray[iq]=(180.0/PI)*delta;
      if(darray[iq]<0.0) darray[iq]+=180.0;
      if(darray[iq]>180.0) darray[iq]-=180.0;
    }
    xgraph.plotline(qarray,darray,nqmax);
  }
	*/
	
  //For Kpi
  double ddeltadq;
  int ell,twoI;
  for(ichannel=0;ichannel<3;ichannel++){
    if(ichannel==0){twoI=1; ell=0;}
    if(ichannel==1){twoI=3; ell=0;}
    if(ichannel==2){twoI=1; ell=1;}
    for(iq=0;iq<nqmax;iq++){
      qarray[iq]=wf->GetQ(iq);
			WaveFunctionRoutines::getphaseshift_kpi(twoI,ell,qarray[iq],&delta,&ddeltadq);
      darray[iq]=(180.0/PI)*delta;
      if(darray[iq]<0.0) darray[iq]+=180.0;
      if(darray[iq]>180.0) darray[iq]-=180.0;
    }
    xgraph.plotline(qarray,darray,nqmax);
  }
	
  Misc::Pause();
  xgraph.closedisplay();
}

