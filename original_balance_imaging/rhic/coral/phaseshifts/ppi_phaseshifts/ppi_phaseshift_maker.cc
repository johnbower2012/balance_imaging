#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
using namespace std;
#include "utilities.h"
const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

int main(){
  int iq,iread,iqsmooth,ichannel;
  const int NREAD=61,nqmax=780;
  double qsmooth,q,tlab,elab,plab,roots;
  double a,reff,tandel;
  double deltaread[NREAD],qread[NREAD];
  
  int ell[3]={0,1,1};
  double eta[nqmax+1];
  double delta[3][nqmax+1]={0.0},ddeltadq[3][nqmax+1]={0.0};
  double delta_1s0,delta_3s1,delta_3p0,delta_1p1,delta_3p1,delta_3p2;
  double eps1,delta_1d2,delta_3d1,delta_3d2,delta_3d3;
  double delta_presmooth,ddeltadq_presmooth;
  double mu=938.28/2.0;
  char dummy[200];
  FILE *fptr;
  double e1,e2,e,MPION=139.58,MPROTON=938.28;
  double delq=0.5;

  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    e2=sqrt(q*q+MPION*MPION);
    e1=sqrt(q*q+MPROTON*MPROTON);
    //eta[iq]=e1*e2/((e1+e2)*137.036*q); // old way is wrong
    e=e1+e2;
    eta[iq]=(pow(e,4)-pow(e1*e1-e2*e2,2))/(4.0*e*e*e*137.036*q);
  }
  
  iqsmooth=2;
  qsmooth=25.0;
  fptr=fopen("s31.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPION*MPION+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MPION));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth) a=tan(deltaread[iread])/qread[iread];
  }
  printf("qmax=%g\n",qread[60]);
  fclose(fptr);
  
  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[0][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[0][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[0][iq]=atan(fabs(a*q))*a/fabs(a);
      ddeltadq[0][iq]=a*pow(cos(delta[0][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  qsmooth=75.0;
  fptr=fopen("p31.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPION*MPION+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MPION));
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth){
      a=tan(deltaread[iread])/pow(qread[iread],3);
    }
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[1][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[1][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[1][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[1][iq]=3*a*q*q*pow(cos(delta[1][iq]),2);
    }
  }

  // ________________________________________________________

  iqsmooth=5;
  qsmooth=75.0;
  fptr=fopen("p33.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<=60;iread++){
    fscanf(fptr,"%lf %lf",&plab,&deltaread[iread]);
    fgets(dummy,200,fptr);
    if(iread>0){
      elab=MPROTON+sqrt(MPION*MPION+plab*plab);
      roots=sqrt(elab*elab-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MPROTON,MPION));
      printf("%g =? %g\n",elab,sqrt(qread[iread]*qread[iread]+MPION*MPION)
	     +sqrt(qread[iread]*qread[iread]+MPROTON*MPROTON));
      printf("plab=%g, q=%g, delta=%g\n",plab,qread[iread],deltaread[iread]);
    }
    else qread[iread]=0.0;
    deltaread[iread]=deltaread[iread]*PI/180.0;
    if(iread==iqsmooth) a=tan(deltaread[iread])/pow(qread[iread],3);
  }
  fclose(fptr);

  iread=0;
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    if(q>qsmooth){
      if(iread>0) iread=iread-1;
      do{
	iread+=1;
      }while(q>qread[iread] && iread<=60);
      delta[2][iq]=((q-qread[iread-1])*deltaread[iread]
		    +((qread[iread]-q)*deltaread[iread-1]))
	/(qread[iread]-qread[iread-1]);
      ddeltadq[2][iq]=(deltaread[iread]-deltaread[iread-1])
	/(qread[iread]-qread[iread-1]);
    }
    else{
      delta[2][iq]=atan(fabs(a*q*q*q))*a/fabs(a);
      ddeltadq[2][iq]=3*a*q*q*pow(cos(delta[2][iq]),2);
    }
  }

  // _______________________________________________________
  
  for(iq=1;iq<nqmax;iq++){
    q=delq*iq;
    for(ichannel=0;ichannel<3;ichannel++)
      CoulWave::phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
					  delta[ichannel][iq],
					  ddeltadq[ichannel][iq]);
  }    

  fptr=fopen("ppi_phaseshiftdat.cc","w");
  fprintf(fptr,"   double delqdata=%g;\n",delq);
  fprintf(fptr,"   int nqdata=%d;\n",nqmax);
  fprintf(fptr,"   double data_delta[%d][%d]={",3,nqmax+1);
  for(ichannel=0;ichannel<3;ichannel++){
    fprintf(fptr,"{0");
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",delta[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=2) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"   double data_ddeltadq[%d][%d]={",3,nqmax+1);
  for(ichannel=0;ichannel<3;ichannel++){
    fprintf(fptr,"{0");
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",ddeltadq[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=2) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");

  fclose(fptr);
}
