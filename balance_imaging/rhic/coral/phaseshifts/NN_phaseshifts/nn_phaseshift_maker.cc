#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
using namespace std;
#include "utilities.h"
const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

double schroedinger(double q);
double vreid(double r);

int main(){
  int iq,iread,iqsmooth,ichannel;
  const int NREAD=35,nqmax=400;
  double qsmooth,q,tlab,elab,plab,roots;
  double a,reff,tandel,delq;
  double deltaread[4][NREAD],qread[NREAD];
  
  int ell[4]={0,1,1,1};
  double eta[nqmax+1];
  double delta[4][nqmax+1]={0.0},ddeltadq[4][nqmax+1]={0.0};
  double delta_1s0,delta_3s1,delta_3p0,delta_1p1,delta_3p1,delta_3p2;
  double eps1,delta_1d2,delta_3d1,delta_3d2,delta_3d3;
  double delta_presmooth,ddeltadq_presmooth;
  double mu;
  char dummy[200];
  FILE *fptr;
  double e1,e2,MNEUTRON=939.6;
  mu=MNEUTRON/2.0;
  delq=0.5;
	
  for(iq=0;iq<=nqmax;iq++){
    q=delq*double(iq);
    e1=e2=sqrt(q*q+MNEUTRON*MNEUTRON);
  }
  
  fptr=fopen("NN_phaseshifts.dat\0","r");
  fgets(dummy,200,fptr);
  for(iread=0;iread<NREAD;iread++){
    fscanf(fptr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&tlab,&delta_1s0,&delta_3s1,&delta_3d1,&eps1,&delta_3p0,
			&delta_1p1,&delta_3p1,&delta_3p2,&delta_1d2,&delta_3d2,&delta_3d3);
    deltaread[0][iread]=delta_1s0;
    deltaread[1][iread]=delta_3p0;
    deltaread[2][iread]=delta_3p1;
    deltaread[3][iread]=delta_3p2;
    if(iread>0){
      elab=MNEUTRON+tlab;
      plab=sqrt(elab*elab-MNEUTRON*MNEUTRON);
      roots=sqrt((elab+MNEUTRON)*(elab+MNEUTRON)-plab*plab);
      qread[iread]=sqrt(Misc::triangle(roots,MNEUTRON,MNEUTRON));
    }
    else qread[iread]=0.0;
    for(ichannel=0;ichannel<4;ichannel++){
      deltaread[ichannel][iread]=deltaread[ichannel][iread]*PI/180.0;
		}
  }
  fclose(fptr);
	
  iqsmooth=6;
  qsmooth=qread[iqsmooth];
  //printf("qsmooth=%g, qread_max=%g\n",qsmooth,qread[NREAD-1]);
	
  for(iq=1;iq<=nqmax;iq++){   
    q=delq*iq;
    if(q<=qread[NREAD-1]){
      iread=0;
      while(q>qread[iread]){
				iread+=1;
      }
      for(ichannel=1;ichannel<4;ichannel++){
				delta[ichannel][iq]=(qread[iread]-q)*deltaread[ichannel][iread-1]
				+(q-qread[iread-1])*deltaread[ichannel][iread];
				delta[ichannel][iq]=delta[ichannel][iq]/(qread[iread]-qread[iread-1]);
				ddeltadq[ichannel][iq]=(deltaread[ichannel][iread]
					-deltaread[ichannel][iread-1])
				/(qread[iread]-qread[iread-1]);
      }
      
    }
    else{
      for(ichannel=0;ichannel<4;ichannel++){
				delta[ichannel][iq]=deltaread[ichannel][NREAD-1];
				ddeltadq[ichannel][iq]=0.0;
      }
      //printf("Warning: qarray goes beyond max. for phaseshift data =%g\n",
      //     qread[NREAD-1]);
      //exit(1);
    }
    
    // Smooth out the p-waves at low q
    if(q<qsmooth){
      for(ichannel=1;ichannel<4;ichannel++){
				delta_presmooth=delta[ichannel][iq];
				ddeltadq_presmooth=ddeltadq[ichannel][iq];
				delta[ichannel][iq]=deltaread[ichannel][iqsmooth]
				*(q*q*q/(qsmooth*qsmooth*qsmooth));
				ddeltadq[ichannel][iq]=3.0*delta[ichannel][iq]/q;
				//delta[ichannel][iq]=ddeltadq[ichannel][iq]=0.0;
				//printf("q=%g, channel=%d : delta=%g=?%g, ddeltadq=%g=?%g\n",
				//     q,ichannel,delta[ichannel][iq],delta_presmooth,
				//     ddeltadq[ichannel][iq],ddeltadq_presmooth);
      }
    }
		
    // Re-do S-waves with good data
    //delta[0][iq]=sdelta[iq];
    //ddeltadq[0][iq]=sddeltadq[iq];
    double deltaq=0.04;
    if(iq>0){
      delta[0][iq]=schroedinger(q);
      ddeltadq[0][iq]=(schroedinger(q+0.5*deltaq)-schroedinger(q-0.5*deltaq))/deltaq;
      printf("iq=%d, q=%g, delta=%g, dd/dq=%g\n",iq,q,(180.0/PI)*delta[0][iq],HBARC*ddeltadq[0][iq]);
    }
  }
  ddeltadq[0][0]=2.0*ddeltadq[0][1]-ddeltadq[0][2];
  fptr=fopen("nn_phaseshiftdat.cc","w");
  fprintf(fptr,"   double delqdata=%g;\n",delq);
  fprintf(fptr,"   int nqdata=%d;\n",nqmax);
  fprintf(fptr,"   double data_delta[%d][%d]={",4,nqmax+1);
  for(ichannel=0;ichannel<4;ichannel++){
    fprintf(fptr,"{0");
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",delta[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=3) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");
  fprintf(fptr,"   double data_ddeltadq[%d][%d]={",4,nqmax+1);
  for(ichannel=0;ichannel<4;ichannel++){
    fprintf(fptr,"{%g",ddeltadq[ichannel][0]);
    for(iq=1;iq<=nqmax;iq++) fprintf(fptr,",%g",ddeltadq[ichannel][iq]);
    fprintf(fptr,"}");
    if(ichannel!=3) fprintf(fptr,",\n      ");
  }
  fprintf(fptr,"};\n");
	
  fclose(fptr);
}

double schroedinger(double q){
  const double PI=4.0*atan(1.0);
  double schr_rmax,schr_delr;
  int schr_nrmax;
  int ir;
  complex<double> psi0,psi1,psi2,psitest,psi00,psi11;
  complex<double> ci(0.0,1.0);
  complex<double> *psi,cg;
  double m1,m2;
  double r0,r1,r2,mu,e,v,phase00,phase0,phase1,phase2,x0,x1,sigma,delta;
  m1=m2=938.28;
  mu=m1*m2/(m1+m2);
  schr_rmax=6.0;
  schr_nrmax=6000;
  psi=new complex<double>[schr_nrmax+1];
  schr_delr=schr_rmax/double(schr_nrmax);
  e=q*q/(2.0*mu);
	
  // First solve without Reid potential
  r1=schr_rmax+2.0*schr_delr;
  r0=schr_rmax+schr_delr;
  x1=r1*q/HBARC;
  x0=r0*q/HBARC;
  psi11=x1*Bessel::h0(x1);
  psi1=psi11;
  psi00=x0*Bessel::h0(x0);
  psi0=psi00;
  for(ir=schr_nrmax;ir>=0;ir--){
    r2=r1;
    psi2=psi1;
    r1=r0;
    psi1=psi0;
    r0=ir*schr_delr;
    x0=r0*q/HBARC;
    v=0.0;
    psi0=2.0*psi1-psi2+(2.0*mu*v-q*q)*schr_delr*schr_delr*psi1/(HBARC*HBARC);
    psi[ir]=psi0;
  }
  phase2=atan2(imag(psi2),real(psi2));
	phase1=atan2(imag(psi1),real(psi1));
  phase0=2.0*phase1-phase2;
  phase00=phase0;
	
  // Now with Reid potential
  r1=schr_rmax+2.0*schr_delr;
  r0=schr_rmax+schr_delr;
  x1=r1*q/HBARC;
  x0=r0*q/HBARC;
  psi1=psi11;
  psi0=psi00;
  for(ir=schr_nrmax;ir>=0;ir--){
    r2=r1;
    psi2=psi1;
    r1=r0;
    psi1=psi0;
    r0=ir*schr_delr;
    v=vreid(r1);
    psi0=2.0*psi1-psi2+(2.0*mu*v-q*q)*schr_delr*schr_delr*psi1/(HBARC*HBARC);
    psi[ir]=psi0;
  }
  phase2=atan2(imag(psi2),real(psi2));
  phase1=atan2(imag(psi1),real(psi1));
  phase0=2.0*phase1-phase2;
  //printf("In Schr_init, q=%g, delta=%g\n",q,(phase00-phase0)*180.0/PI);
  delta=phase00-phase0;
  delete [] psi;
  return delta;
}
double vreid(double r){
  double vr,pmux,f1,f4,f7;
  pmux=r*0.7;
  f1=exp(-pmux);
  f4=(f1*f1*f1*f1);
  f7=f4*(f1*f1*f1);
  vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
  return vr;
}









