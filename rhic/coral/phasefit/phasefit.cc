#include "coral.h"
using namespace std;

const int NWELLS=3;

//#define PIPI
#define KPI
//Use one or the other
#ifdef PIPI
int I=1;
#elif defined KPI
int twoI=3;
#else
#include "phaseshiftdata/read_phase.cc"
string SPECIES="PK";
string CHANNEL="p13";
CReadPhaseShift *phasereader;
#endif

//double m1=139.57,m2=139.57;
//double m1=938.28,m2=139.57;
//double m1=938.28,m2=493.677;
double m1=493.677,m2=139.57;

int ELL=1;

const int NQMAX=32;
const double DELQ=10;

double CalcError(double *x);
double CalcError_Print(double *x,double *delx);
double GetDelta(double q);

double SHRINKSCALE=100.0;
const int NDIM=2*NWELLS;
CGSLMatrix_Real *gslmatrix;
double DELTA[NQMAX];
double BESTX[NDIM];
double MU=m1*m2/(m1+m2);
double SMALLESTERROR=1.0E20;

void GetDelta(){
	const double PI=4.0*atan(1.0);
	double q,ddeltadq;
	int iq;
	for(iq=0;iq<NQMAX;iq++){
		q=(0.5+iq)*DELQ;
		#ifdef PIPI
		WaveFunctionRoutines::getphaseshift_pipi(I,ELL,q,&DELTA[iq],&ddeltadq);
		#elif defined KPI
		WaveFunctionRoutines::getphaseshift_kpi(twoI,ELL,q,&DELTA[iq],&ddeltadq);
		#else
    DELTA[iq]=phasereader->ReadPhaseShift(q);
		#endif
		if(DELTA[iq]<0.0) DELTA[iq]+=PI;
	}
}

void CopyX(double *x0,double *x){
  for(int i=0;i<NDIM;i++) x[i]=x0[i];
}

void ReadPotentialPars(char *filename,double *x,double *delx){
  FILE *fptr=fopen(filename,"r");
  for(int iwell=0;iwell<NWELLS;iwell++){
    fscanf(fptr,"%lf %lf",&x[2*iwell],&delx[2*iwell]);
    fscanf(fptr,"%lf %lf",&x[2*iwell+1],&delx[2*iwell+1]);
  }
  fclose(fptr);
}

bool Newton(int ntries,double *x0,double *delx);
bool SteepestGrad(int ntries, double *x0, double *delx);
bool metropolis(int ntries, double *x0, double *delx);

int main(){
	gslmatrix=new CGSLMatrix_Real(NDIM);
	#ifndef PIPI
	#ifndef KPI
	phasereader=new CReadPhaseShift(SPECIES.c_str(),CHANNEL.c_str());
	#endif
	#endif
  char filename[120];
  double x[NDIM],delx[NDIM];
  int ntries,i,iq;
  bool success=0;
  char mintype;	
	GetDelta();
		
  sprintf(filename,"potentialpars.dat\0");
  ReadPotentialPars(filename,x,delx);
  CalcError_Print(x,delx);
	
  do{
    printf("Enter 's' for Steepest Descent or 'n' for Newton's method, 'd' to change delx, 'x' to change shrinkscale, and '0' to quit: ");
    scanf("%s",&mintype);
    if(mintype=='s'){
      printf("How many tries? ");
      scanf("%d",&ntries);
      success=SteepestGrad(ntries,x,delx);
    }
    if(mintype=='n'){
      printf("How many tries? ");
      scanf("%d",&ntries);
      success=Newton(ntries,x,delx);
    }
    if(mintype=='d'){
      for(i=0;i<NDIM;i++){
				printf("Enter delx[%d] : ",i);
				scanf("%lf",&delx[i]);
      }
    }
    if(mintype=='x'){
      printf("Enter shrinkscale : ");
      scanf("%lf",&SHRINKSCALE);
    }
  }while(mintype=='n' || mintype=='s' || mintype=='d' || mintype=='x');
	
}


double CalcError(double *x){
  double booboo=0.0,q,diff,dd;
	double *V0,*r;
  const double PI=4.0*atan(1.0);
	int iq,i,iwell;
	V0=new double[NWELLS];
	r=new double[NWELLS];
	for(iwell=0;iwell<NWELLS;iwell++){
		r[iwell]=x[2*iwell+1];
		V0[iwell]=x[2*iwell];
	}
  for(iq=0;iq<NQMAX;iq++){
		q=(0.5+iq)*DELQ;
    dd=Misc::CalcDelta_FromSqWells(ELL,MU,NWELLS,q,V0,r);
		
    diff=(180.0/PI)*fabs(DELTA[iq]-dd);
    while(diff>90.0){
      diff=fabs(diff-180.0);
    }
    booboo+=diff*diff*pow(90.0/q,2);
  }
	//booboo=pow((x[0]+35000)/100.0,2)+pow((x[1]-0.3)/0.01,2)+pow((x[2]-600)/100.0,2)+pow((x[3]-1.0)/0.01,2) +pow((x[4]+20)/100.0,2)+pow((x[5]-1.5)/0.01,2);
  if(booboo<SMALLESTERROR){
    SMALLESTERROR=booboo;
    for(i=0;i<NDIM;i++) BESTX[i]=x[i];
  }
	delete[] V0,r;
  return booboo;
}

double CalcError_Print(double *x,double *delx){ 
  const double PI=4.0*atan(1.0);
  double booboo,q,dd;
  double *V0,*r;
	int iwell;
  booboo=CalcError(x);
  int iq,i;
	V0=new double[NWELLS];
	r=new double[NWELLS];
	for(iwell=0;iwell<NWELLS;iwell++){
		r[iwell]=x[2*iwell+1];
		V0[iwell]=x[2*iwell];
	}
  printf("_______________  x,delx = :  __ ____________\n");
  for(i=0;i<NDIM;i++) printf("%g  %g\n",x[i],delx[i]);
  printf("////  error=%g ////\n",booboo);
  printf("____________________________________________\n");
  for(iq=0;iq<NQMAX;iq++){
    q=(iq+0.5)*DELQ;
		dd=Misc::CalcDelta_FromSqWells(ELL,MU,NWELLS,q,V0,r);
    printf("q=%6.1f delta=%7.2f =? %7.2f\n",
			q,dd*180.0/PI,DELTA[iq]*180.0/PI);
  }
  printf("____________________________________________\n");
	delete[] V0,r;
  return booboo;
}

bool SteepestGrad(int ntries, double *x00, double *delx00){
  int success;

  int i,itry;
  double norm;
  double delx[NDIM],du[NDIM],dx[NDIM],dydu[NDIM];
	double y0,y1,y2,*x0,*x1,*x2,*xtemp,del,dx2,dxdotdelx,ycheck;
	x0=new double[NDIM];
	x1=new double[NDIM];
	x2=new double[NDIM];
	
  CopyX(delx00,delx);
	CopyX(x00,x1);
	y1=CalcError(x1);
  for(itry=0;itry<ntries;itry++){
		printf("^^^^^^^^^^^^^^^ itry=%d ^^^^^^^^^^^^^^^^^, y1=%g\n",itry,y1);
		for(i=0;i<NDIM;i++) printf("x1[%d]=%g ",i,x1[i]);
		printf("\n");
		for(i=0;i<NDIM;i++) printf("delx[%d]=%g ",i,delx[i]);
		printf("\n");
		
		// First find direction
    for(i=0;i<NDIM;i++){
			CopyX(x1,x2);
			CopyX(x1,x0);
      x2[i]=x1[i]+0.5*delx[i];
      dydu[i]=CalcError(x2);
      x0[i]=x1[i]-0.5*delx[i];
      dydu[i]-=CalcError(x1);
    }
    norm=0.0;
    for(i=0;i<NDIM;i++){
      du[i]=-dydu[i];
      norm+=du[i]*du[i];
    }
    norm=sqrt(norm);
    for(i=0;i<NDIM;i++){
      du[i]=du[i]/norm;
      dx[i]=delx[i]*du[i];
    }
		printf("Will try this direction\n");
		for(i=0;i<NDIM;i++) printf("dx[%d]=%g ",i,dx[i]);
		printf("\n");
    // Direction dx[i]  is found
		for(i=0;i<NDIM;i++){
			x2[i]=x1[i]+dx[i];
			y2=CalcError(x2);
			x0[i]=x1[i]-dx[i];
			y0=CalcError(x0);
		}		
		
		for(success=0;success<2;success++){
			ycheck=y1;
			printf("__________ success=%d _______________, starting with y=%g\n",success,ycheck);
			while(y2<y1 || y0<y1){
				for(i=0;i<NDIM;i++) dx[i]*=2;		
				if(y2<y0){
					xtemp=x1;
					x1=x2;
					x2=xtemp;
					for(i=0;i<NDIM;i++) x2[i]=x1[i]+dx[i];
					y1=y2;
					y2=CalcError(x2);
				}
				else{
					xtemp=x1;
					x1=x0;
					x0=xtemp;
					for(i=0;i<NDIM;i++) x0[i]=x1[i]-dx[i];	
					y1=y0;
					y0=CalcError(x0);
				}
			}
			if(ycheck<y1){
				printf("OOPS, y1=%g should have been smaller than %g\n",y1,ycheck);
				exit(1);
			}
			ycheck=y1;
			// Now do Newton's method
			del=0.5*(y0-y2)/(y2+y0-2.0*y1);
			for(i=0;i<NDIM;i++) x1[i]+=del*dx[i];
			y1=CalcError(x1);
			if(ycheck<y1){
				printf("Newton failed to help, y1=%g should have been smaller than %g\n",y1,ycheck);
				//exit(1);
				for(i=0;i<NDIM;i++) x1[i]-=del*dx[i];
				y1=ycheck;
			}
			else{
				printf("!!!!!!!!!! NEWTON HELPED !!!!!!!!!!!, y1=%g, ycheck=%g\n",y1,ycheck);
				for(i=0;i<NDIM;i++){
					dx[i]*=0.25;
					x2[i]=x1[i]+dx[i];
					x0[i]=x1[i]-dx[i];
				}
				y0=CalcError(x0);
				y2=CalcError(x2);
			}
		}
		
	}
	for(i=0;i<NDIM;i++){
		delx[i]=0.5*delx[i];
	}
  CopyX(delx,delx00);
  CalcError_Print(BESTX,delx00);
	CopyX(BESTX,x00);
	delete [] x0;
	delete [] x1;
	delete [] x2;
	return 1;
}

bool Newton(int ntries, double *x0, double *delx){
  bool success=1;
  int iq,Nq;
  double shrinktest,shrink[NDIM],mintest;
  char parsfilename[100];
  double x[NDIM],dx[NDIM],**d2ydx2,dydx[NDIM],y,y0;
  int itry,i,j,jgood,jbad;
  d2ydx2=new double *[NDIM];
  for(i=0;i<NDIM;i++) d2ydx2[i]=new double[NDIM];
	
	y0=CalcError(x0);
  for(itry=0;itry<ntries;itry++){
    printf("____________________ itry=%d __________________________\n",itry);
    for(i=0;i<NDIM;i++){
      CopyX(x0,x);
      x[i]+=delx[i];
      dydx[i]=CalcError(x);
      CopyX(x0,x);
      x[i]-=delx[i];
      dydx[i]-=CalcError(x);
      dydx[i]=0.5*dydx[i]/delx[i];
    }
    for(i=1;i<NDIM;i++){
      for(j=0;j<i;j++){
				d2ydx2[i][j]=0.0;
				CopyX(x0,x);
				x[i]+=0.5*delx[i];
				x[j]+=0.5*delx[j];
				d2ydx2[i][j]=CalcError(x);
				CopyX(x0,x);
				x[i]+=0.5*delx[i];
				x[j]-=0.5*delx[j];
				d2ydx2[i][j]-=CalcError(x);
				CopyX(x0,x);
				x[i]-=0.5*delx[i];
				x[j]+=0.5*delx[j];
				d2ydx2[i][j]-=CalcError(x);
				CopyX(x0,x);
				x[i]-=0.5*delx[i];
				x[j]-=0.5*delx[j];
				d2ydx2[i][j]+=CalcError(x);
				d2ydx2[i][j]=d2ydx2[i][j]/(delx[i]*delx[j]);
				d2ydx2[j][i]=d2ydx2[i][j];
      }
    }
    for(i=0;i<NDIM;i++){
      d2ydx2[i][i]=-2.0*y0;
      CopyX(x0,x);
      x[i]+=delx[i];
      d2ydx2[i][i]+=CalcError(x);
      CopyX(x0,x);
      x[i]-=delx[i];
      d2ydx2[i][i]+=CalcError(x);
      d2ydx2[i][i]=d2ydx2[i][i]/(delx[i]*delx[i]);
    }
    gslmatrix->SolveLinearEqs(dydx,d2ydx2,dx);
		for(i=0;i<NDIM;i++){
			for(j=0;j<NDIM;j++){
				printf("M[%d][%d]=%g ",i,j,d2ydx2[i][j]);
			}
			printf("\n");
		}
		for(i=0;i<NDIM;i++) printf("dydx[%d]=%g ",i,dydx[i]);
		printf("\n");
		for(i=0;i<NDIM;i++) printf("dx[%d]=%g ",i,dx[i]);
		printf("\n");
   	for(i=0;i<NDIM;i++) x[i]=x0[i]-dx[i];
		
		y=CalcError(x);
		printf("New y=%g\n",y);
		for(i=0;i<NDIM;i++) printf("x[%d]=%g\n",i,x[i]);
		if(y<y0){
			CopyX(x,x0);
			y0=y;
		}
    
		
		
    printf("----- smallest error=%g ------ \n",SMALLESTERROR);
    //CalcError_Print(BESTX,delx);
    //if(mintest<0.0) break;
  }
	
  for(i=0;i<NDIM;i++) delete [] d2ydx2[i];
  delete [] d2ydx2;
	
  printf("++++++++++ FINISHING NEWTON ++++++++++++\n");
  CalcError_Print(BESTX,delx);
	CopyX(BESTX,x0);
  return success;
}




