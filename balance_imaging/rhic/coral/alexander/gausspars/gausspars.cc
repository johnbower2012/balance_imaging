#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(){
	double x0=0.0,delx,dx,x,dydx,ytarget,y,ystep,*xarray;
	int nsample=30,iy;
	xarray=new double[nsample];
	ystep=1.0/double(nsample);
	int ntrial;
	x=0.0;
	for(iy=0;iy<nsample;iy++){
		ytarget=(0.5+iy)*ystep;
		ntrial=0;
		y=erf(x);
		do{
			ntrial+=1;
			dydx=exp(-x*x);
			delx=(y-ytarget)/dydx;
			if(fabs(delx)>0.1) delx=delx*0.1/fabs(delx);
			x-=delx;	
			y=erf(x);
		}while(fabs(y-ytarget)>1.0E-11 && ntrial<15);
		if(fabs(y-ytarget)>1.0E-11) printf("FAILURE, ntrial=%d, ytarget=%g, y=%g\n",ntrial,ytarget,y);
		printf("y=%12.10f, x=%15.10f\n",y,x);
		xarray[iy]=x;
	}
	
	double q;
	for(q=0.25;q<5;q+=0.5){
		y=0.0;
		for(iy=0;iy<nsample;iy++) y+=cos(q*xarray[iy]);
		y=y/double(nsample);
		printf("q=%g, y=%g =? %g\n",q,y,exp(-0.25*q*q));
	}

  return 0;
}


