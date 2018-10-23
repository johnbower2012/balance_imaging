#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <armadillo>
#include <random>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;

double dumbfunc(double x){
	double x0=0.0,sigma0=0.125;
	//return (0.5/cosh(x/x0))*(1.0-cos(x/sigma0));
	//return 0.5*exp(-0.5*(x-x0)*(x-x0)/(sigma0*sigma0))
		//+0.5*exp(-0.5*(x+x0)*(x+x0)/(sigma0*sigma0));
	//return (x/x0)*(x/x0)*exp(-0.5*x*x/(x0*x0));
	//return 0.75*(x/x0)*(x/x0)/cosh(x/x0);
	//return 0.5*exp(-0.5*x*x/(sigma0*sigma0))+0.5*exp(-0.5*(x-x0)*(x-x0)/(sigma0*sigma0))
	//+0.5*exp(-0.5*(x+x0)*(x+x0)/(sigma0*sigma0));
	//return sigma0*sigma0/((x/x0)*(x/x0)+sigma0*sigma0);
	//return 1.0/cosh(x/sigma0);
	return 1.0/pow(cosh(x),1.0/sigma0);
}

class CCoshExpand{
public:
	int nmax;
	CCoshExpand(int nmaxset);
	double GetG(int n,double x);
	arma::mat a;
	arma::vec I;
	void TestOverlap();
	void CalcOverlap(int n,int nprime);
	void CalcOverlapNumerical(int n,int nprime);
	void CalcOverlapVec(int n);
	void IntegrateCoshInv(double x);
	int seed;
	double GetXFromCoshInvDist();
	mt19937 rng;
	std::uniform_real_distribution<double> uniform_dist;
};

double CCoshExpand::GetG(int n,double x){
	double g,coshx;
	int m;
	coshx=cosh(x);
	g=0.0;
	for(m=0;m<=n;m++){
		g+=a(n,m)*pow(coshx,-m-1);
	}
	return g;
}

CCoshExpand::CCoshExpand(int nmaxset){
	nmax=nmaxset;
	arma::mat A;
	arma::vec x,y,ytest;
	int n,nprime,m,mprime;
	double norm,overlap;
	I.resize(2*nmax+2);
	I(0)=0.5*PI;
	I(1)=1.0;
	for(n=2;n<=2*nmax+1;n++){
		I(n)=I(n-2)*double(n-1)/double(n);
		//printf("I(%d)=%15.10e\n",n,I(n));
	}
	a.resize(nmax+1,nmax+1);
	a.zeros();
	a(0,0)=1.0;
	for(n=1;n<=nmax;n++){
		a(n,n)=1.0;
		A.resize(n,n);
		A.zeros();
		x.resize(n);
		y.resize(n);
		ytest.resize(n);
		for(nprime=0;nprime<n;nprime++){
			for(m=0;m<n;m++){
				A(nprime,m)=0.0;
				for(mprime=0;mprime<=nprime;mprime++){
					A(nprime,m)+=I(m+mprime+1)*a(nprime,mprime);
				}
			}
		}
		for(nprime=0;nprime<n;nprime++){
			y(nprime)=0.0;
			for(mprime=0;mprime<=nprime;mprime++)
				y(nprime)+=I(n+mprime+1)*a(nprime,mprime);
		}
		arma::solve(x,A,y,arma::solve_opts::no_approx);
		for(m=0;m<n;m++)
			a(n,m)=-x(m);
		norm=0.0;
		for(m=0;m<=n;m++){
			for(mprime=0;mprime<=n;mprime++){
				norm+=a(n,m)*a(n,mprime)*I(m+mprime+1);
			}
		}
		if(norm<0.0){
			printf("norm(%d)=%g\n",n,norm);
			printf("negative norm, try making nmax smaller\n");
			exit(1);
		}
		x=x/sqrt(norm);
		y=y/sqrt(norm);
		//ytest=A*x-y;
		//ytest.print("ytest:");
		for(m=0;m<=n;m++)
			a(n,m)=a(n,m)/sqrt(norm);
	}
	seed=123456;
	rng=mt19937(seed);
	uniform_dist=uniform_real_distribution<double>(0,1);
}

void CCoshExpand::TestOverlap(){
	int n,nprime;
	for(n=0;n<=nmax;n++){
		for(nprime=0;nprime<=nmax;nprime++){
			CalcOverlap(n,nprime);
			CalcOverlapNumerical(n,nprime);
			printf("---------------------------\n");
		}
	}
}

void CCoshExpand::CalcOverlapNumerical(int n,int nprime){
	double overlap=0.0,g,gprime;
	double x,dx=0.001;
	for(x=0.5*dx;x<35;x+=dx){
		g=GetG(n,x);
		gprime=GetG(nprime,x);
		overlap+=g*gprime*dx;
	}
	if(abs(overlap)<5.0E-9)
		overlap=0.0;
	printf("numerical: overlap(%d,%d)=%g\n",n,nprime,overlap);
}

void CCoshExpand::CalcOverlap(int n,int nprime){
	double overlap=0.0;
	int m,mprime;
	for(m=0;m<=n;m++){
		for(mprime=0;mprime<=nprime;mprime++){
			overlap+=a(n,m)*I(m+mprime+1)*a(nprime,mprime);
		}
	}
	printf("overlap(%d,%d)=%g\n",n,nprime,overlap);
}

void CCoshExpand::CalcOverlapVec(int n){
	arma::vec overlap,x,y;
	overlap.resize(n);
	overlap.zeros();
	x.resize(n);
	x.zeros();
	y.resize(n);
	y.zeros();
p	arma::mat A;
	A.resize(n,n);
	A.zeros();
	int m,mprime,nprime;
	for(m=0;m<n;m++)
		x(m)=a(n,m);
	for(nprime=0;nprime<n;nprime++){
		for(m=0;m<n;m++){
			A(nprime,m)=0.0;
			for(mprime=0;mprime<=nprime;mprime++){
				A(nprime,m)+=I(m+mprime+1)*a(nprime,mprime);
			}
		}
	}
	for(nprime=0;nprime<n;nprime++){
		for(mprime=0;mprime<=nprime;mprime++){
			y(nprime)+=I(n+mprime+1)*a(nprime,mprime)*a(n,n);
		}
	}
	overlap=A*x+y;
	overlap.print();
	CalcOverlap(n,n);
}

void CCoshExpand::IntegrateCoshInv(double x){
	printf("Testing getting random no. from dist (2/pi)/cosh(x)\n");
	double answer=0,dx=0.0001,xx;
	for(xx=0.5*dx;xx<x;xx+=dx){
		answer+=dx/cosh(xx);
	}
	printf("%g =? %g\n",answer,atan(sinh(x)));
	
	double y,y2bar,dy=0.01;
	int ibin,nbins=1000;
	double *prob=new double[nbins];
	int n,nmax=100000000;
	for(n=0;n<nmax;n++){
		y=GetXFromCoshInvDist();
		ibin=lrint(floor(y/dy));
		if(ibin<nbins)
			prob[ibin]+=1.0/(nmax*dy);
	}
	for(ibin=0;ibin<nbins;ibin++){
		y=(ibin+0.5)*dy;
		printf("%5d %6.3f %10.7f =? %10.7f\n",ibin,y,prob[ibin],(2.0/PI)/cosh(y));
	}
	
	exit(1);
}

double CCoshExpand::GetXFromCoshInvDist(){
	// returns random number from dist (2/pi)/cosh(x)
	return asinh(tan(0.5*PI*uniform_dist(rng)));
}

int main(int argc,char *argv[]){
	int nmax=11;
	double eta0,eta,deta=0.0001;
	int n;
	printf("Enter eta0: ");
	scanf("%lf",&eta0);
	CCoshExpand coshexpand(nmax);
	
	
	coshexpand.IntegrateCoshInv(1.3);
	
	//coshexpand.a.print();
	//coshexpand.TestOverlap();

	arma::vec c(nmax+1);
	for(n=0;n<=nmax;n++){
		c(n)=0.0;
		for(eta=0.5*deta;eta<25;eta+=deta){
			c(n)+=dumbfunc(eta)*coshexpand.GetG(n,eta/eta0)*deta/sqrt(eta0);
		}
		printf("c(%d)=%g\n",n,c(n));
	}
	double dumbo;
	FILE *fptr=fopen("figs/compare.dat","w");
	deta=0.05;
	for(eta=0.5*deta;eta<8;eta+=deta){
		dumbo=0.0;
		for(n=0;n<=nmax;n++){
			dumbo+=c(n)*coshexpand.GetG(n,eta/eta0)/sqrt(eta0);
		}
		printf("%4.2f %g =? %g\n",eta,dumbo,dumbfunc(eta));
		fprintf(fptr,"%4.2f %9.6f %9.6f\n",eta,dumbo,dumbfunc(eta));
	}
	fclose(fptr);
	
	return 0;
}

