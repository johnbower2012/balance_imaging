#include <cstdlib>
#include <cmath>
#include <cstdio>
#include "coral.h"

using namespace std;

int main(){
	FILE *fptr=fopen("xyzt.dat","r");
	FILE *fptrout=fopen("results/xyzt_short.dat","w");
	double x[4],r,tmedian=19;
	double xbar[4]={0.0},sigma[4][4]={0.0};
	CRandom *randy=new CRandom(-1234);
	char dummy[100];
	int alpha,beta;
	long long int norm=0,nabove=0,nbelow=0;
	do{
		fscanf(fptr,"%lf %lf %lf %lf",&x[0],&x[1],&x[2],&x[3]);
		r=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
		if(r<60.0 && x[0]<60.0){
			if(randy->ran()<1.0) fprintf(fptrout,"%g %g %g %g\n",x[0],x[1],x[2],x[3]);
			if(x[0]>tmedian) nabove+=1;
			else nbelow+=1;
			norm+=1;
			for(alpha=0;alpha<4;alpha++){
				xbar[alpha]+=x[alpha];
				for(beta=0;beta<4;beta++){
					sigma[alpha][beta]+=x[alpha]*x[beta];
				}
			}
		}
		fgets(dummy,100,fptr);	
	}while(!feof(fptr));
	printf("tmedian =? %g, nbelow=%lld, nabove=%lld\n",tmedian,nbelow,nabove);
	printf("norm=%lld\n",norm);
	printf("xbar = ");
	for(alpha=0;alpha<4;alpha++){
		xbar[alpha]=xbar[alpha]/double(norm);
		printf("%6.2f ",xbar[alpha]);
	}
	printf("\n");
	printf("sigma =\n");
	for(alpha=0;alpha<4;alpha++){
		for(beta=0;beta<4;beta++){
			sigma[alpha][beta]=sigma[alpha][beta]/double(norm);
			sigma[alpha][beta]-=xbar[alpha]*xbar[beta];
			printf(" %8.2f ",sigma[alpha][beta]);
		}
		printf("\n");
	}
	fclose(fptr);
	fclose(fptrout);
	
}

