#include "coral.h"

using namespace std;

int main(){
	double Rx=3.0,Ry=4.0,Rz=5.0;
	int ievent,nevents=1,ipart,nparts=50000,id=211;
	double r[4],p[4],m=139.57,T=200.0;
	CRandom randy(-1234);
	FILE *fptr=fopen("/Users/scottepratt/data/sdata/oscar_gauss345.tmp","w");
	fprintf(fptr,"OSC1997A\nfinal_id_p_x\nGaussian\n");
	ievent=1; nevents=1;
	fprintf(fptr,"%d %d 0.0 0.0\n",ievent,nparts);
	for(ipart=0;ipart<nparts;ipart++){
		r[0]=1.0;
		r[1]=Rx*randy.gauss();
		r[2]=Ry*randy.gauss();
		r[3]=Rz*randy.gauss();
		randy.generate_boltzmann(m,T,p);
		p[3]=0.0;
		p[1]=sqrt(p[1]*p[1]+p[2]*p[2]);
		p[2]=0.0;
		p[0]=sqrt(m*m+p[1]*p[1]);
		fprintf(fptr,"%d %d %g %g %g %g %g %g %g %g %g\n",ipart,id,0.001*p[1],0.001*p[2],0.001*p[3],0.001*p[0],0.001*m,r[1],r[2],r[3],r[0]);
	}
	fclose(fptr);
  return 0;
}


