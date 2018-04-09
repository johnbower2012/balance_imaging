#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "gslrandom.h"

const double PI=3.14159265358979323844;
const double HBARC=197.3269602;

using namespace std;

int main(){
	FILE *oscarfile;
	double r[4],p[4];
	double MA=139.57;
	double MB=139.57;
	double mass,rdummy1,rdummy2;
	int ident,idummy,i,alpha;
	int ndummy_header=3,ndummy_betweenevents=0;
	int npart,npartmax,ipart,ievent;
	int NEVENTSMAX=5;
	string OSCARfilename="/Users/scottepratt/data/bb/default/oscar.dat";
	char dummy[160];
	int na=0;
	double phi,eta,r0,p0,rx,px,cosheta,sinheta,cphi,sphi;
	CRandom *randy=new CRandom(-1234);
	FILE *joernfile=fopen("joern.dat","w");
	fprintf(joernfile,"# t  x  y  z  px  py\n");
	
	oscarfile=fopen(OSCARfilename.c_str(),"r");
	// Read Header Info
	for(i=0;i<ndummy_header;i++){
		fgets(dummy,120,oscarfile);
		//printf("dummy line=%s\n",dummy);
	}

	//for(ievent=1;ievent<=nevents;ievent++){
	ievent=0;
	do{
		ievent+=1;
		fscanf(oscarfile,"%d %d %lf %lf",&idummy,&npartmax,&rdummy1,&rdummy2);
		fgets(dummy,120,oscarfile);
		//printf("ievent=%d, idummy=%d, npartmax=%d\n",ievent,idummy,npartmax);
		for(npart=0;npart<npartmax;npart++){
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				&ipart,&ident,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
			if(ident==211 || ident==-211 || ident==111){
				for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0;
				r0=r[0];
				p0=p[0];
				eta=atanh(p[3]/p[0]);
				cosheta=cosh(eta);
				sinheta=sinh(eta);
				p[0]=sqrt(MA*MA+p[1]*p[1]+p[2]*p[2]);
				//p[3]=p[3]*cosheta-p0*sinheta;
				p[3]=0.0;
				r[0]=r0*cosheta-r[3]*sinheta;
				r[3]=r[3]*cosheta-r0*sinheta;
				phi=2.0*PI*randy->ran();
				rx=r[1];
				px=p[1];
				cphi=cos(phi);
				sphi=sin(phi);
				r[1]=rx*cphi+r[2]*sphi;
				r[2]=r[2]*cphi-rx*sphi;
				p[1]=px*cphi+p[2]*sphi;
				p[2]=p[2]*cphi-px*sphi;
				fprintf(joernfile,"%10.3e %10.3e %10.3e %10.3e  %10.3e %10.3e\n",r[0],r[1],r[2],r[3],p[1],p[2]);				
			}
		}
		for(i=0;i<ndummy_betweenevents;i++){
			if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
			//printf("dummy line=%s\n",dummy);
		}
	} while(ievent<NEVENTSMAX && !feof(oscarfile));
	fclose(oscarfile);
	fclose(joernfile);
	return 0;
}

#include "gslrandom.cc"