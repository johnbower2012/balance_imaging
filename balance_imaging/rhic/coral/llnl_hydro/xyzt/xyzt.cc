#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "coral.h"

using namespace std;

int main(){
	FILE *oscarfile;
	double r[4],p[4];
	double MA=139.57;
	double MB=139.57;
	double mass,rdummy1,rdummy2,phimax=30.0,cphimax=cos(phimax*PI/180.0),dummy1,dummy2;
	int ident,idummy,i,alpha;
	int ndummy_header=3,ndummy_betweenevents=0;
	int ipart,nparts,ipart_read,ievent,nwrite=0;
	int NEVENTSMAX=2000;
	string OSCARfilename="/Users/scottepratt/data/llnl_hydro/b2.21/bb_oscar.dat";
	char dummy[160];
	int na=0;
	double phi,eta,r0,p0,rx,px,pt,cosheta,sinheta,cphi,sphi;
	CRandom *randy=new CRandom(-1234);
	FILE *xyztfile=fopen("xyzt.dat","w");
	fprintf(xyztfile,"! t  x  y  z  px  py\n");

	oscarfile=fopen(OSCARfilename.c_str(),"r");
	// Read Header Info
	for(i=0;i<ndummy_header;i++){
		fgets(dummy,120,oscarfile);
		//printf("dummy line=%s\n",dummy);
	}

	//for(ievent=1;ievent<=nevents;ievent++){
	ievent=0;
	do{
		if(fscanf(oscarfile,"%d %d %g %g",&ievent,&nparts,&dummy1,&dummy2)){
			for(ipart=0;ipart<nparts;ipart++){
		  fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&ipart_read,&ident,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
			if(ident==211 || ident==-211 || ident==111){
				for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0;
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);
				if(pt>290.0&&pt<310.0){
					cphi=p[1]/pt;
					sphi=p[2]/pt;
					if(cphi>cphimax){
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
						r[1]=r[1]*cphi+r[2]*sphi;
						r[2]=r[2]*cphi-rx*sphi;
						p[1]=pt;
						p[2]=p[2]*cphi-px*sphi;
						fprintf(xyztfile,"%10.3e %10.3e %10.3e %10.3e  %10.3e %10.3e\n",r[0],r[1],r[2],r[3],p[1],p[2]);	
						nwrite+=1;
					}
				}			
			}
		}
		for(i=0;i<ndummy_betweenevents;i++){
			if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
			//printf("dummy line=%s\n",dummy);
		}
		ievent+=1;
	}
	} while(ievent<NEVENTSMAX && !feof(oscarfile));
	fclose(oscarfile);
	fclose(xyztfile);
	printf("nwrite=%d\n",nwrite);
	return 0;
}
