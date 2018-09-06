#include "coral.h"

int main (int argc, char *argv[]){
	FILE *oscarfile;
	string dataroot="/Users/pratt/data/llnl_hydro/";
	//string dataroot="/Users/pratt/data/hydro/";
	//string qualifier[12]={"default","eta0","eta4","INITA0","B0","B4","L0","L1600","c2_0","c2_0.2","collscaling","cgc"};
	string qualifier[12];
	qualifier[0]="1dhydro";
	int iq;
	FILE *fptr=fopen("results/ptbar.dat","w");
	fprintf(fptr,"                pions    kaons  protons\n");
	for(iq=0;iq<1;iq++){
		string OSCARfilename=dataroot+qualifier[iq]+"/bb_oscar.dat";
		int nevents=5000;
		double r[4],p[4];
		double mass,rdummy1,rdummy2,pt,ptbar_pion=0.0,ptbar_proton=0.0,ptbar_kaon=0.0;
		int norm_pion=0,norm_proton=0,norm_kaon=0,eventcount=0;
		int ident,idummy,i,alpha,ibonus;
		int ndummy_header=3,ndummy_betweenevents=0;
		int npart,npartmax,ipart,ievent;

		char dummy[160];
		printf("Opening %s\n",OSCARfilename.c_str());
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
				for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0;
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);
			//printf("ident=%d, pt=%g, p[1]=%g, p[2]=%g\n",ident,pt,p[1],p[2]);
			/* SET KINEMATIC VARIABLES */
				if(ident==211 || ident==-211 || ident==111){
					ptbar_pion+=pt; 
					norm_pion+=1;
				}
				if(ident==311 || ident==-311 || ident==321 || ident==-321){
					ptbar_kaon+=pt; 
					norm_kaon+=1;
				}
				if(ident==2112 || ident==-2112 || ident==2212 || ident==-2212){
					ptbar_proton+=pt; 
					norm_proton+=1;
				}
			}
			for(i=0;i<ndummy_betweenevents;i++){
				if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
			}
			npartmax=0;
			eventcount+=1;
		} while(!feof(oscarfile) && ievent<nevents);
		ptbar_pion=ptbar_pion/double(norm_pion);
		ptbar_kaon=ptbar_kaon/double(norm_kaon);
		ptbar_proton=ptbar_proton/double(norm_proton);
		printf("npions=%d, nkaons=%d, nprotons=%d\n",norm_pion,norm_kaon,norm_proton);
		printf("ptbar_pion=%g, ptbar_kaon=%g, ptbar_proton=%g\n",ptbar_pion,ptbar_kaon,ptbar_proton);
		fprintf(fptr,"%12s  %7.1f  %7.1f  %7.1f\n",qualifier[iq].c_str(),ptbar_pion,ptbar_kaon,ptbar_proton);
		printf("Processed %d events\n",eventcount);
		fclose(oscarfile);
	}
	fclose(fptr);
	return 0;
}
