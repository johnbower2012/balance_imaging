#include "coral.h"

int main (int argc, char *argv[]){
	FILE *oscarfile,*results_file;
	string dataroot="/Users/pratt/data/bb/";
	string qualifier;
	char qread[100];
	printf("Enter qualifier: ");
	scanf("%s",qread);
	qualifier=qread;
	string OSCARfilename=dataroot+qualifier+"/oscar.dat";

	int nevents=3000000;
	const int Nbins=160;
	double r[4],p[4];
	double mt,dmt=0.025,binsize;
	int imt;
	double spectra_pi[Nbins]={0.0},spectra_p[Nbins]={0.0},spectra_k[Nbins]={0.0};
	double mass,rdummy1,rdummy2,pt,ptbar_pion=0.0,ptbar_proton=0.0,ptbar_kaon=0.0;
	double dNdeta=0.0;
	double norm_pion=0,norm_proton=0,norm_kaon=0;
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
		if (!feof(oscarfile)){
			printf("ievent=%d, idummy=%d, npartmax=%d\n",ievent,idummy,npartmax);
			for(npart=0;npart<npartmax;npart++){
				fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&ipart,&ident,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
				pt=sqrt(p[1]*p[1]+p[2]*p[2]);
				if(pt>0.0 && pt<5.0){
					if(ident==211 || ident==-211 || ident==111){
						mass=0.1396;
						mt=sqrt(mass*mass+pt*pt)-mass;
						imt=lrint(floor(mt/dmt));
						ptbar_pion+=pt;
						norm_pion+=1;
						if(imt<Nbins) spectra_pi[imt]+=1.0;
					}
					if(ident==311 || ident==-311 || ident==321 || ident==-321){
						mass=0.494;
						mt=sqrt(mass*mass+pt*pt)-mass;
						imt=lrint(floor(mt/dmt));
						ptbar_kaon+=pt; 
						norm_kaon+=1;
						if(imt<Nbins) spectra_k[imt]+=1.0;
					}
					if(ident==2112 || ident==-2112 || ident==2212 || ident==-2212){
						mass=0.9383;
						mt=sqrt(mass*mass+pt*pt)-mass;
						imt=lrint(floor(mt/dmt));
						ptbar_proton+=pt; 
						norm_proton+=1;
						if(imt<Nbins) spectra_p[imt]+=1.0;
					}
					dNdeta+=pt/sqrt(pt*pt+mass*mass);
				}
			}
		}
		for(i=0;i<ndummy_betweenevents;i++){
			if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
		}
	} while(!feof(oscarfile) && ievent<nevents);
	nevents=ievent;
	printf("dNdeta=%g\n",dNdeta/double(nevents));
	fclose(oscarfile);

	string results_filename="results/spectra_"+qualifier+".dat";
	results_file=fopen(results_filename.c_str(),"w");
	ptbar_pion=ptbar_pion/double(norm_pion);
	ptbar_kaon=ptbar_kaon/double(norm_kaon);
	ptbar_proton=ptbar_proton/double(norm_proton);
	//printf("Nevents=%d, npions/event=%g, nkaons/event=%g, nprotons/event=%g\n",
		//nevents,norm_pion/double(nevents),norm_kaon/double(nevents),norm_proton/double(nevents));
	printf("ptbar_pion=%g, ptbar_kaon=%g, ptbar_proton=%g\n",ptbar_pion,ptbar_kaon,ptbar_proton);
	printf("npions/event=%g, nkaons/event=%g, nbaryons/event=%g\n",norm_pion/double(nevents),norm_kaon/double(nevents),norm_proton/double(nevents));
	fprintf(results_file,"! ptbar_pion=%7.4f  ptbar_kaon=%7.4f  ptbar_proton=%7.4f\n",ptbar_pion,ptbar_kaon,ptbar_proton);
	fprintf(results_file,"!   mt      pions       kaons      protons\n");
	for(imt=0;imt<Nbins;imt++){
		mass=0.1396;
		binsize=double(nevents)*PI*(pow(mass+(imt+1.0)*dmt,2)-pow(mass+double(imt)*dmt,2));
		spectra_pi[imt]=spectra_pi[imt]/(3.0*binsize);
		mass=0.494;
		binsize=double(nevents)*PI*(pow(mass+(imt+1.0)*dmt,2)-pow(mass+double(imt)*dmt,2));
		spectra_k[imt]=spectra_k[imt]/(4.0*binsize);
		mass=0.9383;
		binsize=double(nevents)*PI*(pow(mass+(imt+1.0)*dmt,2)-pow(mass+double(imt)*dmt,2));
		spectra_p[imt]=spectra_p[imt]/(4.0*binsize);
		mt=(0.5+imt)*dmt;
		printf("%6.4f %10.3e %10.3e %10.3e\n",mt,spectra_pi[imt],spectra_k[imt],spectra_p[imt]);
		fprintf(results_file,"%6.4f %10.3e %10.3e %10.3e\n",mt,spectra_pi[imt],spectra_k[imt],spectra_p[imt]);
	}
	fclose(results_file);
	
	norm_proton=0.0;
	ptbar_proton=0.0;
	for(imt=0;imt<Nbins;imt++){
		mass=0.9383;
		binsize=4.0*PI*(pow(mass+(imt+1.0)*dmt,2)-pow(mass+double(imt)*dmt,2));
		norm_proton+=binsize*spectra_p[imt];
		ptbar_proton+=binsize*spectra_p[imt]*sqrt(pow(mass+(imt+0.5)*dmt,2)-mass*mass);
	}
	printf("to check: number baryons=%g, ptbar_proton=%g\n",norm_proton,ptbar_proton/norm_proton);


	return 0;
}
