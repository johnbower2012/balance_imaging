#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>

using namespace std;

void WriteCorrectedInfo(bool modeldata,int NBINS,int NSKIP,string filename,string mixed_filename,string corrected_filename);

int main(int argc, char * const argv[]){

	WriteCorrectedInfo(false,18,0,"../exp_data/star_pipi.dat","model_output/default/mixed_211_211.dat","exp_data_corrected/star_pipi.dat");
	
	WriteCorrectedInfo(false,18,0,"../exp_data/star_KK.dat","model_output/default/mixed_321_321.dat","exp_data_corrected/star_KK.dat");
	
	WriteCorrectedInfo(false,18,0,"../exp_data/star_ppbar.dat","model_output/default/mixed_2212_2212.dat","exp_data_corrected/star_ppbar.dat");

	WriteCorrectedInfo(false,18,0,"../exp_data/star_pK.dat","model_output/default/mixed_321_2212.dat","exp_data_corrected/star_pK.dat");

	//
	const int nruns=12;
	int irun,ij;
	string runName[nruns]={"default","QoverS_0.1","QoverS_0.25","SoverU_0.5","SoverU_1.4",
	"sigmaB_0.1","sigmaB_0.7","onescale_0.0","onescale_0.25","onescale_0.5","onescale_1.0",
	"testrun"};
	string istring[4]={"211","321","2212","321"};
	string jstring[4]={"211","321","2212","2212"};
	string ufilename,cfilename,mfilename,cdirname,command;
	for(irun=0;irun<nruns;irun++){
		for(ij=0;ij<4;ij++){
			ufilename="model_output/"+runName[irun]+"/I"+istring[ij]+"_J"+jstring[ij]+".dat";
			mfilename="model_output/default/mixed_"+istring[ij]+"_"+jstring[ij]+".dat";
			cdirname="output_posterior_corrected/"+runName[irun];
			command="mkdir -p "+cdirname;
			system(command.c_str());
			cfilename=cdirname+"/I"+istring[ij]+"_J"+jstring[ij]+".dat";
			//printf("%s, %s, %s\n",ufilename.c_str(),mfilename.c_str(),cfilename.c_str());
			WriteCorrectedInfo(true,18,6,ufilename,mfilename,cfilename);
		}
	}
	return 0;
}

void WriteCorrectedInfo(bool modeldata,int NBINS,int NSKIP,string uncorrected_filename,string mixed_filename,string corrected_filename){
	double uncertainty,y,B,Bmixed;
	double dummy1,dummy2,dummy3,bmax;
	double sigmaB,sigmaBcorr,sigmaMixed,BQGP,Blocal,Bdecay;
	const double MixedMin=1.0E-2;
	char dummy[300];
	int ibin;
	FILE *correctedfile,*mixedfile,*uncorrectedfile;
	
	uncorrectedfile =fopen(uncorrected_filename.c_str(),"r");
	mixedfile=fopen(mixed_filename.c_str(),"r");
	correctedfile=fopen(corrected_filename.c_str(),"w");
	for(int iskip=0;iskip<NSKIP;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s",dummy);
	}
	bmax=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		if(ibin==0)
			fgets(dummy,300,uncorrectedfile);
		else{
			if(modeldata){
				fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&BQGP,&Blocal,&Bdecay);
			}
			else
				fscanf(uncorrectedfile,"%lf %lf %lf",&y,&B,&sigmaB);
			fgets(dummy,300,uncorrectedfile);
		}
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		fgets(dummy,300,mixedfile);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		if(fabs(B/Bmixed)>bmax)
			bmax=fabs(B/Bmixed);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(modeldata){
			if(fabs(sigmaBcorr)<0.5*bmax && ibin>0)
				fprintf(correctedfile,"%6.2f %g %g %g %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed,
			BQGP/Bmixed,Blocal/Bmixed,Bdecay/Bmixed);
		}
		else{
			if(fabs(sigmaBcorr)<0.5*bmax && ibin>0)
				fprintf(correctedfile,"%6.2f %g %g\n",y,B/Bmixed,fabs(sigmaBcorr));
		}
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
}
	