#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
using namespace std;

void GetBFMoments(string filename,double &mean,double &sigma2,double &norm);
double GetError(double y,double bf,double balnorm){
	double sigma;
	sigma=0.001;
	sigma=sqrt(sigma*sigma+0.12*bf*0.12*bf);
	return sigma;
	//return ERROR_FACTOR;
}

int main(int argc, char * const argv[]){
	bool USE_BALNORM=false;
	const int nbal=4;
	int nruns=40,irun,iread,ibal;
	double y,mean,sigma2,norm,rdummy1,rdummy2,rdummy3,rdummy4,balnorm=1.0,error,meanpt;
	double bf[19],dely[17];
	char dummy[100],wholeline[100];
	char command[100],readfilename[100],writefilename[100],Irun[10];
	string balstring[nbal]={"I211_J211","I321_J321","I321_J2212","I2212_J2212"};
	string balobs[nbal]={"Bpipi","BKK","BpK","Bppbar"};
	char meanptname[100];
	FILE *wfptr,*fptr;
	for(irun=0;irun<nruns;irun++){
		printf("irun=%d\n",irun);
		sprintf(command,"mkdir -p model_output/run%04d",irun);
		system(command);
		sprintf(writefilename,"model_output/run%04d/results.dat",irun);
		wfptr=fopen(writefilename,"w");
		for(ibal=0;ibal<nbal;ibal++){
			sprintf(readfilename,"model_output/run%04d/%s.dat",irun,balstring[ibal].c_str());
			fptr=fopen(readfilename,"r");
			for(iread=0;iread<7;iread++)
				fgets(dummy,100,fptr);
			if(USE_BALNORM)
				balnorm=0.0;
			for(iread=0;iread<17;iread++){
				fscanf(fptr,"%lf %lf %lf %lf %lf %lf",&dely[iread],&bf[iread],&rdummy1,&rdummy2,&rdummy3,&rdummy4);
				balnorm+=bf[iread]*0.1;
			}
			if(USE_BALNORM)
				fprintf(wfptr,"%s_norm %8.5f %8.5f\n",balobs[ibal].c_str(),balnorm,0.04*balnorm+0.001);
			for(iread=0;iread<17;iread++){
				y=(iread+1.5)*0.1;
				error=GetError(y,bf[iread],balnorm);
				if(USE_BALNORM)
					fprintf(wfptr,"%s_ybin%d %8.5f %8.5f\n",balobs[ibal].c_str(),iread+1,bf[iread]/balnorm,error);
				else
					fprintf(wfptr,"%s_ybin%d %8.5f %8.5f\n",balobs[ibal].c_str(),iread+1,bf[iread],error);
			}
			fclose(fptr);
		}
		
		sprintf(readfilename,"model_output/run%04d/meanpt.dat",irun);
		fptr=fopen(readfilename,"r");
		for(iread=0;iread<3;iread++){
			fscanf(fptr,"%s %lf",meanptname,&meanpt);
			fprintf(wfptr,"%s %g %g\n",meanptname,meanpt,0.06*meanpt);
		}
		fclose(fptr);

		fclose(wfptr);
	}
	return 0;
}
