#include <cstdio>
#include <cstdlib>
#include <cmath>

void WriteCorrectedInfo(string filename,string mixed_filename,string corrected_filename);

using namespace std;

int main(int argc, char * const argv[]){
	double uncertainty,y,B,Bmixed;
	double dummy1,dummy2,dummy3;
	double sigmaB,sigmaBcorr,sigmaMixed;
	const double MixedMin=1.0E-2;
	char dummy[300];
	int ibin;
	FILE *starfile,*correctedfile,*mixedfile,*uncorrectedfile;
	
	starfile =fopen("../exp_data/star_pipi.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("exp_data_corrected/star_pipi.dat","w");
	int NBINS=18;
	for(ibin=0;ibin<NBINS;ibin++){
		if(ibin==0)
			fgets(dummy,300,starfile);
		else
			fscanf(starfile,"%lf %lf %lf",&y,&B,&sigmaB);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15 && ibin>0)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(starfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	starfile =fopen("../exp_data/star_KK.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("exp_data_corrected/star_KK.dat","w");
	for(ibin=0;ibin<NBINS;ibin++){
		if(ibin==0)
			fgets(dummy,300,starfile);
		else
			fscanf(starfile,"%lf %lf %lf",&y,&B,&sigmaB);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15 && ibin>0)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(starfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	starfile =fopen("../exp_data/star_ppbar.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("exp_data_corrected/star_ppbar.dat","w");
	for(ibin=0;ibin<NBINS;ibin++){
		if(ibin==0)
			fgets(dummy,300,starfile);
		else
			fscanf(starfile,"%lf %lf %lf",&y,&B,&sigmaB);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15 && ibin>0)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(starfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	double BQGP,BLOCAL,BRES;
	// Now do default model output
	uncorrectedfile =fopen("model_output/default/I211_J211.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/default/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/default/I321_J321.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/default/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/default/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/default/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);

	// Now do onescale (sigma=0.0) model output
	uncorrectedfile =fopen("model_output/onescale_0.0/I211_J211.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.0/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.0/I321_J321.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.0/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.0/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.0/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	// Now do onescale (sigma=0.25) model output
	uncorrectedfile =fopen("model_output/onescale_0.25/I211_J211.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.25/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.25/I321_J321.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.25/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.25/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.25/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	// Now do onescale (sigma=0.5) model output
	uncorrectedfile =fopen("model_output/onescale_0.5/I211_J211.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.5/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.5/I321_J321.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.5/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_0.5/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_0.5/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);

	// Now do onescale (sigma=1.0) model output
	uncorrectedfile =fopen("model_output/onescale_1.0/I211_J211.dat","r");
	mixedfile=fopen("model_output/default/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_1.0/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf %lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_1.0/I321_J321.dat","r");
	mixedfile=fopen("model_output/default/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_1.0/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/onescale_1.0/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/default/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/onescale_1.0/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);

	// Now do perfect acceptance
	
	uncorrectedfile =fopen("model_output/perfect/I211_J211.dat","r");
	mixedfile=fopen("model_output/perfect/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/perfect/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/perfect/I321_J321.dat","r");
	mixedfile=fopen("model_output/perfect/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/perfect/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/perfect/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/perfect/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/perfect/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	// Now do semiperfect acceptance
	
	uncorrectedfile =fopen("model_output/semiperfect/I211_J211.dat","r");
	mixedfile=fopen("model_output/semiperfect/mixed_211_211.dat","r");
	correctedfile=fopen("output_posterior_corrected/semiperfect/I211_J211.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/semiperfect/I321_J321.dat","r");
	mixedfile=fopen("model_output/semiperfect/mixed_321_321.dat","r");
	correctedfile=fopen("output_posterior_corrected/semiperfect/I321_J321.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	uncorrectedfile =fopen("model_output/semiperfect/I2212_J2212.dat","r");
	mixedfile=fopen("model_output/semiperfect/mixed_2212_2212.dat","r");
	correctedfile=fopen("output_posterior_corrected/semiperfect/I2212_J2212.dat","w");
	for(int iskip=0;iskip<6;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		fscanf(uncorrectedfile,"%lf %lf %lf %lf %lf%lf",&y,&B,&sigmaB,&dummy1,&dummy2,&dummy3);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(uncorrectedfile);
	fclose(mixedfile);
	fclose(correctedfile);
	
	return 0;
}

void WriteCorrectedInfo(int NBINS,int NSKIP,string uncorrected_filename,string mixed_filename,string corrected_filename){
	double uncertainty,y,B,Bmixed;
	double dummy1,dummy2,dummy3;
	double sigmaB,sigmaBcorr,sigmaMixed;
	const double MixedMin=1.0E-2;
	char dummy[300];
	int ibin;
	FILE *starfile,*correctedfile,*mixedfile,*uncorrectedfile;
	
	starfile =fopen(uncorrected_filename.c_str(),"r");
	mixedfile=fopen(mixed_filename.c_str(),"r");
	correctedfile=fopen(corrected_filename.c_str(),"w");
	for(int iskip=0;iskip<NSKIP;iskip++){
		fgets(dummy,100,uncorrectedfile);
		fprintf(correctedfile,"%s\n",dummy);
	}
	for(ibin=0;ibin<NBINS;ibin++){
		if(ibin==0)
			fgets(dummy,300,starfile);
		else
			fscanf(starfile,"%lf %lf %lf",&y,&B,&sigmaB);
		fscanf(mixedfile,"%lf %lf %lf",&y,&Bmixed,&sigmaMixed);
		sigmaMixed=sqrt(sigmaMixed*sigmaMixed+MixedMin*MixedMin);
		sigmaBcorr=(B/Bmixed)*sqrt(pow(sigmaB/B,2)+pow(sigmaMixed/Bmixed,2));
		if(fabs(sigmaBcorr)<0.15 && ibin>0)
			fprintf(correctedfile,"%6.2f %g %g %g %g\n",y,B/Bmixed,fabs(sigmaBcorr),B,Bmixed);
	}
	fclose(starfile);
	fclose(mixedfile);
	fclose(correctedfile);
}
	
}
