#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

const double PI=4.0*atan(1.0);
const double HBARC=197.3269602;

using namespace std;
void getsegments(string dumbo,int &nsegments,string *segment);

int main(int argc,char *argv[]){
  FILE *csvfile,*fptr;
	csvfile=fopen("../rhicstat/figs/mcmctrace.csv","r");
	char dirname[100],dummy[400];
	string filename,command;
	string dumbo;
	string segment[50],name[50];
	int ipar,npars,nsegments,isegment,irun=0,iskip,nskip=8333,isample=0,nsamples=40;
	fgets(dummy,400,csvfile);
	dumbo=dummy;
	getsegments(dumbo,npars,name);
	
	iskip=0;
	while(!feof(csvfile)&&isample<nsamples){
		iskip+=1;
		fgets(dummy,400,csvfile);
		if(iskip==nskip){
			iskip=0;
			isample+=1;
			dumbo=dummy;
			getsegments(dumbo,nsegments,segment);
			sprintf(dirname,"model_output_posterior/run%04d",irun);
			command="mkdir -p "+string(dirname);
			system(command.c_str());
			filename=string(dirname)+"/parameters.dat";
			printf("writing to %s\n",filename.c_str());
			fptr=fopen(filename.c_str(),"w");
			for(ipar=0;ipar<npars;ipar++){
				fprintf(fptr,"%s %s\n",name[ipar].c_str(),segment[ipar].c_str());
			}
			fclose(fptr);
			irun+=1;
		}
	}
	printf("isample=%d\n",isample);

  return 0;
}


void getsegments(string dumbo,int &nsegments,string *segment){
  int isegment,length,firstchar,lastchar,quotecheck,ichar;
  length=dumbo.length();
  firstchar=0;
  nsegments=0;
  quotecheck=1;
  for(ichar=0;ichar<=length;ichar++){
    if(dumbo[ichar]=='\"') quotecheck=-quotecheck;
    if((quotecheck==1 && dumbo[ichar]==',') || ichar==length-1){
      lastchar=ichar-1;
      if(dumbo[firstchar]=='\"') firstchar+=1;
      if(dumbo[lastchar]=='\"') lastchar-=1;
      segment[nsegments]=dumbo.substr(firstchar,1+lastchar-firstchar);
      firstchar=ichar+1;
      nsegments+=1;
    }
  }
}