#include "bal.h"
using namespace std;

void CBalance::CalcMixed(){
	if(USE_SAME_ACCEPTANCE){
		acceptance_i=acceptance;
		acceptance_j=acceptance;
	}
	double dcai[4],dcaj[4];
	CPart parti0[10],parti[10],partj0[10],partj[10];
	int ibal,ires,jres,ibin,ni,nj,i,j,k,icode,jcode;
	bool codecheck,idecay,jdecay,accepti,acceptj,findbal;
	double w,effi,effj,phi_j,phi_i,dy,etamax;
	vector<double> norm;
	int isample;
	vector<vector<double> > Mixed;
	vector<vector<long long int> > sigmaN;
	vector<vector<double> > sigmaE; 
	Mixed.resize(NBAL);
	sigmaN.resize(NBAL);
	sigmaE.resize(NBAL);
	norm.resize(NBAL);
	for(ibal=0;ibal<NBAL;ibal++){
		Mixed[ibal].resize(NBINS);
		sigmaN[ibal].resize(NBINS);
		sigmaE[ibal].resize(NBINS);
		norm[ibal]=0.0;
		for(ibin=0;ibin<NBINS;ibin++){
			Mixed[ibal][ibin]=0.0;
			sigmaN[ibal][ibin]=0;
			sigmaE[ibal][ibin]=0.0;
		}
	}

	for(isample=0;isample<NSAMPLE;isample++){
		ires=GetIres();
		jres=GetIres();
		icode=code[ires];
		jcode=code[jres];
		if(abs(icode)==311){
			if(randy->ran()<0.5)
				icode=-icode;
		}
		codecheck=GetProducts(icode,ni,parti0,idecay);
		for(i=0;i<ni;i++){
			icode=parti0[i].resinfo->code;
			if(abs(icode)==311){
				if(randy->ran()<0.5)
					icode=-icode;
			}
			if(icode!=22 && icode!=111 && icode!=-311){
				BoostPart(&parti0[i],&parti[i],idecay);
				CalcDCA(idecay,&parti[i],dcai);
				acceptance_i->CalcAcceptance(accepti,effi,&parti[i],0,dcai);
				if(accepti){
					codecheck=GetProducts(jcode,nj,partj0,jdecay);
					if(abs(jcode)==311){
						if(randy->ran()<0.5)
							jcode=-jcode;
					}
					for(j=0;j<nj;j++){
						jcode=partj0[j].resinfo->code;
						if(jcode!=22 && jcode!=111 && jcode!=-311){
							findbal=false;
							ibal=0;
							while(findbal==false && ibal<NBAL){
								ibal+=1;
								if((abs(icode)==abs(ICODE[ibal])) && (abs(jcode)==abs(JCODE[ibal]))){
									findbal=true;
								}
								if((abs(icode)==abs(JCODE[ibal])) && (abs(jcode)==abs(ICODE[ibal]))){
									findbal=true;
								}
							}
							if(findbal==true){
								BoostPart(&partj0[j],&partj[j],jdecay);
								ibin=lrint(floor(fabs(partj[j].y-parti[i].y)/DELY));
								if(ibin<NBINS)
									norm[ibal]+=effi;
								CalcDCA(jdecay,&partj[j],dcaj);
								acceptance_j->CalcAcceptance(acceptj,effj,&partj[j],0,dcaj);
								if(acceptj){
									ibin=lrint(floor(fabs(partj[j].y-parti[i].y)/DELY));
									if(ibin<NBINS){
										Mixed[ibal][ibin]+=effi*effj;
										sigmaN[ibal][ibin]+=1;
										sigmaE[ibal][ibin]+=effi*effi*effj*effj;
									}
								}
							}
						}
					}
				}
			}
		}
		if(10*(isample+1)%NSAMPLE==0)
			printf("finished %g percent\n",(isample+1)*100.0/double(NSAMPLE));
	}
	string dirname="model_output/"+PARSDIRNAME;
	string command="mkdir -p "+dirname;
	system(command.c_str());
	FILE *fptr;
	double deltaE,deltaN,Ebar,dM,SN;
	etamax=acceptance_i->ETAMAX;
	char filename[100];
	for(ibal=1;ibal<NBAL;ibal++){
		sprintf(filename,"%s/mixed_%d_%d.dat",dirname.c_str(),abs(ICODE[ibal]),abs(JCODE[ibal]));
		fptr=fopen(filename,"w");
		norm[ibal]=norm[ibal]/double(NBINS);
		for(ibin=0;ibin<NBINS;ibin++){
			dy=(ibin+0.5)*DELY;
			if(sigmaN[ibal][ibin]>0){
				SN=double(sigmaN[ibal][ibin]);
				sigmaE[ibal][ibin]=sigmaE[ibal][ibin]/SN;
				Ebar=Mixed[ibal][ibin]/SN;
				sigmaE[ibal][ibin]-=Ebar*Ebar;
				deltaE=sqrt(sigmaE[ibal][ibin]/SN);
				deltaN=sqrt(SN)/SN;
				printf("sigmaN=%g, Mixed=%g, norm=%g\n",SN,Mixed[ibal][ibin],norm[ibal]);
				printf("%5.2f 1.0 =? %g, sigmaN=%lld, Ebar=%g, deltaE=%g\n",dy,Mixed[ibal][ibin]/(norm[ibal]*(2.0*etamax-dy)/etamax),sigmaN[ibal][ibin],Ebar,deltaE);
				Mixed[ibal][ibin]*=(0.5/norm[ibal]);
				dM=Mixed[ibal][ibin]*sqrt(pow(deltaN,2)+pow(deltaE,2));
				fprintf(fptr,"%5.2f %g %g\n",dy,Mixed[ibal][ibin],dM);
			}
		}
		fclose(fptr);
	}
}
