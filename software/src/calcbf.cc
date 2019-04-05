#include "bal.h"
using namespace std;

void CBalance::CalcBF(){
	if(USE_SAME_ACCEPTANCE){
		acceptance_i=acceptance;
		acceptance_j=acceptance;
	}
	double dcai[4],dcaj[4];
	double meanpt[3]={0.0}, mult[3]={0.0};
	CPart parti0[10],parti[10],partj0[10],partj[10];
	int ibal,ires,jres,ibin,iphibin,iphi_j,iplane,ni,nj,i,j,k;
	int ipart,jpart,nparts,mparts,icode,jcode,iicode,jjcode;
	bool codecheck,idecay,jdecay,accepti,acceptj;
	int npiplus=0,npiminus=0,npipairs=0;
	bool paircheckplus,paircheckminus;
	double w,effi,effj,phi_j,phi_i;
	int isample,ab;
	double DELPHI=PI/double(NPHIBINS);
	CResInfo *resinfoi,*resinfoj;
	WriteBFVar();
	Mixed.resize(NBAL);
	sigmaN.resize(NBAL);
	sigmaE.resize(NBAL);
	dB.resize(NBAL);
	for(ibal=0;ibal<NBAL;ibal++){
		Mixed[ibal].resize(NBINS);
		sigmaN[ibal].resize(NBINS);
		sigmaE[ibal].resize(NBINS);
		dB[ibal].resize(NBINS);
		for(ibin=0;ibin<NBINS;ibin++){
			Mixed[ibal][ibin]=0.0;
			sigmaN[ibal][ibin]=0;
			sigmaE[ibal][ibin]=0.0;
		}
	}
	for(ibal=0;ibal<NBAL;ibal++){
		normalization[ibal]=v2[ibal]=v2s[ibal]=v2c[ibal]=cosdphi[ibal]=rhoi[ibal]=rhoj[ibal]=0.0;
		for(ibin=0;ibin<NBINS;ibin++)
			bf[ibal][ibin]=bf_mean[ibal][ibin]=bf_sq[ibal][ibin]=bf_had[ibal][ibin]=bf_resdecay[ibal][ibin]=bf_qgp[ibal][ibin]=0.0;
		for(iphibin=0;iphibin<NPHIBINS;iphibin++)
			bf_phi[ibal][iphibin]=bf_had_phi[ibal][iphibin]=bf_resdecay_phi[ibal][iphibin]=bf_qgp_phi[ibal][iphibin]=0.0;
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
		if(abs(jcode)==311){
			if(randy->ran()<0.5)
				jcode=-jcode;
		}
		codecheck=GetProducts(icode,ni,parti0,idecay);
		codecheck=GetProducts(jcode,nj,partj0,jdecay);
		paircheckplus=paircheckminus=false;
		for(i=0;i<ni;i++){
			icode=parti0[i].resinfo->code;
			if(abs(icode)==211 || abs(icode)==321 || abs(icode)==2212 || abs(icode)==2112){
				if(icode==211){
					paircheckplus=true;
					npiplus+=1;
				}
				if(icode==-211){
					paircheckminus=true;
					npiminus+=1;
				}
				if(paircheckminus && paircheckplus){
					paircheckminus=paircheckplus=false;
					npipairs+=1;
				}
				BoostPart(&parti0[i],&parti[i],idecay);
				if(icode==211 || icode==-211){
					mult[0]+=1.0;
					meanpt[0]+=sqrt(parti[i].p[1]*parti[i].p[1]+parti[i].p[2]*parti[i].p[2]);
				}
				if(icode==321 || icode==-321){
					mult[1]+=1.0;
					meanpt[1]+=sqrt(parti[i].p[1]*parti[i].p[1]+parti[i].p[2]*parti[i].p[2]);
				}
				if(icode==2212 || icode==-2212 || icode==2112 || icode==-2112){
					mult[2]+=1.0;
					meanpt[2]+=sqrt(parti[i].p[1]*parti[i].p[1]+parti[i].p[2]*parti[i].p[2]);
				}
			}
		}
		for(j=0;j<nj;j++){
			jcode=partj0[j].resinfo->code;
			if(abs(jcode)==211 || abs(jcode)==321 || abs(jcode)==2212 || abs(jcode)==2112){
				BoostPart(&partj0[j],&partj[j],jdecay);
				if(jcode==211 || jcode==-211){
					mult[0]+=1.0;
					meanpt[0]+=sqrt(partj[j].p[1]*partj[j].p[1]+partj[j].p[2]*partj[j].p[2]);
				}
				if(jcode==321 || jcode==-321){
					mult[1]+=1.0;
					meanpt[1]+=sqrt(partj[j].p[1]*partj[j].p[1]+partj[j].p[2]*partj[j].p[2]);
				}
				if(jcode==2212 || jcode==-2212 || jcode==2112 || jcode==-2112){
					mult[2]+=1.0;
					meanpt[2]+=sqrt(partj[j].p[1]*partj[j].p[1]+partj[j].p[2]*partj[j].p[2]);
				}
			}
		}
		for(i=0;i<ni;i++){
			if(icode!=22 && icode!=111 && icode!=-311){
				// Do singles and denoms
				BoostPart(&parti0[i],&parti[i],idecay);
				CalcDCA(idecay,&parti[i],dcai);
				acceptance_i->CalcAcceptance(accepti,effi,&parti[i],0,dcai);
				// Note that meanpt is for perfect acceptance, whereas v2 within acceptance and weighted by eff (easily changed)
				for(ibal=0;ibal<NBAL;ibal++){
					if((ibal!=0 && (abs(parti[i].resinfo->code)==abs(ICODE[ibal]))) || (ibal==0 && parti[i].resinfo->charge!=0)){
						// Calculate denominators (and v2 for grins)
						if(accepti){
							denom[ibal]+=effi*rhotot[NRES-1]/double(NSAMPLE);
							rhoi[ibal]+=effi*rhotot[NRES-1]/double(NSAMPLE);
							phi_i=atan2(fabs(parti[i].p[2]),fabs(parti[i].p[1]));
							v2[ibal]+=2.0*DELY*effj*cos(2.0*phi_i)/double(NSAMPLE);
							iplane=lrint(floor(2.0*NPHIBINS_RPLANE*phi_i/PI));
							if(iplane<0 || iplane >=NPHIBINS_RPLANE){
								printf("iplane off, =%d\n",iplane);
								exit(1);
							}
							denom_rplane[ibal][iplane]+=effi*rhotot[NRES-1]/double(NSAMPLE);
						}				
					}
				}
				// Do BF numerator due to resonance decays
				w=rhotot[NRES-1]/double(NSAMPLE);
				if(ni>1){
					for(j=0;j<i;j++){
						if(i!=j){
							IncrementBFArrays(&parti0[i],idecay,&parti0[j],idecay,w,-1);
							IncrementBFArrays(&parti0[j],idecay,&parti0[i],idecay,w,-1);
						}
					}
				}
				// Do BF numerators due to inter-hadron correlations
				for(j=0;j<nj;j++){
					jcode=partj0[j].resinfo->code;
					if(jcode!=22 && jcode!=111 && jcode!=-311){
						w=0.0;
						for(ab=0;ab<4;ab++){
							w=0.25*rhotot[NRES-1]*rhotot[NRES-1]*wB[ires][jres][ab]/double(NSAMPLE);
							IncrementBFArrays(&parti0[i],idecay,&partj0[j],jdecay,w,ab);
							IncrementBFArrays(&partj0[j],jdecay,&parti0[i],idecay,w,ab);
						}
					}
				}
			}
		}
		for(ibal=0;ibal<NBAL;ibal++){
		  double Denom = denom[ibal]*DELY;
		  if(Denom!=0.0){
		    for(ibin=0;ibin<NBINS;ibin++){
		      double BF=bf_had[ibal][ibin]+bf_qgp[ibal][ibin]+bf_resdecay[ibal][ibin];
		      BF = BF/Denom;
		      bf_mean[ibal][ibin] += BF;
		      bf_sq[ibal][ibin] += BF*BF;
		    }
		  }
		}
		if(10*(isample+1)%NSAMPLE==0){
			printf("finished %g percent\n",(isample+1)*100.0/double(NSAMPLE));
			WriteBFVar(isample+1);
		}
		
	}
	// NORMALIZE BFs
	meanpt_pion=meanpt[0]/mult[0];
	meanpt_kaon=meanpt[1]/mult[1];
	meanpt_proton=meanpt[2]/mult[2];
	printf("<pt>= %g for pions, %g for kaons, %g for nucleons\n",meanpt_pion,meanpt_kaon,meanpt_proton);

	double deltaE,deltaN,Ebar,d ,SN;
	for(ibal=0;ibal<NBAL;ibal++){
		denom[ibal]*=DELY;
		normalization[ibal]=0.0;
		for(ibin=0;ibin<NBINS;ibin++){
			//bf_had[ibal][ibin]*=rhotot[NRES-1]/(rhoi[ibal]);
			//bf_qgp[ibal][ibin]*=rhotot[NRES-1]/(rhoi[ibal]);
			bf[ibal][ibin]=bf_had[ibal][ibin]+bf_qgp[ibal][ibin]+bf_resdecay[ibal][ibin];
			bf[ibal][ibin]=bf[ibal][ibin]/denom[ibal];
			bf_had[ibal][ibin]=bf_had[ibal][ibin]/denom[ibal];
			bf_qgp[ibal][ibin]=bf_qgp[ibal][ibin]/denom[ibal];
			bf_resdecay[ibal][ibin]=bf_resdecay[ibal][ibin]/denom[ibal];
			normalization[ibal]+=bf[ibal][ibin]*DELY;
			if(sigmaN[ibal][ibin]>0){
				SN=double(sigmaN[ibal][ibin]);
				sigmaE[ibal][ibin]=sigmaE[ibal][ibin]/SN;
				Ebar=Mixed[ibal][ibin]/SN;
				sigmaE[ibal][ibin]-=Ebar*Ebar;
				deltaE=sqrt(sigmaE[ibal][ibin]/SN);
				deltaN=sqrt(SN)/SN;
				dB[ibal][ibin]=bf[ibal][ibin]*sqrt(pow(deltaN,2)+pow(deltaE,2));
				//dB[ibal][ibin]=bf[ibal][ibin]*deltaN/SN;
			}
		}
		printf("ibal=%d, normalization=%g\n",ibal,normalization[ibal]);
		v2[ibal]=v2[ibal]/denom[ibal];
		cosdphi[ibal]=cosdphi[ibal]*DELY/denom[ibal];
		v2s[ibal]=v2s[ibal]*DELY/denom[ibal];
		v2c[ibal]=(v2c[ibal]*DELY/denom[ibal])-v2[ibal]*cosdphi[ibal];
		gammap[ibal]=v2[ibal]*cosdphi[ibal]+v2c[ibal]-v2s[ibal];
	
		//printf("ICODE=%d, v2=%g, <cosdphi>=%g\nv2*<cosdphi>=%g, v2c=%g, v2s=%g, gammap*mult/2=%g\n",ICODE[ibal], meanpt[ibal],v2[ibal],cosdphi[ibal],v2[ibal]*cosdphi[ibal],v2c[ibal],v2s[ibal],0.5*gammap[ibal]);
		if(ALLCHARGES[ibal])
			printf("dNch/dy=%g\n",mult[ibal]/double(NSAMPLE));

		normalization[ibal]=0.0;
		denom[ibal]=denom[ibal]*DELPHI/DELY;
		for(iphibin=0;iphibin<NPHIBINS;iphibin++){
			bf_phi[ibal][iphibin]=bf_had_phi[ibal][iphibin]+bf_qgp_phi[ibal][iphibin]+bf_resdecay_phi[ibal][iphibin];
			bf_phi[ibal][iphibin]=bf_phi[ibal][iphibin]/denom[ibal];
			bf_had_phi[ibal][iphibin]=bf_had_phi[ibal][iphibin]/denom[ibal];
			bf_qgp_phi[ibal][iphibin]=bf_qgp_phi[ibal][iphibin]/denom[ibal];
			bf_resdecay_phi[ibal][iphibin]=bf_resdecay_phi[ibal][iphibin]/denom[ibal];
			normalization[ibal]+=bf_phi[ibal][iphibin]*DELPHI;
		}
		
		for(iplane=0;iplane<NPHIBINS_RPLANE;iplane++){
			denom_rplane[ibal][iplane]=denom_rplane[ibal][iplane]*DELPHI/DELY;
			for(iphibin=0;iphibin<2*NPHIBINS;iphibin++){
				bf_rplane[ibal][iplane][iphibin]=bf_had_rplane[ibal][iplane][iphibin]+bf_qgp_rplane[ibal][iplane][iphibin]+bf_resdecay_rplane[ibal][iplane][iphibin];
				bf_rplane[ibal][iplane][iphibin]=bf_rplane[ibal][iplane][iphibin]/denom_rplane[ibal][iplane];
			}
		}
	}
	printf("npiplus=%d, npiminus=%d, npipairs=%d\n",npiplus,npiminus,npipairs);
}

void CBalance::IncrementBFArrays(CPart *parti0,bool idecay,CPart *partj0,bool jdecay,double w0,int ab){
	CPart parti,partj;
	int iphibin,ibin,iphi_j,ibal;
	double delyij,phi_i,phi_j,delphi_ij,dcai[4],dcaj[4],effi,effj;
	double scheck,extraweight;
	bool accepti,acceptj;
	double sign,w,pmagi,pmagj,etai,etaj;
	double **b,**b_phi,***b_rplane;
	
	b=bf_resdecay;
	b_phi=bf_resdecay_phi;
	b_rplane=bf_resdecay_rplane;

	BoostParts(parti0,&parti,partj0,&partj,idecay,jdecay,ab,extraweight);
	CalcDCA(jdecay,&partj,dcaj);
	acceptance_j->CalcAcceptance(acceptj,effj,&partj,0,dcaj);
	CalcDCA(idecay,&parti,dcai);
	acceptance_i->CalcAcceptance(accepti,effi,&parti,0,dcai);
	
	if(accepti && acceptj){
		for(ibal=0;ibal<NBAL;ibal++){
			if(IJCheck(ibal,parti0,partj0)){
				if(ALLCHARGES[ibal]){
					pmagi=sqrt(parti.p[1]*parti.p[1]+parti.p[2]*parti.p[2]+parti.p[3]*parti.p[3]);
					pmagj=sqrt(partj.p[1]*partj.p[1]+partj.p[2]*partj.p[2]+partj.p[3]*partj.p[3]);
					etai=0.5*log((pmagi+parti.p[3])/(pmagi-parti.p[3]));
					etaj=0.5*log((pmagj+partj.p[3])/(pmagj-partj.p[3]));
					delyij=etai-etaj;
					sign=-1.0;
					scheck=parti0->resinfo->charge*partj0->resinfo->charge;
					if(scheck<0) sign=1.0;
				}
				else{
					delyij=parti.y-partj.y;
					sign=1.0;
					scheck=double(parti0->resinfo->code)*double(ICODE[ibal])*double(partj0->resinfo->code)*double(JCODE[ibal]);
					if(scheck<0)
						sign=-1.0;
				}
				w=w0*effi*effj*extraweight;

				ibin=lrint(floor(fabs(delyij/DELY)));
				if(ibin<NBINS){
					b[ibal][ibin]+=w*sign;
					sigmaN[ibal][ibin]+=1;
					sigmaE[ibal][ibin]+=w*w;
					Mixed[ibal][ibin]+=w;
				}
				phi_i=atan2(parti.p[2],parti.p[1]);
				phi_j=atan2(partj.p[2],partj.p[1]);
				delphi_ij=phi_i-phi_j;
				if(delphi_ij>PI) delphi_ij-=2.0*PI;
				if(delphi_ij<-PI) delphi_ij+=2.0*PI;
				cosdphi[ibal]+=w*sign*cos(delphi_ij);
				v2s[ibal]+=w*sign*sin(2.0*phi_j)*sin(delphi_ij);
				v2c[ibal]+=w*sign*cos(2.0*phi_j)*cos(delphi_ij);

				iphibin=lrint(floor(NPHIBINS*fabs(delphi_ij)/PI));
				if(iphibin<0 || iphibin>NPHIBINS){
					printf("iphibin=%d???\n",iphibin);
					exit(1);
				}
				b_phi[ibal][iphibin]+=w*sign;
				//now for binning by reaction plane
				
				if(phi_j<0.0){
					phi_i=-phi_i;
					phi_j=-phi_j;
				}
				if(phi_j>0.5*PI){
					phi_i=PI-phi_i;
					phi_j=PI-phi_j;
				}
				delphi_ij=phi_i-phi_j;
				if(delphi_ij>PI) delphi_ij-=2.0*PI;
				if(delphi_ij<-PI) delphi_ij+=2.0*PI;
				iphi_j=lrint(floor(NPHIBINS_RPLANE*(2.0*phi_j/PI)));
				if(iphi_j==0 || iphi_j==NPHIBINS_RPLANE-1){
					if(randy->ran()<0.5) delphi_ij=-delphi_ij; // randomize for bins at 0 and 90 deg
				}
				iphibin=lrint(floor(NPHIBINS*(delphi_ij+PI)/PI));
				if(iphibin<0 || iphibin>=2*NPHIBINS || iphi_j>=NPHIBINS_RPLANE || iphi_j<0){
					printf("iphibin=%d is out of bounds, NPHIBINS=%d and iphi_j=%d\n",iphibin,NPHIBINS,iphi_j);
					exit(1);
				}
				b_rplane[ibal][iphi_j][iphibin]+=w*sign;
			}
		}
	}
}

bool CBalance::IJCheck(int ibal,CPart *parti,CPart *partj){
	bool ijc=false;
	CResInfo *resinfoi=parti->resinfo,*resinfoj=partj->resinfo;
	if(ALLCHARGES[ibal]){
		if(resinfoi->charge!=0 && resinfoj->charge!=0) ijc=true;
	}
	else{
		if(abs(resinfoi->code)==abs(ICODE[ibal]) && abs(resinfoj->code)==abs(JCODE[ibal]))
			ijc=true;
	}
	return ijc;
}

int CBalance::GetIres(){
	double x=rhotot[NRES-1]*randy->ran();
	int ires=0;
	while(x>rhotot[ires]){
		ires+=1;
	}
	return ires;
}
