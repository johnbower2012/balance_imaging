using namespace std;
#include "bal.h"

void CBalance::CalcBF(){
	normalization=0.0;
	if(ALLCHARGES){
		printf("Calculating BF for ALLCHARGES\n");
		CalcBF_ALLCHARGES();
		return;
	}
	else if(GAUSSIAN){
		printf("Calculating BF for normalized GAUSSIAN, ICODE=%d, JCODE=%d\n",ICODE,JCODE);
		CalcBF_Gaussian();
		return;
	}
	else printf("Calculating BF for ICODE=%d, JCODE=%d\n",ICODE,JCODE);
	
	bool icheck[NRES];
	bool jcheck[NRES];
	double dcai[4],dcaj[4];
	CPart parti0[10],partj0[10],parti[10],partj[10];
	int ires,jres,ibin,iphibin,ni,nj,i,j,ievent;
	bool accepti,acceptj,codecheck,idecay,jdecay;
	double efficiency_jbar=0.0,norm_jbar=0.0;
	double effi,effj;
	double rhotot=0.0,wtot,sign,denom=0.0;
	for(ires=0;ires<NRES;ires++) rhotot+=rho_h[ires];
	double wB,w,w_doublecount=1.0,delyij,delphiij,bqgp=0.0,blocal=0.0;
	double ptbari=0.0,ptbarj=0.0,pti,ptj,ptswitch,ptnormi=0.0,ptnormj=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=bf_had[ibin]=bf_resdecay[ibin]=bf_qgp[ibin]=0.0;
	}
	for(iphibin=0;iphibin<NPHIBINS;iphibin++)
		bf_phi[iphibin]=bf_had_phi[iphibin]=bf_resdecay_phi[iphibin]=bf_qgp_phi[iphibin]=0.0;
	if(abs(ICODE)!=abs(JCODE)) w_doublecount=0.5;
	CResInfo *resinfoi;
	for(ires=0;ires<NRES;ires++){
		reslist->GetResInfoptr(code[ires],resinfoi);
		if(resinfoi->code==ICODE || resinfoi->code==JCODE) resinfoi->Print();
		icheck[ires]=resinfoi->CheckForDaughters(ICODE);
		if(icheck[ires]==false) icheck[ires]=resinfoi->CheckForDaughters(-ICODE);
		jcheck[ires]=resinfoi->CheckForDaughters(JCODE);
		if(jcheck[ires]==false) jcheck[ires]=resinfoi->CheckForDaughters(-JCODE);
	}
	
	for(jres=0;jres<NRES;jres++){
		printf("Calculating for jres, ID=%d XXXXXXXX\n",code[jres]);
		if(jcheck[jres] || icheck[jres]){
			for(ievent=0;ievent<NEVENTS;ievent++){
				codecheck=GetProducts(code[jres],nj,partj0,jdecay);
				if(nj>1){
						//printf("-----------------------\n");
						//for(j=0;j<nj;j++) printf("eta=%g\n",partj0[j].eta);
						//printf("-----------------------\n");
				}
				if(codecheck){
					for(j=0;j<nj;j++){
						if(abs(partj0[j].resinfo->code)==abs(JCODE)){
							BoostPart(&partj0[j],&partj[j],jdecay);
							CalcDCA(jdecay,&partj[j],dcaj);
							acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
							denom+=DELY*rho_h[jres]*effj;
						}
						if(abs(partj0[j].resinfo->code)==abs(JCODE) || abs(partj0[j].resinfo->code)==abs(ICODE)){
							/* contribution from same resonance */
							if(icheck[jres] && jcheck[jres]){
								if(nj>1){
									for(i=0;i<nj;i++){
										if(i!=j){
											if((abs(partj0[i].resinfo->code)==abs(ICODE)
												&& abs(partj0[j].resinfo->code)==abs(JCODE))
													|| (abs(partj0[i].resinfo->code)==abs(JCODE)
											&& abs(partj0[j].resinfo->code)==abs(ICODE))){
												BoostParts(&partj0[i],&partj[i],&partj0[j],&partj[j],0.0,jdecay,jdecay);
												CalcDCA(jdecay,&partj[j],dcaj);
												acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
												if(abs(partj0[j].resinfo->code)==abs(JCODE)){
													efficiency_jbar+=effj;
													norm_jbar+=1.0;
												}
												CalcDCA(jdecay,&partj[i],dcai);
												acceptance->CalcAcceptance(accepti,effi,&partj[i],0.0,dcai);
												if(accepti && acceptj && randy->ran()<effi*effj){
													delyij=partj[i].y-partj[j].y;
															 //printf("from decay, deleta=%g\n",partj[i].eta,partj[j].eta);
													sign=-1.0;
													if(partj0[i].resinfo->code*(ICODE/abs(ICODE))*partj0[j].resinfo->code*(JCODE/abs(JCODE))>0) sign=1.0;
													w=w_doublecount*rho_h[jres];
													ibin=lrint(floor(fabs(delyij/DELY)));
													if(ibin<NBINS) bf_resdecay[ibin]+=w*sign;
													delphiij=fabs(atan2(partj[i].p[2],partj[i].p[1])
														-atan2(partj[j].p[2],partj[j].p[1]));
													if(delphiij>PI) delphiij=2.0*PI-delphiij;
													iphibin=lrint(floor(NPHIBINS*delphiij/PI));
													bf_resdecay_phi[iphibin]+=w*sign;
												}
											} 
										}
									}
								}
							}
							/* contribution from other particles */
							for(ires=0;ires<NRES;ires++){
								if((icheck[ires] && jcheck[jres]) || (icheck[jres]&&jcheck[ires])){
									wB=BL_h[ires][jres]+B_h[ires][jres];
									if(fabs(wB)>1.0E-20){
										codecheck=GetProducts(code[ires],ni,parti0,idecay);
										if(codecheck){
											for(i=0;i<ni;i++){
												if((abs(parti0[i].resinfo->code)==abs(ICODE) && abs(partj0[j].resinfo->code)==abs(JCODE)) || (abs(parti0[i].resinfo->code)==abs(JCODE) && abs(partj0[j].resinfo->code)==abs(ICODE))){
													sign=-1.0;
													if(parti0[i].resinfo->code*(ICODE/abs(ICODE))*partj0[j].resinfo->code*(JCODE/abs(JCODE))>0) sign=1.0;
													BoostParts(&parti0[i],&parti[i],&partj0[j],&partj[j],SIGMA_HAD,idecay,jdecay);
													pti=sqrt(parti[i].p[1]*parti[i].p[1]+parti[i].p[2]*parti[i].p[2]);
													ptj=sqrt(partj[j].p[1]*partj[j].p[1]+partj[j].p[2]*partj[j].p[2]);
													if(abs(parti0[i].resinfo->code)==abs(ICODE)){
														ptbari+=pti*rho_h[ires];
														ptnormi+=rho_h[ires];
														ptbarj+=ptj*rho_h[jres];
														ptnormj+=rho_h[jres];
													}
													else{
														ptbarj+=pti*rho_h[ires];
														ptnormj+=rho_h[ires];
														ptbari+=ptj*rho_h[jres];
														ptnormi+=rho_h[jres];
													}
													CalcDCA(jdecay,&partj[j],dcaj);
													acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
													if(abs(partj0[j].resinfo->code)==abs(JCODE)){
														efficiency_jbar+=effj;
														norm_jbar+=1.0;
													}
													
													CalcDCA(idecay,&parti[i],dcai);
													acceptance->CalcAcceptance(accepti,effi,&parti[i],0.0,dcai);
													if(accepti && acceptj && randy->ran()<effi*effj){
														delyij=parti[i].y-partj[j].y;
														w=w_doublecount*0.5*BL_h[ires][jres]*rho_h[jres];
														ibin=lrint(floor(fabs(delyij/DELY)));
														if(ibin<NBINS) bf_had[ibin]+=w*sign;
														delphiij=fabs(atan2(parti[i].p[2],parti[i].p[1])
															-atan2(partj[j].p[2],partj[j].p[1]));
														if(delphiij>PI) delphiij=2.0*PI-delphiij;
														iphibin=lrint(floor(NPHIBINS*delphiij/PI));
														bf_had_phi[iphibin]+=w*sign;
														blocal+=w*sign*DELY;
													}
													
													if(!SINGLEWAVE){
														BoostParts(&parti0[i],&parti[i],&partj0[j],&partj[j],SIGMA_QGP,idecay,jdecay);
														
														CalcDCA(jdecay,&partj[j],dcaj);
														acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
														if(abs(partj0[j].resinfo->code)==abs(JCODE)){
															efficiency_jbar+=effj;
															norm_jbar+=1.0;
														}
														
														CalcDCA(idecay,&parti[i],dcai);
														acceptance->CalcAcceptance(accepti,effi,&parti[i],0.0,dcai);
														if(abs(parti0[i].resinfo->code)==abs(JCODE)){
															efficiency_jbar+=effi;
															norm_jbar+=1.0;
														}
														if(accepti && acceptj && randy->ran()<effi*effj){
															delyij=parti[i].y-partj[j].y;
															w=w_doublecount*0.5*B_h[ires][jres]*rho_h[jres];
															ibin=lrint(floor(fabs(delyij/DELY)));
															if(ibin<NBINS) bf_qgp[ibin]+=w*sign;
															delphiij=fabs(atan2(parti[i].p[2],parti[i].p[1])
																-atan2(partj[j].p[2],partj[j].p[1]));
															if(delphiij>PI) delphiij=2.0*PI-delphiij;
															iphibin=lrint(floor(NPHIBINS*delphiij/PI));
															if(iphibin>=NPHIBINS || iphibin<0){
																printf("iphibin =%d is out of range\n",iphibin);
																exit(1);
															}
															bf_qgp_phi[iphibin]+=w*sign;									
															bqgp+=w*sign*DELY;
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
		//denom=denom*efficiency_jbar/norm_jbar;
		//printf("----------- denom=%g -------------\n",denom);
	printf("<pt_i>=%g, <pt_j>=%g\n",ptbari/ptnormi,ptbarj/ptnormj);
	printf("bqgp=%g, blocal=%g, local/qgp= %g\n",bqgp/denom,blocal/denom,blocal/bqgp);
	normalization=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=bf_had[ibin]+bf_qgp[ibin]+bf_resdecay[ibin];
		bf[ibin]=bf[ibin]/denom;
		bf_had[ibin]=bf_had[ibin]/denom;
		bf_qgp[ibin]=bf_qgp[ibin]/denom;
		bf_resdecay[ibin]=bf_resdecay[ibin]/denom;
		normalization+=bf[ibin]*DELY;
	}
	printf("check 1: normalization=%g\n",normalization);
	normalization=0.0;
	double DELPHI=PI/double(NPHIBINS);
	denom=denom*DELPHI/DELY;
	for(iphibin=0;iphibin<NPHIBINS;iphibin++){
		bf_phi[iphibin]=bf_had_phi[iphibin]+bf_qgp_phi[iphibin]+bf_resdecay_phi[iphibin];
		bf_phi[iphibin]=bf_phi[iphibin]/denom;
		bf_had_phi[iphibin]=bf_had_phi[iphibin]/denom;
		bf_qgp_phi[iphibin]=bf_qgp_phi[iphibin]/denom;
		bf_resdecay_phi[iphibin]=bf_resdecay_phi[iphibin]/denom;
		normalization+=bf_phi[iphibin]*DELPHI;
	}
	printf("check 2: normalization=%g\n",normalization);
	printf("DELPHI=%g, DELY=%g\n",DELPHI,DELY);
}

void CBalance::CalcBF_ALLCHARGES(){
	if(!ALLCHARGES){
		printf("You shouldn't call CBalance::CalcBF_ALLCHARGES unless CBalance::ALLCHARGES==true\n");
		exit(1);
	}
	CPart parti0[10],partj0[10],parti[10],partj[10];
	int ires,jres,ibin,iphibin,ni,nj,i,j,ievent,iq,jq;
	bool accepti,acceptj,codecheck,idecay,jdecay;
	double eta,etai,etaj,delyij,delphiij,dcai[4],dcaj[4];
	double efficiency_jbar=0.0,efficiency_ibar=0.0,norm_jbar=0.0,norm_ibar=0.0;
	double effi,effj;
	double rhotot=0.0,wtot,sign,denom=0.0;
	for(ires=0;ires<NRES;ires++) rhotot+=rho_h[ires];
	double wB,w,w_doublecount=0.5,deleta,deleta0,bqgp=0.0,blocal=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=bf_had[ibin]=bf_resdecay[ibin]=bf_qgp[ibin]=0.0;
	}
	CResInfo *resinfo;
	double qcheck,qtot,ptbar=0.0,pt,ptnorm=0.0,v2bar=0.0,mt,wpt;
	
	for(ievent=0;ievent<NEVENTS;ievent++){
		for(jres=0;jres<NRES;jres++){
			codecheck=GetProducts(code[jres],nj,partj0,jdecay);
			for(j=0;j<nj;j++){
				if(partj0[j].resinfo->charge==0){
					printf("charge of j shouldn't be zero!\n");
					exit(1);
				}
				BoostPart(&partj0[j],&partj[j],jdecay);
				CalcDCA(jdecay,&partj[j],dcaj);
				acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
				denom+=DELY*rho_h[jres]*effj;
				/* Contribution from 2 particles from same mother decay */
				if(nj>1){
					for(i=0;i<nj;i++){
						if(i!=j){
							BoostParts(&partj0[i],&partj[i],&partj0[j],&partj[j],0.0,jdecay,jdecay);
							CalcDCA(jdecay,&partj[j],dcaj);
							acceptance->CalcAcceptance(acceptj,effj,&partj[j],0.0,dcaj);
							efficiency_jbar+=effj*rho_h[jres];
							norm_jbar+=rho_h[jres];
							if(acceptj){
								CalcDCA(jdecay,&partj[i],dcai);
								acceptance->CalcAcceptance(accepti,effi,&partj[i],0.0,dcai);
									//efficiency_ibar+=effi*rho_h[jres]*effj;
									//norm_ibar+=rho_h[jres]*effj;
								
								if(accepti){
									etaj=partj[j].GetPseudoRapidity();
									etai=partj[i].GetPseudoRapidity();
									delyij=etai-etaj;
									sign=-partj0[i].resinfo->charge*partj0[j].resinfo->charge;
									w=rho_h[jres]*effi*effj;
									ibin=lrint(floor(fabs(delyij/DELY)));
									if(ibin<NBINS){
										bf_resdecay[ibin]+=w*sign;
									}
									delphiij=fabs(atan2(partj[i].p[2],partj[i].p[1])
										-atan2(partj[j].p[2],partj[j].p[1]));
									if(delphiij>PI) delphiij=2.0*PI-delphiij;
									iphibin=lrint(floor(NPHIBINS*delphiij/PI));
									bf_resdecay_phi[iphibin]+=w*sign;
								}
							}
						}
					}
				}
				/* Contribution from other particles */
				qtot=0.0;
				for(ires=0;ires<NRES;ires++){
					wB=BL_h[ires][jres]+B_h[ires][jres];
					if(fabs(wB)>1.0E-20){
						codecheck=GetProducts(code[ires],ni,parti0,idecay);
						if(ni>0){
							for(i=0;i<ni;i++){
								sign=-parti0[i].resinfo->charge*partj0[j].resinfo->charge;
								BoostParts(&parti0[i],&parti[i],&partj0[j],&partj[j],SIGMA_HAD,idecay,jdecay);
								pt=sqrt(parti[i].p[1]*parti[i].p[1]+parti[i].p[2]*parti[i].p[2]);
								mt=sqrt(pt*pt+parti[i].mass*parti[i].mass);
								wpt=rho_h[ires]*asinh((pt/mt)*sinh(0.5)); // for average over fixed eta
								//wpt=rho_h[ires]; // for average over fixed y
								
								v2bar+=wpt*(parti[i].p[1]*parti[i].p[1]-parti[i].p[2]-parti[i].p[2]*parti[i].p[2])/(pt*pt);
								ptbar+=wpt*pt;
								ptnorm+=wpt;
								
								pt=sqrt(partj[j].p[1]*partj[j].p[1]+partj[j].p[2]*partj[j].p[2]);
								mt=sqrt(pt*pt+partj[j].mass*partj[j].mass);
								wpt=rho_h[jres]*asinh((pt/mt)*sinh(0.5)); // for average over fixed eta
								//wpt=rho_h[ires]; // for average over fixed y
								v2bar+=wpt*(partj[j].p[1]*partj[j].p[1]-partj[j].p[2]-partj[j].p[2]*partj[j].p[2])/(pt*pt);
								ptbar+=wpt*pt;
								ptnorm+=wpt;
								
								CalcDCA(jdecay,&partj[j],dcaj);
								acceptance->CalcAcceptance(acceptj,effj,&partj[j],0,dcaj);
								efficiency_jbar+=effj*rho_h[jres];
								norm_jbar+=rho_h[jres];
								if(acceptj){								
									CalcDCA(idecay,&parti[i],dcai);
									acceptance->CalcAcceptance(accepti,effi,&parti[i],0,dcai);
									efficiency_ibar+=effi*w_doublecount*B_h[ires][jres]*rho_h[jres]*effj;
									norm_ibar+=w_doublecount*B_h[ires][jres]*rho_h[jres]*effj;									
									if(accepti){
										etaj=partj[j].GetPseudoRapidity();
										etai=parti[i].GetPseudoRapidity();
										delyij=etai-etaj;
										w=w_doublecount*BL_h[ires][jres]*rho_h[jres]*effi*effj;
										qtot-=w*parti0[i].resinfo->charge/rho_h[jres];
										ibin=lrint(floor(fabs(delyij/DELY)));
										if(ibin<NBINS) bf_had[ibin]+=w*sign;
										delphiij=fabs(atan2(parti[i].p[2],parti[i].p[1])
											-atan2(partj[j].p[2],partj[j].p[1]));
										if(delphiij>PI) delphiij=2.0*PI-delphiij;
										iphibin=lrint(floor(NPHIBINS*delphiij/PI));
										bf_had_phi[iphibin]+=w*sign;									
										blocal+=w*sign*DELY;
									}
								}
								
								if(!SINGLEWAVE){
									BoostParts(&parti0[i],&parti[i],&partj0[j],&partj[j],SIGMA_QGP,idecay,jdecay);
									CalcDCA(jdecay,&partj[j],dcaj);
									acceptance->CalcAcceptance(acceptj,effj,&partj[j],0,dcaj);
									efficiency_jbar+=effj*rho_h[jres];
									norm_jbar+=rho_h[jres];
									if(acceptj){
										CalcDCA(idecay,&parti[i],dcai);
										acceptance->CalcAcceptance(accepti,effi,&parti[i],0,dcai);
										efficiency_ibar+=effi*w_doublecount*B_h[ires][jres]*rho_h[jres]*effj;
										norm_ibar+=w_doublecount*B_h[ires][jres]*rho_h[jres]*effj;
										
										if(accepti){
											etaj=partj[j].GetPseudoRapidity();
											etai=parti[i].GetPseudoRapidity();
											delyij=etai-etaj;
											w=w_doublecount*B_h[ires][jres]*rho_h[jres]*effi*effj;
											qtot-=w*parti0[i].resinfo->charge/rho_h[jres];
											ibin=lrint(floor(fabs(delyij/DELY)));
											if(ibin<NBINS) bf_qgp[ibin]+=w*sign;
											delphiij=fabs(atan2(parti[i].p[2],parti[i].p[1])
												-atan2(partj[j].p[2],partj[j].p[1]));
											if(delphiij>PI) delphiij=2.0*PI-delphiij;
											iphibin=lrint(floor(NPHIBINS*delphiij/PI));
											bf_qgp_phi[iphibin]+=w*sign;
											bqgp+=w*sign*DELY;
										}
									}
								}
							}							
						}
					}
				}
					//printf("qtot=%g=?%g\n",0.5*qtot,(2.0*q[jres][0]-q[jres][1]-q[jres][2])/3.0);
			}
		}
	}
	printf("<effj>=%g, <effi>=%g, <pt>=%g <v2>=%g (within fixed eta window)\n",efficiency_jbar/norm_jbar,efficiency_ibar/norm_ibar,ptbar/ptnorm,v2bar/ptnorm);
		//denom=denom*efficiency_jbar/norm_jbar;
		//denom=denom/double(NRES);
		//printf("----------- denom=%g -------------\n",denom);
	printf("bqgp=%g, blocal=%g, local/qgp= %g\n",bqgp/denom,blocal/denom,blocal/bqgp);
	normalization=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=bf_had[ibin]+bf_qgp[ibin]+bf_resdecay[ibin];
		bf[ibin]=bf[ibin]/denom;
		bf_had[ibin]=bf_had[ibin]/denom;
		bf_qgp[ibin]=bf_qgp[ibin]/denom;
		bf_resdecay[ibin]=bf_resdecay[ibin]/denom;
		normalization+=bf[ibin]*DELY;
	}
	double DELPHI=PI/double(NPHIBINS);
	denom=denom*DELPHI/DELY;
	normalization=0.0;
	for(iphibin=0;iphibin<NPHIBINS;iphibin++){
		bf_phi[iphibin]=bf_had_phi[iphibin]+bf_qgp_phi[iphibin]+bf_resdecay_phi[iphibin];
		bf_phi[iphibin]=bf_phi[iphibin]/denom;
		bf_had_phi[iphibin]=bf_had_phi[iphibin]/denom;
		bf_qgp_phi[iphibin]=bf_qgp_phi[iphibin]/denom;
		bf_resdecay_phi[iphibin]=bf_resdecay_phi[iphibin]/denom;
		normalization+=bf_phi[iphibin]*DELPHI;
	}
}

void CBalance::CalcBF_Gaussian(){
	if(ALLCHARGES){
		printf("You shouldn't call CBalance::CalcBF_Gauss unless CBalance::ALLCHARGES==false\n");
		exit(1);
	}
	int ires,jres,ibin,ni,nj,i,j,ievent,iq,jq;
	bool accepti,acceptj,codecheck;
	double eta,etai,etaj,delyij,dcai[4]={0.0},dcaj[4]={0.0};
	double efficiency_jbar=0.0,norm_jbar=0.0,denom=0.0;
	double effi,effj,w,b=0.0;
	double rhotot=0.0,wtot,sign;
	for(ires=0;ires<NRES;ires++) rhotot+=rho_h[ires];
	double u[4]={1.0,0.0,0.0,0.0};
	CPart parti,partj,parti0,partj0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=0.0;
	}
	
	CResInfo *resinfoi,*resinfoj;
	double qcheck,qtot;
	reslist->GetResInfoptr(ICODE,resinfoi);
	reslist->GetResInfoptr(JCODE,resinfoj);
	parti0.resinfo=resinfoi;
	partj0.resinfo=resinfoj;
	
	for(ievent=0;ievent<NEVENTS;ievent++){
		randy->generate_boltzmann(resinfoi->mass,T,parti0.p);
		randy->generate_boltzmann(resinfoj->mass,T,partj0.p);
		BoostParts(&parti0,&parti,&partj0,&partj,SIGMA_GAUSS,false,false);
		acceptance->CalcAcceptance(acceptj,effj,&partj,0,dcaj);
		acceptance->CalcAcceptance(accepti,effi,&parti,0,dcai);
		denom+=effj*DELY;
		
		if(accepti && acceptj){
			delyij=fabs(parti.y-partj.y);;
			ibin=lrint(floor(fabs(delyij/DELY)));
			if(ibin<NBINS) bf[ibin]+=effi*effj;
		}
		
	}
	normalization=0.0;
	for(ibin=0;ibin<NBINS;ibin++){
		bf[ibin]=bf[ibin]/denom;
		normalization+=bf[ibin]*DELY;
	}
}
