#include "blast.h"
using namespace std;

int main(){
	parameterMap parmap;
	double Cnum[80]={0.0},Cden[80]={0.0},r2[80]={0.0};
	double y2bar=0.0;
	double v2=0.0,v2s=0.0,cphibar=0.0,v2c=0.0,cdphi,sdphi;
	long long int ipair,NPAIRS,ntries=0,Nabc=0,Nac=0,Nbc=0;
	int iphi,iphia,iphib,iq,ipt;
	const int NPHI=24;
	double C[NPHI]={0.0},BHBT[NPHI]={0.0};
	double e,gamma,gammav,dely,t,meanpt_pi=0.0,meanpt_k=0.0,meanpt_p=0.0,ptbinsize=25.0;
	int pcent;
	bool accepta,acceptb,acceptc;
	double phia,phib,phic,dphia,dphib,dphi,c2phi;
	double qbc,qac,rac,rbc,ctheta,Pt;
	double psi2_ac_pp,psi2_bc_pp,psi2_ac,psi2_bc_pm,weight,wac_pp,wac_pm,wbc_pp,wbc_pm,sum;
	double c2phia,c2phib,s2phia,s2phib,c2phic,s2phic;
	double Ma=139.58,Mb=139.58,Mc=139.58;
	double rab[4],rc[4],uab[4],uc[4];
	double pa[4],pb[4],pc[4],ya,yb,yc,pta,ptb,ptc,etaab,etac;
	string parsfilename="parameters/blast.dat";
	parameter::ReadParsFromFile(parmap,parsfilename);
	double RHBT=parameter::getD(parmap,"RHBT");
	double etaG=parameter::getD(parmap,"etaG");
	double tau=parameter::getD(parmap,"tau");
	double ymax=parameter::getD(parmap,"ymax");
	NPAIRS=100*llrint(10000*parameter::getD(parmap,"NMILLIONPAIRS"));
	CBlast *blast=new CBlast(tau,RHBT,etaG);
	CRandom *randy=blast->randy;
	printf("Enter percent centrality (0,5,10,20,30,40,50,60,70): ");
	scanf("%d",&pcent);
	blast->SetPars(pcent);
	char filename[80];
	sprintf(filename,"Bhbt_%d.dat",blast->pcent);
	FILE *fptr=fopen(filename,"w");

	CWaveFunction *cw;
	//CWaveFunction_pipluspiplus_sqwell *cw_pp=new CWaveFunction_pipluspiplus_sqwell(string("parameters/wfparameters.dat"));
	//CWaveFunction_pipluspiminus_sqwell *cw_pm=new CWaveFunction_pipluspiminus_sqwell(string("parameters/wfparameters.dat"));

//________________________________________________________________________________//
	for(ipair=0;ipair<NPAIRS;ipair++){
		do{
			Ma=GetM(randy);
			Mb=GetM(randy);
			blast->GetUR(uab,rab);
			blast->GetP(Ma,uab,pa);
			ya=atanh(pa[3]/pa[0]);
			pta=sqrt(pa[1]*pa[1]+pa[2]*pa[2]);
			blast->GetP(Mb,uab,pb);
			yb=atanh(pb[3]/pb[0]);
			ptb=sqrt(pb[1]*pb[1]+pb[2]*pb[2]);
			dely=ymax*(1.0-2.0*randy->ran())-0.5*(ya+yb);
			gamma=cosh(dely);
			gammav=sinh(dely);
			pa[3]=gamma*pa[3]+gammav*pa[0];
			pa[0]=sqrt(Ma*Ma+pa[1]*pa[1]+pa[2]*pa[2]+pa[3]*pa[3]);
			ya=atanh(pa[3]/pa[0]);
			t=rab[0];
			rab[0]=gamma*t+gammav*rab[3];
			rab[3]=gamma*rab[3]+gammav*t;
			pb[3]=gamma*pb[3]+gammav*pb[0];
			pb[0]=sqrt(Mb*Mb+pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);
			yb=atanh(pb[3]/pb[0]);
			ntries+=1;
			accepta=acceptance(pta,ya);
			acceptb=acceptance(ptb,yb);
		}while((!accepta)  && (!acceptb));
		phia=atan2(pa[2],pa[1]);
		phib=atan2(pb[2],pb[1]);
		
		do{
			Mc=GetM(randy);
			blast->GetUR(uc,rc);
			blast->GetP(Mc,uc,pc);
			ptc=sqrt(pc[1]*pc[1]+pc[2]*pc[2]);
			yc=atanh(pc[3]/pc[0]);
			dely=ymax*(1.0-2.0*randy->ran())-yc;
			gamma=cosh(dely);
			gammav=sinh(dely);
			e=pc[0];
			pc[3]=gamma*pc[3]+gammav*pc[0];
			pc[0]=sqrt(Mc*Mc+pc[1]*pc[1]+pc[2]*pc[2]+pc[3]*pc[3]);
			yc=yc+dely;
			t=rc[0];
			rc[0]=gamma*t+gammav*rc[3];
			rc[3]=gamma*rc[3]+gammav*t;
			acceptc=acceptance(ptc,yc);
		} while(!acceptc);
		phic=atan2(pc[2],pc[1]);
		dphia=phia-phic;
		if(fabs(dphia)>PI) dphia=dphia-PI*dphia/fabs(dphia);
		dphib=phib-phic;
		if(fabs(dphib)>PI) dphib=dphib-PI*dphib/fabs(dphib);
		iphia=rint(floor(NPHI*fabs(dphia)/PI));
		iphib=rint(floor(NPHI*fabs(dphib)/PI));
		if(iphia>NPHI || iphib>NPHI) printf("DISASTER, iphia=%d, iphib=%d, NPHI=%d",iphia,iphib,NPHI);
		
		/*
		wac_pm=cw_pm->GetPsiSquared(pa,rab,pc,rc);
		wac_pp=cw_pp->GetPsiSquared(pa,rab,pc,rc);
		wbc_pm=cw_pm->GetPsiSquared(pb,rab,pc,rc);
		wbc_pp=cw_pp->GetPsiSquared(pb,rab,pc,rc);
		*/
		/*
		wac_pm=GetPhiEff2(pa,Ma,rab,pc,Mc,rc,-1,blast);
		wbc_pm=GetPhiEff2(pb,Mb,rab,pc,Mc,rc,-1,blast);
		wac_pp=GetPhiEff2(pa,Ma,rab,pc,Mc,rc,1,blast);
		wbc_pp=GetPhiEff2(pb,Mb,rab,pc,Mc,rc,1,blast);
		*/
		wac_pm=wbc_pm=0.0;
		wac_pp=GetPhiEff2_bosons(pa,Ma,rab,pc,Mc,rc);
		wbc_pp=GetPhiEff2_bosons(pb,Mb,rab,pc,Mc,rc);
				
		weight=wac_pm+wbc_pp-wac_pp-wbc_pm;

		cw->getqrctheta(pa,rab,pc,rc,&qac,&rac,&ctheta);
		iq=lrint(floor(qac/2.0));
		if(iq<40 && fabs(Ma-139.58)<1.0 && fabs(Mc-139.58)<1.0){
			Cnum[iq]+=wac_pp;
			Cden[iq]+=1.0;
			r2[iq]+=rac*rac;
		}
		cw->getqrctheta(pb,rab,pc,rc,&qbc,&rbc,&ctheta);
		iq=lrint(floor(qbc/2.0));
		if(iq<40 && fabs(Mb-139.58)<1.0 && fabs(Mc-139.58)<1.0){
			Cnum[iq]+=wbc_pp;
			Cden[iq]+=1.0;
			r2[iq]+=rbc*rbc;
		}
		c2phic=cos(2.0*phic);
		s2phic=sin(2.0*phic);
		if(accepta){
			C[iphia]+=weight;
			BHBT[iphia]+=weight;
			cdphi=cos(dphia);
			sdphi=sin(dphia);
			c2phia=cos(2.0*phia);
			s2phia=sin(2.0*phia);
			cphibar+=weight*cdphi;
			v2c+=weight*0.5*(c2phic+c2phia)*cos(dphia);
			v2s+=weight*0.5*(s2phic-s2phia)*sin(dphia);
			Nac+=1;
		}
		if(acceptb){
			C[iphib]-=weight;
			BHBT[iphib]-=weight;
			cdphi=cos(dphib);
			sdphi=sin(dphib);
			c2phib=cos(2.0*phib);
			c2phic=cos(2.0*phic);
			s2phib=sin(2.0*phib);
			s2phic=sin(2.0*phic);
			cphibar-=weight*cdphi;
			v2c-=weight*0.5*(c2phic+c2phib)*cdphi;
			v2s-=weight*0.5*(s2phic-s2phib)*sdphi;
			Nbc+=1;
		}
		if((1+ipair)%(NPAIRS/10)==0) printf("finished %g percent\n",100*(1+ipair)/double(NPAIRS));
	}
	for(iphi=0;iphi<NPHI;iphi++)BHBT[iphi]*=0.5*blast->Mult*(double(NPHI)/PI)/double(Nac+Nbc);
	for(iphi=0;iphi<NPHI;iphi++){
		printf("%6.2f %8.4f\n",(0.5+iphi)*180.0/double(NPHI),BHBT[iphi]);
		fprintf(fptr,"%6.2f %8.4f\n",(0.5+iphi)*180.0/double(NPHI),BHBT[iphi]);
	}
	cphibar*=0.5/double(Nac+Nbc);
	v2c*=0.5/double(Nac+Nbc);
	v2s*=0.5/double(Nac+Nbc);
	v2=v2/double(Nabc);
	v2c=v2c-v2*cphibar;
	printf("v2=%g, <cphiB>=%g, v2c=%g, v2s=%g\n",v2,cphibar,v2c,v2s);
	fprintf(fptr,"v2=%g, <cphiB>=%g, v2c=%g, v2s=%g\n",v2,cphibar,v2c,v2s);

	FILE *parsfptr=fopen(parsfilename.c_str(),"r");
	char linechars[120];
	while(!feof(parsfptr)){
		fgets(linechars,120,parsfptr);
		fprintf(fptr,"# %s",linechars);
	}
	fclose(parsfptr);
	fclose(fptr);	
	printf("success rate=%g\n",double(NPAIRS)/double(ntries));
	for(int iq=0;iq<40;iq++){
		printf("C(q=%g)=%g, Cden=%g, Rrms=%g\n",iq*2.0,Cnum[iq]/Cden[iq],Cden[iq],sqrt(0.5*r2[iq]/(3.0*Cden[iq])));
	}
}




