#include "blast.h"
using namespace std;

int main(){
	parameterMap parmap;
	long long int ipair,NPAIRS,NB=0;
	int iphi,iphia,iphib,iphic,iphid;
	const int NPHI=24;
	double BHBT[NPHI]={0.0};
	int pcent;
	bool accepta,acceptb,acceptc,acceptd;
	double phia,phib,phic,phid,dphi;
	double eta,e,g,gv,y;
	double r,q_ac,ctheta_ac,q_ad,ctheta_ad,q_bc,ctheta_bc,q_bd,ctheta_bd;
	double weight;
	double MPI=139.58;
	double rab[4],rcd[4],uab[4],ucd[4];
	double pa[4],pb[4],pc[4],pd[4];
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

	CWaveFunction *cw;
	CWaveFunction_pipluspiplus_nostrong *cw_pp=new CWaveFunction_pipluspiplus_nostrong(string("parameters/wfparameters.dat"));
	CWaveFunction_pipluspiminus_nostrong *cw_pm=new CWaveFunction_pipluspiminus_nostrong(string("parameters/wfparameters.dat"));

//________________________________________________________________________________//
	for(ipair=0;ipair<NPAIRS;ipair++){
		do{
			blast->GetUR(uab,rab);
			blast->GetP(MPI,uab,pa);
			blast->GetP(MPI,uab,pb);
			y=atanh(pa[3]/pa[0]);
			eta=ymax*(1.0-2.0*randy->ran())-y;
			g=cosh(eta);
			gv=sinh(eta);
			rab[0]=g*tau;
			rab[3]=gv*tau;
			e=g*pa[0]+gv*pa[3];
			pa[3]=g*pa[3]+gv*pa[0];
			pa[0]=e;
			e=g*pb[0]+gv*pb[3];
			pb[3]=g*pb[3]+gv*pb[0];
			pb[0]=e;
			accepta=acceptance(pa);
			acceptb=acceptance(pb);
		}while((!accepta)  && (!acceptb));
		phia=atan2(pa[2],pa[1]);
		phib=atan2(pb[2],pb[1]);

		do{
			blast->GetUR(ucd,rcd);
			blast->GetP(MPI,ucd,pc);
			blast->GetP(MPI,ucd,pd);
			y=atanh(pc[3]/pc[0]);
			eta=ymax*(1.0-2.0*randy->ran())-y;
			g=cosh(eta);
			gv=sinh(eta);
			rcd[0]=g*tau;
			rcd[3]=gv*tau;
			e=g*pc[0]+gv*pc[3];
			pc[3]=g*pc[3]+gv*pc[0];
			pc[0]=e;
			e=g*pd[0]+gv*pd[3];
			pd[3]=g*pd[3]+gv*pd[0];
			pd[0]=e;
			acceptc=acceptance(pc);
			acceptd=acceptance(pd);
		}while((!acceptc)  && (!acceptd));
		phic=atan2(pc[2],pc[1]);
		phid=atan2(pd[2],pd[1]);

		if(accepta) NB+=1;
		if(acceptb) NB+=1;
		if(acceptc) NB+=1;
		if(acceptd) NB+=1;

		//weight=GetPhiEff2_bosons(pa,rab,pc,rcd)*GetPhiEff2_bosons(pb,rab,pd,rcd)
			//-GetPhiEff2_bosons(pa,rab,pd,rcd)*GetPhiEff2_bosons(pb,rab,pc,rcd);

		cw->getqrctheta(pa,rab,pc,rcd,&q_ac,&r,&ctheta_ac);
		cw->getqrctheta(pa,rab,pd,rcd,&q_ad,&r,&ctheta_ad);
		cw->getqrctheta(pb,rab,pc,rcd,&q_bc,&r,&ctheta_bc);
		cw->getqrctheta(pb,rab,pd,rcd,&q_bd,&r,&ctheta_bd);
		weight=cw_pp->GetPsiSquared(q_ac,r,ctheta_ac)*cw_pp->GetPsiSquared(q_bd,r,ctheta_bd)
			*cw_pm->GetPsiSquared(q_ad,r,ctheta_ad)*cw_pm->GetPsiSquared(q_bc,r,ctheta_bc);
		weight-=cw_pm->GetPsiSquared(q_ac,r,ctheta_ac)*cw_pm->GetPsiSquared(q_bd,r,ctheta_bd)
			*cw_pp->GetPsiSquared(q_ad,r,ctheta_ad)*cw_pp->GetPsiSquared(q_bc,r,ctheta_bc);

		if(accepta && acceptc){
			dphi=fabs(phia-phic);
			while(dphi>PI){
				dphi-=2.0*PI;
			}
			dphi=fabs(dphi);
			iphi=rint(floor(NPHI*fabs(dphi)/PI));
			if(iphi<NPHI) BHBT[iphi]-=weight;
		}
		if(acceptb && acceptd){
			dphi=fabs(phib-phid);
			while(dphi>PI){
				dphi-=2.0*PI;
			}
			dphi=fabs(dphi);
			iphi=rint(floor(NPHI*fabs(dphi)/PI));
			if(iphi<NPHI) BHBT[iphi]-=weight;
		}
		if(accepta && acceptd){
			dphi=fabs(phia-phid);
			while(dphi>PI){
				dphi-=2.0*PI;
			}
			dphi=fabs(dphi);
			iphi=rint(floor(NPHI*fabs(dphi)/PI));
			if(iphi<NPHI) BHBT[iphi]+=weight;
		}
		if(acceptb && acceptc){
			dphi=fabs(phib-phic);
			while(dphi>PI){
				dphi-=2.0*PI;
			}
			dphi=fabs(dphi);
			iphi=rint(floor(NPHI*fabs(dphi)/PI));
			if(iphi<NPHI) BHBT[iphi]+=weight;
		}	
		if((1+ipair)%(NPAIRS/10)==0) printf("finished %g percent\n",100*(1+ipair)/double(NPAIRS));
	}
	for(iphi=0;iphi<NPHI;iphi++) BHBT[iphi]*=2.0*blast->Mult*(double(NPHI)/PI)/double(NB);

	FILE *fptr=fopen(filename,"a");
	for(iphi=0;iphi<NPHI;iphi++){
		printf("%6.2f %8.4f\n",(0.5+iphi)*180.0/double(NPHI),BHBT[iphi]);
		fprintf(fptr,"%6.2f %8.4f\n",(0.5+iphi)*180.0/double(NPHI),BHBT[iphi]);
	}
	fclose(fptr);
}




