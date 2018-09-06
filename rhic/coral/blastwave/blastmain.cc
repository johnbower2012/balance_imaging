#include "blast.h"
using namespace std;

int main(){
	parameterMap parmap;
	double Cnum[80]={0.0},Cden[80]={0.0},r2[80]={0.0};
	double y2bar=0.0,normy2=0.0;
	double v2_pi[80]={0.0},v2_k[80]={0.0},v2_p[80]={0.0},v2norm_pi[80]={0.0},v2norm_k[80]={0.0},v2norm_p[80]={0.0};
	double v2=0.0,v2s=0.0,cdphibar=0.0,v2c=0.0,v2norm,cdphi,sdphi,s2phia,c2phia,s2phib,c2phib;
	long long int ipair,NPAIRS,ntries=0,Nab=0,NB;
	int iphia,iphib,iphi,iq,ipt;
	const int NPHI=24;
	long double C[NPHI]={0.0},B[NPHI]={0.0};
	double e,gamma,gammav,dely,t,meanpt_pi=0.0,meanpt_k=0.0,meanpt_p=0.0,ptbinsize=25.0;
	double spectranorm_pi,spectranorm_k,spectranorm_p,spectra_pi[80]={0.0},spectra_k[80]={0.0},spectra_p[80]={0.0};
	int pcent;
	bool accepta,acceptb;
	double phia,phib,phic,dphia,dphib,dphi,c2phi;
	double qbc,qac,rac,rbc,ctheta,Pt;
	double psi2_ac_pp,psi2_bc_pp,psi2_ac,psi2_bc_pm,weight,wac_pp,wac_pm,wbc_pp,wbc_pm,sum;
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
	printf("enter percent centrality (0,5,10,20,30,40,50,60,70): ");
	scanf("%d",&pcent);
	blast->SetPars(pcent);
	char filename[80];
	FILE *fptr=fopen("blast.dat","a");

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
			etaab=atanh(rab[3]/rab[0]);
			pb[3]=gamma*pb[3]+gammav*pb[0];
			pb[0]=sqrt(Mb*Mb+pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);
			yb=atanh(pb[3]/pb[0]);
			ntries+=1;		
		} while((!acceptance(pta,ya)) && (!acceptance(ptb,yb)) );
		phia=atan2(pa[2],pa[1]);
		phib=atan2(pb[2],pb[1]);
	
		accepta=acceptance(pta,ya);
		acceptb=acceptance(ptb,yb);
		if(accepta){
			Nab+=1;
			ipt=pta/25.0;
			c2phi=cos(2.0*phia);
			v2+=c2phi;
			if(fabs(Ma-139.57)<1.0 && pta<2000.0){
				meanpt_pi+=pta;
				spectra_pi[ipt]+=1.0;
				spectranorm_pi+=1.0;
				v2_pi[ipt]+=c2phi;
				v2norm_pi[ipt]+=1.0;
			}			
			if(fabs(Ma-493.7)<1.0 && pta<2000.0){
				meanpt_k+=pta;
				spectra_k[ipt]+=1.0;
				spectranorm_k+=1.0;
				v2_k[ipt]+=c2phi;
				v2norm_k[ipt]+=1.0;
			}
			if(fabs(Ma-938.28)<1.0 && pta<2000.0){
				meanpt_p+=pta;
				spectra_p[ipt]+=1.0;
				spectranorm_p+=1.0;
				v2_p[ipt]+=c2phi;
				v2norm_p[ipt]+=1.0;
			}
		}
		if(acceptb){
			Nab+=1;
			c2phi=cos(2.0*phib);
			v2+=c2phi;
			ipt=ptb/25.0;
			if(fabs(Mb-139.58)<1.0 && ptb<2000.0){
				meanpt_pi+=ptb;
				spectra_pi[ipt]+=1.0;
				spectranorm_pi+=1.0;
				v2_pi[ipt]+=c2phi;
				v2norm_pi[ipt]+=1.0;
			}
			if(fabs(Mb-493.7)<1.0 && ptb<2000.0){
				meanpt_k+=ptb;
				spectra_k[ipt]+=1.0;
				spectranorm_k+=1.0;
				v2_k[ipt]+=c2phi;
				v2norm_k[ipt]+=1.0;
			}
			if(fabs(Mb-938.28)<1.0 && ptb<2000.0){
				meanpt_p+=ptb;
				spectra_p[ipt]+=1.0;
				spectranorm_p+=1.0;
				v2_p[ipt]+=c2phi;
				v2norm_p[ipt]+=1.0;
			}
		}

		if(accepta && acceptb){
			dphi=fabs(phia-phib);
			if(dphi>PI) dphi=2.0*PI-dphi;
			iphi=(dphi/PI)*double(NPHI);
			B[iphi]+=1.0;
			cdphi=cos(dphi);
			sdphi=sin(phia-phib);
			c2phia=cos(2.0*phia);
			c2phib=cos(2.0*phib);
			s2phia=sin(2.0*phia);
			s2phib=sin(2.0*phib);
			cdphibar+=cdphi;
			v2c+=cdphi*0.5*(c2phia+c2phib);
			v2s-=sdphi*0.5*(s2phia-s2phib);
			NB+=1;
		}
		if((1+ipair)%(NPAIRS/10)==0) printf("finished %g percent\n",100*(1+ipair)/double(NPAIRS));
	}
	cdphibar=cdphibar/double(NB);
	v2=v2/double(Nab);
	v2c=v2c/double(NB);
	v2c=v2c-cdphibar*v2;
	v2s=v2s/double(NB);

	for(iphi=0;iphi<NPHI;iphi++){
		printf("B(iphi=%2d)=%8.4f\n",iphi,double(B[iphi]/double(NB)));
	}
	//FILE *parsfptr=fopen(parsfilename.c_str(),"r");
	//char linechars[120];
	//while(!feof(parsfptr)){
		//fgets(linechars,120,parsfptr);
		//fprintf(fptr,"# %s",linechars);
	//}
	//fclose(parsfptr);
	printf("success rate=%g\n",double(NPAIRS)/double(ntries));
	printf("v2=%g, cdphibar=%g:  cdphibar*v2=%g, v2c=%g, v2s=%g\n",v2,cdphibar,cdphibar*v2,v2c,v2s);
	fprintf(fptr,"#____________ percent centrality = %d ____________\n",blast->pcent);
	fprintf(fptr,"# NPAIRS=%gx10^6, tau=%g, etaG=%g, ymax=%g\n",double(NPAIRS/1000000),tau,etaG,ymax);
	fprintf(fptr,"# eps=%g, p0=%g, p2=%g, p4=%g, T=%g, Npart=%g, Nch=%g\n",
		blast->eps,blast->p0,blast->p2,blast->p4,blast->T,blast->Npart,blast->Nch);
	fprintf(fptr,"# v2 <cos(phi_B)> v2*<cos(phi_B)>  v2c    v2s\n");
	fprintf(fptr,"%7.4f  %7.4f  %7.4f %7.4f %7.4f\n",v2,cdphibar,cdphibar*v2,v2c,v2s);
	fclose(fptr);


	sprintf(filename,"spectra_%d.dat",pcent);
	fptr=fopen(filename,"w");
	sprintf(filename,"v2_%d.dat",pcent);
	FILE *fptrv2=fopen(filename,"w");
	fprintf(fptr,"#--- mean pt =(%g,%g,%g) ---\n",meanpt_pi/spectranorm_pi,meanpt_k/spectranorm_k,meanpt_p/spectranorm_p);
	printf("----- mean pt =(%g,%g,%g), ----\n",meanpt_pi/spectranorm_pi,meanpt_k/spectranorm_k,meanpt_p/spectranorm_p);
	for(ipt=0;ipt<80;ipt++){
		fprintf(fptr,"%g %g %g %g\n", (ipt+0.5)*25.0,spectra_pi[ipt]/spectranorm_pi,spectra_k[ipt]/spectranorm_k,spectra_p[ipt]/spectranorm_p);
		fprintf(fptrv2,"%g %g %g %g\n", (ipt+0.5)*25.0,v2_pi[ipt]/v2norm_pi[ipt],v2_k[ipt]/v2norm_k[ipt],v2_p[ipt]/v2norm_p[ipt]);
	}
	fclose(fptr);
	fclose(fptrv2);
}




