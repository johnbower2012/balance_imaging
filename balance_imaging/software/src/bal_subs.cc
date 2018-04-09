#include "bal.h"
using namespace std;

CBalance::CBalance(string parsdirname){
	parameter::ReadParsFromFile(parmap,"model_output/"+parsdirname+"/fixed_parameters.dat");
	parameter::ReadParsFromFile(parmap,"model_output/"+parsdirname+"/parameters.dat");
	parameter::PrintPars(parmap);

	GAUSSIAN=parameter::getB(parmap,"GAUSSIAN",false); // simply calc. BF normalized to unity for two species
	OFFDIAG=parameter::getB(parmap,"OFFDIAG",false);
	T=parameter::getD(parmap,"T",100.0); //Blast Wave Temperature
	T=parameter::getD(parmap,"BW_T",T); //Blast Wave Temperature
	LAMBDA_VISC=parameter::getD(parmap,"BW_LAMBDA_VISC",1.0);
	TCHEM=parameter::getD(parmap,"TCHEM",155.0); // Temp. at which chem. ratios are set
	UPERPX=parameter::getD(parmap,"UPERPX",0.65); //Blast Wave umax at surface
	UPERPY=parameter::getD(parmap,"UPERPY",0.65); //Blast Wave umax at surface
	UPERPX=parameter::getD(parmap,"BW_UPERP",UPERPX); //Blast Wave umax at surface
	UPERPY=parameter::getD(parmap,"BW_UPERP",UPERPY); //Blast Wave umax at surface
	SIGMA_QGP=parameter::getD(parmap,"SIGMA_QGP",1.0); //Width of correlation for long-range and local parts
	SIGMA_HAD_RATIO=parameter::getD(parmap,"SIGMA_HAD_RATIO",0.25);
	SIGMAPHI=parameter::getD(parmap,"SIGMAPHI",0.0);
	SIGMA_HAD=SIGMA_QGP*SIGMA_HAD_RATIO;
	BARYON_FUDGE=parameter::getD(parmap,"BARYON_FUDGE",0.7); //Anti-baryons in FS rel. to thermal (Tc=165).
	NSAMPLE=parameter::getI(parmap,"NSAMPLE",1000000); // No. of MC samples for new version
	DELY=parameter::getD(parmap,"DELY",0.1); //bin width for BF
	USE_SAME_ACCEPTANCE=false;
	STAR_ACCEPTANCE_I=parameter::getB(parmap,"STAR_ACCEPTANCE_I",false);
	STAR_ACCEPTANCE_J=parameter::getB(parmap,"STAR_ACCEPTANCE_J",false);
	if(STAR_ACCEPTANCE_I==STAR_ACCEPTANCE_J) USE_SAME_ACCEPTANCE=true;
	SIGMA_GAUSS=parameter::getD(parmap,"SIGMA_GAUSS",0.3);
	CALCRPLANEBF=parameter::getB(parmap,"CALCRPLANEBF",true);
	randy=new CRandom(-1234);
	gslmatrix=new CGSLMatrix_Real(3);
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/resinfo/resonances_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/resinfo/decays_pdg_weak.dat"));
	reslist=new CResList(&parmap);
	b3d=new CB3D(); // this object is used to decay resonances
	b3d->randy=randy;
	bf=NULL;
	PARSDIRNAME=parsdirname;
	if(STAR_ACCEPTANCE_I){
		acceptance_i=new CAcceptance_STAR();
	}
	else{
		acceptance_i=new CAcceptance_CHEAP();
	}
	if(USE_SAME_ACCEPTANCE){
		acceptance_j=acceptance_i;
		acceptance=acceptance_i;
	}
	else{
		if(STAR_ACCEPTANCE_J){
			acceptance_j=new CAcceptance_STAR();
		}
		else{
			acceptance_j=new CAcceptance_CHEAP();
		}
	}

	rapdist=new CDistInvCosh("model_output/"+parsdirname+"/parameters.dat");
}

void CBalance::InitArrays(){
	InitWeightArrays();
	InitBFArrays();
	printf("arrays initialized\n");
	SetRho();
}

void CBalance::InitBFArrays(){
	int iphibin,ibin,ibal;
	if(bf!=NULL){
		for(ibal=0;ibal<NBAL;ibal++){
			delete [] bf[ibal];
			delete [] bf_had[ibal];
			delete [] bf_qgp[ibal];
			delete [] bf_resdecay[ibal];
			delete [] bf_phi[ibal];
			delete [] bf_had_phi[ibal];
			delete [] bf_qgp_phi[ibal];
			delete [] bf_resdecay_phi[ibal];
		}
		delete [] bf;
		delete [] bf_had;
		delete [] bf_qgp;
		delete [] bf_resdecay;
		delete [] bf_phi;
		delete [] bf_had_phi;
		delete [] bf_qgp_phi;
		delete [] bf_resdecay_phi;
	}
	NBAL=11;
	rhoi=new double[NBAL];
	rhoj=new double[NBAL];
	normalization=new double[NBAL];
	v2=new double[NBAL];
	v2c=new double[NBAL];
	v2s=new double[NBAL];
	cosdphi=new double[NBAL];
	gammap=new double[NBAL];
	for(ibal=0;ibal<NBAL;ibal++){
		rhoj[ibal]=normalization[ibal]=v2[ibal]=v2s[ibal]=v2c[ibal]=cosdphi[ibal]=gammap[ibal]=0.0;
	}
	ALLCHARGES=new bool[NBAL];
	for(ibal=0;ibal<NBAL;ibal++)
		ALLCHARGES[ibal]=false;
	ALLCHARGES[0]=true;
	ICODE=new int[NBAL];
	JCODE=new int[NBAL];
	ICODE[1]=211; JCODE[1]=-211;
	ICODE[2]=211; JCODE[2]=-321;
	ICODE[3]=211; JCODE[3]=-2112;
	ICODE[4]=211; JCODE[4]=-2212;
	ICODE[5]=321; JCODE[5]=-321;
	ICODE[6]=321; JCODE[6]=-2112;
	ICODE[7]=321; JCODE[7]=-2212;
	ICODE[8]=2112; JCODE[8]=-2112;
	ICODE[9]=2112; JCODE[9]=-2212;
	ICODE[10]=2212; JCODE[10]=-2212;

	NBINS=lrint((acceptance_i->ETAMAX+acceptance_j->ETAMAX)/DELY);
	denom=new double[NBAL];
	bf=new double*[NBAL];
	bf_had=new double*[NBAL];
	bf_qgp=new double*[NBAL];
	bf_resdecay=new double*[NBAL];
	for(ibal=0;ibal<NBAL;ibal++){
		denom[ibal]=0.0;
		bf[ibal]=new double[NBINS];
		bf_had[ibal]=new double[NBINS];
		bf_qgp[ibal]=new double[NBINS];
		bf_resdecay[ibal]=new double[NBINS];
		for(ibin=0;ibin<NBINS;ibin++){
			bf[ibal][ibin]=bf_had[ibal][ibin]=bf_qgp[ibal][ibin]=bf_resdecay[ibal][ibin]=0.0;
		}
	}
	NPHIBINS=36;
	denom_phi=new double[NBAL];
	bf_phi=new double*[NBAL];
	bf_had_phi=new double*[NBAL];
	bf_qgp_phi=new double*[NBAL];
	bf_resdecay_phi=new double *[NBAL];
	for(ibal=0;ibal<NBAL;ibal++){
		denom_phi[ibal]=0.0;
		bf_phi[ibal]=new double[NPHIBINS];
		bf_had_phi[ibal]=new double[NPHIBINS];
		bf_qgp_phi[ibal]=new double[NPHIBINS];
		bf_resdecay_phi[ibal]=new double[NPHIBINS];
		for(iphibin=0;iphibin<NPHIBINS;iphibin++){
			bf_phi[ibal][iphibin]=bf_qgp_phi[ibal][iphibin]=bf_had_phi[ibal][iphibin]=bf_resdecay_phi[ibal][iphibin]=0.0;
		}
	}
	if(CALCRPLANEBF){
		NPHIBINS_RPLANE=9;
		denom_rplane=new double *[NBAL];
		bf_rplane=new double **[NBAL];
		bf_had_rplane=new double **[NBAL];
		bf_qgp_rplane=new double **[NBAL];
		bf_resdecay_rplane=new double **[NBAL];
		for(ibal=0;ibal<NBAL;ibal++){
			denom_rplane[ibal]=new double[NPHIBINS_RPLANE];
			bf_rplane[ibal]=new double*[NPHIBINS_RPLANE];
			bf_had_rplane[ibal]=new double*[NPHIBINS_RPLANE];
			bf_qgp_rplane[ibal]=new double*[NPHIBINS_RPLANE];
			bf_resdecay_rplane[ibal]=new double*[NPHIBINS_RPLANE];
			for(int iplane=0;iplane<NPHIBINS_RPLANE;iplane++){
				denom_rplane[ibal][iplane]=0.0;
				bf_rplane[ibal][iplane]=new double[2*NPHIBINS];
				bf_qgp_rplane[ibal][iplane]=new double[2*NPHIBINS];
				bf_had_rplane[ibal][iplane]=new double[2*NPHIBINS];
				bf_resdecay_rplane[ibal][iplane]=new double[2*NPHIBINS];
				for(iphibin=0;iphibin<NPHIBINS;iphibin++){
					bf_rplane[ibal][iplane][iphibin]=bf_qgp_rplane[ibal][iplane][iphibin]=bf_had_rplane[ibal][iplane][iphibin]=bf_resdecay_rplane[ibal][iplane][iphibin]=0.0;
				}			
			}
		}
	}
}

void CBalance::InitWeightArrays(){
	int ires;
	CResInfo *resinfo;
	NRES=reslist->resmap.size();
	int a,b,c,d,ab,jres;
	double sign;
	double Q,B,S;
	code=new int[NRES];
	rho_h=new double[NRES];
	rhotot=new double[NRES];
	rho_q=new double[NRES];
	q=new double* [NRES];
	
	for(ires=0;ires<NRES;ires++)
		q[ires]=new double[3];

	CResInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		ires=resinfo->ires;
		code[ires]=resinfo->code;
		Q=resinfo->charge; B=resinfo->baryon; S=resinfo->strange;
		q[ires][0]=Q+B;
		q[ires][1]=2*B-Q+S;
		q[ires][2]=-S;
	}
	
	chi=new double* [3];
	for(a=0;a<3;a++) chi[a]=new double[3];
	chi_inv=new double *[3];
	for(a=0;a<3;a++)
		chi_inv[a]=new double[3];
	wB=new double**[NRES];
	for(ires=0;ires<NRES;ires++){
		wB[ires]=new double *[NRES];
		for(jres=0;jres<NRES;jres++)
			wB[ires][jres]=new double [4];
	}
	g_q=new double* [3];
	for(a=0;a<3;a++)
		g_q[a]=new double[3];
	mu=new double* [3];
	for(a=0;a<3;a++)
		mu[a]=new double[3];
}

void CBalance::SetRho(){
	double mass,Pi,epsiloni,densi,sigma2i,dedti,degen;
	double sdens=0.0,STOT=7000; // STOT = total entropy in one unit of spatial rapidity
	int ires,a,b,ab;
	double baryon,totalhadrons=0.0;
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	totalhadrons=0.0;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		ires=resinfo->ires;
		mass=resinfo->mass;
		degen=2.0*resinfo->spin+1.0;
		if(resinfo->code!=22 && resinfo->code!=333){
			reslist->freegascalc_onespecies(mass,TCHEM,Pi,epsiloni,densi,sigma2i,dedti);
			sdens+=(Pi+epsiloni)*degen/TCHEM;
			rho_h[ires]=densi*degen;
			totalhadrons+=rho_h[ires];
		}
		else
			rho_h[ires]=0.0;
	}
	totalhadrons=0.0;
	for(ires=0;ires<NRES;ires++){
		rho_h[ires]*=STOT/sdens;
		baryon=q[ires][0]+q[ires][1]+q[ires][2];
		if(fabs(baryon)>0.0001)
			rho_h[ires]*=BARYON_FUDGE;
		totalhadrons+=rho_h[ires];
	}
	rhotot[0]=rho_h[0];
	for(ires=1;ires<NRES;ires++){
		rhotot[ires]=rhotot[ires-1]+rho_h[ires];
	}
	for(a=0;a<3;a++){
		for(b=0;b<3;b++)
			g_q[a][b]=0.0;
		g_q[a][a]=rho_q[a];
	}
}

void CBalance::CalcWeights(){
	int a,b,c,d,ab,ires,jres;
	SetRho();
	totalhadrons=0.0;
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			chi[a][b]=0.0;
			for(ires=0;ires<NRES;ires++){
				chi[a][b]+=rho_h[ires]*q[ires][a]*q[ires][b];	
			}
		}
	}
	gslmatrix->Invert(chi,chi_inv);
	
	for(a=0;a<3;a++){
		for(b=0;b<3;b++){
			mu[a][b]=0.0;
			for(c=0;c<3;c++){
				for(d=0;d<3;d++){
					mu[a][b]+=chi_inv[a][c]*g_q[c][d]*chi_inv[d][b];
				}
			}
		}
	}
	/*
	printf("-------  chi   --------\n");
	gslmatrix->Print(chi);
	printf("\n------- chi_inv --------\n");
	gslmatrix->Print(chi_inv);
	*/
	
	/*
	double **A=new double* [3];
	for(a=0;a<3;a++) A[a]=new double[3];
	for(a=0;a<3;a++){
	for(b=0;b<3;b++){
	A[a][b]=0.0;
	for(c=0;c<3;c++){
	A[a][b]+=chi[a][c]*chi_inv[c][b];
	}
	}
	}
	printf("\n------ Unity? ------\n");
	gslmatrix->Print(A);
	 
	for(a=0;a<3;a++) delete [] A[a];
	delete [] A;
	*/

	/*
	printf("\n-------  Mu  --------\n");
	gslmatrix->Print(mu);
	*/
	
	for(jres=0;jres<NRES;jres++){
		totalhadrons+=rho_h[jres];
		for(ires=0;ires<NRES;ires++){
			for(ab=0;ab<4;ab++){
				wB[ires][jres][ab]=0.0;
			}
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					ab=GetAB(a,b);
					wB[ires][jres][ab]-=2.0*q[ires][a]*q[jres][b]*chi_inv[a][b];
				}
			}
		}
	}
	printf("Weights Calculated\n");
}

int CBalance::GetAB(int a,int b){
	int ab;
	if((a==0&&b==0) ||(a==1&&b==1))
		ab=0;
	if((a==0&&b==1) || (a==1&&b==0))
		ab=1;
	if((a==0&&b==2) || (a==2&&b==0) || (a==1&&b==2) || (a==2&&b==1))
		ab=2;
	if(a==2 && b==2)
		ab=3;
	return ab;
}

bool CBalance::GetProducts(int code,int &ni,CPart *parti,bool &decay){
	int motherbaryon=0,daughterbaryons=0;
	CResInfo *resinfoi,*resinfo0;
	resinfo0=reslist->GetResInfoPtr(code);
	motherbaryon=resinfo0->baryon;
	//resinfo0->decay=false; //this turns off decays
	decay=resinfo0->decay;

	double mtot,ucm[4],pprime[4],pt,deltau,mmass;
	bool match;
	int i;
	CPart *mptr;
	array< CPart *, 5> mother;
	array<CPart *, 5> daughter;
	array<CResInfo *, 5> daughterresinfo;
	int alpha,nmothers,ndaughters,idaughter,imother,ntry;
	for(idaughter=0;idaughter<5;idaughter++)
		daughter[idaughter]=new CPart();
	for(imother=0;imother<5;imother++)
		mother[imother]=new CPart();
	ni=0;
	/** Decay the i-particles */
	nmothers=0;
	mmass=resinfo0->mass;
	if(decay){
		mother[0]->resinfo=resinfo0;
		mother[0]->msquared=mmass*mmass;
		randy->generate_boltzmann(mmass,T,mother[0]->p);
		mother[0]->p[3]=mother[0]->p[3]*LAMBDA_VISC;
		mother[0]->p[1]=mother[0]->p[1]*(1.0+0.5*(1.0-LAMBDA_VISC));
		mother[0]->p[2]=mother[0]->p[2]*(1.0+0.5*(1.0-LAMBDA_VISC));
		mother[0]->p[0]=sqrt(mmass*mmass
			+mother[0]->p[1]*mother[0]->p[1]+mother[0]->p[2]*mother[0]->p[2]+mother[0]->p[3]*mother[0]->p[3]);
		mother[0]->y=atanh(mother[0]->p[3]/mother[0]->p[0]);
		for(alpha=1;alpha<4;alpha++)
			mother[0]->r[alpha]=0.0;
		mother[0]->eta=0.0;
		mother[0]->r[0]=mother[0]->tau0=1.0;
		nmothers=1;
	}
	else{
		parti[0].resinfo=resinfo0;
		parti[0].msquared=mmass*mmass;
		randy->generate_boltzmann(mmass,T,parti[0].p);
		parti[0].p[3]=parti[0].p[3]*LAMBDA_VISC;
		parti[0].p[1]=parti[0].p[1]*(1.0+0.5*(1.0-LAMBDA_VISC));
		parti[0].p[2]=parti[0].p[2]*(1.0+0.5*(1.0-LAMBDA_VISC));
		parti[0].p[0]=sqrt(mmass*mmass
					+parti[0].p[1]*parti[0].p[1]+parti[0].p[2]*parti[0].p[2]+parti[0].p[3]*parti[0].p[3]);
		parti[0].y=atanh(parti[0].p[3]/parti[0].p[0]);
		for(alpha=0;alpha<4;alpha++)
			parti[0].r[alpha]=0.0;
		parti[0].eta=0.0;
		parti[0].tau0=parti[0].r[0]=1.0;
		daughterbaryons+=parti[0].resinfo->baryon;
		ni=1;
	}
	imother=0;
	while(imother<nmothers){
		mptr=mother[imother];
		mmass=sqrt(mptr->msquared);
		//if(resinfo0->decay && resinfo0->width<1.0E-5) decay=true;
		deltau=-log(randy->ran())*(HBARC/mptr->resinfo->width);
		for(alpha=0;alpha<4;alpha++)
			mptr->r[alpha]+=deltau*mptr->p[alpha]/mmass;
		mptr->tau0=sqrt(fabs(mptr->r[0]*mptr->r[0]-mptr->r[1]*mptr->r[1]-mptr->r[2]*mptr->r[2]-mptr->r[3]*mptr->r[3]));
		mptr->eta=atanh(mptr->r[3]/mptr->r[0]);
		ntry=0;
		do{
			mtot=0.0;
			if(ntry<25)
				mptr->resinfo->DecayGetResInfoPtr(ndaughters,daughterresinfo);
			else
				mptr->resinfo->DecayGetResInfoPtr_minmass(ndaughters,daughterresinfo);
			for(idaughter=0;idaughter<ndaughters;idaughter++){
				mtot+=daughterresinfo[idaughter]->mass;
				daughter[idaughter]->resinfo=daughterresinfo[idaughter];
			}
			if(ntry>25){
				printf("FATAL: action_perform_decay, ntry too big, mothermass=%g\n",mmass);
				exit(1);
			}
			ntry++;
		}while(mtot>mmass);
		b3d->Decay(mptr,ndaughters,daughter);
		for(idaughter=0;idaughter<ndaughters;idaughter++){
			if(daughter[idaughter]->resinfo->decay){
				mother[nmothers]->Copy(daughter[idaughter]);
				nmothers+=1;
			}
			else{
				
				if(abs(daughter[idaughter]->resinfo->code)!=2212 && abs(daughter[idaughter]->resinfo->code)!=2112 && abs(daughter[idaughter]->resinfo->code)!=211  && abs(daughter[idaughter]->resinfo->code)!=321 && abs(daughter[idaughter]->resinfo->code)!=311 && abs(daughter[idaughter]->resinfo->code)!=111 && abs(daughter[idaughter]->resinfo->code)!=22){
				printf("bizarre decay product, = %d\n", parti[ni].resinfo->code);
				exit(1);
				}
				
				parti[ni].Copy(daughter[idaughter]);
				daughterbaryons+=parti[ni].resinfo->baryon;
				ni+=1;
			}
		}
		imother+=1;
	}

	match=false;
	for(idaughter=0;idaughter<ni;idaughter++){
		if(parti[idaughter].resinfo->code!=22)
			match=true;
	}
	
	for(imother=0;imother<5;imother++)
		delete mother[imother];
	for(idaughter=0;idaughter<5;idaughter++)
		delete daughter[idaughter];
	return match;
}

void CBalance::WriteBF(){
	int ibin,iphi,iplane,ibal;
	string dirname,filename,command;
	FILE *fptr;
	char ijname[60];
	dirname="model_output/"+PARSDIRNAME;
	command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"/meanpt.dat";
	fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"meanpt_pion %g\n",meanpt_pion);
	fprintf(fptr,"meanpt_kaon %g\n",meanpt_kaon);
	fprintf(fptr,"meanpt_proton %g\n",meanpt_proton);
	fclose(fptr);
	for(ibal=0;ibal<NBAL;ibal++){
		if(ALLCHARGES[ibal] || (abs(ICODE[ibal])!=2112 && abs(JCODE[ibal])!=2112)){
			if(!ALLCHARGES[ibal])
				sprintf(ijname,"I%d_J%d.dat",abs(ICODE[ibal]),abs(JCODE[ibal]));
			else
				sprintf(ijname,"allcharges.dat");
			dirname="model_output/"+PARSDIRNAME;
			command="mkdir -p "+dirname;
			system(command.c_str());
			filename=dirname+"/"+string(ijname);
			fptr=fopen(filename.c_str(),"w");
			fprintf(fptr,"# ICODE=%d, JCODE=%d\n",ICODE[ibal],JCODE[ibal]);
			fprintf(fptr,"# Normalized GAUSSIAN BF, SIGMA_GAUSS=%g\n",SIGMA_GAUSS);
			fprintf(fptr,"# NSAMPLE=%d, Tblast=%g, UPERPblast(x,y)=(%g,%g)\n",NSAMPLE,T,UPERPX,UPERPY);
			fprintf(fptr,"# BARYON_FUDGE=%g\n",BARYON_FUDGE);
			fprintf(fptr,"# ------- normalization=%g -----------------\n",normalization[ibal]);
			fprintf(fptr,"#dely    BF\n");
			for(ibin=0;ibin<NBINS;ibin++){
				fprintf(fptr,"%5.3f %7.4f %7.4f\n",DELY*(0.5+ibin),bf[ibal][ibin],dB[ibal][ibin]);
			}
			fclose(fptr);
		}
	}
}

void CBalance::PrintBF(){
	int ibin,iphibin,ibal;
	for(ibal=0;ibal<NBAL;ibal++){
		if(ALLCHARGES[ibal] || (abs(ICODE[ibal])!=2112 && abs(JCODE[ibal])!=2112)){
			printf("# ICODE=%d, JCODE=%d\n",ICODE[ibal],JCODE[ibal]);
			printf("# Normalized GAUSSIAN BF, SIGMA_GAUSS=%g\n",SIGMA_GAUSS);
			printf("# NSAMPLE=%d, Tblast=%g, UPERPblast(x,y)=(%g,%g)\n",NSAMPLE,T,UPERPX,UPERPY);
			printf("# BARYON_FUDGE=%g\n",BARYON_FUDGE);
			printf("# ------- normalization=%g -----------------\n",normalization[ibal]);
			printf("#dely    BF\n");
			for(ibin=0;ibin<NBINS;ibin++){
				printf("%5.3f %7.4f\n",DELY*(0.5+ibin),bf[ibal][ibin]);
			}
		}
	}
}

void CBalance::WriteWeightTables(){
	int ires,jres,icode;
	int ab,ip,jp,ijprint[9]={-999999},ijcode[9]={2112,3122,3222,3112,3322,3312,3334,211,321};
	string iname[9]={"$p$","$\\Lambda$","$\\Sigma^+$","$\\Sigma^-$","$\\Xi^0$","$\\Xi^-$","$\\Omega^-$","$\\pi^+$","$K^+$"};
	string jname[9]={"$\\bar{p}$","$\\bar{\\Lambda}$","$\\bar{\\Sigma}^-$","$\\bar{\\Sigma}^+$","$\\bar{\\Xi}^0$","$\\bar{\\Xi}^+$","$\\bar{\\Omega}^+$","$\\pi^-$","$K^-$"};
	CResInfoMap::iterator rpos;
	CResInfo *resinfo;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		for(ip=0;ip<9;ip++){
			if(ijcode[ip]==resinfo->code){
				ijprint[ip]=resinfo->ires;
				//printf("ip=%2d, rho_h[%d]=%g \t %s\n",ip,resinfo->ires,rho_h[resinfo->ires],iname[ip].c_str());
			}
		}
	}
	string filename,dirname,command;
	dirname="TeXresults/"+PARSDIRNAME;
	command="mkdir -p "+dirname;
	system(command.c_str());
	filename=dirname+"/W.tex";
	FILE *fptr=fopen(filename.c_str(),"w");
	fprintf(fptr,"\\begin{tabular}{|r|c|c|c|c|c|c|c|c|c|c|c|} \\hline\n");
	fprintf(fptr,"                  &");
	for(ip=0;ip<9;ip++){
		fprintf(fptr," %13s",iname[ip].c_str());
		if(ip!=8)
			fprintf(fptr," &");
	}
	fprintf(fptr,"\\\\ \\hline\n"); 
	for(ip=0;ip<9;ip++){
		ires=ijprint[ip];
		fprintf(fptr,"%16s ",jname[ip].c_str());
		for(jp=0;jp<9;jp++){
			jres=ijprint[jp];
			for(ab=0;ab<4;ab++)
				fprintf(fptr," & %6.3f",-wB[ires][jres][ab]);
		}
		fprintf(fptr,"\\\\ \n");
	}
	fprintf(fptr,"\\hline\n\\end{tabular}");
	fclose(fptr);

}

void CBalance::BoostPart(CPart *partj0,CPart *partj,bool jdecay){
	double u[4],pprime[4],yi,yj,mt,phi,oldy,dely,rgx,rgy;
	int alpha;
	partj->Copy(partj0);
	randy->gauss2(&rgx,&rgy);
	u[1]=UPERPX*rgx;
	u[2]=UPERPY*rgy;
	u[3]=0.0;
	u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);

	partj->BoostP(u);
	if(jdecay) partj->BoostR(u);

	partj->y=atanh(partj->p[3]/partj->p[0]);
	if(jdecay) partj->eta=atanh(partj->r[3]/partj->r[0]);

	//partj->y=partj->y-parti->y;
	// boost so j particle has chance of being in acceptance
	oldy=partj->y;
	partj->y=acceptance_j->ETAMIN+randy->ran()*(acceptance_j->ETAMAX-acceptance_j->ETAMIN);
	dely=partj->y-oldy;
	partj->msquared=partj->resinfo->mass*partj->resinfo->mass;
	mt=sqrt(partj->msquared+partj->p[1]*partj->p[1]+partj->p[2]*partj->p[2]);
	partj->p[0]=mt*cosh(partj->y);
	partj->p[3]=mt*sinh(partj->y);
	partj->eta+=dely;
	partj->r[0]=partj->tau0*cosh(partj->eta);
	partj->r[3]=partj->tau0*sinh(partj->eta);
}

void CBalance::CalcDCA(bool decay,CPart *part,double *dca){
	int alpha;
	char nantestc[20];
	string nantests;
	double pdotr,p2,*p=part->p,*r=part->r;
	if(decay){
		p2=p[1]*p[1]+p[2]*p[2]+p[3]*p[3];
		pdotr=(p[1]*r[1]+p[2]*r[2]+p[3]*r[3])/p2;
		for(alpha=1;alpha<4;alpha++) dca[alpha]=(r[alpha]-pdotr*p[alpha])/1.0E13;
		dca[0]=sqrt(dca[1]*dca[1]+dca[2]*dca[2]+dca[3]*dca[3]);
		//if(r[0]!=r[0]){
		sprintf(nantestc,"%g",r[0]);
		nantests=nantestc;
		if(nantests=="NaN" || nantests=="nan" || nantests=="inf" || nantests=="INF"){
			printf("::: dca=(%g,%g,%g,%g)\n",dca[0],dca[1],dca[2],dca[3]);
			printf("::: r=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
			printf("::: p=(%g,%g,%g,%g)\n",r[0],r[1],r[2],r[3]);
			exit(1);
		}
	}
	else{
		for(alpha=0;alpha<4;alpha++)
			dca[alpha]=0.0;
	}
}

void CBalance::BoostParts(CPart *parti0,CPart *parti,CPart *partj0,CPart *partj,bool idecay,bool jdecay,int ab,double &extraweight){
	double u[4],ui[4],uj[4],pprime[4],yi,yj,mt,phi,oldy,dely,rgx,rgy,sigmaux,sigmauy,deleta,oldeta;
	int alpha;
	sigmaux=UPERPX*SIGMAPHI/sqrt(2.0);
	sigmauy=UPERPY*SIGMAPHI/sqrt(2.0);
	if(sigmaux>UPERPX) sigmaux=UPERPX;
	if(sigmauy>UPERPY) sigmauy=UPERPY;
	parti->Copy(parti0);
	partj->Copy(partj0);

	randy->gauss2(&rgx,&rgy);
	u[1]=sqrt(UPERPX*UPERPX-sigmaux*sigmaux)*rgx;
	u[2]=sqrt(UPERPY*UPERPY-sigmauy*sigmauy)*rgy;

	randy->gauss2(&rgx,&rgy);
	ui[1]=u[1]+sigmaux*rgx;
	ui[2]=u[2]+sigmauy*rgy;
	ui[3]=0.0;
	ui[0]=sqrt(1.0+ui[1]*ui[1]+ui[2]*ui[2]+ui[3]*ui[3]);

	randy->gauss2(&rgx,&rgy);
	uj[1]=u[1]+sigmaux*rgx;
	uj[2]=u[2]+sigmauy*rgy;
	uj[3]=0.0;
	uj[0]=sqrt(1.0+uj[1]*uj[1]+uj[2]*uj[2]+uj[3]*uj[3]);

	parti->BoostP(ui);
	if(idecay) parti->BoostR(ui);
	partj->BoostP(uj);
	if(jdecay) partj->BoostR(uj);

	parti->y=atanh(parti->p[3]/parti->p[0]);
	if(idecay)
		parti->eta=atanh(parti->r[3]/parti->r[0]);
	partj->y=atanh(partj->p[3]/partj->p[0]);
	if(jdecay)
		partj->eta=atanh(partj->r[3]/partj->r[0]);

	// boost so j particle has chance of being in acceptance
	oldy=partj->y;
	partj->y=acceptance_j->ETAMIN+randy->ran()*(acceptance_j->ETAMAX-acceptance_j->ETAMIN);
	dely=partj->y-oldy;
	partj->msquared=partj->resinfo->mass*partj->resinfo->mass;
	mt=sqrt(partj->msquared+partj->p[1]*partj->p[1]+partj->p[2]*partj->p[2]);
	partj->p[0]=mt*cosh(partj->y);
	partj->p[3]=mt*sinh(partj->y);
	partj->eta+=dely;
	partj->r[0]=partj->tau0*cosh(partj->eta);
	partj->r[3]=partj->tau0*sinh(partj->eta);

	if(ab>=0){
		// this is for interpair
		rapdist->GenX(ab,deleta,extraweight);
		oldeta=parti->eta;
		parti->eta=partj->eta+deleta;
		parti->y=parti->y+(parti->eta-oldeta);
		parti->msquared=parti->resinfo->mass*parti->resinfo->mass;
		mt=sqrt(parti->msquared+parti->p[1]*parti->p[1]+parti->p[2]*parti->p[2]);
		parti->p[0]=mt*cosh(parti->y);
		parti->p[3]=mt*sinh(parti->y);
		parti->r[0]=parti->tau0*cosh(parti->eta);
		parti->r[3]=parti->tau0*sinh(parti->eta);
	}
	else{
		// this for products from same decay products
		extraweight=1.0;
		parti->eta=parti->eta+dely;
		parti->y=parti->y+dely;
		parti->msquared=parti->resinfo->mass*parti->resinfo->mass;
		mt=sqrt(parti->msquared+parti->p[1]*parti->p[1]+parti->p[2]*parti->p[2]);
		parti->p[0]=mt*cosh(parti->y);
		parti->p[3]=mt*sinh(parti->y);
		parti->r[0]=parti->tau0*cosh(parti->eta);
		parti->r[3]=parti->tau0*sinh(parti->eta);
	}
}
