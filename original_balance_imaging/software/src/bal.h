#ifndef __bal_h__
#define __bal_h__
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <armadillo>
#include "b3d.h"
#include "parametermap.h"
#include "distributions.h"
//#include "TALICEAcceptance.h"

using namespace std;
class CAcceptance;
class CDistribution;

class CBalance{
public:
	int *ICODE,*JCODE,NBAL; // pids of particles for which to calc NBAL balance functions
	/* These are parameters you might change */
	double SIGMA_QGP,SIGMA_HAD,SIGMA_HAD_RATIO,SIGMA_GAUSS,SIGMAPHI;
	// eta widths of long-range, local contributions or for single-gauss width
	double T,UPERPX,UPERPY; // Blast-Wave parameters
	double LAMBDA_VISC; // shear correction to spectra, =lambda_zz. p -> lambda * p
	double TCHEM; // chemical equil. temperature
	int NBINS; // Number of Rapidty Bins for BF
	int NPHIBINS; // Number of Phi Bins for BF
	int NPHIBINS_RPLANE;
	int NEVENTS; // No. of Events for Monte Carlo in Calculating BF
	int NSAMPLE; // No. of Samples for new version MC
	double DELY; // Width of rapidity bins for BF
	double BARYON_FUDGE; // scales baryon popluations 
	double SOVERUD; // strange-quark to up-quark ratio in QGP
	bool *ALLCHARGES,GAUSSIAN,CALCRPLANEBF,ALLPARTS,OFFDIAG,USE_SAME_ACCEPTANCE,STAR_ACCEPTANCE_I,STAR_ACCEPTANCE_J;
	double meanpt_pion,meanpt_kaon,meanpt_proton;
	string PARSDIRNAME;
	parameterMap parmap;
	
	/* These are functions you might call from main program */
	CBalance(string parsdirname);
	void InitArrays();
	int GetAB(int a,int b);
	void SetRho(); // Sets dN/dy for each species
	void WriteBF(); // Writes Balance Functions vs. Dely
	void WriteG(); // Writes Balance Function Numerator vs. Dely
	void PrintBF(); // Prints Balance Functions vs. Dely on Screen
	void PrintG(); // Prints Correlations vs. Dely on Screen
	void WriteWeightTables(); // Writes table of local and QGP weights in tex format
	void WriteGQGPTables(); // Writes table of QGP G_{alpha,beta}
	void PrintGQGPTables(); // Writes table of QGP G_{alpha,beta}
	void CalcWeights(); // Calc.s weights wB[ires][jres][ab] 
	void CalcBF();  // Calculates BF[iy]
	void CalcMixed();
	void CalcBF_ALLCHARGES(); // For all charges, calculates BF in terms of pseudo-rapidity
	// same, except including BF(phi)

	double CalcSoverUD(); // Calc.s s-quark to u-quark ratio in QGP for ms=150, T=250
	void CalcDCA(bool decay,CPart *part,double *dca);
	
	/* You shouldn't change these things from main prog, but you might want to look at them */
	double *rho_h; // rho_h[NRES] Density of hadron species
	double *rhotot;
	double *rho_q; // rho_q[3] Density of quark species
	double **q; // q[ires][iq] charges of hadron species
	double **chi,**chi_inv,**mu,**g_q; // these are correlations and are [3][3] matrices
	double ***wB; // these are weights[ires][jres][ab]
	double **bf,**bf_had,**bf_qgp,**bf_resdecay; // this is the BF[idely] for specific ICODE,JCODE
	double **bf_phi,**bf_had_phi,**bf_qgp_phi,**bf_resdecay_phi; // same for BF[idelphi]
	double ***bf_rplane,***bf_had_rplane,***bf_qgp_rplane,***bf_resdecay_rplane,*denom,*denom_phi,**denom_rplane;
	double *rhoj,*rhoi;
	int NRES; // no. of resonances
	int *code; // pids of resonances
	double *normalization,*v2,*v2s,*v2c,*cosdphi,*gammap;
	
	vector<vector<double> > Mixed,dB,sigmaE;
	vector<vector<long long int> > sigmaN;
	
	CDistribution *rapdist;
	CAcceptance *acceptance,*acceptance_i,*acceptance_j;
	
	CB3D *b3d;
	CResList *reslist;
	CRandom *randy;
	bool GetProducts(int code,int &ni,CPart *parti,bool &decay);
	void BoostParts(CPart *parti0,CPart *parti,CPart *partj0,CPart *partj,bool idecay,bool jdecay,int ab,double &extraweight);
	void BoostPart(CPart *partj0,CPart *partj,bool jdecay);
	double totalhadrons,nquarks;
	CGSLMatrix_Real *gslmatrix;
	void InitWeightArrays();
	void InitBFArrays();
	void IncrementBFArrays(CPart *parti0,bool idecay,CPart *partj0,bool jdecay,double w,int ab);
	bool IJCheck(int ibal,CResInfo *resinfoi,CResInfo *resinfoj);
	bool IJCheck(int ibal,CPart *parti,CPart *partj);
	int GetIres();
};

class CAcceptance{
public:
	CAcceptance();
	double ETAMIN,ETAMAX;
	virtual void CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca);
};

class CAcceptance_CHEAP : public CAcceptance{
public:
	CAcceptance_CHEAP();
	double ptmin,ptmax;
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca);
};

/*
class CAcceptance_ALICE  : public CAcceptance{
public:
	CAcceptance_ALICE();
	string filename;
	double ptmin,ptmax;
	TALICEAcceptance *alice;
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca);
};
*/

class CAcceptance_STAR : public CAcceptance{
public:
	CAcceptance_STAR();
	int centrality;
	double ptmin,ptmax;
	void CalcAcceptance(bool &accept,double &efficiency,CPart *part,int centrality,double *dca);
	void star_acc_eff(int pid,double pt,double eta,double phi,int cen,bool &accept,double &eff);
};

#endif
