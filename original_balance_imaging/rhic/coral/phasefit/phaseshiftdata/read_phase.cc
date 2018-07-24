//read_phase.cc
//1, identify the interaction and partial wave, and read data file
//2, translate the relative momentum into lab momentum and identify its vicinity
//3, interpolation and return phase shift

#include <iostream>
#include "readphaseshift.h"
//#include "lagrange.h"
#include "interpolation.h"

CReadPhaseShift::CReadPhaseShift(const char szInteractionInput[],  const char szPartialWaveInput[])
{
  const double MPI=139.58;
  const double MKAON=493.677;
  const double MPROTON=938.271;
  const double MNEUTRON=939.565;
  const double MLAMBDA=1115.7;
	
  strcpy(szInteraction, szInteractionInput);
  strcpy(szPartialWave, szPartialWaveInput);
  //identify the interaction and partial wave
  if(0==strcmp(szInteraction, "PPi"))
	{
		m1 = MPROTON;
		m2 = MPI;
		//printf("%s\n", szInteraction);
	}
  if(0==strcmp(szInteraction, "PP"))
	{
		m1 = MPROTON;
		m2 = MPROTON;
	}
  if(0==strcmp(szInteraction, "PK"))
	{
		m1 = MPROTON;
		m2 = MKAON;
	}
  if(0==strcmp(szInteraction, "PiPi"))
	{
		m1 = MPI;
		m2 = MPI;
	}
  if(ReadData()==false)
	{
		printf("failed to read from data!\n");
		exit(0);
	}
}

bool CReadPhaseShift::ReadData()
{
  //read data file
  char szFileName[80];
  sprintf(szFileName,"phaseshiftdata/");
  strcat(szFileName, szInteraction);
  strcat(szFileName, szPartialWave);
  strcat(szFileName, ".txt");
  //printf("%s\n", szFileName);
	
  FILE * rfile;
  rfile = fopen(szFileName, "r");
	
  if(NULL != rfile)
	{
		//printf("%s open successfully!\n", szFileName);
		char szTest1[4], szTest2[4], szTest3[4], szTest4[4], szTest5[4], szTest6[4], szTest7[4], szTest8[4], szTest9[4];
		double dum1, dum2, dum3, dum4, dum5, dum6, dum7;
		fscanf(rfile, "%s", szTest1);
		fscanf(rfile, "%s", szTest2);
		fscanf(rfile, "%s", szTest3);
		fscanf(rfile, "%s", szTest4);
		fscanf(rfile, "%s", szTest5);
		fscanf(rfile, "%s", szTest6);
		fscanf(rfile, "%s", szTest7);
		fscanf(rfile, "%s", szTest8);
		fscanf(rfile, "%s", szTest9);
		//printf("%s %s %s %s %s %s %s %s %s is ready!\n", szTest1, szTest2, szTest3, szTest4, szTest5, szTest6, szTest7, szTest8, szTest9);
		for(int ii=0;ii<MAX_DATA;ii++)
		{
			fscanf(rfile, "%lf", &plab[ii]);
			fscanf(rfile, "%lf", &del[ii]);
			fscanf(rfile, "%lf", &dum1);
			fscanf(rfile, "%lf", &dum2);
			fscanf(rfile, "%lf", &dum3);
			fscanf(rfile, "%lf", &dum4);
			fscanf(rfile, "%lf", &dum5);
			fscanf(rfile, "%lf", &dum6);
			fscanf(rfile, "%lf", &dum7);
			//printf("%g %g\n", plab[ii], del[ii]);
		}
	}
  else
	{
		printf("Cannot open file %s!\n", szFileName);
		fclose(rfile);
		return false;
	}
  fclose(rfile);
  return true;
}

double CReadPhaseShift::ReadPhaseShift(double q)
{
  //translate into lab momentum
  double p= sqrt((sqrt(q*q+m1*m1)*sqrt(q*q+m2*m2)+q*q)*(sqrt(q*q+m1*m1)*sqrt(q*q+m2*m2)+q*q)/m1/m1-m2*m2);
  double e1,e2,rootsa,rootsb,s,oldp;
  oldp=p;
  e1=sqrt(m1*m1+q*q);
  e2=sqrt(m2*m2+q*q);
  s=(e1+e2)*(e1+e2);
  p=sqrt(s*s+pow(m1,4)+pow(m2,4)-2.0*m1*m1*s-2.0*m2*m2*s-2.0*m1*m1*m2*m2)/(2.0*m1);
	
  e1=m1; e2=sqrt(m2*m2+p*p);
  //e2=m2; e1=sqrt(m1*m1+p*p);
  rootsa=sqrt(pow(e1+e2,2)-p*p);
  e1=sqrt(m1*m1+q*q);
  e2=sqrt(m2*m2+q*q);
  rootsb=e1+e2;
  if(fabs(rootsa-rootsb)>1.0){
    printf("The translated lab momentum is %g, Li's p=%g, q=%g, roots=%g =? %g\n", p,oldp,q,rootsa,rootsb);
    exit(1);
  }
	
  //search and locate the nearest plab[] to p
  //bi-divid method
  if(p<plab[0])
	{
		printf("Input p=%g smaller than lower bound of plab=%g\n", p, plab[0]);
		return 0.0;
	}
  if(p>plab[MAX_DATA-1])
	{
		printf("Input p=%g greater than upper bound of plab=%g\n", p, plab[MAX_DATA-1]);
		return 0.0;
	}
  int nfloor = 0;
  int ncelling = MAX_DATA - 1;
  int npos = (nfloor + ncelling) / 2;
  do
	{
		if(fabs(p-plab[npos])<1.0e-6)
		{
			nfloor = ncelling = npos;
		}
		if(p<plab[npos])
		{
			ncelling = npos;
		}
		if(p>plab[npos])
		{
			nfloor = npos;
		}
		npos = (ncelling + nfloor) / 2;
	}while(ncelling - nfloor > 1);
	
  //printf("Locate the nearest plab[%d]=%g for p=%g\n", npos, plab[npos], p);
	
  //interpolation
  double InterDel = 0.0;
  double dInterDel = 0.0;
  double * pInterDel = &InterDel;
  double * pdInterDel = &dInterDel;
  //Here we are using 5-point interpolation, which is good enough.
  int npoint = 5;
  if(npos - npoint/2 <= 0)
	{
		polint(&plab[npos], &del[npos], npoint, p, pInterDel, pdInterDel);
	}
  else
	{
		if(npos + npoint/2 > MAX_DATA)
		{
			polint(&plab[npos-npoint], &del[npos-npoint], npoint, p, pInterDel, pdInterDel);
		}
		else
		{
			polint(&plab[npos-npoint/2], &del[npos-npoint/2], npoint, p, pInterDel, pdInterDel);
		}
	}
	
  //printf("The interpolation of p=%g is del=%g\n", p, InterDel);
	
  //return phaseshift
  const double piover180=4.0*atan(1.0)/180.0;
  return *pInterDel*piover180;
}

/*
void main()
{
	CReadPhaseShift rd("PPi", "p11");
	double q;
	do
	{
		scanf("%lf", &q);
		printf("The phase shift of q=%g is del=%g\n", q, rd.ReadPhaseShift(q));
	}while(q>0.0);
	printf("success!\n");
	int a;
	scanf("%d", &a);
}
*/
