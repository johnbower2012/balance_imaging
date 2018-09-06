using namespace std;

#include "utilities.h"
#include "coral.h"

int main(){
  CWaveFunction *wf;
  double q,r,ctheta=1.0,Rx,Ry,Rz,x,y,z;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-time(NULL));
  Rx=Ry=4.0; Rz=4.0;
	
  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);
	
  wf=new CWaveFunction_pipluspiplus_sqwell("parameters/wfparameters.dat");
	//wf->PrintPhaseShifts();
	
  Nq=wf->GetNQMAX();
  c=new double[Nq];
  for(iq=0;iq<Nq;iq++){
    q=wf->GetQ(iq);
    c[iq]=0.0;
    for(imc=0;imc<NMC;imc++){
      x=Rx*sqrt(2.0)*random.gauss();
      y=Ry*sqrt(2.0)*random.gauss();
      z=Rz*sqrt(2.0)*random.gauss();
      r=sqrt(x*x+y*y+z*z);
      ctheta=z/r;
      //c[iq]+=wf->GetPsiSquared(q,r,ctheta);
      c[iq]+=wf->CalcPsiSquared(iq,r,ctheta);
    }
    c[iq]=c[iq]/double(NMC);
    printf("%5.2f : %g\n",q,c[iq]);
  }
	
}

