#include "rhicstat.h"
using namespace std;

int main(int argc, char *argv[]){
	int iz,irun;
	CRHICStat *rhicstat=new CRHICStat("statinfo/statpars.dat");
	//rhicstat->FitExpData();
	//rhicstat->PlotZvsX();
	//rhicstat->randy->reset(-time(NULL));
	rhicstat->Metropolis();
	//rhicstat->CalcCovariance();


	CRunInfo *ri;
/*
	for(irun=0;irun<rhicstat->NRUNS;irun++){
		ri=rhicstat->runinfo[irun];
		printf("%3d   %9.5f   %9.5f   %9.5f\n",irun,ri->netdiff_fit,ri->netdiff_exp,ri->netdiff_fitexp);
	}
*/
/*
	
	double overallcrappiness=0.0,*crappiness=new double[rhicstat->NZ],ctest=0.0,llerror;
	for(iz=0;iz<rhicstat->NZ;iz++)
		crappiness[iz]=0.0;
	printf("irun netdiff_fit netdiff_exp  netdiff_fitexp\n");
	for(irun=0;irun<rhicstat->NTESTRUNS;irun++){
		ri=rhicstat->testinfo[irun];
		llerror=0.0;
		for(iz=0;iz<rhicstat->NZ;iz++){
			crappiness[iz]+=pow(ri->z[iz]-ri->zfit[iz],2);
			llerror+=ri->zfiterror[iz];
		}
		printf("%3d   %9.5f   %9.5f   %9.5f %9.5f\n",irun,ri->netdiff_fit,ri->netdiff_exp,ri->netdiff_fitexp,llerror);
		ctest+=ri->netdiff_fit;
	}
	for(iz=0;iz<rhicstat->NZ;iz++){
		crappiness[iz]=crappiness[iz]/double(rhicstat->NTESTRUNS);
		printf("Crappiness[%d] = %g\n",iz,crappiness[iz]);
		overallcrappiness+=crappiness[iz];
	}
	printf("----- Overall Crappiness = %g=? %g-----\n",overallcrappiness,ctest/double(rhicstat->NTESTRUNS));

*/
	return 0;
}

