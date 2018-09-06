#include "coral.h"
#include "xgraph.h"

int main (int argc, char *argv[]){

	CXGraph xgraph(400,400,50,50);
	xgraph.setaxes(0.0,-8,60.0,-4);
	xgraph.clear();
	xgraph.axesinfo.yintercept=-6.0;
	xgraph.drawaxes();
	//Misc::Pause(1);

	CXGraph thetagraph(400,400,500,50);
	thetagraph.setaxes(0.0,-8,90.0,-4);
	thetagraph.clear();
	thetagraph.axesinfo.yintercept=-6.0;
	thetagraph.drawaxes();

	CCHArray *A=new CCHArray("parameters/pipi/aparsCH_source.dat");
	A->ReadAX("/Users/scottepratt/data/sdata/OSCAR_pionsonly_kt320");
	A->FillRemainderX();

	double *Aprojection[4],*r;
	double phi,*B1,*B2,*B3,*theta;
	int i,ir,itheta,ntheta=40,iphi,nr=A->GetNRADIAL();
	B1=new double[ntheta]; B2=new double[ntheta]; B3=new double[ntheta]; theta=new double[ntheta];

	r=new double[nr];
	for(ir=0;ir<nr;ir++) r[ir]=(0.5+ir)*A->GetRADSTEP();

	for(i=0;i<4;i++) Aprojection[i]=new double[nr];
	A->PrintProjections();
	A->GetProjections(Aprojection);
	for(i=0;i<4;i++)
		for(ir=0;ir<nr;ir++) Aprojection[i][ir]=log10(Aprojection[i][ir]);
	
	xgraph.setcolor("black");
	xgraph.plotpoints(r,Aprojection[0],nr);
	xgraph.setcolor("red");
	xgraph.plotpoints(r,Aprojection[1],nr);
	xgraph.setcolor("green");
	xgraph.plotpoints(r,Aprojection[2],nr);
	xgraph.setcolor("blue");
	xgraph.plotpoints(r,Aprojection[3],nr);

	ir=0;
	while(ir>=0){
		printf("Enter ir for theta/phi plotting (negative to quit): ");
		scanf("%d",&ir);
		for(i=0;i<ntheta;i++){
			theta[i]=(0.5+i)*90.0/double(ntheta);
			phi=0.0;
			B1[i]=log10(A->AExpand(theta[i]*PI/180.0,phi,ir));
			phi=0.5*PI;
			B2[i]=log10(A->AExpand(theta[i]*PI/180.0,phi,ir));
			B3[i]=log10(A->AExpand(0.5*PI,theta[i]*PI/180.0,ir));
		}
		thetagraph.clear();
		thetagraph.drawaxes();
		thetagraph.setcolor("red");
		thetagraph.plotpoints(theta,B1,ntheta);
		thetagraph.setcolor("green");
		thetagraph.plotpoints(theta,B2,ntheta);
		thetagraph.setcolor("blue");
		thetagraph.plotpoints(theta,B3,ntheta);

	}

	return 0;
}
