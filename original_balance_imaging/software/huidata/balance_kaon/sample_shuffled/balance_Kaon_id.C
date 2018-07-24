// Balance function code
// Gary D. Westfall
// Feb, 2010
// use sample tech to generate shuffled events

// Includes
#ifndef __CINT__
#include "TApplication.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>//
#include "string.h"
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TStyle.h>

using namespace std;
#endif

// Globals
//  Arrays for each event
float rap_pos[2000];
float rap_neg[2000];
//  Centrality
int myCentral;
//  Numbers for each event
int num_pos;
int num_neg;
//  Total numbers
double total_num_pos[10];
double total_num_neg[10];
double del_pos_pos[100][10];
double del_pos_neg[100][10];
double del_neg_neg[100][10];
double del_neg_pos[100][10];
// Binning
float binSize=0.04;

//  Canvas and pads
	TCanvas *c;
	TPad *pad0,*pad1,*pad2,*pad3,*pad4,*pad5,*pad6,*pad7,*pad8;	
	
	TH1D *ptPlusHist[10];
	TH1D *etaPlusHist[10];
	TH1D *phiPlusHist[10];
	TH1D *ptMinusHist[10];
	TH1D *etaMinusHist[10];
	TH1D *phiMinusHist[10];
	
	TH1D *Sample_ptHist;
	TH1D *Sample_etaHist;
	TH1D *Sample_phiHist;
	

// Prototypes
void balance_event(void);
void balance_main(void);
void zero_all(void);
void do_balance(void);
float mytheta(float eta);
float myp(float pt,float theta);
float gety(float eta,float pt);
void initCanvas(void);
void initHists(void);
void showItAll();

void balance_main(void){
	FILE *fd;
	int i;
	int totalNumberOfEvents;
	int myDisplay;
	int myLocal;
	int goodMult;
	short pid;
	int Charge;
	float pt;
	float eta;
	float phi;
	float y,p;
	int myindex;
	char name[100];	
	
	myDisplay=200000;
	myLocal=0;

	zero_all();
	
	
//root staff

	TFile *displaySimple = new TFile("../../balance_qinv_kaons_new_bin_TOF_extended/sample_shuffled/displaySimple.root");
//	displaySimple->ls();

for(i=0;i<10;i++) {
		
		sprintf(name,"pt_plus_%d",i);		
		ptPlusHist[i]=(TH1D*)displaySimple->Get(name);
		sprintf(name,"eta_plus_%d",i);
		etaPlusHist[i]=(TH1D*)displaySimple->Get(name);
		sprintf(name,"phi_plus_%d",i);
		phiPlusHist[i]=(TH1D*)displaySimple->Get(name);

		sprintf(name,"pt_minus_%d",i);		
		ptMinusHist[i]=(TH1D*)displaySimple->Get(name);
		sprintf(name,"eta_minus_%d",i);
		etaMinusHist[i]=(TH1D*)displaySimple->Get(name);
		sprintf(name,"phi_minus_%d",i);
		phiMinusHist[i]=(TH1D*)displaySimple->Get(name);
		
	}
	
	initCanvas();
	initHists();
	
// begin shuffle
	totalNumberOfEvents=0;
	
	cout << "Begin file ../../AuAu_identified.bin" << endl;
	fd=fopen("../../AuAu_identified.bin","rb");
	
	fread(&myCentral,sizeof(int),1,fd);
	while(!feof(fd)){
		myLocal++;
		totalNumberOfEvents++;
		if(myLocal == myDisplay){
			cout << "Event " << totalNumberOfEvents << endl;
			myLocal=0;
			showItAll();
		}
		num_pos=0;
		num_neg=0;
		fread(&goodMult,sizeof(int),1,fd);
		myindex=0;
		for(i=0;i<goodMult;i++){ // Loop over tracks in event
			fread(&pid,sizeof(short),1,fd);
			fread(&Charge,sizeof(int),1,fd);
			fread(&pt,sizeof(float),1,fd);
			fread(&eta,sizeof(float),1,fd);
			fread(&phi,sizeof(float),1,fd);
//			cout << "pid = " << pid << ", Charge = " << Charge << ", pt = " << pt << ", eta = " << eta << ", phi = " << phi << endl;
			if((pid == 2)||(pid == 5)){
//count number of pions in each event
				myindex++;
				}
			}

                myindex=myindex+(myindex%2); // make sure my index is even
//          cout<<"My index = "<< myindex <<", goodMult = "<< goodMult <<", my central = " << myCentral << endl;     
			for(i=0;i<myindex;i++){
				Charge=(int)pow(-1,i); // make sure charges are paired
				if(Charge>0){
				pt=ptPlusHist[myCentral]->GetRandom(); // get pt eta phi from sample
				eta=etaPlusHist[myCentral]->GetRandom();
				phi=phiPlusHist[myCentral]->GetRandom();
			}else{
				pt=ptMinusHist[myCentral]->GetRandom(); // get pt eta phi from sample
				eta=etaMinusHist[myCentral]->GetRandom();
				phi=phiMinusHist[myCentral]->GetRandom();
			}
				
				Sample_ptHist->Fill(pt);
				Sample_etaHist->Fill(eta);
				Sample_phiHist->Fill(phi);
				
//			cout<<"My index = "<< myindex << ", charge = " << Charge <<", pt = " << pt <<", eta = " << eta << ", phi = " << phi << endl;

				y=gety(eta,pt);
							 // choose kaons
				if(Charge > 0){
					rap_pos[num_pos]=y;
					num_pos++;
				}
				if(Charge < 0){
					rap_neg[num_neg]=y;
					num_neg++;
				}
			
		}
		balance_event();
		fread(&myCentral,sizeof(int),1,fd);
	}
	fclose(fd);

	cout << "Total number of events read = " << totalNumberOfEvents << endl;
	
	do_balance();
    c->Print("displaySimple.eps");
}

void do_balance(void){
	int i,j;
	float x;
	double my_balance;
	double my_error;
	double factor1,factor2,factor3,factor4;
	double dfactor1,dfactor2,dfactor3,dfactor4;
	double mysum;
	char fileName[80];
	FILE *fb;
	
	for(j=0;j<10;j++){
		sprintf(fileName,"balance_%d.txt",j);
		fb=fopen(fileName,"w");
		mysum=0.0;
		for(i=0;i<26;i++){
			x=0.1*(float)(i)+0.05;	
			factor1=del_pos_neg[i][j]-del_pos_pos[i][j];
			factor2=del_neg_pos[i][j]-del_neg_neg[i][j];
			factor3=factor1/total_num_pos[j];
			factor4=factor2/total_num_neg[j];
			my_balance=(0.5/0.1)*(factor3+factor4);
			mysum+=my_balance*0.1;
			dfactor1=sqrt(del_pos_neg[i][j]+del_pos_pos[i][j]);
			dfactor2=sqrt(del_neg_pos[i][j]+del_neg_neg[i][j]);
			dfactor3=dfactor1/total_num_pos[j];
			dfactor4=dfactor2/total_num_neg[j];
			my_error=(0.5/0.1)*sqrt(dfactor3*dfactor3+dfactor4*dfactor4);
			fprintf(fb,"%f %f %f\n",x,my_balance,my_error);
		}
		fprintf(fb,"%f\n",mysum);
		fclose(fb);
	}
}	

void balance_event(void){
	int i1,i2;
	float rel_rap;
	int index;
	float rappos;
	float rapneg;
	
	total_num_pos[myCentral]+=num_pos; // Increment N+
	total_num_neg[myCentral]+=num_neg; // Increment N-
	
	// Positives and the others
	for(i1=0;i1<num_pos;i1++){
		rappos=rap_pos[i1];
		// N++
		for(i2=0;i2<num_pos;i2++){
			if(i1 != i2){
				rel_rap=fabs(rappos-rap_pos[i2]);
				index=(int)(rel_rap/0.1);
				if((index >=0) && (index < 26)){
					del_pos_pos[index][myCentral]+=1.0;
				}
			}
		}
		// N+-
		for(i2=0;i2<num_neg;i2++){
			rel_rap=fabs(rappos-rap_neg[i2]);
			index=(int)(rel_rap/0.1);
			if((index >=0) && (index < 26)){
				del_pos_neg[index][myCentral]+=1.0;
			}
		}
	}
	
	// Negatives and the others
	for(i1=0;i1<num_neg;i1++){
		rapneg=rap_neg[i1];
		// N--
		for(i2=0;i2<num_neg;i2++){
			if(i1 != i2){
				rel_rap=fabs(rapneg-rap_neg[i2]);
				index=(int)(rel_rap/0.1);
				if((index >=0) && (index < 26)){
					del_neg_neg[index][myCentral]+=1.0;
				}
			}
		}
		// N-+
		for(i2=0;i2<num_pos;i2++){
			rel_rap=fabs(rapneg-rap_pos[i2]);
			index=(int)(rel_rap/0.1);
			if((index >=0) && (index < 26)){
				del_neg_pos[index][myCentral]+=1.0;
			}
		}
	}
}

void zero_all(void){
	int i,j;
	for(j=0;j<10;j++){
		total_num_pos[j]=0.0;
		total_num_neg[j]=0.0;
		for(i=0;i<26;i++){
			del_pos_pos[i][j]=0.0;
			del_pos_neg[i][j]=0.0;
			del_neg_neg[i][j]=0.0;
			del_neg_pos[i][j]=0.0;
		}
	}
}


float mytheta(float eta){
	return 2.0*atan(exp(-eta));
}

float myp(float pt,float theta){
	if(theta > 0.0){
		return pt/sin(theta);
	} else {
		return 0.0;
	}
}

float gety(float eta,float pt){
	float theta;
	float p;
	float pz;
	float E;
	float rap;
	float mK = 0.493677; // Kaon mass from PDG
	theta=2.0*atan(exp(-eta));
	if(theta <= 0.0){
		return 0.0;
	} else {
		p=pt/sin(theta);
		pz=p*cos(theta);
		E=sqrt(p*p+mK*mK);
		rap=0.5*log((E+pz)/(E-pz));
		return rap;
	}
}


void initCanvas(){
	c = new TCanvas("c","STAR Data",10,30,1200,900);
	pad0 = new TPad("pad0","pad0",0.00,0.66,0.33,1.0,21);
	pad1 = new TPad("pad1","pad1",0.33,0.66,0.66,1.0,21);
	pad2 = new TPad("pad2","pad2",0.66,0.66,1.00,1.0,21);
	pad3 = new TPad("pad3","pad3",0.00,0.33,0.33,0.66,21);
	pad4 = new TPad("pad4","pad4",0.33,0.33,0.66,0.66,21);
	pad5 = new TPad("pad5","pad5",0.66,0.33,1.00,0.66,21);
	pad6 = new TPad("pad6","pad6",0.00,0.00,0.33,0.33,21);
	pad7 = new TPad("pad7","pad7",0.33,0.00,0.66,0.33,21);
	pad8 = new TPad("pad8","pad8",0.66,0.00,1.00,0.33,21);

	pad0->Draw();
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();
	pad5->Draw();
	pad6->Draw();
	pad7->Draw();
	pad8->Draw();

	pad0->SetFillColor(0);
	pad1->SetFillColor(0);
	pad2->SetFillColor(0);
	pad3->SetFillColor(0);
	pad4->SetFillColor(0);
	pad5->SetFillColor(0);
	pad6->SetFillColor(0);
	pad7->SetFillColor(0);
	pad8->SetFillColor(0);

	cout << "Canvas initialized" << endl;
}

void showItAll(){

	pad0->cd();
	ptPlusHist[0]->Draw();
	pad1->cd();;
	etaPlusHist[0]->Draw();
	pad2->cd();
	phiPlusHist[0]->Draw();

	pad3->cd();
//	pad3->SetLogy(0);
	Sample_ptHist->SetMinimum(0);
	Sample_ptHist->Draw();
	

 	pad4->cd();
//	pad4->SetLogy(1);
	Sample_etaHist->SetMinimum(0);
	Sample_etaHist->Draw();

	pad5->cd();
//	pad5->SetLogy(1);
	Sample_phiHist->SetMinimum(1);
	Sample_phiHist->Draw();


	c->Update();
}

void initHists(){

	int icolor;
	
	icolor=2;
	Sample_ptHist = new TH1D("Sample_pt","Sample_pt",1000,0,1.8);
	Sample_ptHist->GetXaxis()->SetTitle("Sample_pt");
	Sample_ptHist->GetYaxis()->SetTitle("Counts");
	Sample_ptHist->SetFillColor(icolor);	
	icolor++;
	
	Sample_etaHist = new TH1D("Sample_eta","Sample_eta",1000,-1.2,1.2);
	Sample_etaHist->GetXaxis()->SetTitle("Sample_eta");
	Sample_etaHist->GetYaxis()->SetTitle("Counts");
	Sample_etaHist->SetFillColor(icolor);
	icolor++;
	
	Sample_phiHist = new TH1D("Sample_phi","Sample_phi",1000,-3.2,3.2);
	Sample_phiHist->GetXaxis()->SetTitle("Sample_phi");
	Sample_phiHist->GetYaxis()->SetTitle("Counts");
	Sample_phiHist->SetFillColor(icolor);
	icolor++;
}

#ifndef __CINT__
int main(int argc,char **argv){
	TApplication theApp("App",&argc,argv);
	balance_main();
	theApp.Run();
	return 0;
}
#endif