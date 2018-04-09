// Globals
float xa[40];
float ya[40];
float dya[40];
float avew;
float davew;
// Prototypes
void weightedAverage(int lower,int upper);

void draw_balance(){
	int i;
	float x,dx,y,dy;
	FILE *fp;
	float x0[26],y0[26],dx0[26],dy0[26];
	float x1[26],y1[26],dx1[26],dy1[26];
	float x2[26],y2[26],dx2[26],dy2[26];
	float x3[26],y3[26],dx3[26],dy3[26];
	float x4[26],y4[26],dx4[26],dy4[26];
	float x5[26],y5[26],dx5[26],dy5[26];
	float x6[26],y6[26],dx6[26],dy6[26];
	float x7[26],y7[26],dx7[26],dy7[26];
	float x8[26],y8[26],dx8[26],dy8[26];
	

	float w0,w1,w2,w3,w4,w5,w6,w7,w8;
	float e0,e1,e2,e3,e4,e5,e6,e7,e8;
	
	float xl,yl,xh,yh;
	
	float lowerLimit,upperLimit;
	float myMin,myMax;
	
	int symbolData,symbolMixed,symbolShuffled;
	int colorData,colorMixed,colorShuffled;
	
	float ave_data[9];
	float dave_data[9];
	float ave_shuffled[9];
	float dave_shuffled[9];
	float W;
	float dW;
	float errNum;
	float v1,v2;

	float sum1,sum2;
	float chi2;
	float Wave,dWave;

	// Real data
	fp=fopen("balance_0.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x0[i]=x;
		dx0[i]=0.0;
		y0[i]=y;
		dy0[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("balance_1.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x1[i]=x;
		dx1[i]=0.0;
		y1[i]=y;
		dy1[i]=dy;
	}
	fclose(fp);
	fp=fopen("balance_2.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x2[i]=x;
		dx2[i]=0.0;
		y2[i]=y;
		dy2[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("balance_3.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x3[i]=x;
		dx3[i]=0.0;
		y3[i]=y;
		dy3[i]=dy;
	}
	fclose(fp);
	fp=fopen("balance_4.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x4[i]=x;
		dx4[i]=0.0;
		y4[i]=y;
		dy4[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("balance_5.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x5[i]=x;
		dx5[i]=0.0;
		y5[i]=y;
		dy5[i]=dy;
	}
	fclose(fp);
	fp=fopen("balance_6.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x6[i]=x;
		dx6[i]=0.0;
		y6[i]=y;
		dy6[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("balance_7.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x7[i]=x;
		dx7[i]=0.0;
		y7[i]=y;
		dy7[i]=dy;
	}
	fclose(fp);
	fp=fopen("balance_8.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x8[i]=x;
		dx8[i]=0.0;
		y8[i]=y;
		dy8[i]=dy;
	}
	fclose(fp);
	


	// Set up a dark green
//	color=(TColor*)(gROOT->GetListOfColors()->At(3));
//	color->SetRGB(0.0,0.4,0.0);

	TGraph *gr0 = new TGraphErrors(26,x0,y0,dx0,dy0);
	TGraph *gr1 = new TGraphErrors(26,x1,y1,dx1,dy1);
	TGraph *gr2 = new TGraphErrors(26,x2,y2,dx2,dy2);
	TGraph *gr3 = new TGraphErrors(26,x3,y3,dx3,dy3);
	TGraph *gr4 = new TGraphErrors(26,x4,y4,dx4,dy4);
	TGraph *gr5 = new TGraphErrors(26,x5,y5,dx5,dy5);
	TGraph *gr6 = new TGraphErrors(26,x6,y6,dx6,dy6);
	TGraph *gr7 = new TGraphErrors(26,x7,y7,dx7,dy7);
	TGraph *gr8 = new TGraphErrors(26,x8,y8,dx8,dy8);
		


	TCanvas *c1 = new TCanvas("c1","Balance Function",0,0,1200,800);
	TPad *pad0 = new TPad("pad0","pad0",0.000,0.666,0.333,1.000,21);
	TPad *pad1 = new TPad("pad1","pad1",0.333,0.666,0.666,1.000,21);
	TPad *pad2 = new TPad("pad2","pad2",0.666,0.666,1.000,1.000,21);
	TPad *pad3 = new TPad("pad3","pad3",0.000,0.333,0.333,0.666,21);
	TPad *pad4 = new TPad("pad4","pad4",0.333,0.333,0.666,0.666,21);
	TPad *pad5 = new TPad("pad5","pad5",0.666,0.333,1.000,0.666,21);
	TPad *pad6 = new TPad("pad6","pad6",0.000,0.000,0.333,0.333,21);
	TPad *pad7 = new TPad("pad7","pad7",0.333,0.000,0.666,0.333,21);
	TPad *pad8 = new TPad("pad8","pad8",0.666,0.000,1.000,0.333,21);
	pad0->SetFillColor(0);
	pad1->SetFillColor(0);
	pad2->SetFillColor(0);
	pad3->SetFillColor(0);
	pad4->SetFillColor(0);
	pad5->SetFillColor(0);
	pad6->SetFillColor(0);
	pad7->SetFillColor(0);
	pad8->SetFillColor(0);
	pad0->Draw();
	pad1->Draw();
	pad2->Draw();
	pad3->Draw();
	pad4->Draw();
	pad5->Draw();
	pad6->Draw();
	pad7->Draw();
	pad8->Draw();

	xl=1.0;
	yl=0.52;
	xh=1.95;
	yh=0.58;
	
	TPaveLabel *p0 = new TPaveLabel(xl,yl,xh,yh,"Run7 kaon 0-5%");
	p0->SetBorderSize(0);
	p0->SetFillColor(0);
	TPaveLabel *p1 = new TPaveLabel(xl,yl,xh,yh,"5-10%");
	p1->SetBorderSize(0);
	p1->SetFillColor(0);
	TPaveLabel *p2 = new TPaveLabel(xl,yl,xh,yh,"10-20%");
	p2->SetBorderSize(0);
	p2->SetFillColor(0);
	TPaveLabel *p3 = new TPaveLabel(xl,yl,xh,yh,"20-30%");
	p3->SetBorderSize(0);
	p3->SetFillColor(0);
	TPaveLabel *p4 = new TPaveLabel(xl,yl,xh,yh,"30-40%");
	p4->SetBorderSize(0);
	p4->SetFillColor(0);
	TPaveLabel *p5 = new TPaveLabel(xl,yl,xh,yh,"40-50%");
	p5->SetBorderSize(0);
	p5->SetFillColor(0);
	TPaveLabel *p6 = new TPaveLabel(xl,yl,xh,yh,"50-60%");
	p6->SetBorderSize(0);
	p6->SetFillColor(0);
	TPaveLabel *p7 = new TPaveLabel(xl,yl,xh,yh,"60-70%");
	p7->SetBorderSize(0);
	p7->SetFillColor(0);
	TPaveLabel *p8 = new TPaveLabel(xl,yl,xh,yh,"70-80%");
	p8->SetBorderSize(0);
	p8->SetFillColor(0);
	
	TLine *line0 = new TLine(0,0,2,0);
	TLine *line1 = new TLine(0,0,2,0);
	TLine *line2 = new TLine(0,0,2,0);
	TLine *line3 = new TLine(0,0,2,0);
	TLine *line4 = new TLine(0,0,2,0);
	TLine *line5 = new TLine(0,0,2,0);
	TLine *line6 = new TLine(0,0,2,0);
	TLine *line7 = new TLine(0,0,2,0);
	TLine *line8 = new TLine(0,0,2,0);

	line0->SetLineWidth(1);
	line0->SetLineColor(1);
	line1->SetLineWidth(1);
	line1->SetLineColor(1);
	line2->SetLineWidth(1);
	line2->SetLineColor(1);
	line3->SetLineWidth(1);
	line3->SetLineColor(1);
	line4->SetLineWidth(1);
	line4->SetLineColor(1);
	line5->SetLineWidth(1);
	line5->SetLineColor(1);
	line6->SetLineWidth(1);
	line6->SetLineColor(1);
	line7->SetLineWidth(1);
	line7->SetLineColor(1);
	line8->SetLineWidth(1);
	line8->SetLineColor(1);


	TF1 *f0 = new TF1("f0","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f0->SetParameter(0,1);
	f0->SetParameter(1,0.5);
	TF1 *f1 = new TF1("f1","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f1->SetParameter(0,1);
	f1->SetParameter(1,0.5);
	TF1 *f2 = new TF1("f2","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f2->SetParameter(0,1);
	f2->SetParameter(1,0.5);
	TF1 *f3 = new TF1("f3","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f3->SetParameter(0,1);
	f3->SetParameter(1,0.5);
	TF1 *f4 = new TF1("f4","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f4->SetParameter(0,1);
	f4->SetParameter(1,0.5);
	TF1 *f5 = new TF1("f5","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f5->SetParameter(0,1);
	f5->SetParameter(1,0.5);
	TF1 *f6 = new TF1("f6","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f6->SetParameter(0,1);
	f6->SetParameter(1,0.5);
	TF1 *f7 = new TF1("f7","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f7->SetParameter(0,1);
	f7->SetParameter(1,0.5);
	TF1 *f8 = new TF1("f8","[0]*exp(-(x*x)/(2.0*[1]*[1]))",0.1,2.6);
	f8->SetParameter(0,1);
	f8->SetParameter(1,0.5);
	

	c1->SetFillColor(0);
	c1->SetGrid();
	c1->SetLogy(0);
	c1->GetFrame()->SetFillColor(0);
	c1->GetFrame()->SetBorderSize(12);

	symbolData=20;
	symbolMixed=21;
	symbolShuffled=24;
	
	colorData=2;
	colorMixed=4;
	colorShuffled=3;
	
	// Data
	myMin=-0.01;
	myMax=1.5;

	gr0->SetMinimum(myMin);
	gr0->SetMaximum(myMax);
	gr0->SetMarkerStyle(symbolData);
	gr0->SetMarkerColor(colorData);
	gr0->SetMarkerSize(1.0);
	gr0->SetLineColor(colorData);
	gr0->SetLineWidth(1);
	gr0->SetTitle("");


	gr1->SetMinimum(myMin);
	gr1->SetMaximum(myMax);
	gr1->SetMarkerStyle(symbolData);
	gr1->SetMarkerColor(colorData);
	gr1->SetMarkerSize(1.0);
	gr1->SetLineColor(colorData);
	gr1->SetLineWidth(1);
	gr1->SetTitle("");
	

	gr2->SetMinimum(myMin);
	gr2->SetMaximum(myMax);
	gr2->SetMarkerStyle(symbolData);
	gr2->SetMarkerColor(colorData);
	gr2->SetMarkerSize(1.0);
	gr2->SetLineColor(colorData);
	gr2->SetLineWidth(1);
	gr2->SetTitle("");

	gr3->SetMinimum(myMin);
	gr3->SetMaximum(myMax);
	gr3->SetMarkerStyle(symbolData);
	gr3->SetMarkerColor(colorData);
	gr3->SetMarkerSize(1.0);
	gr3->SetLineColor(colorData);
	gr3->SetLineWidth(1);
	gr3->SetTitle("");
	
	gr4->SetMinimum(myMin);
	gr4->SetMaximum(myMax);
	gr4->SetMarkerStyle(symbolData);
	gr4->SetMarkerColor(colorData);
	gr4->SetMarkerSize(1.0);
	gr4->SetLineColor(colorData);
	gr4->SetLineWidth(1);
	gr4->SetTitle("");
	

	gr5->SetMinimum(myMin);
	gr5->SetMaximum(myMax);
	gr5->SetMarkerStyle(symbolData);
	gr5->SetMarkerColor(colorData);
	gr5->SetMarkerSize(1.0);
	gr5->SetLineColor(colorData);
	gr5->SetLineWidth(1);
	gr5->SetTitle("");
	

	gr6->SetMinimum(myMin);
	gr6->SetMaximum(myMax);
	gr6->SetMarkerStyle(symbolData);
	gr6->SetMarkerColor(colorData);
	gr6->SetMarkerSize(1.0);
	gr6->SetLineColor(colorData);
	gr6->SetLineWidth(1);
	gr6->SetTitle("");
	

	gr7->SetMinimum(myMin);
	gr7->SetMaximum(myMax);
	gr7->SetMarkerStyle(symbolData);
	gr7->SetMarkerColor(colorData);
	gr7->SetMarkerSize(1.0);
	gr7->SetLineColor(colorData);
	gr7->SetLineWidth(1);
	gr7->SetTitle("");


	gr8->SetMinimum(myMin);
	gr8->SetMaximum(myMax);
	gr8->SetMarkerStyle(symbolData);
	gr8->SetMarkerColor(colorData);
	gr8->SetMarkerSize(1.0);
	gr8->SetLineColor(colorData);
	gr8->SetLineWidth(1);
	gr8->SetTitle("");
	

	f0->SetLineColor(colorData);
	f0->SetLineWidth(2);
	f1->SetLineColor(colorData);
	f1->SetLineWidth(2);
	f2->SetLineColor(colorData);
	f2->SetLineWidth(2);
	f3->SetLineColor(colorData);
	f3->SetLineWidth(2);
	f4->SetLineColor(colorData);
	f4->SetLineWidth(2);
	f5->SetLineColor(colorData);
	f5->SetLineWidth(2);
	f6->SetLineColor(colorData);
	f6->SetLineWidth(2);
	f7->SetLineColor(colorData);
	f7->SetLineWidth(2);
	f8->SetLineColor(colorData);
	f8->SetLineWidth(2);

	gr0->Fit("f0","R");
	gr1->Fit("f1","R");
	gr2->Fit("f2","R");
	gr3->Fit("f3","R");
	gr4->Fit("f4","R");
	gr5->Fit("f5","R");
	gr6->Fit("f6","R");
	gr7->Fit("f7","R");
	gr8->Fit("f8","R");

	pad0->cd();
	gr0->Draw("AP");
	gr0->GetXaxis()->SetLimits(0.0,2.6);
	gr0->GetXaxis()->SetTitle("#Delta#eta");
	gr0->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr0->GetXaxis()->CenterTitle();
	gr0->GetYaxis()->CenterTitle();
	gr0->Draw("AP");
	line0->Draw();
	gr0->Draw("P");
	p0->Draw();
	pad1->cd();
	gr1->Draw("AP");
	gr1->Draw("AP");
	gr1->GetXaxis()->SetLimits(0.0,2.6);
	gr1->GetXaxis()->SetTitle("#Delta#eta");
	gr1->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	gr1->Draw("AP");
	line1->Draw();
	gr1->Draw("P");
	p1->Draw();
	pad2->cd();
	gr2->Draw("AP");
	gr2->GetXaxis()->SetLimits(0.0,2.6);
	gr2->GetXaxis()->SetTitle("#Delta#eta");
	gr2->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	gr2->Draw("AP");
	line2->Draw();
	gr2->Draw("P");
	p2->Draw();
	pad3->cd();
	gr3->Draw("AP");
	gr3->GetXaxis()->SetLimits(0.0,2.6);
	gr3->GetXaxis()->SetTitle("#Delta#eta");
	gr3->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr3->GetXaxis()->CenterTitle();
	gr3->GetYaxis()->CenterTitle();
	gr3->Draw("AP");
	line3->Draw();
	gr3->Draw("P");
	p3->Draw();
	pad4->cd();
	gr4->Draw("AP");
	gr4->GetXaxis()->SetLimits(0.0,2.6);
	gr4->GetXaxis()->SetTitle("#Delta#eta");
	gr4->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr4->GetXaxis()->CenterTitle();
	gr4->GetYaxis()->CenterTitle();
	gr4->Draw("AP");
	line4->Draw();
	gr4->Draw("P");
	p4->Draw();
	pad5->cd();
	gr5->Draw("AP");
	gr5->GetXaxis()->SetLimits(0.0,2.6);
	gr5->GetXaxis()->SetTitle("#Delta#eta");
	gr5->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr5->GetXaxis()->CenterTitle();
	gr5->GetYaxis()->CenterTitle();
	gr5->Draw("AP");
	line5->Draw();
	gr5->Draw("P");
	p5->Draw();
	pad6->cd();
	gr6->Draw("AP");
	gr6->GetXaxis()->SetLimits(0.0,2.6);
	gr6->GetXaxis()->SetTitle("#Delta#eta");
	gr6->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr6->GetXaxis()->CenterTitle();
	gr6->GetYaxis()->CenterTitle();
	gr6->Draw("AP");
	line6->Draw();
	gr6->Draw("P");
	p6->Draw();
	pad7->cd();
	gr7->Draw("AP");
	gr7->GetXaxis()->SetLimits(0.0,2.6);
	gr7->GetXaxis()->SetTitle("#Delta#eta");
	gr7->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr7->GetXaxis()->CenterTitle();
	gr7->GetYaxis()->CenterTitle();
	gr7->Draw("AP");
	line7->Draw();
	gr7->Draw("P");
	p7->Draw();
	pad8->cd();
	gr8->Draw("AP");
	gr8->GetXaxis()->SetLimits(0.0,2.6);
	gr8->GetXaxis()->SetTitle("#Delta#eta");
	gr8->GetYaxis()->SetTitle("B(#Delta#eta)");
	gr8->GetXaxis()->CenterTitle();
	gr8->GetYaxis()->CenterTitle();
	gr8->Draw("AP");
	line8->Draw();
	gr8->Draw("P");
	p8->Draw();
	
	c1->Print("balance_all.eps");
	
	w0=fabs(f0->GetParameter(1));
	e0=fabs(f0->GetParError(1));
	w1=fabs(f1->GetParameter(1));
	e1=fabs(f1->GetParError(1));
	w2=fabs(f2->GetParameter(1));
	e2=fabs(f2->GetParError(1));
	w3=fabs(f3->GetParameter(1));
	e3=fabs(f3->GetParError(1));
	w4=fabs(f4->GetParameter(1));
	e4=fabs(f4->GetParError(1));
	w5=fabs(f5->GetParameter(1));
	e5=fabs(f5->GetParError(1));
	w6=fabs(f6->GetParameter(1));
	e6=fabs(f6->GetParError(1));
	w7=fabs(f7->GetParameter(1));
	e7=fabs(f7->GetParError(1));
	w8=fabs(f8->GetParameter(1));
	e8=fabs(f8->GetParError(1));
	
	fp=fopen("width_all.txt","w");
	fprintf(fp,"%f %f %f %f\n",1,0,w0,e0);
	fprintf(fp,"%f %f %f %f\n",2,0,w1,e1);
	fprintf(fp,"%f %f %f %f\n",3,0,w2,e2);
	fprintf(fp,"%f %f %f %f\n",4,0,w3,e3);
	fprintf(fp,"%f %f %f %f\n",5,0,w4,e4);
	fprintf(fp,"%f %f %f %f\n",6,0,w5,e5);
	fprintf(fp,"%f %f %f %f\n",7,0,w6,e6);
	fprintf(fp,"%f %f %f %f\n",8,0,w7,e7);
	fprintf(fp,"%f %f %f %f\n",9,0,w8,e8);
	fclose(fp);

	// Do the weighted averages here in the plot macros instead of later with weightedAverage
	// Real data
	cout << "AMPT" << endl;
	fp=fopen("ave_data_all.txt","w");
		for(i=0;i<26;i++){
			xa[i]=x0[i];
			ya[i]=y0[i];
			dya[i]=dy0[i];
		}
		weightedAverage(2,25);
		cout << "0: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"0 %f %f\n",avew,davew);
		ave_data[0]=avew;
		dave_data[0]=davew;

		for(i=0;i<26;i++){
			xa[i]=x1[i];
			ya[i]=y1[i];
			dya[i]=dy1[i];
		}
		weightedAverage(2,25);
		cout << "1: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"1 %f %f\n",avew,davew);
		ave_data[1]=avew;
		dave_data[1]=davew;

		for(i=0;i<26;i++){
			xa[i]=x2[i];
			ya[i]=y2[i];
			dya[i]=dy2[i];
		}
		weightedAverage(2,25);
		cout << "2: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"2 %f %f\n",avew,davew);
		ave_data[2]=avew;
		dave_data[2]=davew;

		for(i=0;i<26;i++){
			xa[i]=x3[i];
			ya[i]=y3[i];
			dya[i]=dy3[i];
		}
		weightedAverage(2,25);
		cout << "3: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"3 %f %f\n",avew,davew);
		ave_data[3]=avew;
		dave_data[3]=davew;

		for(i=0;i<26;i++){
			xa[i]=x4[i];
			ya[i]=y4[i];
			dya[i]=dy4[i];
		}
		weightedAverage(2,25);
		cout << "4: avew = " << avew << ", davew = " << davew << endl;
		ave_data[4]=avew;
		dave_data[4]=davew;
		fprintf(fp,"4 %f %f\n",avew,davew);

		for(i=0;i<26;i++){
			xa[i]=x5[i];
			ya[i]=y5[i];
			dya[i]=dy5[i];
		}
		weightedAverage(2,25);
		cout << "5: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"5 %f %f\n",avew,davew);
		ave_data[5]=avew;
		dave_data[5]=davew;

		for(i=0;i<26;i++){
			xa[i]=x6[i];
			ya[i]=y6[i];
			dya[i]=dy6[i];
		}
		weightedAverage(2,25);
		cout << "6: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"6 %f %f\n",avew,davew);
		ave_data[6]=avew;
		dave_data[6]=davew;

		for(i=0;i<26;i++){
			xa[i]=x7[i];
			ya[i]=y7[i];
			dya[i]=dy7[i];
		}
		weightedAverage(2,25);
		cout << "7: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"7 %f %f\n",avew,davew);
		ave_data[7]=avew;
		dave_data[7]=davew;

		for(i=0;i<26;i++){
			xa[i]=x8[i];
			ya[i]=y8[i];
			dya[i]=dy8[i];
		}
		weightedAverage(2,25);
		cout << "8: avew = " << avew << ", davew = " << davew << endl;
		fprintf(fp,"8 %f %f\n",avew,davew);
		ave_data[8]=avew;
		dave_data[8]=davew;

	fclose(fp);


	


}

void weightedAverage(int lower,int upper){
	int i;
	float sum1,sum2,sum3,sum4;
	float dsum1,dsum2;
	float f1,f2;
	sum1=0.0;
	sum2=0.0;
	sum3=0.0;
	sum4=0.0;
	for(i=lower;i<upper;i++){
		sum1+=xa[i]*ya[i];
		sum2+=ya[i];
		sum3+=(xa[i]*dya[i])*(xa[i]*dya[i]);
		sum4+=dya[i]*dya[i];
	}
	avew=sum1/sum2;
	dsum1=sqrt(sum3);
	dsum2=sqrt(sum4);
	f1=dsum1/sum1;
	f2=dsum2/sum2;
	davew=avew*sqrt(f1*f1+f2*f2);
}
