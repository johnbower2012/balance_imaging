{
	gROOT->Reset();
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
	
	float x10[26],y10[26],dx10[26],dy10[26];
	float x11[26],y11[26],dx11[26],dy11[26];
	float x12[26],y12[26],dx12[26],dy12[26];
	float x13[26],y13[26],dx13[26],dy13[26];
	float x14[26],y14[26],dx14[26],dy14[26];
	float x15[26],y15[26],dx15[26],dy15[26];
	float x16[26],y16[26],dx16[26],dy16[26];
	float x17[26],y17[26],dx17[26],dy17[26];
	float x18[26],y18[26],dx18[26],dy18[26];

	float x20[26],y20[26],dx20[26],dy20[26];
	float x21[26],y21[26],dx21[26],dy21[26];
	float x22[26],y22[26],dx22[26],dy22[26];
	float x23[26],y23[26],dx23[26],dy23[26];
	float x24[26],y24[26],dx24[26],dy24[26];
	float x25[26],y25[26],dx25[26],dy25[26];
	float x26[26],y26[26],dx26[26],dy26[26];
	float x27[26],y27[26],dx27[26],dy27[26];
	float x28[26],y28[26],dx28[26],dy28[26];

	float w0,w1,w2,w3,w4,w5,w6,w7,w8;
	float e0,e1,e2,e3,e4,e5,e6,e7,e8;
	
	float xl,yl,xh,yh;
	
	float lowerLimit,upperLimit;
	float myMin,myMax;
	
	int symbolData,symboldMixed,symbolShuffled;
	int colorData,colorMixed,colorShuffled;
	
	float sum;
	float integral[9];
	float binSize=0.1;	
	// data-mixed
	fp=fopen("subtract_mixed/balance_0.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x0[i]=x;
		dx0[i]=0.0;
		y0[i]=y;
		dy0[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[0]=sum;
	
	fp=fopen("subtract_mixed/balance_1.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x1[i]=x;
		dx1[i]=0.0;
		y1[i]=y;
		dy1[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[1]=sum;

	fp=fopen("subtract_mixed/balance_2.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x2[i]=x;
		dx2[i]=0.0;
		y2[i]=y;
		dy2[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[2]=sum;
	
	fp=fopen("subtract_mixed/balance_3.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x3[i]=x;
		dx3[i]=0.0;
		y3[i]=y;
		dy3[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[3]=sum;

	fp=fopen("subtract_mixed/balance_4.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x4[i]=x;
		dx4[i]=0.0;
		y4[i]=y;
		dy4[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[4]=sum;
	
	fp=fopen("subtract_mixed/balance_5.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x5[i]=x;
		dx5[i]=0.0;
		y5[i]=y;
		dy5[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[5]=sum;

	fp=fopen("subtract_mixed/balance_6.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x6[i]=x;
		dx6[i]=0.0;
		y6[i]=y;
		dy6[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[6]=sum;
	
	fp=fopen("subtract_mixed/balance_7.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x7[i]=x;
		dx7[i]=0.0;
		y7[i]=y;
		dy7[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[7]=sum;

	fp=fopen("subtract_mixed/balance_8.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x8[i]=x;
		dx8[i]=0.0;
		y8[i]=y;
		dy8[i]=dy;
		sum+=binSize*y;
	}
	fclose(fp);
	integral[8]=sum;
	
	// mixed
	fp=fopen("mixed/balance_0.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x10[i]=x;
		dx10[i]=0.0;
		y10[i]=y;
		dy10[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("mixed/balance_1.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x11[i]=x;
		dx11[i]=0.0;
		y11[i]=y;
		dy11[i]=dy;
	}
	fclose(fp);
	fp=fopen("mixed/balance_2.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x12[i]=x;
		dx12[i]=0.0;
		y12[i]=y;
		dy12[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("mixed/balance_3.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x13[i]=x;
		dx13[i]=0.0;
		y13[i]=y;
		dy13[i]=dy;
	}
	fclose(fp);
	fp=fopen("mixed/balance_4.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x14[i]=x;
		dx14[i]=0.0;
		y14[i]=y;
		dy14[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("mixed/balance_5.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x15[i]=x;
		dx15[i]=0.0;
		y15[i]=y;
		dy15[i]=dy;
	}
	fclose(fp);
	fp=fopen("mixed/balance_6.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x16[i]=x;
		dx16[i]=0.0;
		y16[i]=y;
		dy16[i]=dy;
	}
	fclose(fp);
	
	fp=fopen("mixed/balance_7.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x17[i]=x;
		dx17[i]=0.0;
		y17[i]=y;
		dy17[i]=dy;
	}
	fclose(fp);
	fp=fopen("mixed/balance_8.txt","r");
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x18[i]=x;
		dx18[i]=0.0;
		y18[i]=y;
		dy18[i]=dy;
	}
	fclose(fp);

	// Shuffled
	fp=fopen("sample_shuffled/balance_0.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x20[i]=x;
		dx20[i]=0.0;
		y20[i]=integral[0]*y;
		dy20[i]=integral[0]*dy;
	}
	fclose(fp);
	
	fp=fopen("sample_shuffled/balance_1.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x21[i]=x;
		dx21[i]=0.0;
		y21[i]=integral[1]*y;
		dy21[i]=integral[1]*dy;
	}
	fclose(fp);

	fp=fopen("sample_shuffled/balance_2.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x22[i]=x;
		dx22[i]=0.0;
		y22[i]=integral[2]*y;
		dy22[i]=integral[2]*dy;
	}
	fclose(fp);
	
	fp=fopen("sample_shuffled/balance_3.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x23[i]=x;
		dx23[i]=0.0;
		y23[i]=integral[3]*y;
		dy23[i]=integral[3]*dy;
	}
	fclose(fp);

	fp=fopen("sample_shuffled/balance_4.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x24[i]=x;
		dx24[i]=0.0;
		y24[i]=integral[4]*y;
		dy24[i]=integral[4]*dy;
	}
	fclose(fp);
	
	fp=fopen("sample_shuffled/balance_5.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x25[i]=x;
		dx25[i]=0.0;
		y25[i]=integral[5]*y;
		dy25[i]=integral[5]*dy;
	}
	fclose(fp);

	fp=fopen("sample_shuffled/balance_6.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x26[i]=x;
		dx26[i]=0.0;
		y26[i]=integral[6]*y;
		dy26[i]=integral[6]*dy;
	}
	fclose(fp);
	
	fp=fopen("sample_shuffled/balance_7.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x27[i]=x;
		dx27[i]=0.0;
		y27[i]=integral[7]*y;
		dy27[i]=integral[7]*dy;
	}
	fclose(fp);

	fp=fopen("sample_shuffled/balance_8.txt","r");
	sum=0.0;
	for (i=0;i<26;i++){
		fscanf(fp,"%f %f %f",&x,&y,&dy);
		x28[i]=x;
		dx28[i]=0.0;
		y28[i]=integral[8]*y;
		dy28[i]=integral[8]*dy;
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
		
	TGraph *gr10 = new TGraphErrors(26,x10,y10,dx10,dy10);
	TGraph *gr11 = new TGraphErrors(26,x11,y11,dx11,dy11);
	TGraph *gr12 = new TGraphErrors(26,x12,y12,dx12,dy12);
	TGraph *gr13 = new TGraphErrors(26,x13,y13,dx13,dy13);
	TGraph *gr14 = new TGraphErrors(26,x14,y14,dx14,dy14);
	TGraph *gr15 = new TGraphErrors(26,x15,y15,dx15,dy15);
	TGraph *gr16 = new TGraphErrors(26,x16,y16,dx16,dy16);
	TGraph *gr17 = new TGraphErrors(26,x17,y17,dx17,dy17);
	TGraph *gr18 = new TGraphErrors(26,x18,y18,dx18,dy18);

	TGraph *gr20 = new TGraphErrors(26,x20,y20,dx20,dy20);
	TGraph *gr21 = new TGraphErrors(26,x21,y21,dx21,dy21);
	TGraph *gr22 = new TGraphErrors(26,x22,y22,dx22,dy22);
	TGraph *gr23 = new TGraphErrors(26,x23,y23,dx23,dy23);
	TGraph *gr24 = new TGraphErrors(26,x24,y24,dx24,dy24);
	TGraph *gr25 = new TGraphErrors(26,x25,y25,dx25,dy25);
	TGraph *gr26 = new TGraphErrors(26,x26,y26,dx26,dy26);
	TGraph *gr27 = new TGraphErrors(26,x27,y27,dx27,dy27);
	TGraph *gr28 = new TGraphErrors(26,x28,y28,dx28,dy28);

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
	yl=0.12;
	xh=1.95;
	yh=0.14;
	
	TPaveLabel *p0 = new TPaveLabel(xl-0.5,yl,xh,yh,"kaon 0-5%");
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
	myMax=0.15;
	gr0->SetMinimum(myMin);
	gr0->SetMaximum(myMax);
	gr0->SetMarkerStyle(symbolData);
	gr0->SetMarkerColor(colorData);
	gr0->SetMarkerSize(1.0);
	gr0->SetLineColor(colorData);
	gr0->SetLineWidth(1);
	gr0->SetTitle("");

	gr10->SetMarkerStyle(symbolMixed);
	gr10->SetMarkerColor(colorMixed);
	gr10->SetMarkerSize(1.0);
	gr10->SetLineColor(colorMixed);
	gr10->SetLineWidth(1);
	
	gr20->SetMarkerStyle(symbolShuffled);
	gr20->SetMarkerColor(colorShuffled);
	gr20->SetMarkerSize(1.0);
	gr20->SetLineColor(colorShuffled);
	gr20->SetLineWidth(1);

	gr1->SetMinimum(myMin);
	gr1->SetMaximum(myMax);
	gr1->SetMarkerStyle(symbolData);
	gr1->SetMarkerColor(colorData);
	gr1->SetMarkerSize(1.0);
	gr1->SetLineColor(colorData);
	gr1->SetLineWidth(1);
	gr1->SetTitle("");
	
	gr11->SetMarkerStyle(symbolMixed);
	gr11->SetMarkerColor(colorMixed);
	gr11->SetMarkerSize(1.0);
	gr11->SetLineColor(colorMixed);
	gr11->SetLineWidth(1);

	gr21->SetMarkerStyle(symbolShuffled);
	gr21->SetMarkerColor(colorShuffled);
	gr21->SetMarkerSize(1.0);
	gr21->SetLineColor(colorShuffled);
	gr21->SetLineWidth(1);

	gr2->SetMinimum(myMin);
	gr2->SetMaximum(myMax);
	gr2->SetMarkerStyle(symbolData);
	gr2->SetMarkerColor(colorData);
	gr2->SetMarkerSize(1.0);
	gr2->SetLineColor(colorData);
	gr2->SetLineWidth(1);
	gr2->SetTitle("");
	
	gr12->SetMarkerStyle(symbolMixed);
	gr12->SetMarkerColor(colorMixed);
	gr12->SetMarkerSize(1.0);
	gr12->SetLineColor(colorMixed);
	gr12->SetLineWidth(1);

	gr22->SetMarkerStyle(symbolShuffled);
	gr22->SetMarkerColor(colorShuffled);
	gr22->SetMarkerSize(1.0);
	gr22->SetLineColor(colorShuffled);
	gr22->SetLineWidth(1);

	gr3->SetMinimum(myMin);
	gr3->SetMaximum(myMax);
	gr3->SetMarkerStyle(symbolData);
	gr3->SetMarkerColor(colorData);
	gr3->SetMarkerSize(1.0);
	gr3->SetLineColor(colorData);
	gr3->SetLineWidth(1);
	gr3->SetTitle("");
	
	gr13->SetMarkerStyle(symbolMixed);
	gr13->SetMarkerColor(colorMixed);
	gr13->SetMarkerSize(1.0);
	gr13->SetLineColor(colorMixed);
	gr13->SetLineWidth(1);

	gr23->SetMarkerStyle(symbolShuffled);
	gr23->SetMarkerColor(colorShuffled);
	gr23->SetMarkerSize(1.0);
	gr23->SetLineColor(colorShuffled);
	gr23->SetLineWidth(1);

	gr4->SetMinimum(myMin);
	gr4->SetMaximum(myMax);
	gr4->SetMarkerStyle(symbolData);
	gr4->SetMarkerColor(colorData);
	gr4->SetMarkerSize(1.0);
	gr4->SetLineColor(colorData);
	gr4->SetLineWidth(1);
	gr4->SetTitle("");
	
	gr14->SetMarkerStyle(symbolMixed);
	gr14->SetMarkerColor(colorMixed);
	gr14->SetMarkerSize(1.0);
	gr14->SetLineColor(colorMixed);
	gr14->SetLineWidth(1);

	gr24->SetMarkerStyle(symbolShuffled);
	gr24->SetMarkerColor(colorShuffled);
	gr24->SetMarkerSize(1.0);
	gr24->SetLineColor(colorShuffled);
	gr24->SetLineWidth(1);

	gr5->SetMinimum(myMin);
	gr5->SetMaximum(myMax);
	gr5->SetMarkerStyle(symbolData);
	gr5->SetMarkerColor(colorData);
	gr5->SetMarkerSize(1.0);
	gr5->SetLineColor(colorData);
	gr5->SetLineWidth(1);
	gr5->SetTitle("");
	
	gr15->SetMarkerStyle(symbolMixed);
	gr15->SetMarkerColor(colorMixed);
	gr15->SetMarkerSize(1.0);
	gr15->SetLineColor(colorMixed);
	gr15->SetLineWidth(1);

	gr25->SetMarkerStyle(symbolShuffled);
	gr25->SetMarkerColor(colorShuffled);
	gr25->SetMarkerSize(1.0);
	gr25->SetLineColor(colorShuffled);
	gr25->SetLineWidth(1);

	gr6->SetMinimum(myMin);
	gr6->SetMaximum(myMax);
	gr6->SetMarkerStyle(symbolData);
	gr6->SetMarkerColor(colorData);
	gr6->SetMarkerSize(1.0);
	gr6->SetLineColor(colorData);
	gr6->SetLineWidth(1);
	gr6->SetTitle("");
	
	gr16->SetMarkerStyle(symbolMixed);
	gr16->SetMarkerColor(colorMixed);
	gr16->SetMarkerSize(1.0);
	gr16->SetLineColor(colorMixed);
	gr16->SetLineWidth(1);

	gr26->SetMarkerStyle(symbolShuffled);
	gr26->SetMarkerColor(colorShuffled);
	gr26->SetMarkerSize(1.0);
	gr26->SetLineColor(colorShuffled);
	gr26->SetLineWidth(1);

	gr7->SetMinimum(myMin);
	gr7->SetMaximum(myMax);
	gr7->SetMarkerStyle(symbolData);
	gr7->SetMarkerColor(colorData);
	gr7->SetMarkerSize(1.0);
	gr7->SetLineColor(colorData);
	gr7->SetLineWidth(1);
	gr7->SetTitle("");

	gr17->SetMarkerStyle(symbolMixed);
	gr17->SetMarkerColor(colorMixed);
	gr17->SetMarkerSize(1.0);
	gr17->SetLineColor(colorMixed);
	gr17->SetLineWidth(1);
	
	gr27->SetMarkerStyle(symbolShuffled);
	gr27->SetMarkerColor(colorShuffled);
	gr27->SetMarkerSize(1.0);
	gr27->SetLineColor(colorShuffled);
	gr27->SetLineWidth(1);

	gr8->SetMinimum(myMin);
	gr8->SetMaximum(myMax);
	gr8->SetMarkerStyle(symbolData);
	gr8->SetMarkerColor(colorData);
	gr8->SetMarkerSize(1.0);
	gr8->SetLineColor(colorData);
	gr8->SetLineWidth(1);
	gr8->SetTitle("");
	
	gr18->SetMarkerStyle(symbolMixed);
	gr18->SetMarkerColor(colorMixed);
	gr18->SetMarkerSize(1.0);
	gr18->SetLineColor(colorMixed);
	gr18->SetLineWidth(1);

	gr28->SetMarkerStyle(symbolShuffled);
	gr28->SetMarkerColor(colorShuffled);
	gr28->SetMarkerSize(1.0);
	gr28->SetLineColor(colorShuffled);
	gr28->SetLineWidth(1);

	pad0->cd();
	gr0->Draw("AP");
	gr0->GetXaxis()->SetLimits(0.0,2.0);
	gr0->GetXaxis()->SetTitle("#Deltay");
	gr0->GetYaxis()->SetTitle("B(#Deltay)");
	gr0->GetXaxis()->CenterTitle();
	gr0->GetYaxis()->CenterTitle();
	gr0->Draw("AP");
	line0->Draw();
	gr10->Draw("P");
	gr20->Draw("P");
	gr0->Draw("P");
	p0->Draw();
	pad1->cd();
	gr1->Draw("AP");
	gr1->Draw("AP");
	gr1->GetXaxis()->SetLimits(0.0,2.0);
	gr1->GetXaxis()->SetTitle("#Deltay");
	gr1->GetYaxis()->SetTitle("B(#Deltay)");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	gr1->Draw("AP");
	line1->Draw();
	gr11->Draw("P");
	gr21->Draw("P");
	gr1->Draw("P");
	p1->Draw();
	pad2->cd();
	gr2->Draw("AP");
	gr2->GetXaxis()->SetLimits(0.0,2.0);
	gr2->GetXaxis()->SetTitle("#Deltay");
	gr2->GetYaxis()->SetTitle("B(#Deltay)");
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	gr2->Draw("AP");
	line2->Draw();
	gr12->Draw("P");
	gr22->Draw("P");
	gr2->Draw("P");
	p2->Draw();
	pad3->cd();
	gr3->Draw("AP");
	gr3->GetXaxis()->SetLimits(0.0,2.0);
	gr3->GetXaxis()->SetTitle("#Deltay");
	gr3->GetYaxis()->SetTitle("B(#Deltay)");
	gr3->GetXaxis()->CenterTitle();
	gr3->GetYaxis()->CenterTitle();
	gr3->Draw("AP");
	line3->Draw();
	gr13->Draw("P");
	gr23->Draw("P");
	gr3->Draw("P");
	p3->Draw();
	pad4->cd();
	gr4->Draw("AP");
	gr4->GetXaxis()->SetLimits(0.0,2.0);
	gr4->GetXaxis()->SetTitle("#Deltay");
	gr4->GetYaxis()->SetTitle("B(#Deltay)");
	gr4->GetXaxis()->CenterTitle();
	gr4->GetYaxis()->CenterTitle();
	gr4->Draw("AP");
	line4->Draw();
	gr14->Draw("P");
	gr24->Draw("P");
	gr4->Draw("P");
	p4->Draw();
	pad5->cd();
	gr5->Draw("AP");
	gr5->GetXaxis()->SetLimits(0.0,2.0);
	gr5->GetXaxis()->SetTitle("#Deltay");
	gr5->GetYaxis()->SetTitle("B(#Deltay)");
	gr5->GetXaxis()->CenterTitle();
	gr5->GetYaxis()->CenterTitle();
	gr5->Draw("AP");
	line5->Draw();
	gr15->Draw("P");
	gr25->Draw("P");
	gr5->Draw("P");
	p5->Draw();
	pad6->cd();
	gr6->Draw("AP");
	gr6->GetXaxis()->SetLimits(0.0,2.0);
	gr6->GetXaxis()->SetTitle("#Deltay");
	gr6->GetYaxis()->SetTitle("B(#Deltay)");
	gr6->GetXaxis()->CenterTitle();
	gr6->GetYaxis()->CenterTitle();
	gr6->Draw("AP");
	line6->Draw();
	gr16->Draw("P");
	gr26->Draw("P");
	gr6->Draw("P");
	p6->Draw();
	pad7->cd();
	gr7->Draw("AP");
	gr7->GetXaxis()->SetLimits(0.0,2.0);
	gr7->GetXaxis()->SetTitle("#Deltay");
	gr7->GetYaxis()->SetTitle("B(#Deltay)");
	gr7->GetXaxis()->CenterTitle();
	gr7->GetYaxis()->CenterTitle();
	gr7->Draw("AP");
	line7->Draw();
	gr17->Draw("P");
	gr27->Draw("P");
	gr7->Draw("P");
	p7->Draw();
	pad8->cd();
	gr8->Draw("AP");
	gr8->GetXaxis()->SetLimits(0.0,2.0);
	gr8->GetXaxis()->SetTitle("#Deltay");
	gr8->GetYaxis()->SetTitle("B(#Deltay)");
	gr8->GetXaxis()->CenterTitle();
	gr8->GetYaxis()->CenterTitle();
	gr8->Draw("AP");
	line8->Draw();
	gr18->Draw("P");
	gr28->Draw("P");
	gr8->Draw("P");
	p8->Draw();
	
	c1->Print("balance_kaons_all_sample.eps");
	
}
