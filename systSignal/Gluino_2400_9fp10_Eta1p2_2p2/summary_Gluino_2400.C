#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2400()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 15:56:22 2026) by ROOT version 6.32.13
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(-720.0405,-0.668563,3780.213,3.509956);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx57[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy57[35] = { 100, 0, 100, 100, 100, 100, 100, 24.72208, 100, 100, 100, 100, 100, 0, 100, 100, 0,
   100, 40.44222, 56.50905, 0, 0, 67.69238, 41.82457, 7.06349, 5.361211, 13.27156, 1.510394, 5.546594, 7.38883, 8.184904, 1.58084, 3.331888,
   4.88081, 8.713543 };
   TGraph *graph = new TGraph(35,Graph0_fx57,Graph0_fy57);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph057 = new TH1F("Graph_Graph057","Graph",100,0,3520);
   Graph_Graph057->SetMinimum(1);
   Graph_Graph057->SetMaximum(2000);
   Graph_Graph057->SetDirectory(nullptr);
   Graph_Graph057->SetStats(0);
   Graph_Graph057->SetLineWidth(2);
   Graph_Graph057->SetMarkerStyle(20);
   Graph_Graph057->SetMarkerSize(0.9);
   Graph_Graph057->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph057->GetXaxis()->SetRange(1,101);
   Graph_Graph057->GetXaxis()->SetLabelFont(43);
   Graph_Graph057->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph057->GetXaxis()->SetLabelSize(16);
   Graph_Graph057->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph057->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph057->GetXaxis()->SetTitleFont(42);
   Graph_Graph057->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph057->GetYaxis()->SetLabelFont(43);
   Graph_Graph057->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph057->GetYaxis()->SetLabelSize(16);
   Graph_Graph057->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph057->GetYaxis()->SetTickLength(0.02);
   Graph_Graph057->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph057->GetYaxis()->SetTitleFont(42);
   Graph_Graph057->GetZaxis()->SetLabelFont(42);
   Graph_Graph057->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph057->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph057->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph057->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph057->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph057);
   
   graph->Draw("ap");
   
   TLegend *leg = new TLegend(0.27,0.7,0.5,0.93,NULL,"brNDC");
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","Total","PE1");
   entry->SetLineColor(28);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(28);
   entry->SetMarkerStyle(34);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph0","K","PE1");
   entry->SetLineColor(30);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(30);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","C","PE1");
   entry->SetLineColor(38);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(38);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","PU","PE1");
   entry->SetLineColor(46);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(46);
   entry->SetMarkerStyle(23);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","F^{pixel}","PE1");
   entry->SetLineColor(43);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(43);
   entry->SetMarkerStyle(43);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","Trigger","PE1");
   entry->SetLineColor(40);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(40);
   entry->SetMarkerStyle(39);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   entry=leg->AddEntry("Graph","Jet","PE1");
   entry->SetLineColor(41);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(41);
   entry->SetMarkerStyle(42);
   entry->SetMarkerSize(0.9);
   entry->SetTextFont(62);
   leg->Draw();
   
   Double_t Graph1_fx58[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy58[35] = { 100, 100, 100, 100, 100, 100, 100, 75.27792, 100, 100, 100, 100, 100, 0, 100, 100, 0,
   100, 40.44222, 56.50905, 193.3337, 18.1387, 67.69238, 41.82457, 7.063496, 7.994294, 12.34069, 0.3498316, 3.170323, 3.334844, 2.987695, 0.8113265, 1.665097,
   1.650739, 2.480054 };
   graph = new TGraph(35,Graph1_fx58,Graph1_fy58);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph158 = new TH1F("Graph_Graph158","Graph",100,0,3520);
   Graph_Graph158->SetMinimum(63.47234);
   Graph_Graph158->SetMaximum(205.1393);
   Graph_Graph158->SetDirectory(nullptr);
   Graph_Graph158->SetStats(0);
   Graph_Graph158->SetLineWidth(2);
   Graph_Graph158->SetMarkerStyle(20);
   Graph_Graph158->SetMarkerSize(0.9);
   Graph_Graph158->GetXaxis()->SetLabelFont(42);
   Graph_Graph158->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph158->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph158->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph158->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph158->GetXaxis()->SetTitleFont(42);
   Graph_Graph158->GetYaxis()->SetLabelFont(42);
   Graph_Graph158->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph158->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph158->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph158->GetYaxis()->SetTickLength(0.02);
   Graph_Graph158->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph158->GetYaxis()->SetTitleFont(42);
   Graph_Graph158->GetZaxis()->SetLabelFont(42);
   Graph_Graph158->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph158->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph158->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph158->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph158->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph158);
   
   graph->Draw("p");
   
   Double_t Graph2_fx59[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy59[35] = { 100, 105.7607, 100, 100, 100, 100, 100, 48.74796, 100, 100, 100, 100, 100, 93.1768, 100, 100, 54.13097,
   100, 58.0814, 105.6848, 54.28553, 14.3967, 23.73161, 27.80219, 76.43472, 42.38744, 74.87427, 22.03599, 14.92608, 29.19016, 26.54098, 22.7435, 20.61067,
   17.22265, 23.83257 };
   graph = new TGraph(35,Graph2_fx59,Graph2_fy59);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph259 = new TH1F("Graph_Graph259","Graph",100,0,3520);
   Graph_Graph259->SetMinimum(5.260308);
   Graph_Graph259->SetMaximum(114.8971);
   Graph_Graph259->SetDirectory(nullptr);
   Graph_Graph259->SetStats(0);
   Graph_Graph259->SetLineWidth(2);
   Graph_Graph259->SetMarkerStyle(20);
   Graph_Graph259->SetMarkerSize(0.9);
   Graph_Graph259->GetXaxis()->SetLabelFont(42);
   Graph_Graph259->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph259->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph259->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph259->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph259->GetXaxis()->SetTitleFont(42);
   Graph_Graph259->GetYaxis()->SetLabelFont(42);
   Graph_Graph259->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph259->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph259->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph259->GetYaxis()->SetTickLength(0.02);
   Graph_Graph259->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph259->GetYaxis()->SetTitleFont(42);
   Graph_Graph259->GetZaxis()->SetLabelFont(42);
   Graph_Graph259->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph259->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph259->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph259->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph259->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph259);
   
   graph->Draw("p");
   
   Double_t Graph3_fx60[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy60[35] = { 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 0,
   100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx60,Graph3_fy60);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph360 = new TH1F("Graph_Graph360","Graph",100,0,3520);
   Graph_Graph360->SetMinimum(0.11);
   Graph_Graph360->SetMaximum(110);
   Graph_Graph360->SetDirectory(nullptr);
   Graph_Graph360->SetStats(0);
   Graph_Graph360->SetLineWidth(2);
   Graph_Graph360->SetMarkerStyle(20);
   Graph_Graph360->SetMarkerSize(0.9);
   Graph_Graph360->GetXaxis()->SetLabelFont(42);
   Graph_Graph360->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph360->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph360->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph360->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph360->GetXaxis()->SetTitleFont(42);
   Graph_Graph360->GetYaxis()->SetLabelFont(42);
   Graph_Graph360->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph360->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph360->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph360->GetYaxis()->SetTickLength(0.02);
   Graph_Graph360->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph360->GetYaxis()->SetTitleFont(42);
   Graph_Graph360->GetZaxis()->SetLabelFont(42);
   Graph_Graph360->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph360->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph360->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph360->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph360->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph360);
   
   graph->Draw("p");
   
   Double_t Graph4_fx61[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy61[35] = { 100, 2.507889, 100, 100, 100, 100, 100, 1.549518, 100, 100, 100, 100, 100, 2.507889, 100, 100, 3.647757,
   100, 2.860415, 1.819444, 3.647757, 2.32991, 1.391876, 3.21871, 2.499008, 2.257979, 2.909386, 3.156853, 3.210509, 2.940881, 3.265929, 3.374028, 3.382242,
   3.509986, 3.469431 };
   graph = new TGraph(35,Graph4_fx61,Graph4_fy61);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph461 = new TH1F("Graph_Graph461","Graph",100,0,3520);
   Graph_Graph461->SetMinimum(1.252688);
   Graph_Graph461->SetMaximum(109.8608);
   Graph_Graph461->SetDirectory(nullptr);
   Graph_Graph461->SetStats(0);
   Graph_Graph461->SetLineWidth(2);
   Graph_Graph461->SetMarkerStyle(20);
   Graph_Graph461->SetMarkerSize(0.9);
   Graph_Graph461->GetXaxis()->SetLabelFont(42);
   Graph_Graph461->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph461->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph461->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph461->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph461->GetXaxis()->SetTitleFont(42);
   Graph_Graph461->GetYaxis()->SetLabelFont(42);
   Graph_Graph461->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph461->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph461->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph461->GetYaxis()->SetTickLength(0.02);
   Graph_Graph461->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph461->GetYaxis()->SetTitleFont(42);
   Graph_Graph461->GetZaxis()->SetLabelFont(42);
   Graph_Graph461->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph461->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph461->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph461->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph461->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph461);
   
   graph->Draw("p");
   
   Double_t Graph5_fx62[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy62[35] = { 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 0,
   100, 0, 0, 0, 0, 36.7226, 0, 0, 11.12782, 13.36929, 0, 0, 0, 0.5046427, 0.255692, 0.2676785,
   0.3783047, 0.8916378 };
   graph = new TGraph(35,Graph5_fx62,Graph5_fy62);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph562 = new TH1F("Graph_Graph562","Graph",100,0,3520);
   Graph_Graph562->SetMinimum(0.11);
   Graph_Graph562->SetMaximum(110);
   Graph_Graph562->SetDirectory(nullptr);
   Graph_Graph562->SetStats(0);
   Graph_Graph562->SetLineWidth(2);
   Graph_Graph562->SetMarkerStyle(20);
   Graph_Graph562->SetMarkerSize(0.9);
   Graph_Graph562->GetXaxis()->SetLabelFont(42);
   Graph_Graph562->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph562->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph562->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph562->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph562->GetXaxis()->SetTitleFont(42);
   Graph_Graph562->GetYaxis()->SetLabelFont(42);
   Graph_Graph562->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph562->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph562->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph562->GetYaxis()->SetTickLength(0.02);
   Graph_Graph562->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph562->GetYaxis()->SetTitleFont(42);
   Graph_Graph562->GetZaxis()->SetLabelFont(42);
   Graph_Graph562->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph562->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph562->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph562->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph562->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph562);
   
   graph->Draw("p");
   
   Double_t Graph6_fx63[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy63[35] = { 223.6068, 145.5734, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 93.04145, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 93.21055, 223.6068, 223.6068, 54.25374,
   223.6068, 81.56456, 132.5108, 200.8436, 23.27458, 98.63895, 65.43631, 77.12521, 43.52522, 77.09116, 22.31489, 16.55025, 30.43731, 28.12488, 23.06097, 21.21587,
   18.31631, 25.7314 };
   graph = new TGraph(35,Graph6_fx63,Graph6_fy63);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph663 = new TH1F("Graph_Graph663","Graph",100,0,3520);
   Graph_Graph663->SetMinimum(14.89523);
   Graph_Graph663->SetMaximum(244.3125);
   Graph_Graph663->SetDirectory(nullptr);
   Graph_Graph663->SetStats(0);
   Graph_Graph663->SetLineWidth(2);
   Graph_Graph663->SetMarkerStyle(20);
   Graph_Graph663->SetMarkerSize(0.9);
   Graph_Graph663->GetXaxis()->SetLabelFont(42);
   Graph_Graph663->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph663->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph663->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph663->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph663->GetXaxis()->SetTitleFont(42);
   Graph_Graph663->GetYaxis()->SetLabelFont(42);
   Graph_Graph663->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph663->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph663->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph663->GetYaxis()->SetTickLength(0.02);
   Graph_Graph663->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph663->GetYaxis()->SetTitleFont(42);
   Graph_Graph663->GetZaxis()->SetLabelFont(42);
   Graph_Graph663->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph663->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph663->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph663->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph663->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph663);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.16,0.96,"#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->SetSelected(c1);
}
