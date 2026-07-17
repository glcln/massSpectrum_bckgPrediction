#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2400()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:10:02 2026) by ROOT version 6.32.13
   TCanvas *c2 = new TCanvas("c2", "c2",0,0,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c2->SetHighLightColor(2);
   c2->Range(-418.2588,-0.3883565,3764.329,3.495208);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetLogy();
   c2->SetGridx();
   c2->SetGridy();
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);
   c2->SetFrameLineWidth(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameLineWidth(2);
   c2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx57[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy57[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 1.192093e-05, 0, 100, 0,
   90.97992, 34.89867, 6.903153, 0, 0, 5.222249, 4.829687, 3.619874, 7.212776, 3.406596, 1.303089, 5.12268, 3.894061, 2.391088, 0.5481958, 1.138413,
   2.152216, 0.9474635 };
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
   Graph_Graph057->GetXaxis()->SetLabelSize(22);
   Graph_Graph057->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph057->GetXaxis()->SetTitleOffset(1);
   Graph_Graph057->GetXaxis()->SetTitleFont(42);
   Graph_Graph057->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph057->GetYaxis()->SetLabelFont(43);
   Graph_Graph057->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph057->GetYaxis()->SetLabelSize(22);
   Graph_Graph057->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph057->GetYaxis()->SetTickLength(0.02);
   Graph_Graph057->GetYaxis()->SetTitleOffset(1);
   Graph_Graph057->GetYaxis()->SetTitleFont(42);
   Graph_Graph057->GetZaxis()->SetLabelFont(42);
   Graph_Graph057->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph057->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph057->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph057->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph057->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph057);
   
   graph->Draw("ap");
   
   TLegend *leg = new TLegend(0.12,0.75,0.5,0.93,NULL,"brNDC");
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
   Double_t Graph1_fy58[35] = { 100, 100, 100, 100, 100, 100, 100, 150.0508, 0, 100, 12.52763, 0, 266.8815, 67.90086, 51.24465, 100, 9.706473,
   72.29256, 26.55518, 8.08816, 29.13729, 11.69013, 0, 5.003243, 3.619874, 1.674658, 7.314402, 5.790269, 4.324746, 3.924072, 2.570879, 1.134646, 1.250041,
   2.699816, 1.360822 };
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
   Graph_Graph158->SetMinimum(83.31185);
   Graph_Graph158->SetMaximum(283.5696);
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
   Double_t Graph2_fy59[35] = { 100, 100, 100, 100, 51.08391, 100, 54.53767, 14.24768, 26.91302, 26.50329, 36.25692, 54.53767, 17.32603, 51.741, 30.54602, 100, 8.998895,
   9.911501, 35.25891, 13.77539, 12.25876, 15.07874, 6.410897, 8.499706, 8.69295, 16.6363, 2.471298, 12.12232, 8.784628, 1.170051, 3.824496, 2.090573, 1.139665,
   1.107877, 1.34635 };
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
   Graph_Graph259->SetMinimum(0.9970897);
   Graph_Graph259->SetMaximum(109.8892);
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
   Double_t Graph3_fy60[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0,
   0, 0, 6.903153, 0, 0, 0, 0, 0, 0, 0, 0, 1.352096, 1.357722, 0.8052051, 1.080644, 1.209128,
   1.543343, 1.489186 };
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
   Graph_Graph360->SetMinimum(99.9);
   Graph_Graph360->SetMaximum(101.1);
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
   Double_t Graph4_fy61[35] = { 100, 100, 100, 100, 1.234829, 100, 1.234829, 3.086007, 3.647768, 6.51843, 1.60135, 2.507883, 2.954483, 3.025901, 2.319336, 100, 2.407432,
   2.826846, 1.858914, 2.363837, 2.591968, 2.648139, 2.585244, 2.728426, 3.374922, 2.767694, 3.133917, 2.894354, 3.300071, 3.232229, 3.356528, 3.394115, 3.48376,
   3.466523, 3.546846 };
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
   Graph_Graph461->SetMinimum(1.111346);
   Graph_Graph461->SetMaximum(109.8765);
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
   Double_t Graph5_fy62[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0,
   0, 0, 0, 0, 0, 7.914138, 1.755285, 0, 9.144651, 1.147348, 0, 0, 0, 0.3943563, 0.2723277, 0.1868427,
   0.2573788, 0 };
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
   Graph_Graph562->SetMinimum(99.9);
   Graph_Graph562->SetMaximum(101.1);
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
   Double_t Graph6_fy63[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 112.2991, 223.6068, 113.9117, 150.7573, 27.15911, 103.6577, 38.39362, 54.5953, 267.4596, 85.42139, 59.70304, 223.6068, 13.45331,
   116.661, 56.30045, 18.86994, 31.71714, 19.2624, 8.663427, 11.31581, 10.63788, 18.41888, 9.001895, 13.8041, 11.61179, 6.649946, 6.234431, 4.31813, 4.213787,
   5.248535, 4.399997 };
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
   Graph_Graph663->SetMinimum(3.792408);
   Graph_Graph663->SetMaximum(293.7842);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2400 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
