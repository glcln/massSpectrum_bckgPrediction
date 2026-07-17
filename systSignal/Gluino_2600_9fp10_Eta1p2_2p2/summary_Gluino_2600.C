#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
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
   
   Double_t Graph0_fx64[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy64[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   150.1004, 100, 33.70017, 41.32624, 19.99944, 50.87987, 19.10881, 12.48807, 7.809454, 22.9681, 14.57825, 1.18196, 12.77881, 7.008672, 4.055655, 0.528276,
   5.564618, 3.73081 };
   TGraph *graph = new TGraph(35,Graph0_fx64,Graph0_fy64);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph064 = new TH1F("Graph_Graph064","Graph",100,0,3520);
   Graph_Graph064->SetMinimum(1);
   Graph_Graph064->SetMaximum(2000);
   Graph_Graph064->SetDirectory(nullptr);
   Graph_Graph064->SetStats(0);
   Graph_Graph064->SetLineWidth(2);
   Graph_Graph064->SetMarkerStyle(20);
   Graph_Graph064->SetMarkerSize(0.9);
   Graph_Graph064->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph064->GetXaxis()->SetRange(1,101);
   Graph_Graph064->GetXaxis()->SetLabelFont(43);
   Graph_Graph064->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetXaxis()->SetLabelSize(16);
   Graph_Graph064->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph064->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph064->GetXaxis()->SetTitleFont(42);
   Graph_Graph064->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph064->GetYaxis()->SetLabelFont(43);
   Graph_Graph064->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetYaxis()->SetLabelSize(16);
   Graph_Graph064->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph064->GetYaxis()->SetTickLength(0.02);
   Graph_Graph064->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph064->GetYaxis()->SetTitleFont(42);
   Graph_Graph064->GetZaxis()->SetLabelFont(42);
   Graph_Graph064->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph064->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph064->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph064->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph064);
   
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
   
   Double_t Graph1_fx65[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy65[35] = { 100, 100, 100, 100, 100, 100, 100, 117.269, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   150.1004, 100, 33.70017, 41.32624, 0, 25.08736, 9.421987, 12.48807, 9.51413, 5.998111, 5.149233, 4.477185, 1.691353, 1.736718, 1.865697, 0.376749,
   1.585072, 1.506281 };
   graph = new TGraph(35,Graph1_fx65,Graph1_fy65);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph165 = new TH1F("Graph_Graph165","Graph",100,0,3520);
   Graph_Graph165->SetMinimum(94.98996);
   Graph_Graph165->SetMaximum(155.1105);
   Graph_Graph165->SetDirectory(nullptr);
   Graph_Graph165->SetStats(0);
   Graph_Graph165->SetLineWidth(2);
   Graph_Graph165->SetMarkerStyle(20);
   Graph_Graph165->SetMarkerSize(0.9);
   Graph_Graph165->GetXaxis()->SetLabelFont(42);
   Graph_Graph165->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetXaxis()->SetTitleFont(42);
   Graph_Graph165->GetYaxis()->SetLabelFont(42);
   Graph_Graph165->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetYaxis()->SetTickLength(0.02);
   Graph_Graph165->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetYaxis()->SetTitleFont(42);
   Graph_Graph165->GetZaxis()->SetLabelFont(42);
   Graph_Graph165->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph165);
   
   graph->Draw("p");
   
   Double_t Graph2_fx66[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy66[35] = { 100, 100, 100, 100, 100, 100, 253.2872, 333.206, 100, 100, 100, 100, 100, 100, 100, 100, 111.7553,
   53.79027, 44.45321, 29.17627, 38.76583, 39.99096, 76.84357, 31.28685, 27.46707, 17.73824, 54.61277, 21.24412, 62.31522, 16.79659, 27.20771, 25.43852, 19.32454,
   25.18544, 14.56981 };
   graph = new TGraph(35,Graph2_fx66,Graph2_fy66);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph266 = new TH1F("Graph_Graph266","Graph",100,0,3520);
   Graph_Graph266->SetMinimum(13.11283);
   Graph_Graph266->SetMaximum(365.0696);
   Graph_Graph266->SetDirectory(nullptr);
   Graph_Graph266->SetStats(0);
   Graph_Graph266->SetLineWidth(2);
   Graph_Graph266->SetMarkerStyle(20);
   Graph_Graph266->SetMarkerSize(0.9);
   Graph_Graph266->GetXaxis()->SetLabelFont(42);
   Graph_Graph266->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetXaxis()->SetTitleFont(42);
   Graph_Graph266->GetYaxis()->SetLabelFont(42);
   Graph_Graph266->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetYaxis()->SetTickLength(0.02);
   Graph_Graph266->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetYaxis()->SetTitleFont(42);
   Graph_Graph266->GetZaxis()->SetLabelFont(42);
   Graph_Graph266->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph266);
   
   graph->Draw("p");
   
   Double_t Graph3_fx67[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy67[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx67,Graph3_fy67);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph367 = new TH1F("Graph_Graph367","Graph",100,0,3520);
   Graph_Graph367->SetMinimum(99.9);
   Graph_Graph367->SetMaximum(101.1);
   Graph_Graph367->SetDirectory(nullptr);
   Graph_Graph367->SetStats(0);
   Graph_Graph367->SetLineWidth(2);
   Graph_Graph367->SetMarkerStyle(20);
   Graph_Graph367->SetMarkerSize(0.9);
   Graph_Graph367->GetXaxis()->SetLabelFont(42);
   Graph_Graph367->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetXaxis()->SetTitleFont(42);
   Graph_Graph367->GetYaxis()->SetLabelFont(42);
   Graph_Graph367->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetYaxis()->SetTickLength(0.02);
   Graph_Graph367->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetYaxis()->SetTitleFont(42);
   Graph_Graph367->GetZaxis()->SetLabelFont(42);
   Graph_Graph367->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph367);
   
   graph->Draw("p");
   
   Double_t Graph4_fx68[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy68[35] = { 100, 100, 100, 100, 100, 100, 4.550165, 4.550159, 100, 100, 100, 100, 100, 100, 100, 100, 1.60315,
   2.507889, 3.111076, 3.503931, 3.824914, 3.74648, 2.271569, 2.955568, 2.518785, 3.149438, 2.68873, 3.023463, 3.391111, 3.103769, 3.349578, 3.234613, 3.494143,
   3.491759, 3.508627 };
   graph = new TGraph(35,Graph4_fx68,Graph4_fy68);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph468 = new TH1F("Graph_Graph468","Graph",100,0,3520);
   Graph_Graph468->SetMinimum(1.442835);
   Graph_Graph468->SetMaximum(109.8397);
   Graph_Graph468->SetDirectory(nullptr);
   Graph_Graph468->SetStats(0);
   Graph_Graph468->SetLineWidth(2);
   Graph_Graph468->SetMarkerStyle(20);
   Graph_Graph468->SetMarkerSize(0.9);
   Graph_Graph468->GetXaxis()->SetLabelFont(42);
   Graph_Graph468->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetXaxis()->SetTitleFont(42);
   Graph_Graph468->GetYaxis()->SetLabelFont(42);
   Graph_Graph468->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetYaxis()->SetTickLength(0.02);
   Graph_Graph468->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetYaxis()->SetTitleFont(42);
   Graph_Graph468->GetZaxis()->SetLabelFont(42);
   Graph_Graph468->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph468);
   
   graph->Draw("p");
   
   Double_t Graph5_fx69[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy69[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.033944, 0.7043719,
   0.2831936, 0.4048884 };
   graph = new TGraph(35,Graph5_fx69,Graph5_fy69);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph569 = new TH1F("Graph_Graph569","Graph",100,0,3520);
   Graph_Graph569->SetMinimum(99.9);
   Graph_Graph569->SetMaximum(101.1);
   Graph_Graph569->SetDirectory(nullptr);
   Graph_Graph569->SetStats(0);
   Graph_Graph569->SetLineWidth(2);
   Graph_Graph569->SetMarkerStyle(20);
   Graph_Graph569->SetMarkerSize(0.9);
   Graph_Graph569->GetXaxis()->SetLabelFont(42);
   Graph_Graph569->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetXaxis()->SetTitleFont(42);
   Graph_Graph569->GetYaxis()->SetLabelFont(42);
   Graph_Graph569->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetYaxis()->SetTickLength(0.02);
   Graph_Graph569->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetYaxis()->SetTitleFont(42);
   Graph_Graph569->GetZaxis()->SetLabelFont(42);
   Graph_Graph569->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph569);
   
   graph->Draw("p");
   
   Double_t Graph6_fx70[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy70[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 272.3511, 353.269, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 111.7668,
   218.9976, 148.276, 55.99049, 70.23628, 44.86971, 95.54178, 37.96739, 32.75192, 21.81904, 59.60952, 26.44795, 62.57898, 21.399, 28.34813, 26.02903, 19.64861,
   26.07635, 15.51701 };
   graph = new TGraph(35,Graph6_fx70,Graph6_fy70);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph670 = new TH1F("Graph_Graph670","Graph",100,0,3520);
   Graph_Graph670->SetMinimum(13.96531);
   Graph_Graph670->SetMaximum(387.0442);
   Graph_Graph670->SetDirectory(nullptr);
   Graph_Graph670->SetStats(0);
   Graph_Graph670->SetLineWidth(2);
   Graph_Graph670->SetMarkerStyle(20);
   Graph_Graph670->SetMarkerSize(0.9);
   Graph_Graph670->GetXaxis()->SetLabelFont(42);
   Graph_Graph670->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetXaxis()->SetTitleFont(42);
   Graph_Graph670->GetYaxis()->SetLabelFont(42);
   Graph_Graph670->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetYaxis()->SetTickLength(0.02);
   Graph_Graph670->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetYaxis()->SetTitleFont(42);
   Graph_Graph670->GetZaxis()->SetLabelFont(42);
   Graph_Graph670->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph670);
   
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
