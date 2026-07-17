#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:54:12 2026) by ROOT version 6.32.13
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
   69.58598, 100, 33.70017, 41.32624, 0, 50.87987, 19.10881, 20.59215, 12.49506, 21.65165, 9.199572, 3.392792, 13.48184, 5.020147, 4.939633, 1.044106,
   4.867697, 3.669608 };
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
   69.58598, 100, 33.70017, 41.32624, 0, 25.08736, 9.421987, 14.88259, 9.030557, 5.654329, 4.706335, 3.756881, 2.615643, 1.241243, 1.778114, 0.3416896,
   0.9781957, 2.310097 };
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
   Graph_Graph165->SetMinimum(98.2731);
   Graph_Graph165->SetMaximum(118.9959);
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
   54.41906, 44.45321, 29.17627, 38.76583, 35.18309, 76.84357, 31.28685, 37.60712, 14.64928, 65.30868, 26.93002, 51.705, 15.93082, 21.51392, 28.55684, 20.24896,
   24.19811, 15.19542 };
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
   Graph_Graph266->SetMinimum(13.18435);
   Graph_Graph266->SetMaximum(365.0617);
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
   4.659158, 3.111076, 3.503931, 3.824914, 4.05612, 2.271569, 2.955568, 2.405214, 3.17477, 2.64715, 3.017533, 3.264022, 3.192437, 3.360379, 3.219306, 3.460348,
   3.469872, 3.516769 };
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
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.8942187, 0.5827904,
   0.2331376, 0.309056 };
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
   112.5502, 148.276, 55.99049, 70.23628, 35.41613, 95.54178, 37.96739, 45.44897, 21.50251, 69.08688, 29.00196, 52.05465, 21.27403, 22.38042, 29.21333, 20.57186,
   24.94474, 16.18861 };
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
   Graph_Graph670->SetMinimum(14.56975);
   Graph_Graph670->SetMaximum(386.977);
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
