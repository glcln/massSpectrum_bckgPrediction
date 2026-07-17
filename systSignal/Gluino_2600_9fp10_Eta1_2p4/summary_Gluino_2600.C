#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:15:05 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx64[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy64[35] = { 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   1.192093e-05, 100, 0, 0, 0, 0, 0, 1.192093e-05, 14.01694, 17.79807, 15.37266, 3.879344, 6.408387, 2.823257, 1.578939, 0.9381175,
   2.754408, 1.235175 };
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
   Graph_Graph064->GetXaxis()->SetLabelSize(22);
   Graph_Graph064->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph064->GetXaxis()->SetTitleOffset(1);
   Graph_Graph064->GetXaxis()->SetTitleFont(42);
   Graph_Graph064->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph064->GetYaxis()->SetLabelFont(43);
   Graph_Graph064->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetYaxis()->SetLabelSize(22);
   Graph_Graph064->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph064->GetYaxis()->SetTickLength(0.02);
   Graph_Graph064->GetYaxis()->SetTitleOffset(1);
   Graph_Graph064->GetYaxis()->SetTitleFont(42);
   Graph_Graph064->GetZaxis()->SetLabelFont(42);
   Graph_Graph064->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph064->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph064->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph064->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph064);
   
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
   
   Double_t Graph1_fx65[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy65[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0.09607077,
   31.63019, 100, 25.20577, 100, 36.42242, 17.17554, 11.97512, 17.07587, 24.71802, 17.79807, 15.37266, 7.33999, 7.011956, 4.111135, 2.318692, 0.7572412,
   4.005199, 2.484238 };
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
   Graph_Graph165->SetMinimum(0.08646369);
   Graph_Graph165->SetMaximum(109.9904);
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
   Double_t Graph2_fy66[35] = { 100, 100, 100, 100, 100, 27.74895, 33.25831, 100, 100, 100, 100, 100, 100, 100, 100, 100, 51.08392,
   23.84821, 100, 17.3983, 20.74047, 29.95122, 33.17901, 14.53331, 18.33145, 33.28388, 1.250219, 19.22644, 16.45044, 3.352749, 4.944533, 6.831896, 0.9085417,
   5.83756, 1.88123 };
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
   Graph_Graph266->SetMinimum(0.8176875);
   Graph_Graph266->SetMaximum(109.9091);
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
   Double_t Graph3_fy67[35] = { 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   0, 100, 0, 0, 0, 0, 0, 0, 0, 8.89411, 3.756821, 0, 2.529883, 2.033341, 2.376652, 1.203048,
   1.751238, 3.373682 };
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
   Double_t Graph4_fy68[35] = { 100, 100, 100, 100, 100, 4.550165, 4.550159, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1.60315,
   4.023891, 100, 3.024828, 5.389786, 4.05612, 2.282321, 3.109622, 2.372968, 3.31161, 2.62233, 3.252065, 3.163505, 3.148651, 3.358233, 3.234804, 3.469026,
   3.458095, 3.520203 };
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
   Double_t Graph5_fy69[35] = { 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0,
   0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.394398, 0.6673396, 0.584352,
   0.2433836, 0.323981 };
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
   Double_t Graph6_fy70[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 103.8783, 105.4837, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 51.10916,
   39.81706, 223.6068, 30.77631, 102.2703, 47.32991, 37.43066, 19.08639, 25.16462, 43.88894, 26.85308, 29.44457, 18.69625, 10.85312, 8.045671, 8.405749, 3.969936,
   8.528233, 5.916898 };
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
   Graph_Graph670->SetMinimum(3.572943);
   Graph_Graph670->SetMaximum(245.5705);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2600 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
