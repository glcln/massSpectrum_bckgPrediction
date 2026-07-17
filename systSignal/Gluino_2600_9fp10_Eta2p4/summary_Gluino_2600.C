#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
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
   
   Double_t Graph0_fx64[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy64[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 21.3025, 36.83226, 16.78014, 55.4904, 100, 100, 46.08107, 0,
   0, 0, 0, 16.99353, 20.63361, 0.1066923, 18.3874, 5.906171, 5.618811, 2.797431, 5.909157, 2.103221, 3.594804, 3.995144, 0.6121159, 0.7209778,
   2.136511, 1.791799 };
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
   Double_t Graph1_fy65[35] = { 100, 100, 0, 100, 35.99792, 100, 94.83073, 100, 100, 57.61015, 31.86309, 47.38583, 55.4904, 100, 100, 46.20805, 0.09607077,
   0.5160809, 219.2467, 25.77857, 30.65702, 7.60138, 5.773556, 18.3874, 9.212006, 3.47684, 8.242571, 8.299327, 6.215119, 1.396453, 4.182148, 1.004684, 0.4945755,
   2.312452, 2.347875 };
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
   Graph_Graph165->SetMinimum(88.07533);
   Graph_Graph165->SetMaximum(231.1714);
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
   Double_t Graph2_fy66[35] = { 100, 100, 47.33641, 100, 31.28117, 25.98869, 34.35631, 100, 100, 66.15778, 42.67295, 19.41261, 29.70948, 100, 14.12714, 40.3645, 51.08392,
   24.14898, 25.33882, 17.67333, 18.29686, 26.46173, 26.18648, 13.07487, 8.070898, 15.46028, 7.270825, 2.337843, 10.16535, 3.867728, 0.6226063, 3.192669, 2.288014,
   2.50029, 1.398492 };
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
   Graph_Graph266->SetMinimum(0.5603456);
   Graph_Graph266->SetMaximum(109.9377);
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
   Double_t Graph3_fy67[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1.989448, 2.403808, 2.768791, 1.522255, 1.362276, 0.9653926, 0.9922504,
   1.301265, 1.667821 };
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
   Double_t Graph4_fy68[35] = { 100, 100, 0.9231091, 100, 3.647768, 3.866875, 4.550153, 100, 100, 3.647768, 3.134429, 3.647768, 1.061857, 100, 3.647768, 5.199045, 1.60315,
   3.899753, 3.647768, 2.582586, 3.564632, 3.375453, 3.181958, 2.807069, 2.413988, 3.367698, 2.885044, 3.116655, 3.072083, 3.230822, 3.364897, 3.411508, 3.459406,
   3.512537, 3.521109 };
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
   Graph_Graph468->SetMinimum(0.8307981);
   Graph_Graph468->SetMaximum(109.9077);
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
   Double_t Graph5_fy69[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 3.643978, 0, 0, 0.5878091, 0.9470582, 0.1521409, 0.4908323,
   0.2862513, 0.1496136 };
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
   Double_t Graph6_fy70[35] = { 223.6068, 223.6068, 47.34541, 223.6068, 47.82957, 103.3942, 100.965, 223.6068, 223.6068, 90.34869, 64.82806, 54.0106, 83.91752, 223.6068, 142.172, 76.90885, 51.10916,
   24.46727, 220.7363, 31.36162, 39.70031, 34.57085, 27.00374, 29.24081, 13.8098, 17.14704, 11.87063, 11.16935, 12.7863, 6.525908, 6.856931, 4.914019, 4.353329,
   5.494977, 5.085142 };
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
   Graph_Graph670->SetMinimum(3.917996);
   Graph_Graph670->SetMaximum(245.5321);
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
