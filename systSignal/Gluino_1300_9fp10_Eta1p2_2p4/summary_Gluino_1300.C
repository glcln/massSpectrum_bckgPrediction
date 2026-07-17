#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1300()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 16:05:59 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx15[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy15[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 100, 49.43507, 39.40901, 100, 0, 0, 58.90679, 28.38901, 16.77189,
   40.22198, 36.47787, 24.10053, 43.9172, 22.10586, 12.55963, 22.39437, 8.351588, 9.604764, 8.013659, 8.816492, 1.213324, 6.000119, 4.421282, 11.02319, 4.185963,
   7.830966, 1.192093e-05 };
   TGraph *graph = new TGraph(35,Graph0_fx15,Graph0_fy15);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph015 = new TH1F("Graph_Graph015","Graph",100,0,3520);
   Graph_Graph015->SetMinimum(1);
   Graph_Graph015->SetMaximum(2000);
   Graph_Graph015->SetDirectory(nullptr);
   Graph_Graph015->SetStats(0);
   Graph_Graph015->SetLineWidth(2);
   Graph_Graph015->SetMarkerStyle(20);
   Graph_Graph015->SetMarkerSize(0.9);
   Graph_Graph015->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph015->GetXaxis()->SetRange(1,101);
   Graph_Graph015->GetXaxis()->SetLabelFont(43);
   Graph_Graph015->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetXaxis()->SetLabelSize(16);
   Graph_Graph015->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph015->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph015->GetXaxis()->SetTitleFont(42);
   Graph_Graph015->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph015->GetYaxis()->SetLabelFont(43);
   Graph_Graph015->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetYaxis()->SetLabelSize(16);
   Graph_Graph015->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph015->GetYaxis()->SetTickLength(0.02);
   Graph_Graph015->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph015->GetYaxis()->SetTitleFont(42);
   Graph_Graph015->GetZaxis()->SetLabelFont(42);
   Graph_Graph015->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph015->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph015->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph015->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph015);
   
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
   
   Double_t Graph1_fx16[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy16[35] = { 100, 100, 0, 100, 239.6609, 100, 100, 0, 100, 2.739704, 39.40901, 100, 0, 59.4857, 0, 19.73249, 16.44656,
   17.15362, 15.55685, 36.59022, 43.9172, 0, 13.45115, 20.50989, 11.40249, 4.496896, 5.42829, 5.093014, 0.9733319, 2.633417, 1.578546, 10.26736, 2.034366,
   4.373264, 1.192093e-05 };
   graph = new TGraph(35,Graph1_fx16,Graph1_fy16);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph116 = new TH1F("Graph_Graph116","Graph",100,0,3520);
   Graph_Graph116->SetMinimum(86.03391);
   Graph_Graph116->SetMaximum(253.627);
   Graph_Graph116->SetDirectory(nullptr);
   Graph_Graph116->SetStats(0);
   Graph_Graph116->SetLineWidth(2);
   Graph_Graph116->SetMarkerStyle(20);
   Graph_Graph116->SetMarkerSize(0.9);
   Graph_Graph116->GetXaxis()->SetLabelFont(42);
   Graph_Graph116->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetXaxis()->SetTitleFont(42);
   Graph_Graph116->GetYaxis()->SetLabelFont(42);
   Graph_Graph116->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetYaxis()->SetTickLength(0.02);
   Graph_Graph116->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetYaxis()->SetTitleFont(42);
   Graph_Graph116->GetZaxis()->SetLabelFont(42);
   Graph_Graph116->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph116);
   
   graph->Draw("p");
   
   Double_t Graph2_fx17[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy17[35] = { 100, 100, 95.056, 100, 341.3144, 100, 100, 46.15802, 70.84902, 93.50281, 22.12948, 100, 62.18934, 37.95864, 50.01815, 28.6711, 104.4808,
   31.43879, 60.05592, 10.73374, 42.09863, 59.21656, 52.14, 32.24055, 39.5532, 30.15809, 21.3887, 19.38431, 24.41587, 20.35459, 23.06091, 25.08667, 33.17645,
   10.84372, 20.70045 };
   graph = new TGraph(35,Graph2_fx17,Graph2_fy17);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph217 = new TH1F("Graph_Graph217","Graph",100,0,3520);
   Graph_Graph217->SetMinimum(9.660362);
   Graph_Graph217->SetMaximum(374.3724);
   Graph_Graph217->SetDirectory(nullptr);
   Graph_Graph217->SetStats(0);
   Graph_Graph217->SetLineWidth(2);
   Graph_Graph217->SetMarkerStyle(20);
   Graph_Graph217->SetMarkerSize(0.9);
   Graph_Graph217->GetXaxis()->SetLabelFont(42);
   Graph_Graph217->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetXaxis()->SetTitleFont(42);
   Graph_Graph217->GetYaxis()->SetLabelFont(42);
   Graph_Graph217->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetYaxis()->SetTickLength(0.02);
   Graph_Graph217->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetYaxis()->SetTitleFont(42);
   Graph_Graph217->GetZaxis()->SetLabelFont(42);
   Graph_Graph217->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph217);
   
   graph->Draw("p");
   
   Double_t Graph3_fx18[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy18[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx18,Graph3_fy18);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph318 = new TH1F("Graph_Graph318","Graph",100,0,3520);
   Graph_Graph318->SetMinimum(99.9);
   Graph_Graph318->SetMaximum(101.1);
   Graph_Graph318->SetDirectory(nullptr);
   Graph_Graph318->SetStats(0);
   Graph_Graph318->SetLineWidth(2);
   Graph_Graph318->SetMarkerStyle(20);
   Graph_Graph318->SetMarkerSize(0.9);
   Graph_Graph318->GetXaxis()->SetLabelFont(42);
   Graph_Graph318->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetXaxis()->SetTitleFont(42);
   Graph_Graph318->GetYaxis()->SetLabelFont(42);
   Graph_Graph318->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetYaxis()->SetTickLength(0.02);
   Graph_Graph318->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetYaxis()->SetTitleFont(42);
   Graph_Graph318->GetZaxis()->SetLabelFont(42);
   Graph_Graph318->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph318);
   
   graph->Draw("p");
   
   Double_t Graph4_fx19[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy19[35] = { 100, 100, 3.722632, 100, 3.647768, 100, 100, 1.60315, 0.9231031, 5.064124, 2.54131, 100, 3.120363, 1.234829, 3.682387, 3.82185, 2.471733,
   2.480364, 2.032423, 2.367902, 2.301538, 2.642179, 3.063571, 1.804096, 3.149343, 3.029573, 2.817023, 3.189898, 3.280115, 3.346062, 3.263509, 3.356004, 3.068817,
   3.180349, 2.977443 };
   graph = new TGraph(35,Graph4_fx19,Graph4_fy19);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph419 = new TH1F("Graph_Graph419","Graph",100,0,3520);
   Graph_Graph419->SetMinimum(0.8307928);
   Graph_Graph419->SetMaximum(109.9077);
   Graph_Graph419->SetDirectory(nullptr);
   Graph_Graph419->SetStats(0);
   Graph_Graph419->SetLineWidth(2);
   Graph_Graph419->SetMarkerStyle(20);
   Graph_Graph419->SetMarkerSize(0.9);
   Graph_Graph419->GetXaxis()->SetLabelFont(42);
   Graph_Graph419->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetXaxis()->SetTitleFont(42);
   Graph_Graph419->GetYaxis()->SetLabelFont(42);
   Graph_Graph419->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetYaxis()->SetTickLength(0.02);
   Graph_Graph419->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetYaxis()->SetTitleFont(42);
   Graph_Graph419->GetZaxis()->SetLabelFont(42);
   Graph_Graph419->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph419);
   
   graph->Draw("p");
   
   Double_t Graph5_fx20[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy20[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 15.19713,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1801908, 0.8548021, 0.7603705, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph5_fx20,Graph5_fy20);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph520 = new TH1F("Graph_Graph520","Graph",100,0,3520);
   Graph_Graph520->SetMinimum(99.9);
   Graph_Graph520->SetMaximum(101.1);
   Graph_Graph520->SetDirectory(nullptr);
   Graph_Graph520->SetStats(0);
   Graph_Graph520->SetLineWidth(2);
   Graph_Graph520->SetMarkerStyle(20);
   Graph_Graph520->SetMarkerSize(0.9);
   Graph_Graph520->GetXaxis()->SetLabelFont(42);
   Graph_Graph520->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetXaxis()->SetTitleFont(42);
   Graph_Graph520->GetYaxis()->SetLabelFont(42);
   Graph_Graph520->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetYaxis()->SetTickLength(0.02);
   Graph_Graph520->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetYaxis()->SetTitleFont(42);
   Graph_Graph520->GetZaxis()->SetLabelFont(42);
   Graph_Graph520->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph520);
   
   graph->Draw("p");
   
   Double_t Graph6_fx21[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy21[35] = { 223.6068, 223.6068, 95.12887, 223.6068, 417.0685, 223.6068, 223.6068, 46.18585, 158.1785, 105.9233, 60.01927, 223.6068, 62.26757, 70.57571, 77.36527, 45.07711, 107.1173,
   53.91293, 71.9965, 45.17189, 75.06686, 63.26335, 55.37728, 44.32687, 42.12054, 32.11172, 23.64524, 22.12682, 24.68428, 21.64352, 23.75911, 29.45391, 33.64157,
   14.42741, 20.91349 };
   graph = new TGraph(35,Graph6_fx21,Graph6_fy21);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph621 = new TH1F("Graph_Graph621","Graph",100,0,3520);
   Graph_Graph621->SetMinimum(12.98467);
   Graph_Graph621->SetMaximum(457.3326);
   Graph_Graph621->SetDirectory(nullptr);
   Graph_Graph621->SetStats(0);
   Graph_Graph621->SetLineWidth(2);
   Graph_Graph621->SetMarkerStyle(20);
   Graph_Graph621->SetMarkerSize(0.9);
   Graph_Graph621->GetXaxis()->SetLabelFont(42);
   Graph_Graph621->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetXaxis()->SetTitleFont(42);
   Graph_Graph621->GetYaxis()->SetLabelFont(42);
   Graph_Graph621->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetYaxis()->SetTickLength(0.02);
   Graph_Graph621->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetYaxis()->SetTitleFont(42);
   Graph_Graph621->GetZaxis()->SetLabelFont(42);
   Graph_Graph621->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph621);
   
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
