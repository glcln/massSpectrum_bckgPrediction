#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:55 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx37[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy37[35] = { 100, 100, 100, 100, 100, 100, 100, 49.55408, 100, 50.79404, 0, 1.192093e-05, 0, 0, 43.5695, 45.27968, 35.61391,
   22.25349, 43.34724, 15.74968, 8.202707, 15.74835, 20.61926, 15.98032, 6.886125, 8.595938, 7.528031, 13.45668, 7.687884, 7.662308, 1.758361, 3.328443, 5.18589,
   6.370139, 2.743888 };
   TGraph *graph = new TGraph(35,Graph0_fx37,Graph0_fy37);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph037 = new TH1F("Graph_Graph037","Graph",100,0,3520);
   Graph_Graph037->SetMinimum(1);
   Graph_Graph037->SetMaximum(2000);
   Graph_Graph037->SetDirectory(nullptr);
   Graph_Graph037->SetStats(0);
   Graph_Graph037->SetLineWidth(2);
   Graph_Graph037->SetMarkerStyle(20);
   Graph_Graph037->SetMarkerSize(0.9);
   Graph_Graph037->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph037->GetXaxis()->SetRange(1,101);
   Graph_Graph037->GetXaxis()->SetLabelFont(43);
   Graph_Graph037->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph037->GetXaxis()->SetLabelSize(16);
   Graph_Graph037->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph037->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph037->GetXaxis()->SetTitleFont(42);
   Graph_Graph037->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph037->GetYaxis()->SetLabelFont(43);
   Graph_Graph037->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph037->GetYaxis()->SetLabelSize(16);
   Graph_Graph037->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph037->GetYaxis()->SetTickLength(0.02);
   Graph_Graph037->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph037->GetYaxis()->SetTitleFont(42);
   Graph_Graph037->GetZaxis()->SetLabelFont(42);
   Graph_Graph037->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph037->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph037->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph037->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph037->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph037);
   
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
   leg->Draw();
   
   Double_t Graph1_fx38[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy38[35] = { 100, 100, 100, 100, 100, 100, 100, 49.55408, 100, 0, 0, 1.192093e-05, 0, 0, 43.5695, 28.03841, 0,
   33.80045, 43.34724, 0, 8.945489, 6.997424, 6.268322, 4.65818, 5.194235, 2.749753, 4.850185, 2.750534, 0.9443998, 1.981807, 0.5667508, 0.8205175, 1.251996,
   1.977527, 1.394129 };
   graph = new TGraph(35,Graph1_fx38,Graph1_fy38);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph138 = new TH1F("Graph_Graph138","Graph",100,0,3520);
   Graph_Graph138->SetMinimum(44.50949);
   Graph_Graph138->SetMaximum(105.0446);
   Graph_Graph138->SetDirectory(nullptr);
   Graph_Graph138->SetStats(0);
   Graph_Graph138->SetLineWidth(2);
   Graph_Graph138->SetMarkerStyle(20);
   Graph_Graph138->SetMarkerSize(0.9);
   Graph_Graph138->GetXaxis()->SetLabelFont(42);
   Graph_Graph138->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph138->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph138->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph138->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph138->GetXaxis()->SetTitleFont(42);
   Graph_Graph138->GetYaxis()->SetLabelFont(42);
   Graph_Graph138->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph138->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph138->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph138->GetYaxis()->SetTickLength(0.02);
   Graph_Graph138->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph138->GetYaxis()->SetTitleFont(42);
   Graph_Graph138->GetZaxis()->SetLabelFont(42);
   Graph_Graph138->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph138->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph138->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph138->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph138->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph138);
   
   graph->Draw("p");
   
   Double_t Graph2_fx39[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy39[35] = { 100, 100, 100, 100, 100, 100, 16.41645, 24.72953, 100, 28.88285, 35.26473, 77.00236, 23.1173, 25.14429, 17.13395, 14.66124, 33.68246,
   54.94324, 20.74479, 29.6043, 23.53987, 113.0076, 27.66666, 75.9308, 24.85277, 23.56844, 75.96805, 66.2976, 31.12402, 43.43454, 40.15603, 36.53541, 35.52992,
   32.79652, 46.4936 };
   graph = new TGraph(35,Graph2_fx39,Graph2_fy39);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph239 = new TH1F("Graph_Graph239","Graph",100,0,3520);
   Graph_Graph239->SetMinimum(4.826606);
   Graph_Graph239->SetMaximum(122.8422);
   Graph_Graph239->SetDirectory(nullptr);
   Graph_Graph239->SetStats(0);
   Graph_Graph239->SetLineWidth(2);
   Graph_Graph239->SetMarkerStyle(20);
   Graph_Graph239->SetMarkerSize(0.9);
   Graph_Graph239->GetXaxis()->SetLabelFont(42);
   Graph_Graph239->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph239->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph239->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph239->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph239->GetXaxis()->SetTitleFont(42);
   Graph_Graph239->GetYaxis()->SetLabelFont(42);
   Graph_Graph239->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph239->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph239->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph239->GetYaxis()->SetTickLength(0.02);
   Graph_Graph239->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph239->GetYaxis()->SetTitleFont(42);
   Graph_Graph239->GetZaxis()->SetLabelFont(42);
   Graph_Graph239->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph239->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph239->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph239->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph239->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph239);
   
   graph->Draw("p");
   
   Double_t Graph3_fx40[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy40[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx40,Graph3_fy40);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph340 = new TH1F("Graph_Graph340","Graph",100,0,3520);
   Graph_Graph340->SetMinimum(99.9);
   Graph_Graph340->SetMaximum(101.1);
   Graph_Graph340->SetDirectory(nullptr);
   Graph_Graph340->SetStats(0);
   Graph_Graph340->SetLineWidth(2);
   Graph_Graph340->SetMarkerStyle(20);
   Graph_Graph340->SetMarkerSize(0.9);
   Graph_Graph340->GetXaxis()->SetLabelFont(42);
   Graph_Graph340->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph340->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph340->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph340->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph340->GetXaxis()->SetTitleFont(42);
   Graph_Graph340->GetYaxis()->SetLabelFont(42);
   Graph_Graph340->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph340->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph340->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph340->GetYaxis()->SetTickLength(0.02);
   Graph_Graph340->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph340->GetYaxis()->SetTitleFont(42);
   Graph_Graph340->GetZaxis()->SetLabelFont(42);
   Graph_Graph340->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph340->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph340->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph340->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph340->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph340);
   
   graph->Draw("p");
   
   Double_t Graph4_fx41[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy41[35] = { 100, 100, 100, 100, 100, 100, 2.803659, 1.359153, 100, 2.835542, 2.281326, 3.163213, 0.9149969, 2.110279, 0.8152604, 1.042372, 2.556038,
   1.4642, 1.709199, 1.983315, 1.676643, 1.641387, 1.980478, 1.952934, 1.810837, 1.835269, 1.910603, 1.939541, 1.987988, 1.957232, 1.984751, 2.004683, 1.974934,
   1.990676, 2.078325 };
   graph = new TGraph(35,Graph4_fx41,Graph4_fy41);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph441 = new TH1F("Graph_Graph441","Graph",100,0,3520);
   Graph_Graph441->SetMinimum(0.7337344);
   Graph_Graph441->SetMaximum(109.9185);
   Graph_Graph441->SetDirectory(nullptr);
   Graph_Graph441->SetStats(0);
   Graph_Graph441->SetLineWidth(2);
   Graph_Graph441->SetMarkerStyle(20);
   Graph_Graph441->SetMarkerSize(0.9);
   Graph_Graph441->GetXaxis()->SetLabelFont(42);
   Graph_Graph441->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph441->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph441->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph441->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph441->GetXaxis()->SetTitleFont(42);
   Graph_Graph441->GetYaxis()->SetLabelFont(42);
   Graph_Graph441->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph441->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph441->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph441->GetYaxis()->SetTickLength(0.02);
   Graph_Graph441->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph441->GetYaxis()->SetTitleFont(42);
   Graph_Graph441->GetZaxis()->SetLabelFont(42);
   Graph_Graph441->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph441->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph441->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph441->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph441->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph441);
   
   graph->Draw("p");
   
   Double_t Graph5_fx42[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy42[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 142.3986, 74.32773, 223.6068, 58.50037, 35.33845, 77.0673, 23.1354, 25.23269, 63.95967, 55.24889, 49.08556,
   68.25388, 64.73974, 33.59167, 26.53757, 114.3258, 35.12566, 77.75841, 26.36926, 25.30397, 76.51791, 67.73316, 32.13491, 44.19308, 40.24747, 36.7506, 35.98245,
   33.52707, 46.64169 };
   graph = new TGraph(35,Graph5_fx42,Graph5_fy42);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph542 = new TH1F("Graph_Graph542","Graph",100,0,3520);
   Graph_Graph542->SetMinimum(3.088259);
   Graph_Graph542->SetMaximum(243.6539);
   Graph_Graph542->SetDirectory(nullptr);
   Graph_Graph542->SetStats(0);
   Graph_Graph542->SetLineWidth(2);
   Graph_Graph542->SetMarkerStyle(20);
   Graph_Graph542->SetMarkerSize(0.9);
   Graph_Graph542->GetXaxis()->SetLabelFont(42);
   Graph_Graph542->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph542->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph542->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph542->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph542->GetXaxis()->SetTitleFont(42);
   Graph_Graph542->GetYaxis()->SetLabelFont(42);
   Graph_Graph542->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph542->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph542->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph542->GetYaxis()->SetTickLength(0.02);
   Graph_Graph542->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph542->GetYaxis()->SetTitleFont(42);
   Graph_Graph542->GetZaxis()->SetLabelFont(42);
   Graph_Graph542->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph542->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph542->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph542->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph542->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph542);
   
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
