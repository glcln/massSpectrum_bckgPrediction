#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:21 2026) by ROOT version 6.32.13
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
   Double_t Graph0_fy37[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 100, 50.79404, 47.5344, 16.37311, 0, 0, 100, 33.33947, 35.61391,
   69.78165, 65.12205, 23.32451, 15.9013, 35.21984, 16.48809, 30.5739, 12.22681, 4.012132, 3.526175, 16.74457, 6.838405, 8.506304, 1.625228, 3.245211, 5.142826,
   6.424642, 4.456866 };
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
   Double_t Graph1_fy38[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 100, 0, 0,
   24.97638, 23.3086, 0, 0, 16.70763, 7.821643, 5.960464e-06, 0, 5.960464e-06, 0, 4.258454, 0.952363, 1.665574, 0.4566073, 0.7503867, 0.7407546,
   1.326144, 2.264476 };
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
   Graph_Graph138->SetMinimum(99.9);
   Graph_Graph138->SetMaximum(101.1);
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
   Double_t Graph2_fy39[35] = { 100, 100, 100, 100, 100, 100, 100, 24.72953, 100, 28.88285, 22.06556, 26.62921, 23.1173, 25.14429, 100, 14.54455, 33.68246,
   94.56793, 125.3133, 25.9168, 24.43519, 138.1641, 23.80236, 106.4152, 25.18825, 17.70713, 67.83899, 70.52579, 30.01871, 42.61378, 39.83301, 37.22956, 35.93593,
   33.16821, 44.10337 };
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
   Graph_Graph239->SetMinimum(2.182587);
   Graph_Graph239->SetMaximum(150.5261);
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
   Double_t Graph3_fy40[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 100, 0, 0,
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
   Double_t Graph4_fy41[35] = { 100, 100, 100, 100, 100, 100, 100, 1.359153, 100, 2.835542, 3.17356, 2.272475, 0.9149969, 2.110279, 100, 1.06492, 2.556038,
   1.659453, 2.29519, 1.987028, 1.630598, 1.85892, 1.841807, 1.809943, 2.211857, 1.86547, 1.824498, 1.947826, 1.96777, 1.970416, 1.984775, 2.023995, 1.989776,
   2.05667, 2.096176 };
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
   Graph_Graph441->SetMinimum(0.8234972);
   Graph_Graph441->SetMaximum(109.9085);
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
   Double_t Graph5_fy42[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 24.76685, 223.6068, 58.50037, 52.50219, 31.34258, 23.1354, 25.23269, 223.6068, 36.38953, 49.08556,
   120.163, 143.1533, 34.92366, 29.19912, 143.5701, 30.04963, 110.735, 28.08621, 18.25157, 67.95507, 72.63743, 30.86529, 43.531, 39.91814, 37.43302, 36.3641,
   33.87321, 44.43526 };
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
   Graph_Graph542->SetMinimum(16.42641);
   Graph_Graph542->SetMaximum(244.1423);
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
