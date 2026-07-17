#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1800()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:55:04 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx36[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy36[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 139.294, 100, 100, 23.72633, 65.14432, 100,
   100, 64.42568, 69.43282, 0, 124.2456, 18.21475, 10.87024, 4.470694, 2.242911, 3.399038, 10.90248, 5.512416, 5.83818, 1.204956, 6.354463, 5.682731,
   3.737473, 6.47577 };
   TGraph *graph = new TGraph(35,Graph0_fx36,Graph0_fy36);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph036 = new TH1F("Graph_Graph036","Graph",100,0,3520);
   Graph_Graph036->SetMinimum(1);
   Graph_Graph036->SetMaximum(2000);
   Graph_Graph036->SetDirectory(nullptr);
   Graph_Graph036->SetStats(0);
   Graph_Graph036->SetLineWidth(2);
   Graph_Graph036->SetMarkerStyle(20);
   Graph_Graph036->SetMarkerSize(0.9);
   Graph_Graph036->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph036->GetXaxis()->SetRange(1,101);
   Graph_Graph036->GetXaxis()->SetLabelFont(43);
   Graph_Graph036->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetXaxis()->SetLabelSize(16);
   Graph_Graph036->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph036->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph036->GetXaxis()->SetTitleFont(42);
   Graph_Graph036->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph036->GetYaxis()->SetLabelFont(43);
   Graph_Graph036->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetYaxis()->SetLabelSize(16);
   Graph_Graph036->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph036->GetYaxis()->SetTickLength(0.02);
   Graph_Graph036->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph036->GetYaxis()->SetTitleFont(42);
   Graph_Graph036->GetZaxis()->SetLabelFont(42);
   Graph_Graph036->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph036->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph036->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph036->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph036);
   
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
   
   Double_t Graph1_fx37[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy37[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 139.294, 100, 100, 0, 65.14432, 100,
   49.94995, 44.87506, 53.22217, 0, 0, 5.811632, 25.53954, 11.04062, 11.04808, 7.176107, 3.293759, 2.987289, 1.632875, 0.7915974, 2.657419, 2.800661,
   1.807404, 1.47351 };
   graph = new TGraph(35,Graph1_fx37,Graph1_fy37);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph137 = new TH1F("Graph_Graph137","Graph",100,0,3520);
   Graph_Graph137->SetMinimum(96.0706);
   Graph_Graph137->SetMaximum(143.2234);
   Graph_Graph137->SetDirectory(nullptr);
   Graph_Graph137->SetStats(0);
   Graph_Graph137->SetLineWidth(2);
   Graph_Graph137->SetMarkerStyle(20);
   Graph_Graph137->SetMarkerSize(0.9);
   Graph_Graph137->GetXaxis()->SetLabelFont(42);
   Graph_Graph137->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetXaxis()->SetTitleFont(42);
   Graph_Graph137->GetYaxis()->SetLabelFont(42);
   Graph_Graph137->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetYaxis()->SetTickLength(0.02);
   Graph_Graph137->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetYaxis()->SetTitleFont(42);
   Graph_Graph137->GetZaxis()->SetLabelFont(42);
   Graph_Graph137->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph137);
   
   graph->Draw("p");
   
   Double_t Graph2_fx38[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy38[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 54.75391, 100, 100, 69.78245, 33.50092, 100, 74.11238, 61.54933, 100,
   68.54632, 12.86974, 24.4803, 49.13321, 30.78724, 61.00661, 16.26924, 26.0963, 73.52899, 48.51663, 34.2016, 13.89034, 24.80097, 21.79881, 19.68087, 17.60989,
   26.96909, 16.92555 };
   graph = new TGraph(35,Graph2_fx38,Graph2_fy38);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph238 = new TH1F("Graph_Graph238","Graph",100,0,3520);
   Graph_Graph238->SetMinimum(4.156713);
   Graph_Graph238->SetMaximum(108.713);
   Graph_Graph238->SetDirectory(nullptr);
   Graph_Graph238->SetStats(0);
   Graph_Graph238->SetLineWidth(2);
   Graph_Graph238->SetMarkerStyle(20);
   Graph_Graph238->SetMarkerSize(0.9);
   Graph_Graph238->GetXaxis()->SetLabelFont(42);
   Graph_Graph238->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetXaxis()->SetTitleFont(42);
   Graph_Graph238->GetYaxis()->SetLabelFont(42);
   Graph_Graph238->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetYaxis()->SetTickLength(0.02);
   Graph_Graph238->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetYaxis()->SetTitleFont(42);
   Graph_Graph238->GetZaxis()->SetLabelFont(42);
   Graph_Graph238->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph238);
   
   graph->Draw("p");
   
   Double_t Graph3_fx39[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy39[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 0, 0, 100, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx39,Graph3_fy39);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph339 = new TH1F("Graph_Graph339","Graph",100,0,3520);
   Graph_Graph339->SetMinimum(99.9);
   Graph_Graph339->SetMaximum(101.1);
   Graph_Graph339->SetDirectory(nullptr);
   Graph_Graph339->SetStats(0);
   Graph_Graph339->SetLineWidth(2);
   Graph_Graph339->SetMarkerStyle(20);
   Graph_Graph339->SetMarkerSize(0.9);
   Graph_Graph339->GetXaxis()->SetLabelFont(42);
   Graph_Graph339->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetXaxis()->SetTitleFont(42);
   Graph_Graph339->GetYaxis()->SetLabelFont(42);
   Graph_Graph339->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetYaxis()->SetTickLength(0.02);
   Graph_Graph339->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetYaxis()->SetTitleFont(42);
   Graph_Graph339->GetZaxis()->SetLabelFont(42);
   Graph_Graph339->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph339);
   
   graph->Draw("p");
   
   Double_t Graph4_fx40[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy40[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 3.111076, 100, 100, 4.550165, 0.9231091, 100, 4.738003, 2.256489, 100,
   3.297698, 2.327281, 2.691734, 2.530158, 2.042806, 2.276754, 1.772487, 2.398944, 2.567494, 2.795374, 2.659941, 3.201938, 3.292239, 3.279352, 3.356063, 3.401518,
   3.47116, 3.693223 };
   graph = new TGraph(35,Graph4_fx40,Graph4_fy40);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph440 = new TH1F("Graph_Graph440","Graph",100,0,3520);
   Graph_Graph440->SetMinimum(0.8307981);
   Graph_Graph440->SetMaximum(109.9077);
   Graph_Graph440->SetDirectory(nullptr);
   Graph_Graph440->SetStats(0);
   Graph_Graph440->SetLineWidth(2);
   Graph_Graph440->SetMarkerStyle(20);
   Graph_Graph440->SetMarkerSize(0.9);
   Graph_Graph440->GetXaxis()->SetLabelFont(42);
   Graph_Graph440->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetXaxis()->SetTitleFont(42);
   Graph_Graph440->GetYaxis()->SetLabelFont(42);
   Graph_Graph440->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetYaxis()->SetTickLength(0.02);
   Graph_Graph440->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetYaxis()->SetTitleFont(42);
   Graph_Graph440->GetZaxis()->SetLabelFont(42);
   Graph_Graph440->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph440);
   
   graph->Draw("p");
   
   Double_t Graph5_fx41[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy41[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 0, 0, 100, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2067506, 0.3705263, 0.6354034, 0.1601219, 0.2904236,
   0, 0.8969367 };
   graph = new TGraph(35,Graph5_fx41,Graph5_fy41);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph541 = new TH1F("Graph_Graph541","Graph",100,0,3520);
   Graph_Graph541->SetMinimum(99.9);
   Graph_Graph541->SetMaximum(101.1);
   Graph_Graph541->SetDirectory(nullptr);
   Graph_Graph541->SetStats(0);
   Graph_Graph541->SetLineWidth(2);
   Graph_Graph541->SetMarkerStyle(20);
   Graph_Graph541->SetMarkerSize(0.9);
   Graph_Graph541->GetXaxis()->SetLabelFont(42);
   Graph_Graph541->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetXaxis()->SetTitleFont(42);
   Graph_Graph541->GetYaxis()->SetLabelFont(42);
   Graph_Graph541->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetYaxis()->SetTickLength(0.02);
   Graph_Graph541->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetYaxis()->SetTitleFont(42);
   Graph_Graph541->GetZaxis()->SetLabelFont(42);
   Graph_Graph541->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph541);
   
   graph->Draw("p");
   
   Double_t Graph6_fx42[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy42[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 54.84222, 223.6068, 223.6068, 209.0358, 145.3381, 223.6068, 77.96174, 110.8196, 223.6068,
   131.1658, 79.59576, 90.88479, 49.19831, 128.0195, 63.97298, 32.22204, 28.78635, 74.43249, 49.24152, 36.14606, 15.57256, 25.74252, 22.0912, 21.11968, 19.02146,
   27.50665, 18.55319 };
   graph = new TGraph(35,Graph6_fx42,Graph6_fy42);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph642 = new TH1F("Graph_Graph642","Graph",100,0,3520);
   Graph_Graph642->SetMinimum(14.0153);
   Graph_Graph642->SetMaximum(244.4102);
   Graph_Graph642->SetDirectory(nullptr);
   Graph_Graph642->SetStats(0);
   Graph_Graph642->SetLineWidth(2);
   Graph_Graph642->SetMarkerStyle(20);
   Graph_Graph642->SetMarkerSize(0.9);
   Graph_Graph642->GetXaxis()->SetLabelFont(42);
   Graph_Graph642->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetXaxis()->SetTitleFont(42);
   Graph_Graph642->GetYaxis()->SetLabelFont(42);
   Graph_Graph642->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetYaxis()->SetTickLength(0.02);
   Graph_Graph642->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetYaxis()->SetTitleFont(42);
   Graph_Graph642->GetZaxis()->SetLabelFont(42);
   Graph_Graph642->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph642);
   
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
