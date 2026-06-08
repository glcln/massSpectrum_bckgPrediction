#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1400()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:54 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx19[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy19[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 27.31032, 45.85263, 31.7054, 127.2741, 21.48189, 15.53789, 18.4504,
   31.82392, 19.10694, 34.74795, 16.51501, 0.3931165, 1.063257, 22.16705, 7.535982, 11.04661, 10.03911, 9.12354, 4.225028, 3.018057, 8.351099, 6.326222, 4.306716,
   7.177806, 5.584931 };
   TGraph *graph = new TGraph(35,Graph0_fx19,Graph0_fy19);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph019 = new TH1F("Graph_Graph019","Graph",100,0,3520);
   Graph_Graph019->SetMinimum(1);
   Graph_Graph019->SetMaximum(2000);
   Graph_Graph019->SetDirectory(nullptr);
   Graph_Graph019->SetStats(0);
   Graph_Graph019->SetLineWidth(2);
   Graph_Graph019->SetMarkerStyle(20);
   Graph_Graph019->SetMarkerSize(0.9);
   Graph_Graph019->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph019->GetXaxis()->SetRange(1,101);
   Graph_Graph019->GetXaxis()->SetLabelFont(43);
   Graph_Graph019->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph019->GetXaxis()->SetLabelSize(16);
   Graph_Graph019->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph019->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph019->GetXaxis()->SetTitleFont(42);
   Graph_Graph019->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph019->GetYaxis()->SetLabelFont(43);
   Graph_Graph019->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph019->GetYaxis()->SetLabelSize(16);
   Graph_Graph019->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph019->GetYaxis()->SetTickLength(0.02);
   Graph_Graph019->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph019->GetYaxis()->SetTitleFont(42);
   Graph_Graph019->GetZaxis()->SetLabelFont(42);
   Graph_Graph019->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph019->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph019->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph019->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph019->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph019);
   
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
   
   Double_t Graph1_fx20[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy20[35] = { 100, 100, 100, 32.42584, 0, 100, 0, 0, 0, 25.22911, 32.37425, 37.23769, 5.155593, 100, 0, 15.53789, 23.28366,
   25.50403, 11.08288, 20.15537, 16.51501, 11.18983, 10.87582, 13.50038, 3.417975, 5.857432, 3.399777, 2.219594, 0.7285714, 0.5513132, 2.893126, 1.988131, 0.6225467,
   1.58841, 1.66657 };
   graph = new TGraph(35,Graph1_fx20,Graph1_fy20);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph120 = new TH1F("Graph_Graph120","Graph",100,0,3520);
   Graph_Graph120->SetMinimum(25.66842);
   Graph_Graph120->SetMaximum(106.7574);
   Graph_Graph120->SetDirectory(nullptr);
   Graph_Graph120->SetStats(0);
   Graph_Graph120->SetLineWidth(2);
   Graph_Graph120->SetMarkerStyle(20);
   Graph_Graph120->SetMarkerSize(0.9);
   Graph_Graph120->GetXaxis()->SetLabelFont(42);
   Graph_Graph120->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph120->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph120->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph120->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph120->GetXaxis()->SetTitleFont(42);
   Graph_Graph120->GetYaxis()->SetLabelFont(42);
   Graph_Graph120->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph120->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph120->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph120->GetYaxis()->SetTickLength(0.02);
   Graph_Graph120->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph120->GetYaxis()->SetTitleFont(42);
   Graph_Graph120->GetZaxis()->SetLabelFont(42);
   Graph_Graph120->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph120->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph120->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph120->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph120->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph120);
   
   graph->Draw("p");
   
   Double_t Graph2_fx21[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy21[35] = { 100, 100, 100, 149.2289, 25.3091, 100, 29.50337, 30.32486, 24.55111, 213.9225, 30.23605, 124.8454, 26.78555, 33.93324, 83.5133, 27.51592, 19.2637,
   123.3517, 26.26516, 25.95093, 24.13324, 41.26838, 23.03386, 41.53763, 15.04145, 43.44556, 47.96166, 35.53811, 38.49498, 35.99485, 36.56956, 47.46777, 31.75687,
   17.54479, 32.85889 };
   graph = new TGraph(35,Graph2_fx21,Graph2_fy21);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph221 = new TH1F("Graph_Graph221","Graph",100,0,3520);
   Graph_Graph221->SetMinimum(13.5373);
   Graph_Graph221->SetMaximum(233.8106);
   Graph_Graph221->SetDirectory(nullptr);
   Graph_Graph221->SetStats(0);
   Graph_Graph221->SetLineWidth(2);
   Graph_Graph221->SetMarkerStyle(20);
   Graph_Graph221->SetMarkerSize(0.9);
   Graph_Graph221->GetXaxis()->SetLabelFont(42);
   Graph_Graph221->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph221->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph221->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph221->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph221->GetXaxis()->SetTitleFont(42);
   Graph_Graph221->GetYaxis()->SetLabelFont(42);
   Graph_Graph221->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph221->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph221->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph221->GetYaxis()->SetTickLength(0.02);
   Graph_Graph221->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph221->GetYaxis()->SetTitleFont(42);
   Graph_Graph221->GetZaxis()->SetLabelFont(42);
   Graph_Graph221->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph221->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph221->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph221->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph221->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph221);
   
   graph->Draw("p");
   
   Double_t Graph3_fx22[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy22[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx22,Graph3_fy22);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph322 = new TH1F("Graph_Graph322","Graph",100,0,3520);
   Graph_Graph322->SetMinimum(99.9);
   Graph_Graph322->SetMaximum(101.1);
   Graph_Graph322->SetDirectory(nullptr);
   Graph_Graph322->SetStats(0);
   Graph_Graph322->SetLineWidth(2);
   Graph_Graph322->SetMarkerStyle(20);
   Graph_Graph322->SetMarkerSize(0.9);
   Graph_Graph322->GetXaxis()->SetLabelFont(42);
   Graph_Graph322->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph322->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph322->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph322->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph322->GetXaxis()->SetTitleFont(42);
   Graph_Graph322->GetYaxis()->SetLabelFont(42);
   Graph_Graph322->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph322->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph322->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph322->GetYaxis()->SetTickLength(0.02);
   Graph_Graph322->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph322->GetYaxis()->SetTitleFont(42);
   Graph_Graph322->GetZaxis()->SetLabelFont(42);
   Graph_Graph322->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph322->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph322->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph322->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph322->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph322);
   
   graph->Draw("p");
   
   Double_t Graph4_fx23[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy23[35] = { 100, 100, 100, 2.292418, 2.803659, 100, 0.84185, 0.7267118, 2.248728, 1.372117, 1.66949, 2.350748, 2.132672, 0.6862044, 2.527964, 1.744366, 1.82004,
   1.074624, 1.826352, 2.113408, 1.935589, 1.577759, 2.386886, 2.02812, 2.193898, 2.022386, 1.852459, 1.89805, 1.925272, 1.915145, 1.959032, 1.891667, 2.0446,
   2.038169, 1.841056 };
   graph = new TGraph(35,Graph4_fx23,Graph4_fy23);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph423 = new TH1F("Graph_Graph423","Graph",100,0,3520);
   Graph_Graph423->SetMinimum(0.617584);
   Graph_Graph423->SetMaximum(109.9314);
   Graph_Graph423->SetDirectory(nullptr);
   Graph_Graph423->SetStats(0);
   Graph_Graph423->SetLineWidth(2);
   Graph_Graph423->SetMarkerStyle(20);
   Graph_Graph423->SetMarkerSize(0.9);
   Graph_Graph423->GetXaxis()->SetLabelFont(42);
   Graph_Graph423->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph423->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph423->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph423->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph423->GetXaxis()->SetTitleFont(42);
   Graph_Graph423->GetYaxis()->SetLabelFont(42);
   Graph_Graph423->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph423->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph423->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph423->GetYaxis()->SetTickLength(0.02);
   Graph_Graph423->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph423->GetYaxis()->SetTitleFont(42);
   Graph_Graph423->GetZaxis()->SetLabelFont(42);
   Graph_Graph423->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph423->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph423->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph423->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph423->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph423);
   
   graph->Draw("p");
   
   Double_t Graph5_fx24[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy24[35] = { 223.6068, 223.6068, 223.6068, 152.7283, 25.46391, 223.6068, 29.51538, 30.33357, 24.65388, 215.4094, 52.0668, 138.134, 41.87871, 165.3803, 86.26896, 35.2565, 35.45347,
   129.9231, 34.36713, 47.87041, 33.64001, 42.78943, 25.60606, 49.02169, 17.30699, 45.25422, 49.15379, 36.8066, 38.78082, 36.17609, 37.67336, 47.96604, 32.11876,
   19.13159, 33.42253 };
   graph = new TGraph(35,Graph5_fx24,Graph5_fy24);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph524 = new TH1F("Graph_Graph524","Graph",100,0,3520);
   Graph_Graph524->SetMinimum(15.57629);
   Graph_Graph524->SetMaximum(244.2368);
   Graph_Graph524->SetDirectory(nullptr);
   Graph_Graph524->SetStats(0);
   Graph_Graph524->SetLineWidth(2);
   Graph_Graph524->SetMarkerStyle(20);
   Graph_Graph524->SetMarkerSize(0.9);
   Graph_Graph524->GetXaxis()->SetLabelFont(42);
   Graph_Graph524->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph524->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph524->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph524->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph524->GetXaxis()->SetTitleFont(42);
   Graph_Graph524->GetYaxis()->SetLabelFont(42);
   Graph_Graph524->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph524->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph524->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph524->GetYaxis()->SetTickLength(0.02);
   Graph_Graph524->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph524->GetYaxis()->SetTitleFont(42);
   Graph_Graph524->GetZaxis()->SetLabelFont(42);
   Graph_Graph524->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph524->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph524->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph524->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph524->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph524);
   
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
