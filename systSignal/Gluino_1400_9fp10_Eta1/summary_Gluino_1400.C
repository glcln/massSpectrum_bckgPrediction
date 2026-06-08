#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1400()
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
   
   Double_t Graph0_fx19[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy19[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 228.5594, 66.17436, 100, 45.20653, 0, 47.65553,
   52.61921, 35.5909, 63.56396, 20.33415, 16.58136, 1.192093e-05, 43.67359, 1.004195, 10.55649, 7.881415, 9.541297, 4.544079, 3.050244, 8.849848, 5.512023, 5.083734,
   12.51917, 4.167008 };
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
   Double_t Graph1_fy20[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 57.66947, 32.8795, 114.0335, 33.01589, 100, 0, 0, 0,
   0, 20.64432, 36.86994, 20.33415, 16.58136, 0, 0, 5.130672, 2.833849, 0.07197261, 2.030855, 0.6624699, 0.3459036, 2.497232, 1.471961, 1.471663,
   2.739036, 1.192093e-05 };
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
   Graph_Graph120->SetMinimum(98.59665);
   Graph_Graph120->SetMaximum(115.4368);
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
   Double_t Graph2_fy21[35] = { 100, 100, 100, 218.7126, 25.3091, 100, 31.63878, 27.44938, 24.55111, 26.33525, 29.09412, 23.04946, 30.47314, 100, 188.752, 24.03646, 24.39324,
   210.8921, 21.84011, 18.66727, 22.39756, 67.5807, 23.57038, 23.84483, 21.7206, 57.13304, 51.68179, 32.18755, 40.14076, 35.59308, 33.71355, 45.73558, 32.23565,
   15.54607, 37.75015 };
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
   Graph_Graph221->SetMinimum(13.99146);
   Graph_Graph221->SetMaximum(239.0292);
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
   Double_t Graph3_fy22[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
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
   Double_t Graph4_fy23[35] = { 100, 100, 100, 2.953422, 2.803659, 100, 0.7658124, 0.7658124, 2.248728, 0.7297277, 1.120788, 4.264009, 2.472895, 100, 2.517456, 1.395035, 2.769458,
   1.278871, 2.08599, 1.540977, 1.734912, 1.906019, 2.293688, 2.398777, 2.220047, 2.166051, 1.871145, 1.93271, 1.934409, 1.922393, 1.982027, 1.913971, 2.104652,
   2.111804, 1.869512 };
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
   Graph_Graph423->SetMinimum(0.656755);
   Graph_Graph423->SetMaximum(109.927);
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
   Double_t Graph5_fy24[35] = { 223.6068, 223.6068, 223.6068, 218.7325, 25.46391, 223.6068, 31.64804, 27.46006, 24.65388, 63.40225, 43.91795, 256.5005, 80.02389, 223.6068, 194.1064, 24.0769, 53.60736,
   217.3612, 46.62877, 75.83279, 36.49131, 71.55884, 23.68172, 49.81678, 22.45095, 58.2095, 52.31282, 33.68879, 40.44886, 35.7769, 35.00126, 46.12977, 32.73494,
   20.25764, 38.02542 };
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
   Graph_Graph524->SetMinimum(18.23188);
   Graph_Graph524->SetMaximum(280.1248);
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
