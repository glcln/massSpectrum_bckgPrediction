#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
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
   
   Double_t Graph0_fx7[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy7[35] = { 100, 100, 0, 100, 100, 100, 46.42441, 52.17081, 0, 13.00907, 1.205385, 116.1401, 19.75127, 95.32021, 32.19221, 2.756488, 33.96618,
   94.20033, 80.66568, 16.0188, 16.89047, 23.71116, 14.41557, 9.256495, 4.469681, 9.27986, 10.27529, 4.888219, 3.591287, 7.436407, 8.184165, 4.964852, 18.13331,
   0.0346303, 3.954995 };
   TGraph *graph = new TGraph(35,Graph0_fx7,Graph0_fy7);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph07 = new TH1F("Graph_Graph07","Graph",100,0,3520);
   Graph_Graph07->SetMinimum(1);
   Graph_Graph07->SetMaximum(2000);
   Graph_Graph07->SetDirectory(nullptr);
   Graph_Graph07->SetStats(0);
   Graph_Graph07->SetLineWidth(2);
   Graph_Graph07->SetMarkerStyle(20);
   Graph_Graph07->SetMarkerSize(0.9);
   Graph_Graph07->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph07->GetXaxis()->SetRange(1,101);
   Graph_Graph07->GetXaxis()->SetLabelFont(43);
   Graph_Graph07->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph07->GetXaxis()->SetLabelSize(16);
   Graph_Graph07->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph07->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph07->GetXaxis()->SetTitleFont(42);
   Graph_Graph07->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph07->GetYaxis()->SetLabelFont(43);
   Graph_Graph07->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph07->GetYaxis()->SetLabelSize(16);
   Graph_Graph07->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph07->GetYaxis()->SetTickLength(0.02);
   Graph_Graph07->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph07->GetYaxis()->SetTitleFont(42);
   Graph_Graph07->GetZaxis()->SetLabelFont(42);
   Graph_Graph07->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph07->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph07->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph07->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph07->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph07);
   
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
   
   Double_t Graph1_fx8[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy8[35] = { 100, 100, 0, 100, 100, 100, 46.42441, 52.17081, 55.21128, 16.72954, 20.96955, 68.15215, 11.45847, 42.39326, 28.3689, 49.69147, 0,
   34.90425, 28.44741, 1.718658, 8.159041, 23.98763, 8.343286, 6.75174, 2.577507, 6.361085, 2.076471, 1.449955, 0.8623719, 2.671498, 2.213496, 3.648949, 6.128001,
   3.000736, 2.084446 };
   graph = new TGraph(35,Graph1_fx8,Graph1_fy8);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph18 = new TH1F("Graph_Graph18","Graph",100,0,3520);
   Graph_Graph18->SetMinimum(99.9);
   Graph_Graph18->SetMaximum(101.1);
   Graph_Graph18->SetDirectory(nullptr);
   Graph_Graph18->SetStats(0);
   Graph_Graph18->SetLineWidth(2);
   Graph_Graph18->SetMarkerStyle(20);
   Graph_Graph18->SetMarkerSize(0.9);
   Graph_Graph18->GetXaxis()->SetLabelFont(42);
   Graph_Graph18->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph18->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph18->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph18->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph18->GetXaxis()->SetTitleFont(42);
   Graph_Graph18->GetYaxis()->SetLabelFont(42);
   Graph_Graph18->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph18->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph18->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph18->GetYaxis()->SetTickLength(0.02);
   Graph_Graph18->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph18->GetYaxis()->SetTitleFont(42);
   Graph_Graph18->GetZaxis()->SetLabelFont(42);
   Graph_Graph18->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph18->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph18->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph18->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph18->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph18);
   
   graph->Draw("p");
   
   Double_t Graph2_fx9[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy9[35] = { 100, 100, 421.6929, 100, 100, 100, 32.14455, 23.26024, 29.75043, 62.08229, 20.02852, 19.90974, 22.3869, 21.07467, 20.11895, 32.47726, 64.52052,
   30.87213, 27.57621, 32.97998, 27.25786, 32.14111, 37.60148, 20.77674, 40.41966, 30.43981, 36.08125, 38.94223, 39.70607, 32.95885, 28.96262, 38.34497, 39.25953,
   13.52791, 48.00491 };
   graph = new TGraph(35,Graph2_fx9,Graph2_fy9);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph29 = new TH1F("Graph_Graph29","Graph",100,0,3520);
   Graph_Graph29->SetMinimum(12.17512);
   Graph_Graph29->SetMaximum(462.5094);
   Graph_Graph29->SetDirectory(nullptr);
   Graph_Graph29->SetStats(0);
   Graph_Graph29->SetLineWidth(2);
   Graph_Graph29->SetMarkerStyle(20);
   Graph_Graph29->SetMarkerSize(0.9);
   Graph_Graph29->GetXaxis()->SetLabelFont(42);
   Graph_Graph29->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph29->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph29->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph29->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph29->GetXaxis()->SetTitleFont(42);
   Graph_Graph29->GetYaxis()->SetLabelFont(42);
   Graph_Graph29->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph29->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph29->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph29->GetYaxis()->SetTickLength(0.02);
   Graph_Graph29->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph29->GetYaxis()->SetTitleFont(42);
   Graph_Graph29->GetZaxis()->SetLabelFont(42);
   Graph_Graph29->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph29->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph29->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph29->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph29->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph29);
   
   graph->Draw("p");
   
   Double_t Graph3_fx10[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy10[35] = { 100, 100, 0, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx10,Graph3_fy10);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph310 = new TH1F("Graph_Graph310","Graph",100,0,3520);
   Graph_Graph310->SetMinimum(99.9);
   Graph_Graph310->SetMaximum(101.1);
   Graph_Graph310->SetDirectory(nullptr);
   Graph_Graph310->SetStats(0);
   Graph_Graph310->SetLineWidth(2);
   Graph_Graph310->SetMarkerStyle(20);
   Graph_Graph310->SetMarkerSize(0.9);
   Graph_Graph310->GetXaxis()->SetLabelFont(42);
   Graph_Graph310->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph310->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph310->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph310->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph310->GetXaxis()->SetTitleFont(42);
   Graph_Graph310->GetYaxis()->SetLabelFont(42);
   Graph_Graph310->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph310->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph310->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph310->GetYaxis()->SetTickLength(0.02);
   Graph_Graph310->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph310->GetYaxis()->SetTitleFont(42);
   Graph_Graph310->GetZaxis()->SetLabelFont(42);
   Graph_Graph310->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph310->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph310->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph310->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph310->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph310);
   
   graph->Draw("p");
   
   Double_t Graph4_fx11[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy11[35] = { 100, 100, 0.6861925, 100, 100, 100, 1.857603, 0.9448171, 1.88511, 1.720887, 0.9384513, 0.7290006, 1.860726, 2.501905, 2.226472, 0.6861925, 2.132004,
   1.736772, 1.378757, 2.213824, 1.382339, 1.904792, 1.68156, 1.862413, 2.023786, 1.888394, 1.957238, 1.890534, 1.924783, 2.01059, 1.874787, 2.112126, 1.8206,
   2.341115, 1.968253 };
   graph = new TGraph(35,Graph4_fx11,Graph4_fy11);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph411 = new TH1F("Graph_Graph411","Graph",100,0,3520);
   Graph_Graph411->SetMinimum(0.6175733);
   Graph_Graph411->SetMaximum(109.9314);
   Graph_Graph411->SetDirectory(nullptr);
   Graph_Graph411->SetStats(0);
   Graph_Graph411->SetLineWidth(2);
   Graph_Graph411->SetMarkerStyle(20);
   Graph_Graph411->SetMarkerSize(0.9);
   Graph_Graph411->GetXaxis()->SetLabelFont(42);
   Graph_Graph411->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph411->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph411->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph411->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph411->GetXaxis()->SetTitleFont(42);
   Graph_Graph411->GetYaxis()->SetLabelFont(42);
   Graph_Graph411->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph411->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph411->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph411->GetYaxis()->SetTickLength(0.02);
   Graph_Graph411->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph411->GetYaxis()->SetTitleFont(42);
   Graph_Graph411->GetZaxis()->SetLabelFont(42);
   Graph_Graph411->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph411->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph411->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph411->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph411->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph411);
   
   graph->Draw("p");
   
   Double_t Graph5_fx12[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy12[35] = { 223.6068, 223.6068, 421.6935, 223.6068, 223.6068, 223.6068, 73.12438, 77.36613, 62.74493, 65.62229, 29.03786, 136.1256, 32.03193, 106.4591, 47.44325, 59.43134, 72.94617,
   105.11, 89.88078, 36.7714, 33.11736, 46.62944, 41.15967, 23.79937, 40.79787, 32.50734, 37.62421, 39.32008, 39.9239, 33.9524, 30.23621, 38.89424, 43.71492,
   14.05314, 48.2528 };
   graph = new TGraph(35,Graph5_fx12,Graph5_fy12);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph512 = new TH1F("Graph_Graph512","Graph",100,0,3520);
   Graph_Graph512->SetMinimum(12.64782);
   Graph_Graph512->SetMaximum(462.4575);
   Graph_Graph512->SetDirectory(nullptr);
   Graph_Graph512->SetStats(0);
   Graph_Graph512->SetLineWidth(2);
   Graph_Graph512->SetMarkerStyle(20);
   Graph_Graph512->SetMarkerSize(0.9);
   Graph_Graph512->GetXaxis()->SetLabelFont(42);
   Graph_Graph512->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph512->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph512->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph512->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph512->GetXaxis()->SetTitleFont(42);
   Graph_Graph512->GetYaxis()->SetLabelFont(42);
   Graph_Graph512->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph512->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph512->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph512->GetYaxis()->SetTickLength(0.02);
   Graph_Graph512->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph512->GetYaxis()->SetTitleFont(42);
   Graph_Graph512->GetZaxis()->SetLabelFont(42);
   Graph_Graph512->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph512->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph512->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph512->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph512->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph512);
   
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
