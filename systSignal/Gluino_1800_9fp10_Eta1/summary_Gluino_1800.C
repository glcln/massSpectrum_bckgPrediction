#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1800()
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
   
   Double_t Graph0_fx31[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy31[35] = { 100, 100, 100, 100, 100, 57.01979, 130.9829, 0, 0, 34.40939, 100, 37.10669, 51.81791, 50.83652, 265.9648, 5.960464e-06, 58.3087,
   36.7024, 135.4266, 58.01213, 0, 0, 31.56491, 0.2924025, 13.71965, 4.453588, 3.86712, 7.186043, 14.23861, 4.70252, 1.882386, 5.41935, 7.126045,
   7.591343, 6.071615 };
   TGraph *graph = new TGraph(35,Graph0_fx31,Graph0_fy31);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph031 = new TH1F("Graph_Graph031","Graph",100,0,3520);
   Graph_Graph031->SetMinimum(1);
   Graph_Graph031->SetMaximum(2000);
   Graph_Graph031->SetDirectory(nullptr);
   Graph_Graph031->SetStats(0);
   Graph_Graph031->SetLineWidth(2);
   Graph_Graph031->SetMarkerStyle(20);
   Graph_Graph031->SetMarkerSize(0.9);
   Graph_Graph031->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph031->GetXaxis()->SetRange(1,101);
   Graph_Graph031->GetXaxis()->SetLabelFont(43);
   Graph_Graph031->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph031->GetXaxis()->SetLabelSize(16);
   Graph_Graph031->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph031->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph031->GetXaxis()->SetTitleFont(42);
   Graph_Graph031->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph031->GetYaxis()->SetLabelFont(43);
   Graph_Graph031->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph031->GetYaxis()->SetLabelSize(16);
   Graph_Graph031->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph031->GetYaxis()->SetTickLength(0.02);
   Graph_Graph031->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph031->GetYaxis()->SetTitleFont(42);
   Graph_Graph031->GetZaxis()->SetLabelFont(42);
   Graph_Graph031->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph031->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph031->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph031->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph031->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph031);
   
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
   
   Double_t Graph1_fx32[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy32[35] = { 100, 100, 100, 100, 100, 57.01979, 100, 100.1288, 0, 34.40939, 100, 0, 0, 24.83575, 129.9349, 0, 0,
   0, 0, 0, 0, 15.09281, 16.29871, 0, 0, 3.825247, 3.321517, 3.191555, 2.02148, 0.8958757, 0.1705527, 1.315689, 1.673579,
   1.71541, 2.366638 };
   graph = new TGraph(35,Graph1_fx32,Graph1_fy32);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph132 = new TH1F("Graph_Graph132","Graph",100,0,3520);
   Graph_Graph132->SetMinimum(49.72828);
   Graph_Graph132->SetMaximum(137.2264);
   Graph_Graph132->SetDirectory(nullptr);
   Graph_Graph132->SetStats(0);
   Graph_Graph132->SetLineWidth(2);
   Graph_Graph132->SetMarkerStyle(20);
   Graph_Graph132->SetMarkerSize(0.9);
   Graph_Graph132->GetXaxis()->SetLabelFont(42);
   Graph_Graph132->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph132->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph132->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph132->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph132->GetXaxis()->SetTitleFont(42);
   Graph_Graph132->GetYaxis()->SetLabelFont(42);
   Graph_Graph132->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph132->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph132->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph132->GetYaxis()->SetTickLength(0.02);
   Graph_Graph132->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph132->GetYaxis()->SetTitleFont(42);
   Graph_Graph132->GetZaxis()->SetLabelFont(42);
   Graph_Graph132->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph132->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph132->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph132->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph132->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph132);
   
   graph->Draw("p");
   
   Double_t Graph2_fx33[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy33[35] = { 100, 100, 100, 100, 100, 30.49747, 19.08832, 30.59461, 16.60683, 30.49694, 100, 32.95188, 29.79082, 85.06232, 421.6929, 29.36035, 193.9571,
   137.8025, 25.3091, 23.38433, 30.9604, 22.72325, 26.08899, 37.61002, 29.82622, 28.21718, 80.78674, 54.69376, 41.87063, 36.95549, 39.37392, 47.85608, 39.79882,
   36.00775, 51.59035 };
   graph = new TGraph(35,Graph2_fx33,Graph2_fy33);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph233 = new TH1F("Graph_Graph233","Graph",100,0,3520);
   Graph_Graph233->SetMinimum(14.94615);
   Graph_Graph233->SetMaximum(462.2015);
   Graph_Graph233->SetDirectory(nullptr);
   Graph_Graph233->SetStats(0);
   Graph_Graph233->SetLineWidth(2);
   Graph_Graph233->SetMarkerStyle(20);
   Graph_Graph233->SetMarkerSize(0.9);
   Graph_Graph233->GetXaxis()->SetLabelFont(42);
   Graph_Graph233->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph233->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph233->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph233->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph233->GetXaxis()->SetTitleFont(42);
   Graph_Graph233->GetYaxis()->SetLabelFont(42);
   Graph_Graph233->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph233->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph233->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph233->GetYaxis()->SetTickLength(0.02);
   Graph_Graph233->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph233->GetYaxis()->SetTitleFont(42);
   Graph_Graph233->GetZaxis()->SetLabelFont(42);
   Graph_Graph233->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph233->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph233->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph233->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph233->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph233);
   
   graph->Draw("p");
   
   Double_t Graph3_fx34[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy34[35] = { 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx34,Graph3_fy34);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph334 = new TH1F("Graph_Graph334","Graph",100,0,3520);
   Graph_Graph334->SetMinimum(99.9);
   Graph_Graph334->SetMaximum(101.1);
   Graph_Graph334->SetDirectory(nullptr);
   Graph_Graph334->SetStats(0);
   Graph_Graph334->SetLineWidth(2);
   Graph_Graph334->SetMarkerStyle(20);
   Graph_Graph334->SetMarkerSize(0.9);
   Graph_Graph334->GetXaxis()->SetLabelFont(42);
   Graph_Graph334->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph334->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph334->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph334->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph334->GetXaxis()->SetTitleFont(42);
   Graph_Graph334->GetYaxis()->SetLabelFont(42);
   Graph_Graph334->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph334->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph334->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph334->GetYaxis()->SetTickLength(0.02);
   Graph_Graph334->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph334->GetYaxis()->SetTitleFont(42);
   Graph_Graph334->GetZaxis()->SetLabelFont(42);
   Graph_Graph334->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph334->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph334->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph334->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph334->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph334);
   
   graph->Draw("p");
   
   Double_t Graph4_fx35[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy35[35] = { 100, 100, 100, 100, 100, 1.93069, 0.7658124, 0.6861925, 1.018858, 2.080703, 100, 1.595163, 3.035951, 2.054685, 0.6861925, 1.471412, 1.176935,
   1.938158, 0.6861925, 1.141596, 0.7658124, 2.082121, 1.454198, 1.618856, 2.09052, 1.767403, 1.987302, 2.013385, 1.983297, 2.001059, 1.975119, 2.029467, 1.974368,
   2.017313, 1.897061 };
   graph = new TGraph(35,Graph4_fx35,Graph4_fy35);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph435 = new TH1F("Graph_Graph435","Graph",100,0,3520);
   Graph_Graph435->SetMinimum(0.6175733);
   Graph_Graph435->SetMaximum(109.9314);
   Graph_Graph435->SetDirectory(nullptr);
   Graph_Graph435->SetStats(0);
   Graph_Graph435->SetLineWidth(2);
   Graph_Graph435->SetMarkerStyle(20);
   Graph_Graph435->SetMarkerSize(0.9);
   Graph_Graph435->GetXaxis()->SetLabelFont(42);
   Graph_Graph435->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph435->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph435->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph435->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph435->GetXaxis()->SetTitleFont(42);
   Graph_Graph435->GetYaxis()->SetLabelFont(42);
   Graph_Graph435->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph435->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph435->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph435->GetYaxis()->SetTickLength(0.02);
   Graph_Graph435->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph435->GetYaxis()->SetTitleFont(42);
   Graph_Graph435->GetZaxis()->SetLabelFont(42);
   Graph_Graph435->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph435->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph435->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph435->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph435->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph435);
   
   graph->Draw("p");
   
   Double_t Graph5_fx36[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy36[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 86.23419, 165.896, 104.7009, 16.63806, 57.46655, 223.6068, 49.65156, 59.84819, 102.1811, 515.2142, 29.3972, 202.5355,
   142.6196, 137.773, 62.55827, 30.96987, 27.35826, 44.09923, 37.64598, 32.89685, 28.8756, 80.97181, 55.29273, 44.31599, 37.31794, 39.46872, 48.22264, 40.51451,
   36.89443, 52.03487 };
   graph = new TGraph(35,Graph5_fx36,Graph5_fy36);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph536 = new TH1F("Graph_Graph536","Graph",100,0,3520);
   Graph_Graph536->SetMinimum(14.97425);
   Graph_Graph536->SetMaximum(565.0719);
   Graph_Graph536->SetDirectory(nullptr);
   Graph_Graph536->SetStats(0);
   Graph_Graph536->SetLineWidth(2);
   Graph_Graph536->SetMarkerStyle(20);
   Graph_Graph536->SetMarkerSize(0.9);
   Graph_Graph536->GetXaxis()->SetLabelFont(42);
   Graph_Graph536->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph536->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph536->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph536->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph536->GetXaxis()->SetTitleFont(42);
   Graph_Graph536->GetYaxis()->SetLabelFont(42);
   Graph_Graph536->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph536->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph536->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph536->GetYaxis()->SetTickLength(0.02);
   Graph_Graph536->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph536->GetYaxis()->SetTitleFont(42);
   Graph_Graph536->GetZaxis()->SetLabelFont(42);
   Graph_Graph536->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph536->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph536->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph536->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph536->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph536);
   
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
