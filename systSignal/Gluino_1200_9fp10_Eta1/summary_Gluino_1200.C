#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
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
   
   Double_t Graph0_fx7[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy7[35] = { 100, 100, 0, 100, 100, 100, 31.95414, 52.17081, 0, 15.62265, 1.205385, 101.5408, 15.35309, 100, 0, 104.2526, 49.31997,
   0, 52.21827, 27.38346, 18.99071, 0.7394433, 8.930075, 12.21943, 6.498223, 8.174801, 11.41056, 4.717243, 3.979528, 7.604152, 7.723993, 7.114685, 16.94376,
   4.537153, 0 };
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
   Double_t Graph1_fy8[35] = { 100, 100, 0, 100, 100, 100, 31.95414, 52.17081, 0, 0, 24.99062, 101.5408, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 9.088439, 6.626821, 12.19109, 1.929015, 3.688282, 2.110797, 1.128531, 0.5622625, 2.613139, 2.469879, 2.091372, 8.187294,
   0, 0 };
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
   Graph_Graph18->SetMinimum(99.84592);
   Graph_Graph18->SetMaximum(101.6949);
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
   Double_t Graph2_fy9[35] = { 100, 100, 421.6929, 100, 100, 100, 25.37447, 23.26024, 27.44938, 76.67009, 20.02852, 25.3091, 20.31058, 100, 23.04946, 33.30547, 221.9099,
   30.87213, 27.57621, 22.59401, 27.24866, 41.75229, 34.97542, 26.1955, 40.72503, 23.03231, 33.77388, 40.33451, 40.48673, 36.91146, 25.58342, 33.2958, 43.01637,
   18.27147, 60.69012 };
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
   Graph_Graph29->SetMinimum(16.44432);
   Graph_Graph29->SetMaximum(462.035);
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
   Double_t Graph3_fy10[35] = { 100, 100, 0, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
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
   Double_t Graph4_fy11[35] = { 100, 100, 0.6861925, 100, 100, 100, 1.630151, 0.9448171, 0.7658124, 1.503354, 0.9384513, 0.7658124, 2.204525, 100, 2.803665, 0.6861925, 1.977164,
   1.736772, 1.378757, 2.337343, 1.418364, 2.049315, 1.680219, 1.741618, 2.007043, 1.940101, 1.943111, 1.912665, 1.952308, 2.011716, 1.916194, 2.250713, 1.822621,
   2.386904, 2.020293 };
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
   Double_t Graph5_fy12[35] = { 223.6068, 223.6068, 421.6935, 223.6068, 223.6068, 223.6068, 51.85224, 77.36613, 27.46006, 78.26002, 32.06254, 145.8157, 25.55576, 223.6068, 23.21935, 109.4456, 227.3332,
   30.92095, 59.06857, 35.57818, 33.24377, 42.78551, 36.73914, 31.41934, 41.33406, 24.79278, 35.7646, 40.6701, 40.73254, 37.83061, 26.9062, 34.18579, 46.9878,
   18.97708, 60.72374 };
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
   Graph_Graph512->SetMinimum(17.07937);
   Graph_Graph512->SetMaximum(461.9651);
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
