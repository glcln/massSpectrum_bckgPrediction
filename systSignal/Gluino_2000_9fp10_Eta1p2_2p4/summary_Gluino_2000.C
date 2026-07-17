#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 16:06:00 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx43[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy43[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 100, 100,
   0, 25.83271, 14.05017, 28.61946, 3.995663, 207.3519, 9.943104, 7.217711, 1.558352, 10.72342, 1.174915, 9.215916, 9.580982, 3.046113, 4.020119, 6.635463,
   8.451092, 1.467955 };
   TGraph *graph = new TGraph(35,Graph0_fx43,Graph0_fy43);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph043 = new TH1F("Graph_Graph043","Graph",100,0,3520);
   Graph_Graph043->SetMinimum(1);
   Graph_Graph043->SetMaximum(2000);
   Graph_Graph043->SetDirectory(nullptr);
   Graph_Graph043->SetStats(0);
   Graph_Graph043->SetLineWidth(2);
   Graph_Graph043->SetMarkerStyle(20);
   Graph_Graph043->SetMarkerSize(0.9);
   Graph_Graph043->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph043->GetXaxis()->SetRange(1,101);
   Graph_Graph043->GetXaxis()->SetLabelFont(43);
   Graph_Graph043->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph043->GetXaxis()->SetLabelSize(16);
   Graph_Graph043->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph043->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph043->GetXaxis()->SetTitleFont(42);
   Graph_Graph043->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph043->GetYaxis()->SetLabelFont(43);
   Graph_Graph043->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph043->GetYaxis()->SetLabelSize(16);
   Graph_Graph043->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph043->GetYaxis()->SetTickLength(0.02);
   Graph_Graph043->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph043->GetYaxis()->SetTitleFont(42);
   Graph_Graph043->GetZaxis()->SetLabelFont(42);
   Graph_Graph043->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph043->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph043->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph043->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph043->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph043);
   
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
   
   Double_t Graph1_fx44[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy44[35] = { 100, 100, 100, 100, 100, 100, 0, 100, 100, 100, 100, 0, 100, 0, 0, 100, 100,
   85.04453, 26.34826, 14.33058, 28.61946, 3.995663, 207.3519, 0, 7.089668, 5.171323, 5.27122, 2.011478, 1.885128, 2.260828, 0.903672, 1.181245, 2.031887,
   0.4418015, 0.8272171 };
   graph = new TGraph(35,Graph1_fx44,Graph1_fy44);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph144 = new TH1F("Graph_Graph144","Graph",100,0,3520);
   Graph_Graph144->SetMinimum(89.26481);
   Graph_Graph144->SetMaximum(218.0871);
   Graph_Graph144->SetDirectory(nullptr);
   Graph_Graph144->SetStats(0);
   Graph_Graph144->SetLineWidth(2);
   Graph_Graph144->SetMarkerStyle(20);
   Graph_Graph144->SetMarkerSize(0.9);
   Graph_Graph144->GetXaxis()->SetLabelFont(42);
   Graph_Graph144->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph144->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph144->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph144->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph144->GetXaxis()->SetTitleFont(42);
   Graph_Graph144->GetYaxis()->SetLabelFont(42);
   Graph_Graph144->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph144->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph144->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph144->GetYaxis()->SetTickLength(0.02);
   Graph_Graph144->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph144->GetYaxis()->SetTitleFont(42);
   Graph_Graph144->GetZaxis()->SetLabelFont(42);
   Graph_Graph144->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph144->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph144->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph144->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph144->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph144);
   
   graph->Draw("p");
   
   Double_t Graph2_fx45[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy45[35] = { 100, 100, 100, 100, 100, 100, 54.13097, 65.51646, 100, 100, 100, 300.4514, 100, 92.76472, 124.4167, 100, 100,
   73.68394, 60.38054, 53.83511, 76.8599, 100.4092, 297.1491, 28.17622, 60.81168, 20.84193, 37.6928, 14.26826, 15.16874, 24.5892, 23.51752, 20.77532, 23.23298,
   16.82072, 12.8428 };
   graph = new TGraph(35,Graph2_fx45,Graph2_fy45);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph245 = new TH1F("Graph_Graph245","Graph",100,0,3520);
   Graph_Graph245->SetMinimum(11.55852);
   Graph_Graph245->SetMaximum(329.2123);
   Graph_Graph245->SetDirectory(nullptr);
   Graph_Graph245->SetStats(0);
   Graph_Graph245->SetLineWidth(2);
   Graph_Graph245->SetMarkerStyle(20);
   Graph_Graph245->SetMarkerSize(0.9);
   Graph_Graph245->GetXaxis()->SetLabelFont(42);
   Graph_Graph245->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph245->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph245->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph245->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph245->GetXaxis()->SetTitleFont(42);
   Graph_Graph245->GetYaxis()->SetLabelFont(42);
   Graph_Graph245->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph245->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph245->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph245->GetYaxis()->SetTickLength(0.02);
   Graph_Graph245->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph245->GetYaxis()->SetTitleFont(42);
   Graph_Graph245->GetZaxis()->SetLabelFont(42);
   Graph_Graph245->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph245->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph245->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph245->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph245->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph245);
   
   graph->Draw("p");
   
   Double_t Graph3_fx46[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy46[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 100, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx46,Graph3_fy46);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph346 = new TH1F("Graph_Graph346","Graph",100,0,3520);
   Graph_Graph346->SetMinimum(99.9);
   Graph_Graph346->SetMaximum(101.1);
   Graph_Graph346->SetDirectory(nullptr);
   Graph_Graph346->SetStats(0);
   Graph_Graph346->SetLineWidth(2);
   Graph_Graph346->SetMarkerStyle(20);
   Graph_Graph346->SetMarkerSize(0.9);
   Graph_Graph346->GetXaxis()->SetLabelFont(42);
   Graph_Graph346->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph346->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph346->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph346->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph346->GetXaxis()->SetTitleFont(42);
   Graph_Graph346->GetYaxis()->SetLabelFont(42);
   Graph_Graph346->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph346->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph346->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph346->GetYaxis()->SetTickLength(0.02);
   Graph_Graph346->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph346->GetYaxis()->SetTitleFont(42);
   Graph_Graph346->GetZaxis()->SetLabelFont(42);
   Graph_Graph346->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph346->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph346->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph346->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph346->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph346);
   
   graph->Draw("p");
   
   Double_t Graph4_fx47[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy47[35] = { 100, 100, 100, 100, 100, 100, 3.647757, 5.389774, 100, 100, 100, 6.518417, 100, 1.60315, 1.601553, 100, 100,
   6.518417, 3.003311, 2.921629, 3.020394, 2.769327, 4.174704, 3.280032, 3.07833, 2.594495, 2.696204, 2.511513, 2.980602, 3.118634, 3.321075, 3.325963, 3.418303,
   3.520346, 3.382123 };
   graph = new TGraph(35,Graph4_fx47,Graph4_fy47);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph447 = new TH1F("Graph_Graph447","Graph",100,0,3520);
   Graph_Graph447->SetMinimum(1.441398);
   Graph_Graph447->SetMaximum(109.8398);
   Graph_Graph447->SetDirectory(nullptr);
   Graph_Graph447->SetStats(0);
   Graph_Graph447->SetLineWidth(2);
   Graph_Graph447->SetMarkerStyle(20);
   Graph_Graph447->SetMarkerSize(0.9);
   Graph_Graph447->GetXaxis()->SetLabelFont(42);
   Graph_Graph447->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph447->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph447->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph447->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph447->GetXaxis()->SetTitleFont(42);
   Graph_Graph447->GetYaxis()->SetLabelFont(42);
   Graph_Graph447->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph447->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph447->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph447->GetYaxis()->SetTickLength(0.02);
   Graph_Graph447->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph447->GetYaxis()->SetTitleFont(42);
   Graph_Graph447->GetZaxis()->SetLabelFont(42);
   Graph_Graph447->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph447->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph447->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph447->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph447->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph447);
   
   graph->Draw("p");
   
   Double_t Graph5_fx48[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy48[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 100, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4507899, 0.6691337, 0.7289529, 0.2602696,
   0.5494237, 1.178157 };
   graph = new TGraph(35,Graph5_fx48,Graph5_fy48);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph548 = new TH1F("Graph_Graph548","Graph",100,0,3520);
   Graph_Graph548->SetMinimum(99.9);
   Graph_Graph548->SetMaximum(101.1);
   Graph_Graph548->SetDirectory(nullptr);
   Graph_Graph548->SetStats(0);
   Graph_Graph548->SetLineWidth(2);
   Graph_Graph548->SetMarkerStyle(20);
   Graph_Graph548->SetMarkerSize(0.9);
   Graph_Graph548->GetXaxis()->SetLabelFont(42);
   Graph_Graph548->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph548->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph548->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph548->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph548->GetXaxis()->SetTitleFont(42);
   Graph_Graph548->GetYaxis()->SetLabelFont(42);
   Graph_Graph548->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph548->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph548->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph548->GetYaxis()->SetTickLength(0.02);
   Graph_Graph548->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph548->GetYaxis()->SetTitleFont(42);
   Graph_Graph548->GetZaxis()->SetLabelFont(42);
   Graph_Graph548->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph548->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph548->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph548->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph548->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph548);
   
   graph->Draw("p");
   
   Double_t Graph6_fx49[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy49[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 54.25374, 119.6723, 223.6068, 223.6068, 223.6068, 300.5221, 223.6068, 92.77858, 124.427, 223.6068, 223.6068,
   112.7137, 70.82647, 57.5285, 86.91786, 100.6062, 417.4981, 30.05867, 61.72435, 21.68613, 39.63325, 14.67369, 18.0959, 26.66949, 23.96245, 21.45303, 24.48702,
   19.15582, 13.38714 };
   graph = new TGraph(35,Graph6_fx49,Graph6_fy49);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph649 = new TH1F("Graph_Graph649","Graph",100,0,3520);
   Graph_Graph649->SetMinimum(12.04843);
   Graph_Graph649->SetMaximum(457.9092);
   Graph_Graph649->SetDirectory(nullptr);
   Graph_Graph649->SetStats(0);
   Graph_Graph649->SetLineWidth(2);
   Graph_Graph649->SetMarkerStyle(20);
   Graph_Graph649->SetMarkerSize(0.9);
   Graph_Graph649->GetXaxis()->SetLabelFont(42);
   Graph_Graph649->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph649->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph649->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph649->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph649->GetXaxis()->SetTitleFont(42);
   Graph_Graph649->GetYaxis()->SetLabelFont(42);
   Graph_Graph649->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph649->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph649->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph649->GetYaxis()->SetTickLength(0.02);
   Graph_Graph649->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph649->GetYaxis()->SetTitleFont(42);
   Graph_Graph649->GetZaxis()->SetLabelFont(42);
   Graph_Graph649->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph649->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph649->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph649->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph649->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph649);
   
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
