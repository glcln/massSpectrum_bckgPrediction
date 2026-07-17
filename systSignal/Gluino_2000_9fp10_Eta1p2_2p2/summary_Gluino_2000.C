#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 15:56:22 2026) by ROOT version 6.32.13
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
   0, 25.83271, 14.05017, 28.61946, 3.995663, 207.3519, 9.943104, 7.565403, 4.883921, 11.00954, 1.235759, 9.239328, 9.439158, 3.193927, 4.300952, 6.413752,
   9.094917, 1.576996 };
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
   85.04453, 26.34826, 14.33058, 28.61946, 3.995663, 207.3519, 0, 2.613986, 5.171323, 5.411858, 2.115643, 1.570553, 2.387935, 0.7734835, 1.236057, 1.872981,
   0.4754543, 0.8886695 };
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
   73.68394, 60.38054, 53.83511, 76.8599, 100.4092, 297.1491, 28.17622, 57.03634, 20.84193, 32.19904, 16.14411, 14.88223, 23.83621, 22.63818, 20.73389, 22.89636,
   19.72129, 12.28118 };
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
   Graph_Graph245->SetMinimum(11.05306);
   Graph_Graph245->SetMaximum(329.2684);
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
   6.518417, 3.003311, 2.921629, 3.020394, 2.769327, 4.174704, 3.280032, 2.966976, 2.594495, 2.701223, 2.452672, 2.903056, 3.12959, 3.33488, 3.353512, 3.397572,
   3.487396, 3.355396 };
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
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.4761398, 0.6091297, 0.7662773, 0.2726436,
   0.5912721, 1.265687 };
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
   112.7137, 70.82647, 57.5285, 86.91786, 100.6062, 417.4981, 30.05867, 57.67162, 22.17459, 34.56259, 16.51214, 17.82527, 25.9376, 23.11727, 21.47478, 24.09214,
   22.0008, 12.85934 };
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
   Graph_Graph649->SetMinimum(11.57341);
   Graph_Graph649->SetMaximum(457.962);
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
