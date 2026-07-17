#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:10:02 2026) by ROOT version 6.32.13
   TCanvas *c2 = new TCanvas("c2", "c2",0,0,800,600);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c2->SetHighLightColor(2);
   c2->Range(-418.2588,-0.3883565,3764.329,3.495208);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetLogy();
   c2->SetGridx();
   c2->SetGridy();
   c2->SetRightMargin(0.05);
   c2->SetTopMargin(0.05);
   c2->SetFrameLineWidth(2);
   c2->SetFrameBorderMode(0);
   c2->SetFrameLineWidth(2);
   c2->SetFrameBorderMode(0);
   
   Double_t Graph0_fx43[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy43[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 100, 13.45271, 46.67209, 0, 0, 12.40454, 14.23863,
   48.28161, 21.19664, 8.117038, 10.77319, 2.987784, 7.288551, 5.967599, 3.138518, 1.973623, 6.044436, 7.29422, 2.487862, 3.954828, 0.3438532, 1.237023, 2.172232,
   2.102029, 1.90028 };
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
   Graph_Graph043->GetXaxis()->SetLabelSize(22);
   Graph_Graph043->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph043->GetXaxis()->SetTitleOffset(1);
   Graph_Graph043->GetXaxis()->SetTitleFont(42);
   Graph_Graph043->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph043->GetYaxis()->SetLabelFont(43);
   Graph_Graph043->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph043->GetYaxis()->SetLabelSize(22);
   Graph_Graph043->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph043->GetYaxis()->SetTickLength(0.02);
   Graph_Graph043->GetYaxis()->SetTitleOffset(1);
   Graph_Graph043->GetYaxis()->SetTitleFont(42);
   Graph_Graph043->GetZaxis()->SetLabelFont(42);
   Graph_Graph043->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph043->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph043->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph043->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph043->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph043);
   
   graph->Draw("ap");
   
   TLegend *leg = new TLegend(0.12,0.75,0.5,0.93,NULL,"brNDC");
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
   Double_t Graph1_fy44[35] = { 100, 100, 100, 100, 100, 100, 58.67825, 145.0913, 0, 0, 100, 13.45271, 93.75925, 31.95499, 0, 12.40454, 39.36345,
   73.70476, 44.39882, 14.31395, 18.99794, 5.303204, 14.0043, 5.754673, 6.246102, 4.205, 1.7241, 9.094453, 3.544426, 3.439748, 0.6122351, 0.8095503, 2.405894,
   4.097617, 2.201647 };
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
   Graph_Graph144->SetMinimum(50.03695);
   Graph_Graph144->SetMaximum(153.7326);
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
   Double_t Graph2_fy45[35] = { 100, 100, 100, 100, 100, 100, 23.12768, 17.80117, 27.03399, 19.66195, 100, 27.84765, 48.57019, 33.45624, 33.74184, 56.52951, 24.63786,
   18.49145, 22.29505, 7.146239, 23.41332, 19.09516, 20.00351, 19.57861, 5.60354, 6.439328, 4.891312, 3.53936, 3.379619, 1.167953, 3.411764, 3.070009, 4.85518,
   0.5729675, 0.4724145 };
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
   Graph_Graph245->SetMinimum(0.425173);
   Graph_Graph245->SetMaximum(109.9528);
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
   Double_t Graph3_fy46[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 3.219008, 1.559162, 1.972556, 1.070833, 0.8127809, 1.530147, 1.482379,
   3.593898, 3.484374 };
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
   Double_t Graph4_fy47[35] = { 100, 100, 100, 100, 100, 100, 4.177231, 1.60315, 3.647757, 3.647757, 100, 2.66149, 0.9231091, 3.059781, 2.56598, 3.072083, 3.423321,
   3.260857, 3.049409, 2.706003, 2.812874, 2.779841, 3.06977, 3.356576, 2.901912, 2.729976, 3.320146, 2.899957, 3.16987, 3.332436, 3.415895, 3.427017, 3.465414,
   3.400147, 3.353786 };
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
   Graph_Graph447->SetMinimum(0.8307981);
   Graph_Graph447->SetMaximum(109.9077);
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
   Double_t Graph5_fy48[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 100, 0, 100, 0, 0, 11.73369, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 1.92737, 0, 0.9748459, 1.52958, 0.3517568, 0.6263018, 0.2931237, 0.2495468,
   0.4661083, 0.5084217 };
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
   Double_t Graph6_fy49[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 63.20978, 146.188, 27.27898, 19.99746, 223.6068, 33.83084, 115.4513, 46.36598, 33.83927, 59.26862, 48.69253,
   90.08925, 54.10102, 18.14295, 32.14156, 20.23372, 25.66724, 21.5248, 9.417263, 8.396126, 9.209671, 12.62072, 6.641784, 6.410028, 4.945922, 5.069156, 6.948829,
   6.78339, 5.663059 };
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
   Graph_Graph649->SetMinimum(4.45133);
   Graph_Graph649->SetMaximum(245.4729);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2000 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
