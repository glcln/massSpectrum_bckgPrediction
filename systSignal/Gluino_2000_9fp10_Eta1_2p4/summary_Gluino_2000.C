#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2000()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:15:05 2026) by ROOT version 6.32.13
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
   Double_t Graph0_fy43[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 100, 85.04453,
   100, 30.99655, 13.24387, 33.40651, 0, 0, 9.612846, 6.89683, 0, 15.14854, 5.476069, 0.6486058, 3.374881, 0.4763126, 0.8629441, 4.315102,
   0.8086145, 0.6863236 };
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
   Double_t Graph1_fy44[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 100, 100,
   100, 64.92587, 23.35485, 58.91057, 4.485452, 0, 6.291425, 13.72567, 2.267015, 16.83658, 5.779219, 1.985395, 5.665487, 0.9583592, 1.337659, 4.311812,
   6.373727, 1.108372 };
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
   Graph_Graph144->SetMinimum(99.9);
   Graph_Graph144->SetMaximum(101.1);
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
   Double_t Graph2_fy45[35] = { 100, 100, 100, 100, 100, 100, 11.54944, 38.8168, 100, 100, 100, 29.92773, 100, 52.13418, 43.56629, 34.58839, 45.27677,
   100, 17.10592, 17.92077, 21.91409, 24.98558, 32.05928, 16.74592, 21.87612, 8.903158, 12.80614, 2.582431, 6.042194, 1.905966, 4.069453, 4.214746, 6.06547,
   3.519785, 6.340623 };
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
   Graph_Graph245->SetMinimum(1.715369);
   Graph_Graph245->SetMaximum(109.8094);
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
   Double_t Graph3_fy46[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 0, 0,
   100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.106015, 4.965293, 1.88874, 1.0185, 2.252841, 2.868521,
   4.13025, 8.279782 };
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
   Double_t Graph4_fy47[35] = { 100, 100, 100, 100, 100, 100, 3.647757, 5.389774, 100, 100, 100, 6.518417, 100, 1.60315, 1.601553, 3.647768, 6.518417,
   100, 2.772772, 3.048432, 2.758741, 2.769327, 4.174704, 3.292251, 3.103638, 2.711165, 2.744567, 2.54004, 3.063333, 3.172314, 3.32576, 3.367984, 3.384197,
   3.518963, 3.404856 };
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
   Double_t Graph5_fy48[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 100, 0, 100, 0, 0, 0, 0,
   100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3681898, 0.8721828, 0.5570471, 0.2476573,
   0.5009592, 1.106244 };
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
   Double_t Graph6_fy49[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 12.1118, 39.18921, 223.6068, 223.6068, 223.6068, 30.62938, 223.6068, 52.15882, 43.59571, 145.6354, 139.0146,
   223.6068, 74.00308, 32.4237, 71.23404, 25.53561, 32.32994, 20.57313, 26.91018, 9.578935, 26.16255, 9.662651, 8.654984, 7.794305, 5.459286, 6.05943, 9.678995,
   9.116426, 11.04766 };
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
   Graph_Graph649->SetMinimum(4.913357);
   Graph_Graph649->SetMaximum(245.4215);
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
