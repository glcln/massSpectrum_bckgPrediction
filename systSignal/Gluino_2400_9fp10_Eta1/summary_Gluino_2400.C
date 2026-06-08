#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2400()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:22 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx49[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy49[35] = { 100, 100, 100, 100, 100, 100, 0, 100, 105.7256, 0, 0, 100, 100, 59.71383, 30.79019, 50.97281, 100,
   179.9476, 35.95353, 0, 100, 43.62598, 12.14254, 31.61217, 7.828069, 17.45448, 3.965575, 7.03339, 7.545435, 8.14302, 7.774282, 1.279265, 2.279103,
   5.136824, 4.811823 };
   TGraph *graph = new TGraph(35,Graph0_fx49,Graph0_fy49);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph049 = new TH1F("Graph_Graph049","Graph",100,0,3520);
   Graph_Graph049->SetMinimum(1);
   Graph_Graph049->SetMaximum(2000);
   Graph_Graph049->SetDirectory(nullptr);
   Graph_Graph049->SetStats(0);
   Graph_Graph049->SetLineWidth(2);
   Graph_Graph049->SetMarkerStyle(20);
   Graph_Graph049->SetMarkerSize(0.9);
   Graph_Graph049->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph049->GetXaxis()->SetRange(1,101);
   Graph_Graph049->GetXaxis()->SetLabelFont(43);
   Graph_Graph049->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetXaxis()->SetLabelSize(16);
   Graph_Graph049->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph049->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetXaxis()->SetTitleFont(42);
   Graph_Graph049->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph049->GetYaxis()->SetLabelFont(43);
   Graph_Graph049->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetYaxis()->SetLabelSize(16);
   Graph_Graph049->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph049->GetYaxis()->SetTickLength(0.02);
   Graph_Graph049->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetYaxis()->SetTitleFont(42);
   Graph_Graph049->GetZaxis()->SetLabelFont(42);
   Graph_Graph049->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph049->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph049->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph049);
   
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
   
   Double_t Graph1_fx50[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy50[35] = { 100, 100, 100, 100, 100, 100, 0, 100, 105.7256, 100, 0, 100, 100, 5.960464e-06, 0, 49.02719, 1.757115,
   87.07249, 17.77995, 26.03616, 53.10911, 20.75535, 12.08375, 1.192093e-05, 0, 6.647968, 4.475933, 0, 3.867561, 1.107156, 0.9936929, 0.2384484, 0.3165126,
   0.7883191, 2.065492 };
   graph = new TGraph(35,Graph1_fx50,Graph1_fy50);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph150 = new TH1F("Graph_Graph150","Graph",100,0,3520);
   Graph_Graph150->SetMinimum(99.42744);
   Graph_Graph150->SetMaximum(106.2982);
   Graph_Graph150->SetDirectory(nullptr);
   Graph_Graph150->SetStats(0);
   Graph_Graph150->SetLineWidth(2);
   Graph_Graph150->SetMarkerStyle(20);
   Graph_Graph150->SetMarkerSize(0.9);
   Graph_Graph150->GetXaxis()->SetLabelFont(42);
   Graph_Graph150->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetXaxis()->SetTitleFont(42);
   Graph_Graph150->GetYaxis()->SetLabelFont(42);
   Graph_Graph150->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetYaxis()->SetTickLength(0.02);
   Graph_Graph150->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetYaxis()->SetTitleFont(42);
   Graph_Graph150->GetZaxis()->SetLabelFont(42);
   Graph_Graph150->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph150);
   
   graph->Draw("p");
   
   Double_t Graph2_fx51[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy51[35] = { 100, 100, 100, 100, 100, 100, 27.44938, 16.41645, 23.95937, 23.59124, 20.92748, 100, 100, 22.55688, 19.7679, 22.81093, 20.65356,
   34.19884, 24.19728, 27.15669, 32.39196, 25.55862, 33.61632, 27.44234, 53.32633, 22.1107, 33.51657, 41.34746, 15.42134, 30.00105, 34.68246, 39.06856, 37.22148,
   40.44384, 33.99102 };
   graph = new TGraph(35,Graph2_fx51,Graph2_fy51);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph251 = new TH1F("Graph_Graph251","Graph",100,0,3520);
   Graph_Graph251->SetMinimum(6.963477);
   Graph_Graph251->SetMaximum(108.4579);
   Graph_Graph251->SetDirectory(nullptr);
   Graph_Graph251->SetStats(0);
   Graph_Graph251->SetLineWidth(2);
   Graph_Graph251->SetMarkerStyle(20);
   Graph_Graph251->SetMarkerSize(0.9);
   Graph_Graph251->GetXaxis()->SetLabelFont(42);
   Graph_Graph251->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetXaxis()->SetTitleFont(42);
   Graph_Graph251->GetYaxis()->SetLabelFont(42);
   Graph_Graph251->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetYaxis()->SetTickLength(0.02);
   Graph_Graph251->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetYaxis()->SetTitleFont(42);
   Graph_Graph251->GetZaxis()->SetLabelFont(42);
   Graph_Graph251->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph251);
   
   graph->Draw("p");
   
   Double_t Graph3_fx52[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy52[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 100, 100, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx52,Graph3_fy52);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph352 = new TH1F("Graph_Graph352","Graph",100,0,3520);
   Graph_Graph352->SetMinimum(99.9);
   Graph_Graph352->SetMaximum(101.1);
   Graph_Graph352->SetDirectory(nullptr);
   Graph_Graph352->SetStats(0);
   Graph_Graph352->SetLineWidth(2);
   Graph_Graph352->SetMarkerStyle(20);
   Graph_Graph352->SetMarkerSize(0.9);
   Graph_Graph352->GetXaxis()->SetLabelFont(42);
   Graph_Graph352->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetXaxis()->SetTitleFont(42);
   Graph_Graph352->GetYaxis()->SetLabelFont(42);
   Graph_Graph352->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetYaxis()->SetTickLength(0.02);
   Graph_Graph352->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetYaxis()->SetTitleFont(42);
   Graph_Graph352->GetZaxis()->SetLabelFont(42);
   Graph_Graph352->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph352);
   
   graph->Draw("p");
   
   Double_t Graph4_fx53[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy53[35] = { 100, 100, 100, 100, 100, 100, 3.598678, 2.28132, 2.803665, 0.9150028, 1.133335, 100, 100, 1.921612, 1.537561, 0.7252336, 1.779002,
   0.9149969, 0.9474993, 1.889336, 0.9282708, 1.960611, 2.761793, 2.197134, 2.330357, 1.840973, 1.862407, 1.716542, 1.868939, 1.897043, 2.049547, 1.981139, 2.007389,
   1.942319, 1.994169 };
   graph = new TGraph(35,Graph4_fx53,Graph4_fy53);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph453 = new TH1F("Graph_Graph453","Graph",100,0,3520);
   Graph_Graph453->SetMinimum(0.6527102);
   Graph_Graph453->SetMaximum(109.9275);
   Graph_Graph453->SetDirectory(nullptr);
   Graph_Graph453->SetStats(0);
   Graph_Graph453->SetLineWidth(2);
   Graph_Graph453->SetMarkerStyle(20);
   Graph_Graph453->SetMarkerSize(0.9);
   Graph_Graph453->GetXaxis()->SetLabelFont(42);
   Graph_Graph453->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetXaxis()->SetTitleFont(42);
   Graph_Graph453->GetYaxis()->SetLabelFont(42);
   Graph_Graph453->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetYaxis()->SetTickLength(0.02);
   Graph_Graph453->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetYaxis()->SetTitleFont(42);
   Graph_Graph453->GetZaxis()->SetLabelFont(42);
   Graph_Graph453->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph453);
   
   graph->Draw("p");
   
   Double_t Graph5_fx54[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy54[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 27.68427, 142.3893, 151.4521, 102.7491, 20.95814, 223.6068, 223.6068, 63.86115, 36.62198, 74.31525, 102.1412,
   202.8131, 46.85284, 37.66877, 117.7738, 54.69093, 37.83044, 41.91944, 53.94819, 29.00218, 34.09676, 41.97651, 17.69753, 31.16402, 35.61602, 39.1404, 37.34652,
   40.8226, 34.44976 };
   graph = new TGraph(35,Graph5_fx54,Graph5_fy54);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph554 = new TH1F("Graph_Graph554","Graph",100,0,3520);
   Graph_Graph554->SetMinimum(15.92777);
   Graph_Graph554->SetMaximum(244.1977);
   Graph_Graph554->SetDirectory(nullptr);
   Graph_Graph554->SetStats(0);
   Graph_Graph554->SetLineWidth(2);
   Graph_Graph554->SetMarkerStyle(20);
   Graph_Graph554->SetMarkerSize(0.9);
   Graph_Graph554->GetXaxis()->SetLabelFont(42);
   Graph_Graph554->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetXaxis()->SetTitleFont(42);
   Graph_Graph554->GetYaxis()->SetLabelFont(42);
   Graph_Graph554->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetYaxis()->SetTickLength(0.02);
   Graph_Graph554->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetYaxis()->SetTitleFont(42);
   Graph_Graph554->GetZaxis()->SetLabelFont(42);
   Graph_Graph554->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph554);
   
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
