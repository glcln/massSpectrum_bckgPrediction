#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:55:15 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx8[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy8[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 52.46599, 59.05666, 50.1147, 57.30361, 33.84479, 36.71448,
   246.3426, 240.4939, 15.45682, 49.30299, 20.69393, 14.0201, 15.74587, 24.5326, 11.72053, 12.37498, 5.844766, 3.786266, 6.947428, 4.042804, 3.695714, 5.859923,
   11.07467, 2.022791 };
   TGraph *graph = new TGraph(35,Graph0_fx8,Graph0_fy8);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph08 = new TH1F("Graph_Graph08","Graph",100,0,3520);
   Graph_Graph08->SetMinimum(1);
   Graph_Graph08->SetMaximum(2000);
   Graph_Graph08->SetDirectory(nullptr);
   Graph_Graph08->SetStats(0);
   Graph_Graph08->SetLineWidth(2);
   Graph_Graph08->SetMarkerStyle(20);
   Graph_Graph08->SetMarkerSize(0.9);
   Graph_Graph08->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph08->GetXaxis()->SetRange(1,101);
   Graph_Graph08->GetXaxis()->SetLabelFont(43);
   Graph_Graph08->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetXaxis()->SetLabelSize(16);
   Graph_Graph08->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph08->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph08->GetXaxis()->SetTitleFont(42);
   Graph_Graph08->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph08->GetYaxis()->SetLabelFont(43);
   Graph_Graph08->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetYaxis()->SetLabelSize(16);
   Graph_Graph08->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph08->GetYaxis()->SetTickLength(0.02);
   Graph_Graph08->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph08->GetYaxis()->SetTitleFont(42);
   Graph_Graph08->GetZaxis()->SetLabelFont(42);
   Graph_Graph08->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph08->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph08->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph08->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph08);
   
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
   
   Double_t Graph1_fx9[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy9[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 52.46599, 105.4664, 48.83487, 0, 0, 0,
   0, 0, 22.14085, 17.33272, 20.69392, 14.0201, 0, 1.909816, 4.637999, 9.960771, 3.362703, 1.80428, 3.390622, 3.69103, 2.002919, 4.590487,
   6.894541, 2.022791 };
   graph = new TGraph(35,Graph1_fx9,Graph1_fy9);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph19 = new TH1F("Graph_Graph19","Graph",100,0,3520);
   Graph_Graph19->SetMinimum(99.45336);
   Graph_Graph19->SetMaximum(106.013);
   Graph_Graph19->SetDirectory(nullptr);
   Graph_Graph19->SetStats(0);
   Graph_Graph19->SetLineWidth(2);
   Graph_Graph19->SetMarkerStyle(20);
   Graph_Graph19->SetMarkerSize(0.9);
   Graph_Graph19->GetXaxis()->SetLabelFont(42);
   Graph_Graph19->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetXaxis()->SetTitleFont(42);
   Graph_Graph19->GetYaxis()->SetLabelFont(42);
   Graph_Graph19->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetYaxis()->SetTickLength(0.02);
   Graph_Graph19->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetYaxis()->SetTitleFont(42);
   Graph_Graph19->GetZaxis()->SetLabelFont(42);
   Graph_Graph19->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph19);
   
   graph->Draw("p");
   
   Double_t Graph2_fx10[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy10[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 61.44222, 100, 100, 34.07614, 16.765, 41.02626, 39.17492, 12.43657, 26.76177,
   54.63989, 32.79752, 19.14533, 15.01803, 22.20265, 18.57994, 12.37167, 19.50402, 15.84203, 35.47984, 26.37888, 22.24848, 22.34533, 16.15682, 12.34866, 22.36709,
   28.84579, 36.24296 };
   graph = new TGraph(35,Graph2_fx10,Graph2_fy10);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph210 = new TH1F("Graph_Graph210","Graph",100,0,3520);
   Graph_Graph210->SetMinimum(3.58353);
   Graph_Graph210->SetMaximum(108.7651);
   Graph_Graph210->SetDirectory(nullptr);
   Graph_Graph210->SetStats(0);
   Graph_Graph210->SetLineWidth(2);
   Graph_Graph210->SetMarkerStyle(20);
   Graph_Graph210->SetMarkerSize(0.9);
   Graph_Graph210->GetXaxis()->SetLabelFont(42);
   Graph_Graph210->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetXaxis()->SetTitleFont(42);
   Graph_Graph210->GetYaxis()->SetLabelFont(42);
   Graph_Graph210->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetYaxis()->SetTickLength(0.02);
   Graph_Graph210->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetYaxis()->SetTitleFont(42);
   Graph_Graph210->GetZaxis()->SetLabelFont(42);
   Graph_Graph210->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph210);
   
   graph->Draw("p");
   
   Double_t Graph3_fx11[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy11[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx11,Graph3_fy11);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph311 = new TH1F("Graph_Graph311","Graph",100,0,3520);
   Graph_Graph311->SetMinimum(99.9);
   Graph_Graph311->SetMaximum(101.1);
   Graph_Graph311->SetDirectory(nullptr);
   Graph_Graph311->SetStats(0);
   Graph_Graph311->SetLineWidth(2);
   Graph_Graph311->SetMarkerStyle(20);
   Graph_Graph311->SetMarkerSize(0.9);
   Graph_Graph311->GetXaxis()->SetLabelFont(42);
   Graph_Graph311->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetXaxis()->SetTitleFont(42);
   Graph_Graph311->GetYaxis()->SetLabelFont(42);
   Graph_Graph311->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetYaxis()->SetTickLength(0.02);
   Graph_Graph311->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetYaxis()->SetTitleFont(42);
   Graph_Graph311->GetZaxis()->SetLabelFont(42);
   Graph_Graph311->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph311);
   
   graph->Draw("p");
   
   Double_t Graph4_fx12[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy12[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 2.645302, 100, 100, 1.338339, 2.225983, 1.77803, 2.830219, 1.881361, 3.112805,
   1.60315, 1.234829, 1.507318, 2.514625, 1.900434, 2.711868, 2.677894, 3.001988, 2.73037, 3.021538, 3.184497, 3.377211, 3.336453, 3.302419, 3.068519, 3.60291,
   4.001307, 3.315735 };
   graph = new TGraph(35,Graph4_fx12,Graph4_fy12);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph412 = new TH1F("Graph_Graph412","Graph",100,0,3520);
   Graph_Graph412->SetMinimum(1.111346);
   Graph_Graph412->SetMaximum(109.8765);
   Graph_Graph412->SetDirectory(nullptr);
   Graph_Graph412->SetStats(0);
   Graph_Graph412->SetLineWidth(2);
   Graph_Graph412->SetMarkerStyle(20);
   Graph_Graph412->SetMarkerSize(0.9);
   Graph_Graph412->GetXaxis()->SetLabelFont(42);
   Graph_Graph412->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetXaxis()->SetTitleFont(42);
   Graph_Graph412->GetYaxis()->SetLabelFont(42);
   Graph_Graph412->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetYaxis()->SetTickLength(0.02);
   Graph_Graph412->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetYaxis()->SetTitleFont(42);
   Graph_Graph412->GetZaxis()->SetLabelFont(42);
   Graph_Graph412->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph412);
   
   graph->Draw("p");
   
   Double_t Graph5_fx13[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy13[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 0, 0, 0, 0, 18.78004, 0,
   0, 0, 0, 27.44176, 0, 0, 8.20663, 0, 1.263231, 0.7664502, 0.7620811, 0.5866349, 0.2463341, 0.5689502, 0, 0,
   0, 5.480587 };
   graph = new TGraph(35,Graph5_fx13,Graph5_fy13);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph513 = new TH1F("Graph_Graph513","Graph",100,0,3520);
   Graph_Graph513->SetMinimum(99.9);
   Graph_Graph513->SetMaximum(101.1);
   Graph_Graph513->SetDirectory(nullptr);
   Graph_Graph513->SetStats(0);
   Graph_Graph513->SetLineWidth(2);
   Graph_Graph513->SetMarkerStyle(20);
   Graph_Graph513->SetMarkerSize(0.9);
   Graph_Graph513->GetXaxis()->SetLabelFont(42);
   Graph_Graph513->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetXaxis()->SetTitleFont(42);
   Graph_Graph513->GetYaxis()->SetLabelFont(42);
   Graph_Graph513->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetYaxis()->SetTickLength(0.02);
   Graph_Graph513->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetYaxis()->SetTitleFont(42);
   Graph_Graph513->GetZaxis()->SetLabelFont(42);
   Graph_Graph513->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph513);
   
   graph->Draw("p");
   
   Double_t Graph6_fx14[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy14[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 61.49914, 223.6068, 223.6068, 81.65987, 122.0527, 81.13349, 69.47221, 36.10648, 45.53938,
   252.3347, 242.7231, 33.13527, 54.43411, 36.78378, 27.30742, 20.20301, 31.54227, 20.42809, 38.9911, 27.41269, 22.89087, 23.87904, 17.37575, 13.40057, 23.84699,
   31.91039, 36.50657 };
   graph = new TGraph(35,Graph6_fx14,Graph6_fy14);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph614 = new TH1F("Graph_Graph614","Graph",100,0,3520);
   Graph_Graph614->SetMinimum(12.06051);
   Graph_Graph614->SetMaximum(276.2281);
   Graph_Graph614->SetDirectory(nullptr);
   Graph_Graph614->SetStats(0);
   Graph_Graph614->SetLineWidth(2);
   Graph_Graph614->SetMarkerStyle(20);
   Graph_Graph614->SetMarkerSize(0.9);
   Graph_Graph614->GetXaxis()->SetLabelFont(42);
   Graph_Graph614->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetXaxis()->SetTitleFont(42);
   Graph_Graph614->GetYaxis()->SetLabelFont(42);
   Graph_Graph614->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetYaxis()->SetTickLength(0.02);
   Graph_Graph614->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetYaxis()->SetTitleFont(42);
   Graph_Graph614->GetZaxis()->SetLabelFont(42);
   Graph_Graph614->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph614);
   
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
