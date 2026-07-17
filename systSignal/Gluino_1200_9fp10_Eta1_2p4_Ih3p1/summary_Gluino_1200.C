#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:54:11 2026) by ROOT version 6.32.13
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
   246.3426, 240.4939, 13.67091, 49.30299, 29.95057, 12.06601, 9.746957, 17.49303, 5.772018, 11.02198, 5.780005, 3.795528, 7.095689, 3.790474, 3.599215, 7.746351,
   10.60201, 2.022791 };
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
   Double_t Graph1_fy9[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 39.6133, 100, 100, 52.46599, 105.4664, 48.83487, 0, 0, 0,
   0, 66.10919, 11.55415, 17.33272, 21.35723, 12.06601, 9.914022, 19.53444, 5.653268, 9.158527, 3.196883, 2.025712, 2.982771, 4.696834, 1.096535, 2.704072,
   2.332306, 4.499567 };
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
   Graph_Graph19->SetMinimum(33.02799);
   Graph_Graph19->SetMaximum(112.0517);
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
   Double_t Graph2_fy10[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 35.65569, 100, 100, 34.07614, 16.765, 41.02626, 39.17492, 12.43657, 26.76177,
   54.63989, 32.79752, 12.82049, 15.01803, 34.97789, 19.01549, 11.50043, 23.91474, 15.38808, 36.40426, 26.53611, 22.3756, 21.34161, 14.54406, 12.61861, 22.36709,
   32.17374, 36.24296 };
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
   Graph_Graph210->SetMinimum(2.650473);
   Graph_Graph210->SetMaximum(108.85);
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
   Double_t Graph4_fy12[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 3.042412, 100, 100, 1.338339, 2.225983, 1.77803, 2.830219, 1.881361, 3.112805,
   1.60315, 1.234829, 1.75463, 2.514625, 2.157664, 2.842307, 3.082871, 2.940285, 2.799284, 3.033137, 3.186285, 3.382194, 3.352535, 3.256011, 3.083646, 3.60291,
   3.963315, 3.315735 };
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
   0, 0, 0, 27.44176, 0, 0, 6.531734, 0, 0.9921372, 0.7078648, 0.7442892, 0.5780399, 0.2381563, 0.5334437, 0, 0,
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
   Double_t Graph6_fy14[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 53.3835, 223.6068, 223.6068, 81.65987, 122.0527, 81.13349, 69.47221, 36.10648, 45.53938,
   252.3347, 251.565, 22.08701, 54.43411, 50.80626, 25.70691, 18.30451, 35.61125, 17.60411, 39.24071, 27.53082, 23.0351, 22.93359, 16.07978, 13.52387, 24.09534,
   34.18625, 36.72716 };
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
   Graph_Graph614->SetMinimum(12.17148);
   Graph_Graph614->SetMaximum(276.2157);
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
