#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 15:56:21 2026) by ROOT version 6.32.13
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
   246.3426, 240.4939, 10.70721, 49.30299, 27.11611, 9.318453, 8.51034, 17.36291, 5.101001, 7.983696, 5.265927, 3.271616, 5.856371, 4.788423, 6.053758, 6.562841,
   16.53971, 1.058227 };
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
   277.9222, 66.10919, 33.84409, 17.33272, 51.27766, 13.50777, 28.08836, 19.72927, 3.826493, 6.405413, 2.785814, 1.565266, 2.530539, 5.451429, 2.276158, 5.141133,
   6.830633, 2.611566 };
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
   Graph_Graph19->SetMinimum(82.20778);
   Graph_Graph19->SetMaximum(295.7144);
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
   Double_t Graph2_fy10[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 117.3527, 100, 100, 34.07614, 16.765, 41.02626, 39.17492, 12.43657, 26.76177,
   54.63989, 32.79752, 35.94542, 15.01803, 30.61811, 12.39731, 22.55538, 17.21539, 14.79421, 35.93156, 29.17218, 20.78644, 19.47889, 15.97743, 16.37242, 23.70424,
   21.9975, 34.61355 };
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
   Graph_Graph210->SetMinimum(1.901766);
   Graph_Graph210->SetMaximum(127.8483);
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
   Double_t Graph4_fy12[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 1.234829, 100, 100, 1.338339, 2.225983, 1.77803, 2.830219, 1.881361, 3.112805,
   1.60315, 1.234829, 2.249932, 2.514625, 1.856244, 2.561271, 3.225696, 2.879691, 2.94137, 3.029561, 3.188014, 3.387237, 3.34903, 3.22268, 3.008211, 3.597534,
   3.762925, 3.225946 };
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
   0, 0, 0, 27.44176, 0, 0, 5.703026, 0, 0.4288495, 0.8215725, 0.8998156, 0.4499793, 0.2711654, 0.6191373, 0, 0,
   0, 5.778873 };
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
   Double_t Graph6_fy14[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 117.3592, 223.6068, 223.6068, 81.65987, 122.0527, 81.13349, 69.47221, 36.10648, 45.53938,
   375.3845, 251.565, 50.56879, 54.43411, 65.61704, 20.72551, 37.15551, 31.5496, 16.37628, 37.48365, 29.94445, 21.37061, 20.76882, 17.84127, 17.85873, 25.38376,
   28.6054, 34.87757 };
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
   Graph_Graph614->SetMinimum(14.73866);
   Graph_Graph614->SetMaximum(411.2853);
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
