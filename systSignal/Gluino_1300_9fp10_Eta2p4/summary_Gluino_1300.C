#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1300()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:54 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx13[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy13[35] = { 100, 100, 0, 0, 55.84391, 100, 0, 0, 0, 0, 0, 78.26144, 21.45518, 32.96331, 21.32539, 26.56662, 20.60658,
   59.81266, 16.94951, 1.587713, 34.60788, 25.98128, 0.5107284, 11.8367, 21.17282, 6.08707, 11.0462, 7.59095, 0.8882344, 6.322992, 7.838649, 7.997775, 5.489361,
   5.438662, 2.545696 };
   TGraph *graph = new TGraph(35,Graph0_fx13,Graph0_fy13);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph013 = new TH1F("Graph_Graph013","Graph",100,0,3520);
   Graph_Graph013->SetMinimum(1);
   Graph_Graph013->SetMaximum(2000);
   Graph_Graph013->SetDirectory(nullptr);
   Graph_Graph013->SetStats(0);
   Graph_Graph013->SetLineWidth(2);
   Graph_Graph013->SetMarkerStyle(20);
   Graph_Graph013->SetMarkerSize(0.9);
   Graph_Graph013->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph013->GetXaxis()->SetRange(1,101);
   Graph_Graph013->GetXaxis()->SetLabelFont(43);
   Graph_Graph013->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetXaxis()->SetLabelSize(16);
   Graph_Graph013->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph013->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetXaxis()->SetTitleFont(42);
   Graph_Graph013->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph013->GetYaxis()->SetLabelFont(43);
   Graph_Graph013->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetYaxis()->SetLabelSize(16);
   Graph_Graph013->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph013->GetYaxis()->SetTickLength(0.02);
   Graph_Graph013->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetYaxis()->SetTitleFont(42);
   Graph_Graph013->GetZaxis()->SetLabelFont(42);
   Graph_Graph013->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph013->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph013->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph013);
   
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
   
   Double_t Graph1_fx14[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy14[35] = { 100, 100, 0, 100, 55.84391, 100, 48.46771, 0, 100, 25.97064, 98.60365, 171.414, 35.99737, 17.38531, 18.77546, 26.08832, 13.03653,
   51.43906, 13.89275, 1.915896, 1.344907, 16.23132, 5.408704, 5.634689, 13.52805, 2.664852, 3.289366, 1.493615, 0.6036401, 1.740474, 2.997851, 1.115304, 2.477223,
   5.279458, 0 };
   graph = new TGraph(35,Graph1_fx14,Graph1_fy14);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph114 = new TH1F("Graph_Graph114","Graph",100,0,3520);
   Graph_Graph114->SetMinimum(92.8586);
   Graph_Graph114->SetMaximum(178.5555);
   Graph_Graph114->SetDirectory(nullptr);
   Graph_Graph114->SetStats(0);
   Graph_Graph114->SetLineWidth(2);
   Graph_Graph114->SetMarkerStyle(20);
   Graph_Graph114->SetMarkerSize(0.9);
   Graph_Graph114->GetXaxis()->SetLabelFont(42);
   Graph_Graph114->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetXaxis()->SetTitleFont(42);
   Graph_Graph114->GetYaxis()->SetLabelFont(42);
   Graph_Graph114->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetYaxis()->SetTickLength(0.02);
   Graph_Graph114->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetYaxis()->SetTitleFont(42);
   Graph_Graph114->GetZaxis()->SetLabelFont(42);
   Graph_Graph114->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph114);
   
   graph->Draw("p");
   
   Double_t Graph2_fx15[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy15[35] = { 100, 100, 27.44938, 27.44938, 249.2065, 100, 30.59461, 25.18971, 23.59124, 103.8535, 20.16381, 30.59461, 74.18942, 68.55835, 158.6791, 26.91041, 57.4389,
   25.59432, 29.73971, 33.16046, 15.77933, 72.57782, 20.83315, 23.36727, 38.94018, 57.70451, 39.5405, 38.31015, 40.54686, 34.463, 39.76637, 32.21499, 24.68441,
   17.57034, 17.60458 };
   graph = new TGraph(35,Graph2_fx15,Graph2_fy15);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph215 = new TH1F("Graph_Graph215","Graph",100,0,3520);
   Graph_Graph215->SetMinimum(14.2014);
   Graph_Graph215->SetMaximum(272.5492);
   Graph_Graph215->SetDirectory(nullptr);
   Graph_Graph215->SetStats(0);
   Graph_Graph215->SetLineWidth(2);
   Graph_Graph215->SetMarkerStyle(20);
   Graph_Graph215->SetMarkerSize(0.9);
   Graph_Graph215->GetXaxis()->SetLabelFont(42);
   Graph_Graph215->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetXaxis()->SetTitleFont(42);
   Graph_Graph215->GetYaxis()->SetLabelFont(42);
   Graph_Graph215->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetYaxis()->SetTickLength(0.02);
   Graph_Graph215->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetYaxis()->SetTitleFont(42);
   Graph_Graph215->GetZaxis()->SetLabelFont(42);
   Graph_Graph215->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph215);
   
   graph->Draw("p");
   
   Double_t Graph3_fx16[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy16[35] = { 100, 100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx16,Graph3_fy16);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph316 = new TH1F("Graph_Graph316","Graph",100,0,3520);
   Graph_Graph316->SetMinimum(99.9);
   Graph_Graph316->SetMaximum(101.1);
   Graph_Graph316->SetDirectory(nullptr);
   Graph_Graph316->SetStats(0);
   Graph_Graph316->SetLineWidth(2);
   Graph_Graph316->SetMarkerStyle(20);
   Graph_Graph316->SetMarkerSize(0.9);
   Graph_Graph316->GetXaxis()->SetLabelFont(42);
   Graph_Graph316->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetXaxis()->SetTitleFont(42);
   Graph_Graph316->GetYaxis()->SetLabelFont(42);
   Graph_Graph316->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetYaxis()->SetTickLength(0.02);
   Graph_Graph316->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetYaxis()->SetTitleFont(42);
   Graph_Graph316->GetZaxis()->SetLabelFont(42);
   Graph_Graph316->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph316);
   
   graph->Draw("p");
   
   Double_t Graph4_fx17[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy17[35] = { 100, 100, 3.598678, 2.803653, 2.684188, 100, 2.298272, 3.117466, 2.803665, 2.820194, 2.803653, 0.9149969, 1.653767, 2.384239, 2.20083, 2.423215, 2.159798,
   0.9355307, 1.523256, 2.369523, 2.252066, 2.212524, 1.858038, 2.014667, 1.828051, 1.90112, 1.86249, 1.961488, 1.944184, 1.932955, 1.894611, 1.898593, 2.122265,
   1.997524, 1.645178 };
   graph = new TGraph(35,Graph4_fx17,Graph4_fy17);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph417 = new TH1F("Graph_Graph417","Graph",100,0,3520);
   Graph_Graph417->SetMinimum(0.8234972);
   Graph_Graph417->SetMaximum(109.9085);
   Graph_Graph417->SetDirectory(nullptr);
   Graph_Graph417->SetStats(0);
   Graph_Graph417->SetLineWidth(2);
   Graph_Graph417->SetMarkerStyle(20);
   Graph_Graph417->SetMarkerSize(0.9);
   Graph_Graph417->GetXaxis()->SetLabelFont(42);
   Graph_Graph417->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetXaxis()->SetTitleFont(42);
   Graph_Graph417->GetYaxis()->SetLabelFont(42);
   Graph_Graph417->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetYaxis()->SetTickLength(0.02);
   Graph_Graph417->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetYaxis()->SetTitleFont(42);
   Graph_Graph417->GetZaxis()->SetLabelFont(42);
   Graph_Graph417->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph417);
   
   graph->Draw("p");
   
   Double_t Graph5_fx18[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy18[35] = { 223.6068, 223.6068, 27.68427, 103.7368, 261.4348, 223.6068, 57.36228, 25.38188, 102.7833, 107.0886, 100.6833, 190.9044, 85.22289, 78.06895, 161.2178, 46.00465, 62.43776,
   82.9426, 36.97385, 33.338, 38.12575, 78.80938, 21.60989, 26.86904, 46.37861, 58.11694, 41.22813, 39.1327, 40.60765, 35.13466, 40.68643, 33.26588, 25.49693,
   19.23951, 17.8636 };
   graph = new TGraph(35,Graph5_fx18,Graph5_fy18);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph518 = new TH1F("Graph_Graph518","Graph",100,0,3520);
   Graph_Graph518->SetMinimum(16.07724);
   Graph_Graph518->SetMaximum(285.7919);
   Graph_Graph518->SetDirectory(nullptr);
   Graph_Graph518->SetStats(0);
   Graph_Graph518->SetLineWidth(2);
   Graph_Graph518->SetMarkerStyle(20);
   Graph_Graph518->SetMarkerSize(0.9);
   Graph_Graph518->GetXaxis()->SetLabelFont(42);
   Graph_Graph518->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetXaxis()->SetTitleFont(42);
   Graph_Graph518->GetYaxis()->SetLabelFont(42);
   Graph_Graph518->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetYaxis()->SetTickLength(0.02);
   Graph_Graph518->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetYaxis()->SetTitleFont(42);
   Graph_Graph518->GetZaxis()->SetLabelFont(42);
   Graph_Graph518->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph518);
   
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
