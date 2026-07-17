#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy29[35] = { 100, 100, 100, 100, 100, 0, 100, 100, 100, 12.51066, 100, 0, 0, 0, 100, 70.32389, 243.2832,
   380.3062, 55.77448, 25.24694, 27.70723, 0, 27.63655, 78.61607, 8.280969, 14.48209, 4.119051, 6.791997, 8.955359, 0.6964624, 5.411565, 4.014146, 5.999083,
   2.735835, 2.612448 };
   TGraph *graph = new TGraph(35,Graph0_fx29,Graph0_fy29);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph029 = new TH1F("Graph_Graph029","Graph",100,0,3520);
   Graph_Graph029->SetMinimum(1);
   Graph_Graph029->SetMaximum(2000);
   Graph_Graph029->SetDirectory(nullptr);
   Graph_Graph029->SetStats(0);
   Graph_Graph029->SetLineWidth(2);
   Graph_Graph029->SetMarkerStyle(20);
   Graph_Graph029->SetMarkerSize(0.9);
   Graph_Graph029->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph029->GetXaxis()->SetRange(1,101);
   Graph_Graph029->GetXaxis()->SetLabelFont(43);
   Graph_Graph029->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph029->GetXaxis()->SetLabelSize(16);
   Graph_Graph029->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph029->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph029->GetXaxis()->SetTitleFont(42);
   Graph_Graph029->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph029->GetYaxis()->SetLabelFont(43);
   Graph_Graph029->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph029->GetYaxis()->SetLabelSize(16);
   Graph_Graph029->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph029->GetYaxis()->SetTickLength(0.02);
   Graph_Graph029->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph029->GetYaxis()->SetTitleFont(42);
   Graph_Graph029->GetZaxis()->SetLabelFont(42);
   Graph_Graph029->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph029->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph029->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph029->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph029->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph029);
   
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
   
   Double_t Graph1_fx30[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy30[35] = { 100, 100, 100, 100, 100, 0, 100, 100, 100, 100, 100, 0, 0, 0, 100, 37.36437, 129.2608,
   0, 0, 0, 0, 0, 0, 0, 12.08262, 5.384541, 4.231519, 3.999269, 2.467573, 0.3240228, 2.18153, 1.304543, 2.847373,
   1.489007, 1.325691 };
   graph = new TGraph(35,Graph1_fx30,Graph1_fy30);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph130 = new TH1F("Graph_Graph130","Graph",100,0,3520);
   Graph_Graph130->SetMinimum(97.07392);
   Graph_Graph130->SetMaximum(132.1869);
   Graph_Graph130->SetDirectory(nullptr);
   Graph_Graph130->SetStats(0);
   Graph_Graph130->SetLineWidth(2);
   Graph_Graph130->SetMarkerStyle(20);
   Graph_Graph130->SetMarkerSize(0.9);
   Graph_Graph130->GetXaxis()->SetLabelFont(42);
   Graph_Graph130->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph130->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph130->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph130->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph130->GetXaxis()->SetTitleFont(42);
   Graph_Graph130->GetYaxis()->SetLabelFont(42);
   Graph_Graph130->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph130->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph130->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph130->GetYaxis()->SetTickLength(0.02);
   Graph_Graph130->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph130->GetYaxis()->SetTitleFont(42);
   Graph_Graph130->GetZaxis()->SetLabelFont(42);
   Graph_Graph130->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph130->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph130->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph130->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph130->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph130);
   
   graph->Draw("p");
   
   Double_t Graph2_fx31[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy31[35] = { 100, 100, 100, 100, 100, 219.3418, 100, 100, 100, 15.68012, 100, 66.62581, 54.29229, 112.2731, 100, 46.37568, 44.56636,
   72.88062, 16.99251, 86.63124, 73.03217, 37.22375, 15.68468, 28.51801, 18.75391, 42.48054, 14.07829, 20.03422, 16.96368, 22.17284, 23.47646, 25.88817, 26.90433,
   15.38306, 29.04555 };
   graph = new TGraph(35,Graph2_fx31,Graph2_fy31);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph231 = new TH1F("Graph_Graph231","Graph",100,0,3520);
   Graph_Graph231->SetMinimum(12.67047);
   Graph_Graph231->SetMaximum(239.8682);
   Graph_Graph231->SetDirectory(nullptr);
   Graph_Graph231->SetStats(0);
   Graph_Graph231->SetLineWidth(2);
   Graph_Graph231->SetMarkerStyle(20);
   Graph_Graph231->SetMarkerSize(0.9);
   Graph_Graph231->GetXaxis()->SetLabelFont(42);
   Graph_Graph231->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph231->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph231->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph231->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph231->GetXaxis()->SetTitleFont(42);
   Graph_Graph231->GetYaxis()->SetLabelFont(42);
   Graph_Graph231->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph231->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph231->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph231->GetYaxis()->SetTickLength(0.02);
   Graph_Graph231->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph231->GetYaxis()->SetTitleFont(42);
   Graph_Graph231->GetZaxis()->SetLabelFont(42);
   Graph_Graph231->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph231->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph231->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph231->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph231->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph231);
   
   graph->Draw("p");
   
   Double_t Graph3_fx32[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy32[35] = { 100, 100, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 0, 0, 100, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx32,Graph3_fy32);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph332 = new TH1F("Graph_Graph332","Graph",100,0,3520);
   Graph_Graph332->SetMinimum(99.9);
   Graph_Graph332->SetMaximum(101.1);
   Graph_Graph332->SetDirectory(nullptr);
   Graph_Graph332->SetStats(0);
   Graph_Graph332->SetLineWidth(2);
   Graph_Graph332->SetMarkerStyle(20);
   Graph_Graph332->SetMarkerSize(0.9);
   Graph_Graph332->GetXaxis()->SetLabelFont(42);
   Graph_Graph332->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph332->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph332->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph332->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph332->GetXaxis()->SetTitleFont(42);
   Graph_Graph332->GetYaxis()->SetLabelFont(42);
   Graph_Graph332->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph332->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph332->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph332->GetYaxis()->SetTickLength(0.02);
   Graph_Graph332->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph332->GetYaxis()->SetTitleFont(42);
   Graph_Graph332->GetZaxis()->SetLabelFont(42);
   Graph_Graph332->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph332->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph332->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph332->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph332->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph332);
   
   graph->Draw("p");
   
   Double_t Graph4_fx33[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy33[35] = { 100, 100, 100, 100, 100, 1.234829, 100, 100, 100, 3.396153, 100, 4.187942, 3.647768, 3.495669, 100, 1.239753, 0.9231031,
   2.042806, 2.672142, 2.387643, 2.248597, 1.535571, 2.71616, 0.9230971, 2.601242, 2.36249, 2.564597, 2.907431, 3.132427, 3.350723, 3.333926, 3.395915, 3.224194,
   3.340137, 3.692949 };
   graph = new TGraph(35,Graph4_fx33,Graph4_fy33);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph433 = new TH1F("Graph_Graph433","Graph",100,0,3520);
   Graph_Graph433->SetMinimum(0.8307874);
   Graph_Graph433->SetMaximum(109.9077);
   Graph_Graph433->SetDirectory(nullptr);
   Graph_Graph433->SetStats(0);
   Graph_Graph433->SetLineWidth(2);
   Graph_Graph433->SetMarkerStyle(20);
   Graph_Graph433->SetMarkerSize(0.9);
   Graph_Graph433->GetXaxis()->SetLabelFont(42);
   Graph_Graph433->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph433->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph433->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph433->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph433->GetXaxis()->SetTitleFont(42);
   Graph_Graph433->GetYaxis()->SetLabelFont(42);
   Graph_Graph433->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph433->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph433->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph433->GetYaxis()->SetTickLength(0.02);
   Graph_Graph433->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph433->GetYaxis()->SetTitleFont(42);
   Graph_Graph433->GetZaxis()->SetLabelFont(42);
   Graph_Graph433->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph433->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph433->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph433->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph433->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph433);
   
   graph->Draw("p");
   
   Double_t Graph5_fx34[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy34[35] = { 100, 100, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 0, 0, 100, 0, 0,
   0, 14.16489, 0, 0, 0, 0, 0, 0, 0, 2.386403, 0.605619, 0.5776346, 0.4109383, 0.8392811, 0.8486092, 0.8183241,
   0, 0 };
   graph = new TGraph(35,Graph5_fx34,Graph5_fy34);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph534 = new TH1F("Graph_Graph534","Graph",100,0,3520);
   Graph_Graph534->SetMinimum(99.9);
   Graph_Graph534->SetMaximum(101.1);
   Graph_Graph534->SetDirectory(nullptr);
   Graph_Graph534->SetStats(0);
   Graph_Graph534->SetLineWidth(2);
   Graph_Graph534->SetMarkerStyle(20);
   Graph_Graph534->SetMarkerSize(0.9);
   Graph_Graph534->GetXaxis()->SetLabelFont(42);
   Graph_Graph534->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph534->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph534->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph534->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph534->GetXaxis()->SetTitleFont(42);
   Graph_Graph534->GetYaxis()->SetLabelFont(42);
   Graph_Graph534->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph534->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph534->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph534->GetYaxis()->SetTickLength(0.02);
   Graph_Graph534->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph534->GetYaxis()->SetTitleFont(42);
   Graph_Graph534->GetZaxis()->SetLabelFont(42);
   Graph_Graph534->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph534->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph534->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph534->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph534->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph534);
   
   graph->Draw("p");
   
   Double_t Graph6_fx35[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy35[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 219.3453, 223.6068, 223.6068, 223.6068, 102.0486, 223.6068, 66.7573, 54.4147, 112.3275, 223.6068, 92.16174, 279.0736,
   387.2319, 58.36676, 90.26672, 78.14375, 37.25541, 31.89303, 83.63382, 23.93825, 45.2648, 15.48056, 21.72437, 19.5925, 22.43774, 24.41932, 26.44891, 27.89866,
   16.04672, 29.42557 };
   graph = new TGraph(35,Graph6_fx35,Graph6_fy35);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph635 = new TH1F("Graph_Graph635","Graph",100,0,3520);
   Graph_Graph635->SetMinimum(13.93251);
   Graph_Graph635->SetMaximum(424.407);
   Graph_Graph635->SetDirectory(nullptr);
   Graph_Graph635->SetStats(0);
   Graph_Graph635->SetLineWidth(2);
   Graph_Graph635->SetMarkerStyle(20);
   Graph_Graph635->SetMarkerSize(0.9);
   Graph_Graph635->GetXaxis()->SetLabelFont(42);
   Graph_Graph635->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph635->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph635->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph635->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph635->GetXaxis()->SetTitleFont(42);
   Graph_Graph635->GetYaxis()->SetLabelFont(42);
   Graph_Graph635->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph635->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph635->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph635->GetYaxis()->SetTickLength(0.02);
   Graph_Graph635->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph635->GetYaxis()->SetTitleFont(42);
   Graph_Graph635->GetZaxis()->SetLabelFont(42);
   Graph_Graph635->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph635->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph635->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph635->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph635->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph635);
   
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
