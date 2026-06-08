#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
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
   
   Double_t Graph0_fx55[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy55[35] = { 100, 100, 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 100, 52.89886, 106.7396,
   98.19983, 47.69672, 84.89902, 28.23936, 42.81411, 13.98923, 19.52841, 1.721025, 19.17278, 13.24539, 10.09866, 9.273624, 7.510239, 10.98349, 3.042746, 2.145082,
   3.704768, 4.571211 };
   TGraph *graph = new TGraph(35,Graph0_fx55,Graph0_fy55);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph055 = new TH1F("Graph_Graph055","Graph",100,0,3520);
   Graph_Graph055->SetMinimum(1);
   Graph_Graph055->SetMaximum(2000);
   Graph_Graph055->SetDirectory(nullptr);
   Graph_Graph055->SetStats(0);
   Graph_Graph055->SetLineWidth(2);
   Graph_Graph055->SetMarkerStyle(20);
   Graph_Graph055->SetMarkerSize(0.9);
   Graph_Graph055->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph055->GetXaxis()->SetRange(1,101);
   Graph_Graph055->GetXaxis()->SetLabelFont(43);
   Graph_Graph055->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph055->GetXaxis()->SetLabelSize(16);
   Graph_Graph055->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph055->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph055->GetXaxis()->SetTitleFont(42);
   Graph_Graph055->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph055->GetYaxis()->SetLabelFont(43);
   Graph_Graph055->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph055->GetYaxis()->SetLabelSize(16);
   Graph_Graph055->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph055->GetYaxis()->SetTickLength(0.02);
   Graph_Graph055->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph055->GetYaxis()->SetTitleFont(42);
   Graph_Graph055->GetZaxis()->SetLabelFont(42);
   Graph_Graph055->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph055->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph055->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph055->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph055->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph055);
   
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
   
   Double_t Graph1_fx56[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy56[35] = { 100, 100, 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 100, 52.89886, 106.7396,
   45.40884, 22.05557, 100, 28.23936, 0, 0, 0, 5.960464e-06, 6.445265, 4.458916, 3.294551, 2.23242, 0.829643, 2.072048, 0.4469991, 0.3690183,
   0.7661462, 0.848186 };
   graph = new TGraph(35,Graph1_fx56,Graph1_fy56);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph156 = new TH1F("Graph_Graph156","Graph",100,0,3520);
   Graph_Graph156->SetMinimum(99.32604);
   Graph_Graph156->SetMaximum(107.4136);
   Graph_Graph156->SetDirectory(nullptr);
   Graph_Graph156->SetStats(0);
   Graph_Graph156->SetLineWidth(2);
   Graph_Graph156->SetMarkerStyle(20);
   Graph_Graph156->SetMarkerSize(0.9);
   Graph_Graph156->GetXaxis()->SetLabelFont(42);
   Graph_Graph156->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph156->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph156->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph156->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph156->GetXaxis()->SetTitleFont(42);
   Graph_Graph156->GetYaxis()->SetLabelFont(42);
   Graph_Graph156->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph156->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph156->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph156->GetYaxis()->SetTickLength(0.02);
   Graph_Graph156->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph156->GetYaxis()->SetTitleFont(42);
   Graph_Graph156->GetZaxis()->SetLabelFont(42);
   Graph_Graph156->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph156->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph156->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph156->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph156->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph156);
   
   graph->Draw("p");
   
   Double_t Graph2_fx57[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy57[35] = { 100, 100, 100, 100, 31.63878, 100, 31.46752, 100, 100, 32.87451, 23.15911, 30.76717, 428.0302, 100, 100, 20.24186, 19.08832,
   16.73798, 18.52283, 33.93324, 96.36926, 21.53542, 71.81394, 63.36883, 18.83631, 86.97018, 45.06075, 25.78949, 15.8388, 20.36961, 36.18748, 40.84605, 37.8088,
   35.89821, 38.90873 };
   graph = new TGraph(35,Graph2_fx57,Graph2_fy57);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph257 = new TH1F("Graph_Graph257","Graph",100,0,3520);
   Graph_Graph257->SetMinimum(14.25492);
   Graph_Graph257->SetMaximum(469.2494);
   Graph_Graph257->SetDirectory(nullptr);
   Graph_Graph257->SetStats(0);
   Graph_Graph257->SetLineWidth(2);
   Graph_Graph257->SetMarkerStyle(20);
   Graph_Graph257->SetMarkerSize(0.9);
   Graph_Graph257->GetXaxis()->SetLabelFont(42);
   Graph_Graph257->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph257->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph257->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph257->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph257->GetXaxis()->SetTitleFont(42);
   Graph_Graph257->GetYaxis()->SetLabelFont(42);
   Graph_Graph257->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph257->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph257->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph257->GetYaxis()->SetTickLength(0.02);
   Graph_Graph257->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph257->GetYaxis()->SetTitleFont(42);
   Graph_Graph257->GetZaxis()->SetLabelFont(42);
   Graph_Graph257->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph257->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph257->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph257->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph257->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph257);
   
   graph->Draw("p");
   
   Double_t Graph3_fx58[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy58[35] = { 100, 100, 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 100, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx58,Graph3_fy58);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph358 = new TH1F("Graph_Graph358","Graph",100,0,3520);
   Graph_Graph358->SetMinimum(99.9);
   Graph_Graph358->SetMaximum(101.1);
   Graph_Graph358->SetDirectory(nullptr);
   Graph_Graph358->SetStats(0);
   Graph_Graph358->SetLineWidth(2);
   Graph_Graph358->SetMarkerStyle(20);
   Graph_Graph358->SetMarkerSize(0.9);
   Graph_Graph358->GetXaxis()->SetLabelFont(42);
   Graph_Graph358->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph358->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph358->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph358->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph358->GetXaxis()->SetTitleFont(42);
   Graph_Graph358->GetYaxis()->SetLabelFont(42);
   Graph_Graph358->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph358->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph358->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph358->GetYaxis()->SetTickLength(0.02);
   Graph_Graph358->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph358->GetYaxis()->SetTitleFont(42);
   Graph_Graph358->GetZaxis()->SetLabelFont(42);
   Graph_Graph358->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph358->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph358->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph358->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph358->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph358);
   
   graph->Draw("p");
   
   Double_t Graph4_fx59[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy59[35] = { 100, 100, 100, 100, 2.281326, 100, 2.301472, 100, 100, 2.803653, 1.971108, 2.803653, 0.6861925, 100, 100, 2.803659, 2.803653,
   1.891363, 1.936394, 1.413596, 1.780581, 1.817054, 1.927984, 2.085721, 1.737118, 2.091223, 1.974171, 1.895225, 2.051723, 1.888388, 1.986456, 1.993006, 2.04832,
   1.995522, 1.98409 };
   graph = new TGraph(35,Graph4_fx59,Graph4_fy59);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph459 = new TH1F("Graph_Graph459","Graph",100,0,3520);
   Graph_Graph459->SetMinimum(0.6175733);
   Graph_Graph459->SetMaximum(109.9314);
   Graph_Graph459->SetDirectory(nullptr);
   Graph_Graph459->SetStats(0);
   Graph_Graph459->SetLineWidth(2);
   Graph_Graph459->SetMarkerStyle(20);
   Graph_Graph459->SetMarkerSize(0.9);
   Graph_Graph459->GetXaxis()->SetLabelFont(42);
   Graph_Graph459->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph459->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph459->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph459->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph459->GetXaxis()->SetTitleFont(42);
   Graph_Graph459->GetYaxis()->SetLabelFont(42);
   Graph_Graph459->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph459->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph459->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph459->GetYaxis()->SetTickLength(0.02);
   Graph_Graph459->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph459->GetYaxis()->SetTitleFont(42);
   Graph_Graph459->GetZaxis()->SetLabelFont(42);
   Graph_Graph459->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph459->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph459->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph459->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph459->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph459);
   
   graph->Draw("p");
   
   Double_t Graph5_fx60[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy60[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 31.72092, 223.6068, 31.55157, 223.6068, 223.6068, 32.99384, 23.24284, 30.89465, 428.0308, 223.6068, 223.6068, 77.5511, 152.1805,
   109.4939, 55.75186, 135.5039, 104.3318, 47.95961, 73.18919, 66.34243, 18.99437, 89.31585, 47.2196, 27.9558, 18.60271, 21.80778, 37.92638, 41.01012, 37.92675,
   36.15212, 39.23572 };
   graph = new TGraph(35,Graph5_fx60,Graph5_fy60);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph560 = new TH1F("Graph_Graph560","Graph",100,0,3520);
   Graph_Graph560->SetMinimum(16.74244);
   Graph_Graph560->SetMaximum(468.9736);
   Graph_Graph560->SetDirectory(nullptr);
   Graph_Graph560->SetStats(0);
   Graph_Graph560->SetLineWidth(2);
   Graph_Graph560->SetMarkerStyle(20);
   Graph_Graph560->SetMarkerSize(0.9);
   Graph_Graph560->GetXaxis()->SetLabelFont(42);
   Graph_Graph560->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph560->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph560->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph560->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph560->GetXaxis()->SetTitleFont(42);
   Graph_Graph560->GetYaxis()->SetLabelFont(42);
   Graph_Graph560->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph560->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph560->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph560->GetYaxis()->SetTickLength(0.02);
   Graph_Graph560->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph560->GetYaxis()->SetTitleFont(42);
   Graph_Graph560->GetZaxis()->SetLabelFont(42);
   Graph_Graph560->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph560->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph560->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph560->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph560->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph560);
   
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
