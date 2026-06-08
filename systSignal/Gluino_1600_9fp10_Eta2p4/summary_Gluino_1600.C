#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx25[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy25[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 100, 30.33572, 25.13007, 26.07455, 5.753839, 56.92727, 29.81358, 50.21745,
   0, 53.80991, 20.3342, 0, 14.9345, 5.962479, 9.997833, 3.718829, 8.74486, 7.540727, 6.585455, 8.988833, 1.928836, 5.105031, 5.477917, 6.415105,
   3.726828, 2.866697 };
   TGraph *graph = new TGraph(35,Graph0_fx25,Graph0_fy25);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph025 = new TH1F("Graph_Graph025","Graph",100,0,3520);
   Graph_Graph025->SetMinimum(1);
   Graph_Graph025->SetMaximum(2000);
   Graph_Graph025->SetDirectory(nullptr);
   Graph_Graph025->SetStats(0);
   Graph_Graph025->SetLineWidth(2);
   Graph_Graph025->SetMarkerStyle(20);
   Graph_Graph025->SetMarkerSize(0.9);
   Graph_Graph025->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph025->GetXaxis()->SetRange(1,101);
   Graph_Graph025->GetXaxis()->SetLabelFont(43);
   Graph_Graph025->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetXaxis()->SetLabelSize(16);
   Graph_Graph025->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph025->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetXaxis()->SetTitleFont(42);
   Graph_Graph025->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph025->GetYaxis()->SetLabelFont(43);
   Graph_Graph025->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetYaxis()->SetLabelSize(16);
   Graph_Graph025->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph025->GetYaxis()->SetTickLength(0.02);
   Graph_Graph025->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetYaxis()->SetTitleFont(42);
   Graph_Graph025->GetZaxis()->SetLabelFont(42);
   Graph_Graph025->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph025->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph025->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph025);
   
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
   
   Double_t Graph1_fx26[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy26[35] = { 100, 100, 100, 100, 0, 46.98205, 0, 100, 100, 100, 30.33572, 26.04973, 26.07455, 5.753839, 56.92727, 0, 0,
   1.192093e-05, 0, 0, 0, 0, 6.417519, 5.610728, 0, 5.824447, 3.722018, 1.800489, 1.503444, 0.4695177, 1.170695, 1.324761, 1.90177,
   1.29354, 0.8566141 };
   graph = new TGraph(35,Graph1_fx26,Graph1_fy26);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph126 = new TH1F("Graph_Graph126","Graph",100,0,3520);
   Graph_Graph126->SetMinimum(99.9);
   Graph_Graph126->SetMaximum(101.1);
   Graph_Graph126->SetDirectory(nullptr);
   Graph_Graph126->SetStats(0);
   Graph_Graph126->SetLineWidth(2);
   Graph_Graph126->SetMarkerStyle(20);
   Graph_Graph126->SetMarkerSize(0.9);
   Graph_Graph126->GetXaxis()->SetLabelFont(42);
   Graph_Graph126->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetXaxis()->SetTitleFont(42);
   Graph_Graph126->GetYaxis()->SetLabelFont(42);
   Graph_Graph126->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetYaxis()->SetTickLength(0.02);
   Graph_Graph126->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetYaxis()->SetTitleFont(42);
   Graph_Graph126->GetZaxis()->SetLabelFont(42);
   Graph_Graph126->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph126);
   
   graph->Draw("p");
   
   Double_t Graph2_fx27[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy27[35] = { 100, 100, 25.3091, 100, 19.08832, 19.72026, 36.23116, 20.9009, 100, 100, 120.0113, 26.73502, 79.04928, 133.3809, 19.71855, 26.7654, 24.054,
   23.15606, 23.6035, 22.3958, 79.21188, 100.886, 26.92083, 53.83383, 81.8122, 24.05903, 56.57682, 26.06591, 34.21172, 38.92252, 38.67911, 38.08707, 46.05551,
   44.3951, 39.56836 };
   graph = new TGraph(35,Graph2_fx27,Graph2_fy27);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph227 = new TH1F("Graph_Graph227","Graph",100,0,3520);
   Graph_Graph227->SetMinimum(7.659054);
   Graph_Graph227->SetMaximum(144.8102);
   Graph_Graph227->SetDirectory(nullptr);
   Graph_Graph227->SetStats(0);
   Graph_Graph227->SetLineWidth(2);
   Graph_Graph227->SetMarkerStyle(20);
   Graph_Graph227->SetMarkerSize(0.9);
   Graph_Graph227->GetXaxis()->SetLabelFont(42);
   Graph_Graph227->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetXaxis()->SetTitleFont(42);
   Graph_Graph227->GetYaxis()->SetLabelFont(42);
   Graph_Graph227->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetYaxis()->SetTickLength(0.02);
   Graph_Graph227->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetYaxis()->SetTitleFont(42);
   Graph_Graph227->GetZaxis()->SetLabelFont(42);
   Graph_Graph227->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph227);
   
   graph->Draw("p");
   
   Double_t Graph3_fx28[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy28[35] = { 100, 100, 0, 100, 0, 0, 0, 0, 100, 100, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx28,Graph3_fy28);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph328 = new TH1F("Graph_Graph328","Graph",100,0,3520);
   Graph_Graph328->SetMinimum(99.9);
   Graph_Graph328->SetMaximum(101.1);
   Graph_Graph328->SetDirectory(nullptr);
   Graph_Graph328->SetStats(0);
   Graph_Graph328->SetLineWidth(2);
   Graph_Graph328->SetMarkerStyle(20);
   Graph_Graph328->SetMarkerSize(0.9);
   Graph_Graph328->GetXaxis()->SetLabelFont(42);
   Graph_Graph328->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetXaxis()->SetTitleFont(42);
   Graph_Graph328->GetYaxis()->SetLabelFont(42);
   Graph_Graph328->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetYaxis()->SetTickLength(0.02);
   Graph_Graph328->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetYaxis()->SetTitleFont(42);
   Graph_Graph328->GetZaxis()->SetLabelFont(42);
   Graph_Graph328->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph328);
   
   graph->Draw("p");
   
   Double_t Graph4_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy29[35] = { 100, 100, 1.413608, 100, 2.803653, 1.60706, 0.9150028, 1.284808, 100, 100, 2.483177, 2.473193, 2.424741, 1.622707, 0.8163929, 1.533878, 0.7258177,
   1.667655, 1.926917, 1.986784, 1.72264, 1.595861, 2.173954, 1.595509, 1.772642, 1.999146, 1.789516, 2.100056, 1.922148, 1.949906, 1.972973, 1.955563, 1.877099,
   1.956928, 1.935953 };
   graph = new TGraph(35,Graph4_fx29,Graph4_fy29);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph429 = new TH1F("Graph_Graph429","Graph",100,0,3520);
   Graph_Graph429->SetMinimum(0.6532359);
   Graph_Graph429->SetMaximum(109.9274);
   Graph_Graph429->SetDirectory(nullptr);
   Graph_Graph429->SetStats(0);
   Graph_Graph429->SetLineWidth(2);
   Graph_Graph429->SetMarkerStyle(20);
   Graph_Graph429->SetMarkerSize(0.9);
   Graph_Graph429->GetXaxis()->SetLabelFont(42);
   Graph_Graph429->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetXaxis()->SetTitleFont(42);
   Graph_Graph429->GetYaxis()->SetLabelFont(42);
   Graph_Graph429->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetYaxis()->SetTickLength(0.02);
   Graph_Graph429->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetYaxis()->SetTitleFont(42);
   Graph_Graph429->GetZaxis()->SetLabelFont(42);
   Graph_Graph429->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph429);
   
   graph->Draw("p");
   
   Double_t Graph5_fx30[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy30[35] = { 223.6068, 223.6068, 103.1627, 223.6068, 19.29312, 50.97827, 36.24271, 142.9633, 223.6068, 223.6068, 127.4731, 45.06647, 87.26071, 133.6388, 82.89098, 40.09475, 55.68585,
   23.21603, 58.79068, 30.315, 79.23061, 101.9979, 28.39354, 55.06418, 81.91586, 26.32927, 57.22635, 27.02687, 35.45696, 39.02186, 39.08194, 38.55142, 46.57686,
   44.61297, 39.72852 };
   graph = new TGraph(35,Graph5_fx30,Graph5_fy30);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph530 = new TH1F("Graph_Graph530","Graph",100,0,3520);
   Graph_Graph530->SetMinimum(17.3638);
   Graph_Graph530->SetMaximum(244.0382);
   Graph_Graph530->SetDirectory(nullptr);
   Graph_Graph530->SetStats(0);
   Graph_Graph530->SetLineWidth(2);
   Graph_Graph530->SetMarkerStyle(20);
   Graph_Graph530->SetMarkerSize(0.9);
   Graph_Graph530->GetXaxis()->SetLabelFont(42);
   Graph_Graph530->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetXaxis()->SetTitleFont(42);
   Graph_Graph530->GetYaxis()->SetLabelFont(42);
   Graph_Graph530->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetYaxis()->SetTickLength(0.02);
   Graph_Graph530->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetYaxis()->SetTitleFont(42);
   Graph_Graph530->GetZaxis()->SetLabelFont(42);
   Graph_Graph530->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph530);
   
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
