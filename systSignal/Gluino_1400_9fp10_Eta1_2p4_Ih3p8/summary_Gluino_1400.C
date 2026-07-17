#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1400()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:55:18 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx22[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy22[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 0, 32.80989, 29.49659, 65.90914, 100, 53.09324, 29.47217, 28.80791,
   41.86884, 7.858026, 0, 20.01693, 59.72072, 26.95533, 10.34371, 1.540804, 7.800257, 12.39359, 6.396341, 3.33963, 3.338677, 6.047976, 10.08471, 3.835082,
   6.443608, 6.174707 };
   TGraph *graph = new TGraph(35,Graph0_fx22,Graph0_fy22);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph022 = new TH1F("Graph_Graph022","Graph",100,0,3520);
   Graph_Graph022->SetMinimum(1);
   Graph_Graph022->SetMaximum(2000);
   Graph_Graph022->SetDirectory(nullptr);
   Graph_Graph022->SetStats(0);
   Graph_Graph022->SetLineWidth(2);
   Graph_Graph022->SetMarkerStyle(20);
   Graph_Graph022->SetMarkerSize(0.9);
   Graph_Graph022->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph022->GetXaxis()->SetRange(1,101);
   Graph_Graph022->GetXaxis()->SetLabelFont(43);
   Graph_Graph022->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph022->GetXaxis()->SetLabelSize(16);
   Graph_Graph022->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph022->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph022->GetXaxis()->SetTitleFont(42);
   Graph_Graph022->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph022->GetYaxis()->SetLabelFont(43);
   Graph_Graph022->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph022->GetYaxis()->SetLabelSize(16);
   Graph_Graph022->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph022->GetYaxis()->SetTickLength(0.02);
   Graph_Graph022->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph022->GetYaxis()->SetTitleFont(42);
   Graph_Graph022->GetZaxis()->SetLabelFont(42);
   Graph_Graph022->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph022->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph022->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph022->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph022->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph022);
   
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
   
   Double_t Graph1_fx23[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy23[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 0, 32.80989, 29.49659, 65.90914, 100, 0, 0, 5.960464e-06,
   0, 0, 0, 0, 0, 8.276617, 4.898304, 3.849518, 2.584004, 4.056227, 3.06555, 2.217335, 1.063466, 2.764255, 6.240678, 2.08202,
   2.96191, 2.838308 };
   graph = new TGraph(35,Graph1_fx23,Graph1_fy23);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph123 = new TH1F("Graph_Graph123","Graph",100,0,3520);
   Graph_Graph123->SetMinimum(99.9);
   Graph_Graph123->SetMaximum(101.1);
   Graph_Graph123->SetDirectory(nullptr);
   Graph_Graph123->SetStats(0);
   Graph_Graph123->SetLineWidth(2);
   Graph_Graph123->SetMarkerStyle(20);
   Graph_Graph123->SetMarkerSize(0.9);
   Graph_Graph123->GetXaxis()->SetLabelFont(42);
   Graph_Graph123->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph123->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph123->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph123->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph123->GetXaxis()->SetTitleFont(42);
   Graph_Graph123->GetYaxis()->SetLabelFont(42);
   Graph_Graph123->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph123->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph123->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph123->GetYaxis()->SetTickLength(0.02);
   Graph_Graph123->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph123->GetYaxis()->SetTitleFont(42);
   Graph_Graph123->GetZaxis()->SetLabelFont(42);
   Graph_Graph123->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph123->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph123->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph123->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph123->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph123);
   
   graph->Draw("p");
   
   Double_t Graph2_fx24[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy24[35] = { 100, 100, 100, 100, 100, 100, 68.75135, 65.38313, 100, 298.7019, 110.633, 34.51376, 17.61767, 100, 36.47068, 42.41244, 46.57565,
   175.0072, 30.45465, 52.62339, 53.38706, 37.93631, 42.76069, 32.93658, 13.43863, 37.78808, 20.18869, 23.15565, 26.61533, 23.81127, 28.70398, 27.7131, 31.66667,
   18.87351, 17.87533 };
   graph = new TGraph(35,Graph2_fx24,Graph2_fy24);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph224 = new TH1F("Graph_Graph224","Graph",100,0,3520);
   Graph_Graph224->SetMinimum(12.09477);
   Graph_Graph224->SetMaximum(327.2282);
   Graph_Graph224->SetDirectory(nullptr);
   Graph_Graph224->SetStats(0);
   Graph_Graph224->SetLineWidth(2);
   Graph_Graph224->SetMarkerStyle(20);
   Graph_Graph224->SetMarkerSize(0.9);
   Graph_Graph224->GetXaxis()->SetLabelFont(42);
   Graph_Graph224->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph224->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph224->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph224->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph224->GetXaxis()->SetTitleFont(42);
   Graph_Graph224->GetYaxis()->SetLabelFont(42);
   Graph_Graph224->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph224->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph224->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph224->GetYaxis()->SetTickLength(0.02);
   Graph_Graph224->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph224->GetYaxis()->SetTitleFont(42);
   Graph_Graph224->GetZaxis()->SetLabelFont(42);
   Graph_Graph224->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph224->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph224->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph224->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph224->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph224);
   
   graph->Draw("p");
   
   Double_t Graph3_fx25[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy25[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx25,Graph3_fy25);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph325 = new TH1F("Graph_Graph325","Graph",100,0,3520);
   Graph_Graph325->SetMinimum(99.9);
   Graph_Graph325->SetMaximum(101.1);
   Graph_Graph325->SetDirectory(nullptr);
   Graph_Graph325->SetStats(0);
   Graph_Graph325->SetLineWidth(2);
   Graph_Graph325->SetMarkerStyle(20);
   Graph_Graph325->SetMarkerSize(0.9);
   Graph_Graph325->GetXaxis()->SetLabelFont(42);
   Graph_Graph325->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph325->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph325->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph325->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph325->GetXaxis()->SetTitleFont(42);
   Graph_Graph325->GetYaxis()->SetLabelFont(42);
   Graph_Graph325->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph325->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph325->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph325->GetYaxis()->SetTickLength(0.02);
   Graph_Graph325->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph325->GetYaxis()->SetTitleFont(42);
   Graph_Graph325->GetZaxis()->SetLabelFont(42);
   Graph_Graph325->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph325->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph325->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph325->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph325->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph325);
   
   graph->Draw("p");
   
   Double_t Graph4_fx26[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy26[35] = { 100, 100, 100, 100, 100, 100, 4.550165, 1.234829, 100, 3.647768, 1.443058, 2.042806, 1.234829, 100, 3.362823, 3.669822, 1.432204,
   1.207829, 3.011644, 2.023339, 2.554774, 2.73034, 3.512108, 2.370024, 3.409523, 2.648747, 2.895927, 3.100145, 3.250444, 3.383219, 3.429985, 3.365648, 3.32917,
   3.302121, 3.339517 };
   graph = new TGraph(35,Graph4_fx26,Graph4_fy26);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph426 = new TH1F("Graph_Graph426","Graph",100,0,3520);
   Graph_Graph426->SetMinimum(1.087046);
   Graph_Graph426->SetMaximum(109.8792);
   Graph_Graph426->SetDirectory(nullptr);
   Graph_Graph426->SetStats(0);
   Graph_Graph426->SetLineWidth(2);
   Graph_Graph426->SetMarkerStyle(20);
   Graph_Graph426->SetMarkerSize(0.9);
   Graph_Graph426->GetXaxis()->SetLabelFont(42);
   Graph_Graph426->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph426->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph426->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph426->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph426->GetXaxis()->SetTitleFont(42);
   Graph_Graph426->GetYaxis()->SetLabelFont(42);
   Graph_Graph426->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph426->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph426->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph426->GetYaxis()->SetTickLength(0.02);
   Graph_Graph426->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph426->GetYaxis()->SetTitleFont(42);
   Graph_Graph426->GetZaxis()->SetLabelFont(42);
   Graph_Graph426->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph426->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph426->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph426->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph426->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph426);
   
   graph->Draw("p");
   
   Double_t Graph5_fx27[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy27[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 59.61923, 40.11477, 0, 0, 0, 0, 0, 1.414526, 0.3943443, 0.4651606, 0.6243825, 0.8776069, 1.440668, 0,
   2.039236, 1.789707 };
   graph = new TGraph(35,Graph5_fx27,Graph5_fy27);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph527 = new TH1F("Graph_Graph527","Graph",100,0,3520);
   Graph_Graph527->SetMinimum(99.9);
   Graph_Graph527->SetMaximum(101.1);
   Graph_Graph527->SetDirectory(nullptr);
   Graph_Graph527->SetStats(0);
   Graph_Graph527->SetLineWidth(2);
   Graph_Graph527->SetMarkerStyle(20);
   Graph_Graph527->SetMarkerSize(0.9);
   Graph_Graph527->GetXaxis()->SetLabelFont(42);
   Graph_Graph527->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph527->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph527->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph527->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph527->GetXaxis()->SetTitleFont(42);
   Graph_Graph527->GetYaxis()->SetLabelFont(42);
   Graph_Graph527->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph527->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph527->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph527->GetYaxis()->SetTickLength(0.02);
   Graph_Graph527->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph527->GetYaxis()->SetTitleFont(42);
   Graph_Graph527->GetZaxis()->SetLabelFont(42);
   Graph_Graph527->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph527->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph527->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph527->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph527->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph527);
   
   graph->Draw("p");
   
   Double_t Graph6_fx28[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy28[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 68.90176, 65.39479, 223.6068, 298.7242, 119.978, 54.17998, 94.868, 223.6068, 64.50047, 51.77732, 54.78356,
   179.9499, 31.59595, 52.66228, 57.07348, 70.80383, 51.34106, 34.94884, 14.47116, 38.76179, 24.20793, 24.41528, 27.11109, 24.30433, 29.66315, 30.33135, 32.13882,
   20.43052, 19.41295 };
   graph = new TGraph(35,Graph6_fx28,Graph6_fy28);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph628 = new TH1F("Graph_Graph628","Graph",100,0,3520);
   Graph_Graph628->SetMinimum(13.02404);
   Graph_Graph628->SetMaximum(327.1495);
   Graph_Graph628->SetDirectory(nullptr);
   Graph_Graph628->SetStats(0);
   Graph_Graph628->SetLineWidth(2);
   Graph_Graph628->SetMarkerStyle(20);
   Graph_Graph628->SetMarkerSize(0.9);
   Graph_Graph628->GetXaxis()->SetLabelFont(42);
   Graph_Graph628->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph628->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph628->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph628->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph628->GetXaxis()->SetTitleFont(42);
   Graph_Graph628->GetYaxis()->SetLabelFont(42);
   Graph_Graph628->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph628->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph628->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph628->GetYaxis()->SetTickLength(0.02);
   Graph_Graph628->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph628->GetYaxis()->SetTitleFont(42);
   Graph_Graph628->GetZaxis()->SetLabelFont(42);
   Graph_Graph628->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph628->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph628->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph628->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph628->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph628);
   
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
