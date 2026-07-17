#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1800()
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
   
   Double_t Graph0_fx36[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy36[35] = { 100, 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 0, 100, 0, 31.49567, 12.4284, 100,
   0, 1.192093e-05, 260.1717, 71.23677, 48.88728, 7.29419, 5.222201, 12.50305, 4.351997, 10.21537, 10.37201, 2.624571, 3.14672, 1.747531, 3.431714, 1.572382,
   0.5474806, 4.69451 };
   TGraph *graph = new TGraph(35,Graph0_fx36,Graph0_fy36);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph036 = new TH1F("Graph_Graph036","Graph",100,0,3520);
   Graph_Graph036->SetMinimum(1);
   Graph_Graph036->SetMaximum(2000);
   Graph_Graph036->SetDirectory(nullptr);
   Graph_Graph036->SetStats(0);
   Graph_Graph036->SetLineWidth(2);
   Graph_Graph036->SetMarkerStyle(20);
   Graph_Graph036->SetMarkerSize(0.9);
   Graph_Graph036->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph036->GetXaxis()->SetRange(1,101);
   Graph_Graph036->GetXaxis()->SetLabelFont(43);
   Graph_Graph036->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetXaxis()->SetLabelSize(22);
   Graph_Graph036->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph036->GetXaxis()->SetTitleOffset(1);
   Graph_Graph036->GetXaxis()->SetTitleFont(42);
   Graph_Graph036->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph036->GetYaxis()->SetLabelFont(43);
   Graph_Graph036->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetYaxis()->SetLabelSize(22);
   Graph_Graph036->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph036->GetYaxis()->SetTickLength(0.02);
   Graph_Graph036->GetYaxis()->SetTitleOffset(1);
   Graph_Graph036->GetYaxis()->SetTitleFont(42);
   Graph_Graph036->GetZaxis()->SetLabelFont(42);
   Graph_Graph036->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph036->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph036->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph036->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph036);
   
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
   
   Double_t Graph1_fx37[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy37[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 100, 0, 100, 100, 95.77658, 12.4284, 100,
   33.56739, 36.10858, 234.4341, 65.23187, 48.88728, 26.59252, 40.92343, 7.058001, 26.93885, 24.44131, 11.54083, 3.116131, 4.917932, 2.096581, 4.20382, 3.110033,
   1.712728, 6.942731 };
   graph = new TGraph(35,Graph1_fx37,Graph1_fy37);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph137 = new TH1F("Graph_Graph137","Graph",100,0,3520);
   Graph_Graph137->SetMinimum(86.55659);
   Graph_Graph137->SetMaximum(247.8775);
   Graph_Graph137->SetDirectory(nullptr);
   Graph_Graph137->SetStats(0);
   Graph_Graph137->SetLineWidth(2);
   Graph_Graph137->SetMarkerStyle(20);
   Graph_Graph137->SetMarkerSize(0.9);
   Graph_Graph137->GetXaxis()->SetLabelFont(42);
   Graph_Graph137->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetXaxis()->SetTitleFont(42);
   Graph_Graph137->GetYaxis()->SetLabelFont(42);
   Graph_Graph137->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetYaxis()->SetTickLength(0.02);
   Graph_Graph137->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetYaxis()->SetTitleFont(42);
   Graph_Graph137->GetZaxis()->SetLabelFont(42);
   Graph_Graph137->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph137);
   
   graph->Draw("p");
   
   Double_t Graph2_fx38[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy38[35] = { 100, 100, 54.12139, 100, 100, 100, 100, 100, 28.89853, 100, 100, 26.18185, 100, 19.20915, 26.78412, 16.44696, 100,
   25.57343, 14.44529, 29.0803, 125.548, 27.53997, 13.1533, 22.99118, 19.64142, 22.36671, 20.10467, 3.703415, 5.340743, 3.217328, 5.130434, 3.832471, 4.420328,
   6.614077, 13.53555 };
   graph = new TGraph(35,Graph2_fx38,Graph2_fy38);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph238 = new TH1F("Graph_Graph238","Graph",100,0,3520);
   Graph_Graph238->SetMinimum(2.895595);
   Graph_Graph238->SetMaximum(137.7811);
   Graph_Graph238->SetDirectory(nullptr);
   Graph_Graph238->SetStats(0);
   Graph_Graph238->SetLineWidth(2);
   Graph_Graph238->SetMarkerStyle(20);
   Graph_Graph238->SetMarkerSize(0.9);
   Graph_Graph238->GetXaxis()->SetLabelFont(42);
   Graph_Graph238->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetXaxis()->SetTitleFont(42);
   Graph_Graph238->GetYaxis()->SetLabelFont(42);
   Graph_Graph238->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetYaxis()->SetTickLength(0.02);
   Graph_Graph238->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetYaxis()->SetTitleFont(42);
   Graph_Graph238->GetZaxis()->SetLabelFont(42);
   Graph_Graph238->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph238);
   
   graph->Draw("p");
   
   Double_t Graph3_fx39[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy39[35] = { 100, 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 0, 100, 0, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 3.647971, 2.334797, 2.144241, 1.756752, 2.08326, 1.910609, 5.922556,
   5.031866, 6.125373 };
   graph = new TGraph(35,Graph3_fx39,Graph3_fy39);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph339 = new TH1F("Graph_Graph339","Graph",100,0,3520);
   Graph_Graph339->SetMinimum(99.9);
   Graph_Graph339->SetMaximum(101.1);
   Graph_Graph339->SetDirectory(nullptr);
   Graph_Graph339->SetStats(0);
   Graph_Graph339->SetLineWidth(2);
   Graph_Graph339->SetMarkerStyle(20);
   Graph_Graph339->SetMarkerSize(0.9);
   Graph_Graph339->GetXaxis()->SetLabelFont(42);
   Graph_Graph339->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetXaxis()->SetTitleFont(42);
   Graph_Graph339->GetYaxis()->SetLabelFont(42);
   Graph_Graph339->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetYaxis()->SetTickLength(0.02);
   Graph_Graph339->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetYaxis()->SetTitleFont(42);
   Graph_Graph339->GetZaxis()->SetLabelFont(42);
   Graph_Graph339->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph339);
   
   graph->Draw("p");
   
   Double_t Graph4_fx40[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy40[35] = { 100, 100, 3.111076, 100, 100, 100, 100, 100, 3.111076, 100, 100, 2.438831, 100, 6.518424, 3.111076, 2.150273, 100,
   3.437066, 2.504969, 2.83612, 2.530158, 2.582157, 1.930928, 2.421832, 2.480173, 2.700233, 2.850604, 2.772307, 3.214777, 3.285384, 3.290665, 3.361976, 3.371513,
   3.469491, 3.694642 };
   graph = new TGraph(35,Graph4_fx40,Graph4_fy40);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph440 = new TH1F("Graph_Graph440","Graph",100,0,3520);
   Graph_Graph440->SetMinimum(1.737835);
   Graph_Graph440->SetMaximum(109.8069);
   Graph_Graph440->SetDirectory(nullptr);
   Graph_Graph440->SetStats(0);
   Graph_Graph440->SetLineWidth(2);
   Graph_Graph440->SetMarkerStyle(20);
   Graph_Graph440->SetMarkerSize(0.9);
   Graph_Graph440->GetXaxis()->SetLabelFont(42);
   Graph_Graph440->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetXaxis()->SetTitleFont(42);
   Graph_Graph440->GetYaxis()->SetLabelFont(42);
   Graph_Graph440->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetYaxis()->SetTickLength(0.02);
   Graph_Graph440->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetYaxis()->SetTitleFont(42);
   Graph_Graph440->GetZaxis()->SetLabelFont(42);
   Graph_Graph440->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph440);
   
   graph->Draw("p");
   
   Double_t Graph5_fx41[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy41[35] = { 100, 100, 0, 100, 100, 100, 100, 100, 0, 100, 100, 0, 100, 0, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1920164, 0.7777691, 0.3464937, 0.1710773, 0.3031433,
   0, 0.9248257 };
   graph = new TGraph(35,Graph5_fx41,Graph5_fy41);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph541 = new TH1F("Graph_Graph541","Graph",100,0,3520);
   Graph_Graph541->SetMinimum(99.9);
   Graph_Graph541->SetMaximum(101.1);
   Graph_Graph541->SetDirectory(nullptr);
   Graph_Graph541->SetStats(0);
   Graph_Graph541->SetLineWidth(2);
   Graph_Graph541->SetMarkerStyle(20);
   Graph_Graph541->SetMarkerSize(0.9);
   Graph_Graph541->GetXaxis()->SetLabelFont(42);
   Graph_Graph541->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetXaxis()->SetTitleFont(42);
   Graph_Graph541->GetYaxis()->SetLabelFont(42);
   Graph_Graph541->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetYaxis()->SetTickLength(0.02);
   Graph_Graph541->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetYaxis()->SetTitleFont(42);
   Graph_Graph541->GetZaxis()->SetLabelFont(42);
   Graph_Graph541->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph541);
   
   graph->Draw("p");
   
   Double_t Graph6_fx42[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy42[35] = { 223.6068, 223.6068, 113.7489, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 29.06551, 223.6068, 223.6068, 26.29519, 223.6068, 102.0367, 104.3657, 24.16726, 223.6068,
   42.33891, 38.97141, 351.429, 158.4253, 74.46509, 30.61218, 47.2912, 24.45563, 35.38647, 33.57623, 16.35914, 7.749499, 7.636689, 6.995664, 7.686996, 8.838991,
   9.183475, 17.45341 };
   graph = new TGraph(35,Graph6_fx42,Graph6_fy42);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph642 = new TH1F("Graph_Graph642","Graph",100,0,3520);
   Graph_Graph642->SetMinimum(6.296098);
   Graph_Graph642->SetMaximum(385.8724);
   Graph_Graph642->SetDirectory(nullptr);
   Graph_Graph642->SetStats(0);
   Graph_Graph642->SetLineWidth(2);
   Graph_Graph642->SetMarkerStyle(20);
   Graph_Graph642->SetMarkerSize(0.9);
   Graph_Graph642->GetXaxis()->SetLabelFont(42);
   Graph_Graph642->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetXaxis()->SetTitleFont(42);
   Graph_Graph642->GetYaxis()->SetLabelFont(42);
   Graph_Graph642->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetYaxis()->SetTickLength(0.02);
   Graph_Graph642->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetYaxis()->SetTitleFont(42);
   Graph_Graph642->GetZaxis()->SetLabelFont(42);
   Graph_Graph642->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph642);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1800 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
