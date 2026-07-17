#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2200()
{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jul  8 13:55:25 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx50[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy50[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 51.39296, 100, 100, 98.98171, 75.67915, 100,
   100, 0, 0, 0, 1.192093e-05, 8.589983, 25.34635, 16.81004, 11.9267, 0.4103899, 9.443987, 4.312706, 7.297111, 4.939246, 0.7892728, 4.320347,
   3.239238, 6.245065 };
   TGraph *graph = new TGraph(35,Graph0_fx50,Graph0_fy50);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph050 = new TH1F("Graph_Graph050","Graph",100,0,3520);
   Graph_Graph050->SetMinimum(1);
   Graph_Graph050->SetMaximum(2000);
   Graph_Graph050->SetDirectory(nullptr);
   Graph_Graph050->SetStats(0);
   Graph_Graph050->SetLineWidth(2);
   Graph_Graph050->SetMarkerStyle(20);
   Graph_Graph050->SetMarkerSize(0.9);
   Graph_Graph050->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph050->GetXaxis()->SetRange(1,101);
   Graph_Graph050->GetXaxis()->SetLabelFont(43);
   Graph_Graph050->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetXaxis()->SetLabelSize(16);
   Graph_Graph050->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph050->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph050->GetXaxis()->SetTitleFont(42);
   Graph_Graph050->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph050->GetYaxis()->SetLabelFont(43);
   Graph_Graph050->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetYaxis()->SetLabelSize(16);
   Graph_Graph050->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph050->GetYaxis()->SetTickLength(0.02);
   Graph_Graph050->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph050->GetYaxis()->SetTitleFont(42);
   Graph_Graph050->GetZaxis()->SetLabelFont(42);
   Graph_Graph050->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph050->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph050->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph050->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph050);
   
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
   
   Double_t Graph1_fx51[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy51[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 51.39296, 100, 100, 98.98171, 75.67915, 100,
   100, 0, 0, 0, 0, 8.589994, 20.2109, 2.365756, 2.2421, 5.395091, 2.110636, 1.599443, 1.186919, 1.598048, 0.116545, 1.474625,
   1.660573, 2.147758 };
   graph = new TGraph(35,Graph1_fx51,Graph1_fy51);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph151 = new TH1F("Graph_Graph151","Graph",100,0,3520);
   Graph_Graph151->SetMinimum(46.53226);
   Graph_Graph151->SetMaximum(104.8607);
   Graph_Graph151->SetDirectory(nullptr);
   Graph_Graph151->SetStats(0);
   Graph_Graph151->SetLineWidth(2);
   Graph_Graph151->SetMarkerStyle(20);
   Graph_Graph151->SetMarkerSize(0.9);
   Graph_Graph151->GetXaxis()->SetLabelFont(42);
   Graph_Graph151->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetXaxis()->SetTitleFont(42);
   Graph_Graph151->GetYaxis()->SetLabelFont(42);
   Graph_Graph151->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetYaxis()->SetTickLength(0.02);
   Graph_Graph151->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetYaxis()->SetTitleFont(42);
   Graph_Graph151->GetZaxis()->SetLabelFont(42);
   Graph_Graph151->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph151);
   
   graph->Draw("p");
   
   Double_t Graph2_fx52[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy52[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 261.4552, 62.81124, 100, 100, 43.30014, 61.05094, 100,
   100, 202.3992, 90.07781, 81.80632, 26.5945, 12.14176, 46.16048, 15.48322, 22.65323, 40.92112, 19.81486, 21.7417, 23.11814, 23.34439, 25.83896, 21.1136,
   27.01958, 23.21839 };
   graph = new TGraph(35,Graph2_fx52,Graph2_fy52);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph252 = new TH1F("Graph_Graph252","Graph",100,0,3520);
   Graph_Graph252->SetMinimum(10.92758);
   Graph_Graph252->SetMaximum(286.3866);
   Graph_Graph252->SetDirectory(nullptr);
   Graph_Graph252->SetStats(0);
   Graph_Graph252->SetLineWidth(2);
   Graph_Graph252->SetMarkerStyle(20);
   Graph_Graph252->SetMarkerSize(0.9);
   Graph_Graph252->GetXaxis()->SetLabelFont(42);
   Graph_Graph252->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetXaxis()->SetTitleFont(42);
   Graph_Graph252->GetYaxis()->SetLabelFont(42);
   Graph_Graph252->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetYaxis()->SetTickLength(0.02);
   Graph_Graph252->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetYaxis()->SetTitleFont(42);
   Graph_Graph252->GetZaxis()->SetLabelFont(42);
   Graph_Graph252->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph252);
   
   graph->Draw("p");
   
   Double_t Graph3_fx53[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy53[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 0, 0, 100,
   100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx53,Graph3_fy53);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph353 = new TH1F("Graph_Graph353","Graph",100,0,3520);
   Graph_Graph353->SetMinimum(99.9);
   Graph_Graph353->SetMaximum(101.1);
   Graph_Graph353->SetDirectory(nullptr);
   Graph_Graph353->SetStats(0);
   Graph_Graph353->SetLineWidth(2);
   Graph_Graph353->SetMarkerStyle(20);
   Graph_Graph353->SetMarkerSize(0.9);
   Graph_Graph353->GetXaxis()->SetLabelFont(42);
   Graph_Graph353->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetXaxis()->SetTitleFont(42);
   Graph_Graph353->GetYaxis()->SetLabelFont(42);
   Graph_Graph353->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetYaxis()->SetTickLength(0.02);
   Graph_Graph353->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetYaxis()->SetTitleFont(42);
   Graph_Graph353->GetZaxis()->SetLabelFont(42);
   Graph_Graph353->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph353);
   
   graph->Draw("p");
   
   Double_t Graph4_fx54[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy54[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 3.647768, 1.234829, 100, 100, 3.722632, 5.820227, 100,
   100, 3.647768, 3.647768, 3.625548, 1.529455, 2.959657, 3.444386, 2.623045, 2.557504, 2.608836, 2.872002, 3.09056, 2.872729, 3.260779, 3.331649, 3.422678,
   3.633428, 3.325665 };
   graph = new TGraph(35,Graph4_fx54,Graph4_fy54);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph454 = new TH1F("Graph_Graph454","Graph",100,0,3520);
   Graph_Graph454->SetMinimum(1.111346);
   Graph_Graph454->SetMaximum(109.8765);
   Graph_Graph454->SetDirectory(nullptr);
   Graph_Graph454->SetStats(0);
   Graph_Graph454->SetLineWidth(2);
   Graph_Graph454->SetMarkerStyle(20);
   Graph_Graph454->SetMarkerSize(0.9);
   Graph_Graph454->GetXaxis()->SetLabelFont(42);
   Graph_Graph454->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetXaxis()->SetTitleFont(42);
   Graph_Graph454->GetYaxis()->SetLabelFont(42);
   Graph_Graph454->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetYaxis()->SetTickLength(0.02);
   Graph_Graph454->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetYaxis()->SetTitleFont(42);
   Graph_Graph454->GetZaxis()->SetLabelFont(42);
   Graph_Graph454->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph454);
   
   graph->Draw("p");
   
   Double_t Graph5_fx55[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy55[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 0, 0, 100,
   100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.683611, 1.38765, 1.163161, 0.2996683, 0.7035732, 0.2682149,
   0, 0 };
   graph = new TGraph(35,Graph5_fx55,Graph5_fy55);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph555 = new TH1F("Graph_Graph555","Graph",100,0,3520);
   Graph_Graph555->SetMinimum(99.9);
   Graph_Graph555->SetMaximum(101.1);
   Graph_Graph555->SetDirectory(nullptr);
   Graph_Graph555->SetStats(0);
   Graph_Graph555->SetLineWidth(2);
   Graph_Graph555->SetMarkerStyle(20);
   Graph_Graph555->SetMarkerSize(0.9);
   Graph_Graph555->GetXaxis()->SetLabelFont(42);
   Graph_Graph555->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetXaxis()->SetTitleFont(42);
   Graph_Graph555->GetYaxis()->SetLabelFont(42);
   Graph_Graph555->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetYaxis()->SetTickLength(0.02);
   Graph_Graph555->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetYaxis()->SetTitleFont(42);
   Graph_Graph555->GetZaxis()->SetLabelFont(42);
   Graph_Graph555->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph555);
   
   graph->Draw("p");
   
   Double_t Graph6_fx56[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy56[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 297.2745, 96.06898, 223.6068, 223.6068, 146.5726, 123.3522, 223.6068,
   223.6068, 202.432, 90.15164, 81.88662, 26.63844, 17.42864, 56.5117, 23.12541, 25.82601, 41.35963, 22.23782, 22.43682, 24.4409, 24.13593, 26.06508, 21.87095,
   27.50472, 24.36735 };
   graph = new TGraph(35,Graph6_fx56,Graph6_fy56);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph656 = new TH1F("Graph_Graph656","Graph",100,0,3520);
   Graph_Graph656->SetMinimum(15.68578);
   Graph_Graph656->SetMaximum(325.2591);
   Graph_Graph656->SetDirectory(nullptr);
   Graph_Graph656->SetStats(0);
   Graph_Graph656->SetLineWidth(2);
   Graph_Graph656->SetMarkerStyle(20);
   Graph_Graph656->SetMarkerSize(0.9);
   Graph_Graph656->GetXaxis()->SetLabelFont(42);
   Graph_Graph656->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetXaxis()->SetTitleFont(42);
   Graph_Graph656->GetYaxis()->SetLabelFont(42);
   Graph_Graph656->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetYaxis()->SetTickLength(0.02);
   Graph_Graph656->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetYaxis()->SetTitleFont(42);
   Graph_Graph656->GetZaxis()->SetLabelFont(42);
   Graph_Graph656->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph656);
   
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
