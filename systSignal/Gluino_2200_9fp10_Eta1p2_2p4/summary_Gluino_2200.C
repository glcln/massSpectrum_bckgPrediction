#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2200()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Jul  6 16:06:00 2026) by ROOT version 6.32.13
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
   Double_t Graph0_fy50[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 32.05665, 100, 100, 54.36042, 75.67915, 0,
   100, 0, 100, 0, 0, 8.589983, 25.34635, 15.24467, 5.49463, 5.445206, 9.754419, 2.320433, 7.234573, 4.58231, 0.8940578, 3.971261,
   2.54035, 6.029892 };
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
   Double_t Graph1_fy51[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 100, 32.05665, 100, 100, 9.280098, 75.67915, 0,
   100, 0, 100, 0, 0, 8.589994, 20.2109, 2.145445, 5.746675, 5.445206, 2.18001, 0.7430077, 1.149261, 1.724756, 0.3016055, 1.620054,
   1.77747, 2.288455 };
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
   Graph_Graph151->SetMinimum(99.9);
   Graph_Graph151->SetMaximum(101.1);
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
   Double_t Graph2_fy52[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 38.32729, 261.4552, 63.87022, 100, 100, 37.81927, 61.05094, 514.413,
   100, 202.3992, 100, 81.80632, 23.42031, 12.14176, 46.16048, 18.15161, 24.68246, 40.92112, 21.08856, 26.5551, 23.48018, 22.04509, 25.27494, 21.12971,
   26.51122, 25.3442 };
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
   Graph_Graph252->SetMaximum(564.6401);
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
   Double_t Graph3_fy53[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 0, 100, 100, 0, 0, 0,
   100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
   Double_t Graph4_fy54[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 4.550165, 3.647768, 1.936507, 100, 100, 3.446949, 5.820227, 2.042818,
   100, 3.647768, 100, 3.625548, 1.490265, 2.959657, 3.444386, 2.802479, 2.620161, 2.608836, 2.931321, 3.015363, 2.906466, 3.235722, 3.353083, 3.408325,
   3.612959, 3.275061 };
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
   Graph_Graph454->SetMinimum(1.341239);
   Graph_Graph454->SetMaximum(109.851);
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
   Double_t Graph5_fy55[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 0, 100, 100, 0, 0, 0,
   100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 1.738954, 1.47078, 1.154447, 0.3395855, 0.6693482, 0.2975762,
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
   Double_t Graph6_fy56[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 38.59644, 297.2745, 78.34802, 223.6068, 223.6068, 66.95786, 123.3522, 514.417,
   223.6068, 202.432, 223.6068, 81.88662, 23.46768, 17.42864, 56.5117, 23.96534, 26.06347, 41.72103, 23.52066, 26.83658, 24.76744, 22.81289, 25.51385, 21.82835,
   26.93531, 26.35624 };
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
   Graph_Graph656->SetMaximum(564.1158);
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
