#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2400()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:55 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx49[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy49[35] = { 100, 100, 0, 100, 100, 100, 100, 100, 58.92458, 0, 0, 100, 100, 43.69671, 0, 0, 67.54655,
   179.9476, 33.5132, 0, 66.28979, 20.81307, 12.61127, 11.74359, 5.985618, 2.988327, 2.574229, 4.792303, 7.381856, 7.587189, 8.208013, 1.88244, 2.217865,
   5.274427, 4.846895 };
   TGraph *graph = new TGraph(35,Graph0_fx49,Graph0_fy49);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph049 = new TH1F("Graph_Graph049","Graph",100,0,3520);
   Graph_Graph049->SetMinimum(1);
   Graph_Graph049->SetMaximum(2000);
   Graph_Graph049->SetDirectory(nullptr);
   Graph_Graph049->SetStats(0);
   Graph_Graph049->SetLineWidth(2);
   Graph_Graph049->SetMarkerStyle(20);
   Graph_Graph049->SetMarkerSize(0.9);
   Graph_Graph049->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph049->GetXaxis()->SetRange(1,101);
   Graph_Graph049->GetXaxis()->SetLabelFont(43);
   Graph_Graph049->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetXaxis()->SetLabelSize(16);
   Graph_Graph049->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph049->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetXaxis()->SetTitleFont(42);
   Graph_Graph049->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph049->GetYaxis()->SetLabelFont(43);
   Graph_Graph049->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetYaxis()->SetLabelSize(16);
   Graph_Graph049->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph049->GetYaxis()->SetTickLength(0.02);
   Graph_Graph049->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetYaxis()->SetTitleFont(42);
   Graph_Graph049->GetZaxis()->SetLabelFont(42);
   Graph_Graph049->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph049->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph049->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph049->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph049->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph049);
   
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
   
   Double_t Graph1_fx50[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy50[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 58.92458, 100, 0, 100, 100, 1.192093e-05, 0, 100, 1.186872,
   87.07249, 13.4258, 17.4288, 35.20592, 9.901964, 12.55021, 0, 0, 0, 0, 2.375281, 4.688621, 1.154101, 1.527572, 0.2285004, 0.3132701,
   0.9796321, 2.261102 };
   graph = new TGraph(35,Graph1_fx50,Graph1_fy50);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph150 = new TH1F("Graph_Graph150","Graph",100,0,3520);
   Graph_Graph150->SetMinimum(54.81704);
   Graph_Graph150->SetMaximum(104.1075);
   Graph_Graph150->SetDirectory(nullptr);
   Graph_Graph150->SetStats(0);
   Graph_Graph150->SetLineWidth(2);
   Graph_Graph150->SetMarkerStyle(20);
   Graph_Graph150->SetMarkerSize(0.9);
   Graph_Graph150->GetXaxis()->SetLabelFont(42);
   Graph_Graph150->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetXaxis()->SetTitleFont(42);
   Graph_Graph150->GetYaxis()->SetLabelFont(42);
   Graph_Graph150->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetYaxis()->SetTickLength(0.02);
   Graph_Graph150->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetYaxis()->SetTitleFont(42);
   Graph_Graph150->GetZaxis()->SetLabelFont(42);
   Graph_Graph150->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph150->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph150->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph150->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph150->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph150);
   
   graph->Draw("p");
   
   Double_t Graph2_fx51[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy51[35] = { 100, 100, 23.1173, 100, 100, 100, 100, 16.41645, 22.27921, 23.59124, 19.96554, 100, 100, 23.86921, 18.07167, 14.21868, 17.83574,
   34.19884, 30.05189, 24.06104, 25.50798, 22.37631, 30.9772, 40.95699, 38.36626, 32.99677, 33.72689, 51.06199, 23.80619, 30.06325, 35.55865, 39.70108, 38.16103,
   39.37952, 32.18884 };
   graph = new TGraph(35,Graph2_fx51,Graph2_fy51);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph251 = new TH1F("Graph_Graph251","Graph",100,0,3520);
   Graph_Graph251->SetMinimum(5.640543);
   Graph_Graph251->SetMaximum(108.5781);
   Graph_Graph251->SetDirectory(nullptr);
   Graph_Graph251->SetStats(0);
   Graph_Graph251->SetLineWidth(2);
   Graph_Graph251->SetMarkerStyle(20);
   Graph_Graph251->SetMarkerSize(0.9);
   Graph_Graph251->GetXaxis()->SetLabelFont(42);
   Graph_Graph251->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetXaxis()->SetTitleFont(42);
   Graph_Graph251->GetYaxis()->SetLabelFont(42);
   Graph_Graph251->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetYaxis()->SetTickLength(0.02);
   Graph_Graph251->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetYaxis()->SetTitleFont(42);
   Graph_Graph251->GetZaxis()->SetLabelFont(42);
   Graph_Graph251->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph251->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph251->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph251->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph251->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph251);
   
   graph->Draw("p");
   
   Double_t Graph3_fx52[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy52[35] = { 100, 100, 0, 100, 100, 100, 100, 0, 0, 0, 0, 100, 100, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx52,Graph3_fy52);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph352 = new TH1F("Graph_Graph352","Graph",100,0,3520);
   Graph_Graph352->SetMinimum(99.9);
   Graph_Graph352->SetMaximum(101.1);
   Graph_Graph352->SetDirectory(nullptr);
   Graph_Graph352->SetStats(0);
   Graph_Graph352->SetLineWidth(2);
   Graph_Graph352->SetMarkerStyle(20);
   Graph_Graph352->SetMarkerSize(0.9);
   Graph_Graph352->GetXaxis()->SetLabelFont(42);
   Graph_Graph352->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetXaxis()->SetTitleFont(42);
   Graph_Graph352->GetYaxis()->SetLabelFont(42);
   Graph_Graph352->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetYaxis()->SetTickLength(0.02);
   Graph_Graph352->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetYaxis()->SetTitleFont(42);
   Graph_Graph352->GetZaxis()->SetLabelFont(42);
   Graph_Graph352->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph352->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph352->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph352->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph352->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph352);
   
   graph->Draw("p");
   
   Double_t Graph4_fx53[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy53[35] = { 100, 100, 1.802301, 100, 100, 100, 100, 2.28132, 1.866275, 0.9150028, 1.231527, 100, 100, 1.785344, 1.852858, 0.7658124, 1.49861,
   0.9149969, 1.400441, 2.023935, 2.052748, 1.91648, 2.56902, 2.137649, 2.109307, 1.877654, 1.918101, 1.673514, 1.961201, 1.837355, 2.006692, 1.978451, 2.003092,
   1.957411, 2.029461 };
   graph = new TGraph(35,Graph4_fx53,Graph4_fy53);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph453 = new TH1F("Graph_Graph453","Graph",100,0,3520);
   Graph_Graph453->SetMinimum(0.6892312);
   Graph_Graph453->SetMaximum(109.9234);
   Graph_Graph453->SetDirectory(nullptr);
   Graph_Graph453->SetStats(0);
   Graph_Graph453->SetLineWidth(2);
   Graph_Graph453->SetMarkerStyle(20);
   Graph_Graph453->SetMarkerSize(0.9);
   Graph_Graph453->GetXaxis()->SetLabelFont(42);
   Graph_Graph453->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetXaxis()->SetTitleFont(42);
   Graph_Graph453->GetYaxis()->SetLabelFont(42);
   Graph_Graph453->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetYaxis()->SetTickLength(0.02);
   Graph_Graph453->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetYaxis()->SetTitleFont(42);
   Graph_Graph453->GetZaxis()->SetLabelFont(42);
   Graph_Graph453->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph453->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph453->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph453->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph453->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph453);
   
   graph->Draw("p");
   
   Double_t Graph5_fx54[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy54[35] = { 223.6068, 223.6068, 102.6531, 223.6068, 223.6068, 223.6068, 223.6068, 142.3893, 86.27896, 102.7491, 20.00348, 223.6068, 223.6068, 49.82298, 18.16641, 101.0087, 69.8878,
   202.8131, 46.99429, 29.77907, 79.3011, 32.18082, 35.81534, 42.66095, 38.88761, 33.18498, 33.87933, 51.36862, 25.4373, 31.0817, 36.58073, 39.79555, 38.27915,
   39.79143, 32.6932 };
   graph = new TGraph(35,Graph5_fx54,Graph5_fy54);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph554 = new TH1F("Graph_Graph554","Graph",100,0,3520);
   Graph_Graph554->SetMinimum(16.34977);
   Graph_Graph554->SetMaximum(244.1508);
   Graph_Graph554->SetDirectory(nullptr);
   Graph_Graph554->SetStats(0);
   Graph_Graph554->SetLineWidth(2);
   Graph_Graph554->SetMarkerStyle(20);
   Graph_Graph554->SetMarkerSize(0.9);
   Graph_Graph554->GetXaxis()->SetLabelFont(42);
   Graph_Graph554->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetXaxis()->SetTitleFont(42);
   Graph_Graph554->GetYaxis()->SetLabelFont(42);
   Graph_Graph554->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetYaxis()->SetTickLength(0.02);
   Graph_Graph554->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetYaxis()->SetTitleFont(42);
   Graph_Graph554->GetZaxis()->SetLabelFont(42);
   Graph_Graph554->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph554->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph554->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph554->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph554->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph554);
   
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
