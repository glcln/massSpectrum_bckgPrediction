#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2200()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:10:24 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx50[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy50[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 1.192093e-05, 7.312131, 4.270744, 1.192093e-05, 19.99837, 5.076921, 4.096746, 3.49499, 2.185023, 5.05296, 2.820086, 1.872998, 0.7868767, 2.054995,
   2.447116, 2.026761 };
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
   Graph_Graph050->GetXaxis()->SetLabelSize(22);
   Graph_Graph050->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph050->GetXaxis()->SetTitleOffset(1);
   Graph_Graph050->GetXaxis()->SetTitleFont(42);
   Graph_Graph050->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph050->GetYaxis()->SetLabelFont(43);
   Graph_Graph050->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetYaxis()->SetLabelSize(22);
   Graph_Graph050->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph050->GetYaxis()->SetTickLength(0.02);
   Graph_Graph050->GetYaxis()->SetTitleOffset(1);
   Graph_Graph050->GetYaxis()->SetTitleFont(42);
   Graph_Graph050->GetZaxis()->SetLabelFont(42);
   Graph_Graph050->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph050->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph050->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph050->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph050);
   
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
   
   Double_t Graph1_fx51[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy51[35] = { 100, 100, 100, 100, 100, 100, 0, 30.09958, 412.6682, 0, 0, 0, 0, 0, 0, 17.69736, 21.12389,
   12.36625, 0, 1.192093e-05, 0, 1.192093e-05, 8.548201, 12.08581, 0.7494092, 4.096746, 6.126368, 0.7652164, 4.796708, 1.871753, 1.630378, 0.5928755, 1.642281,
   2.277768, 1.293159 };
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
   Graph_Graph151->SetMinimum(68.73318);
   Graph_Graph151->SetMaximum(443.935);
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
   Double_t Graph2_fy52[35] = { 100, 100, 100, 100, 100, 100, 91.52737, 27.84832, 26.5033, 35.13874, 67.16209, 22.62713, 24.57522, 24.46201, 28.89853, 15.60764, 29.88533,
   27.62736, 31.75251, 22.75421, 15.64739, 5.525923, 25.056, 25.65176, 5.159593, 25.55894, 9.418535, 6.079042, 4.372108, 2.602994, 1.110315, 5.251723, 2.130151,
   0.4184902, 5.541152 };
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
   Graph_Graph252->SetMinimum(0.3766412);
   Graph_Graph252->SetMaximum(109.9582);
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
   Double_t Graph3_fy53[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 9.315443, 0, 0, 0, 3.844881, 0, 1.361936, 1.498717, 1.625961, 0.6852269, 0.7820547, 0.9374976,
   1.579475, 1.124609 };
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
   Double_t Graph4_fy54[35] = { 100, 100, 100, 100, 100, 100, 3.647757, 2.761465, 3.111076, 2.499926, 1.234829, 3.647768, 2.881527, 2.946353, 6.518424, 1.234829, 3.54861,
   2.998042, 1.393676, 3.213477, 4.364234, 3.892887, 2.624536, 3.461051, 2.548707, 2.794194, 3.011405, 3.029907, 3.12413, 3.368902, 3.424883, 3.401363, 3.476489,
   3.549683, 3.553891 };
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
   Double_t Graph5_fy55[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 4.609644, 0, 0, 0, 0, 0, 0, 1.997107, 1.601666, 0, 0.2204895, 0.3911614, 0.237149, 0.1681149,
   0.3000855, 0 };
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
   Double_t Graph6_fy56[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 91.60003, 41.09914, 413.5301, 35.22755, 67.17345, 22.91927, 24.74357, 24.63881, 29.62456, 23.62879, 36.7688,
   30.41682, 31.78308, 22.98, 17.81445, 12.27628, 26.60382, 34.8711, 7.710645, 26.63487, 12.14598, 7.304097, 8.925384, 5.676723, 4.42701, 6.382174, 4.941914,
   5.142639, 7.097836 };
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
   Graph_Graph656->SetMinimum(3.984309);
   Graph_Graph656->SetMaximum(454.4404);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2200 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
