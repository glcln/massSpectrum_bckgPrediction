#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy29[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 27.51509, 19.36631,
   13.48532, 0, 13.30897, 35.061, 27.17528, 0, 12.2113, 1.522791, 3.771043, 4.474926, 2.893722, 3.29597, 0.1722336, 1.995069, 2.657104, 1.150501,
   3.022337, 0.9047449 };
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
   Graph_Graph029->GetXaxis()->SetLabelSize(22);
   Graph_Graph029->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph029->GetXaxis()->SetTitleOffset(1);
   Graph_Graph029->GetXaxis()->SetTitleFont(42);
   Graph_Graph029->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph029->GetYaxis()->SetLabelFont(43);
   Graph_Graph029->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph029->GetYaxis()->SetLabelSize(22);
   Graph_Graph029->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph029->GetYaxis()->SetTickLength(0.02);
   Graph_Graph029->GetYaxis()->SetTitleOffset(1);
   Graph_Graph029->GetYaxis()->SetTitleFont(42);
   Graph_Graph029->GetZaxis()->SetLabelFont(42);
   Graph_Graph029->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph029->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph029->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph029->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph029->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph029);
   
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
   
   Double_t Graph1_fx30[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy30[35] = { 100, 100, 100, 100, 0, 75.46519, 320.2316, 14.4192, 0, 100, 0, 27.14648, 37.71579, 100, 100, 30.40602, 47.21861,
   32.87966, 44.19194, 0.1112461, 26.58557, 28.29906, 10.97984, 19.52612, 16.0478, 4.581064, 9.1465, 3.475809, 2.50296, 0.7405877, 1.701391, 3.664768, 0.9668469,
   4.461277, 0.9047449 };
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
   Graph_Graph130->SetMinimum(77.97684);
   Graph_Graph130->SetMaximum(342.2547);
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
   Double_t Graph2_fy31[35] = { 100, 100, 100, 100, 30.75011, 38.50652, 43.9416, 17.94869, 24.00274, 17.60072, 47.97527, 9.90755, 26.3146, 100, 100, 28.18118, 43.31357,
   33.28905, 24.14929, 13.84249, 35.60442, 21.81557, 14.26194, 15.25032, 9.132248, 4.823589, 11.79178, 5.263048, 3.673446, 1.07919, 0.9123325, 4.358613, 0.1232147,
   11.60578, 10.62864 };
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
   Graph_Graph231->SetMinimum(0.1108932);
   Graph_Graph231->SetMaximum(109.9877);
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
   Double_t Graph3_fy32[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 100, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 4.24273, 0.8445501, 0.9999812, 1.334852, 1.436973, 1.786113, 2.704567,
   3.258562, 4.036963 };
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
   Double_t Graph4_fy33[35] = { 100, 100, 100, 100, 4.904199, 2.805918, 1.234829, 2.548844, 1.287436, 0.9231091, 3.118026, 2.44565, 2.405941, 100, 100, 2.373648, 2.800655,
   2.895236, 3.41059, 2.507925, 1.267648, 3.317964, 2.70611, 2.980065, 2.387971, 3.358912, 3.242576, 3.240061, 3.329265, 3.417456, 3.443074, 3.409874, 3.480291,
   3.389895, 3.472304 };
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
   Graph_Graph433->SetMinimum(0.8307981);
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
   Double_t Graph5_fy34[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 22.83404, 0, 0, 0, 21.63039, 100, 100, 0, 0,
   0, 0, 0, 0, 0, 20.30948, 0, 0, 10.96953, 0, 0.7340252, 0.2218306, 0.252974, 0.243336, 0.6789684, 0.7375062,
   0, 0.5482376 };
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
   Double_t Graph6_fy35[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 31.13873, 84.76804, 323.2346, 23.16388, 24.03724, 101.5413, 48.07648, 29.00125, 46.05136, 223.6068, 223.6068, 49.81385, 66.99672,
   48.7798, 50.47522, 19.36609, 56.61579, 45.01401, 18.20118, 27.78197, 18.68023, 8.352019, 16.46951, 7.704921, 6.535014, 3.899182, 4.650487, 7.369238, 4.658405,
   13.63236, 11.95655 };
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
   Graph_Graph635->SetMinimum(3.509264);
   Graph_Graph635->SetMaximum(355.1682);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1600 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
