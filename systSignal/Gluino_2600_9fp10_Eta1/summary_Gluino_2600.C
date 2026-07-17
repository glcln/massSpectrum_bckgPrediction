#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2600()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:10:25 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx64[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy64[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 21.3025, 36.83226, 16.78014, 55.4904, 100, 100, 46.08107, 100,
   0, 0, 0, 22.549, 39.07511, 0.1355171, 31.3894, 9.035235, 7.899302, 8.496511, 1.656926, 2.606189, 3.106213, 4.713929, 0.6076396, 0.6009817,
   1.949036, 1.877558 };
   TGraph *graph = new TGraph(35,Graph0_fx64,Graph0_fy64);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph064 = new TH1F("Graph_Graph064","Graph",100,0,3520);
   Graph_Graph064->SetMinimum(1);
   Graph_Graph064->SetMaximum(2000);
   Graph_Graph064->SetDirectory(nullptr);
   Graph_Graph064->SetStats(0);
   Graph_Graph064->SetLineWidth(2);
   Graph_Graph064->SetMarkerStyle(20);
   Graph_Graph064->SetMarkerSize(0.9);
   Graph_Graph064->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph064->GetXaxis()->SetRange(1,101);
   Graph_Graph064->GetXaxis()->SetLabelFont(43);
   Graph_Graph064->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetXaxis()->SetLabelSize(22);
   Graph_Graph064->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph064->GetXaxis()->SetTitleOffset(1);
   Graph_Graph064->GetXaxis()->SetTitleFont(42);
   Graph_Graph064->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph064->GetYaxis()->SetLabelFont(43);
   Graph_Graph064->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetYaxis()->SetLabelSize(22);
   Graph_Graph064->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph064->GetYaxis()->SetTickLength(0.02);
   Graph_Graph064->GetYaxis()->SetTitleOffset(1);
   Graph_Graph064->GetYaxis()->SetTitleFont(42);
   Graph_Graph064->GetZaxis()->SetLabelFont(42);
   Graph_Graph064->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph064->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph064->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph064->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph064->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph064);
   
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
   
   Double_t Graph1_fx65[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy65[35] = { 100, 100, 0, 100, 0, 100, 120.3678, 100, 100, 57.61015, 31.86309, 47.38583, 55.4904, 100, 100, 46.08107, 100,
   27.70239, 0, 26.56601, 22.549, 39.07511, 3.243208, 31.3894, 7.622159, 14.92022, 4.92456, 3.536785, 6.760204, 1.277542, 4.211307, 0.5222321, 0.289619,
   1.920378, 2.348351 };
   graph = new TGraph(35,Graph1_fx65,Graph1_fy65);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph165 = new TH1F("Graph_Graph165","Graph",100,0,3520);
   Graph_Graph165->SetMinimum(97.96322);
   Graph_Graph165->SetMaximum(122.4045);
   Graph_Graph165->SetDirectory(nullptr);
   Graph_Graph165->SetStats(0);
   Graph_Graph165->SetLineWidth(2);
   Graph_Graph165->SetMarkerStyle(20);
   Graph_Graph165->SetMarkerSize(0.9);
   Graph_Graph165->GetXaxis()->SetLabelFont(42);
   Graph_Graph165->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetXaxis()->SetTitleFont(42);
   Graph_Graph165->GetYaxis()->SetLabelFont(42);
   Graph_Graph165->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetYaxis()->SetTickLength(0.02);
   Graph_Graph165->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetYaxis()->SetTitleFont(42);
   Graph_Graph165->GetZaxis()->SetLabelFont(42);
   Graph_Graph165->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph165->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph165->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph165->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph165->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph165);
   
   graph->Draw("p");
   
   Double_t Graph2_fx66[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy66[35] = { 100, 100, 47.33641, 100, 31.28117, 27.03399, 34.71781, 100, 100, 66.15778, 42.67295, 19.41261, 29.70948, 100, 14.12714, 40.3645, 100,
   21.64087, 21.96645, 18.05142, 22.62325, 26.29202, 23.32653, 12.04359, 26.97548, 9.410334, 13.91073, 3.049886, 9.046721, 6.185931, 1.008797, 2.405882, 2.270752,
   1.743037, 1.596004 };
   graph = new TGraph(35,Graph2_fx66,Graph2_fy66);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph266 = new TH1F("Graph_Graph266","Graph",100,0,3520);
   Graph_Graph266->SetMinimum(0.907917);
   Graph_Graph266->SetMaximum(109.8991);
   Graph_Graph266->SetDirectory(nullptr);
   Graph_Graph266->SetStats(0);
   Graph_Graph266->SetLineWidth(2);
   Graph_Graph266->SetMarkerStyle(20);
   Graph_Graph266->SetMarkerSize(0.9);
   Graph_Graph266->GetXaxis()->SetLabelFont(42);
   Graph_Graph266->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetXaxis()->SetTitleFont(42);
   Graph_Graph266->GetYaxis()->SetLabelFont(42);
   Graph_Graph266->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetYaxis()->SetTickLength(0.02);
   Graph_Graph266->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetYaxis()->SetTitleFont(42);
   Graph_Graph266->GetZaxis()->SetLabelFont(42);
   Graph_Graph266->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph266->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph266->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph266->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph266->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph266);
   
   graph->Draw("p");
   
   Double_t Graph3_fx67[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy67[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 3.441584, 1.759863, 3.622925, 1.414752, 1.249957, 0.6706238, 1.036853,
   1.188374, 1.355445 };
   graph = new TGraph(35,Graph3_fx67,Graph3_fy67);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph367 = new TH1F("Graph_Graph367","Graph",100,0,3520);
   Graph_Graph367->SetMinimum(99.9);
   Graph_Graph367->SetMaximum(101.1);
   Graph_Graph367->SetDirectory(nullptr);
   Graph_Graph367->SetStats(0);
   Graph_Graph367->SetLineWidth(2);
   Graph_Graph367->SetMarkerStyle(20);
   Graph_Graph367->SetMarkerSize(0.9);
   Graph_Graph367->GetXaxis()->SetLabelFont(42);
   Graph_Graph367->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetXaxis()->SetTitleFont(42);
   Graph_Graph367->GetYaxis()->SetLabelFont(42);
   Graph_Graph367->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetYaxis()->SetTickLength(0.02);
   Graph_Graph367->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetYaxis()->SetTitleFont(42);
   Graph_Graph367->GetZaxis()->SetLabelFont(42);
   Graph_Graph367->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph367->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph367->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph367->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph367->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph367);
   
   graph->Draw("p");
   
   Double_t Graph4_fx68[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy68[35] = { 100, 100, 0.9231091, 100, 3.647768, 3.647757, 4.550165, 100, 100, 3.647768, 3.134429, 3.647768, 1.061857, 100, 3.647768, 5.199045, 100,
   3.567016, 3.647768, 1.974607, 2.967942, 2.767122, 3.260028, 2.593124, 2.735531, 3.390467, 3.049934, 3.12624, 3.029108, 3.225148, 3.369284, 3.445506, 3.443825,
   3.517365, 3.523147 };
   graph = new TGraph(35,Graph4_fx68,Graph4_fy68);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph468 = new TH1F("Graph_Graph468","Graph",100,0,3520);
   Graph_Graph468->SetMinimum(0.8307981);
   Graph_Graph468->SetMaximum(109.9077);
   Graph_Graph468->SetDirectory(nullptr);
   Graph_Graph468->SetStats(0);
   Graph_Graph468->SetLineWidth(2);
   Graph_Graph468->SetMarkerStyle(20);
   Graph_Graph468->SetMarkerSize(0.9);
   Graph_Graph468->GetXaxis()->SetLabelFont(42);
   Graph_Graph468->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetXaxis()->SetTitleFont(42);
   Graph_Graph468->GetYaxis()->SetLabelFont(42);
   Graph_Graph468->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetYaxis()->SetTickLength(0.02);
   Graph_Graph468->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetYaxis()->SetTitleFont(42);
   Graph_Graph468->GetZaxis()->SetLabelFont(42);
   Graph_Graph468->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph468->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph468->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph468->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph468->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph468);
   
   graph->Draw("p");
   
   Double_t Graph5_fx69[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy69[35] = { 100, 100, 0, 100, 0, 0, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 6.303841, 0, 0, 0.817728, 1.096243, 0.07704496, 0.445348,
   0.2879024, 0.1924813 };
   graph = new TGraph(35,Graph5_fx69,Graph5_fy69);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph569 = new TH1F("Graph_Graph569","Graph",100,0,3520);
   Graph_Graph569->SetMinimum(99.9);
   Graph_Graph569->SetMaximum(101.1);
   Graph_Graph569->SetDirectory(nullptr);
   Graph_Graph569->SetStats(0);
   Graph_Graph569->SetLineWidth(2);
   Graph_Graph569->SetMarkerStyle(20);
   Graph_Graph569->SetMarkerSize(0.9);
   Graph_Graph569->GetXaxis()->SetLabelFont(42);
   Graph_Graph569->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetXaxis()->SetTitleFont(42);
   Graph_Graph569->GetYaxis()->SetLabelFont(42);
   Graph_Graph569->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetYaxis()->SetTickLength(0.02);
   Graph_Graph569->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetYaxis()->SetTitleFont(42);
   Graph_Graph569->GetZaxis()->SetLabelFont(42);
   Graph_Graph569->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph569->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph569->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph569->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph569->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph569);
   
   graph->Draw("p");
   
   Double_t Graph6_fx70[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy70[35] = { 223.6068, 223.6068, 47.34541, 223.6068, 31.49314, 103.654, 125.3572, 223.6068, 223.6068, 90.34869, 64.82806, 54.0106, 83.91752, 223.6068, 142.172, 76.83263, 223.6068,
   35.33374, 22.26726, 32.17927, 39.21141, 61.25892, 23.77587, 46.0691, 29.57858, 19.62299, 17.63794, 6.117724, 12.51546, 7.870801, 7.340869, 4.330297, 4.305388,
   4.930402, 5.083012 };
   graph = new TGraph(35,Graph6_fx70,Graph6_fy70);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph670 = new TH1F("Graph_Graph670","Graph",100,0,3520);
   Graph_Graph670->SetMinimum(3.87485);
   Graph_Graph670->SetMaximum(245.5369);
   Graph_Graph670->SetDirectory(nullptr);
   Graph_Graph670->SetStats(0);
   Graph_Graph670->SetLineWidth(2);
   Graph_Graph670->SetMarkerStyle(20);
   Graph_Graph670->SetMarkerSize(0.9);
   Graph_Graph670->GetXaxis()->SetLabelFont(42);
   Graph_Graph670->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetXaxis()->SetTitleFont(42);
   Graph_Graph670->GetYaxis()->SetLabelFont(42);
   Graph_Graph670->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetYaxis()->SetTickLength(0.02);
   Graph_Graph670->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetYaxis()->SetTitleFont(42);
   Graph_Graph670->GetZaxis()->SetLabelFont(42);
   Graph_Graph670->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph670->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph670->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph670->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph670->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph670);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2600 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
