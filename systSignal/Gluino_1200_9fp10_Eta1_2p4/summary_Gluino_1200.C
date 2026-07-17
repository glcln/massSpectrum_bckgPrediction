#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
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
   
   Double_t Graph0_fx8[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy8[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 30.67296, 52.40382, 0, 10.29285, 2.118099, 1.960748, 3.377378, 3.121591, 1.954687, 2.265567, 1.42827, 1.793432, 5.244553, 1.465178,
   15.11131, 10.2593 };
   TGraph *graph = new TGraph(35,Graph0_fx8,Graph0_fy8);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph08 = new TH1F("Graph_Graph08","Graph",100,0,3520);
   Graph_Graph08->SetMinimum(1);
   Graph_Graph08->SetMaximum(2000);
   Graph_Graph08->SetDirectory(nullptr);
   Graph_Graph08->SetStats(0);
   Graph_Graph08->SetLineWidth(2);
   Graph_Graph08->SetMarkerStyle(20);
   Graph_Graph08->SetMarkerSize(0.9);
   Graph_Graph08->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph08->GetXaxis()->SetRange(1,101);
   Graph_Graph08->GetXaxis()->SetLabelFont(43);
   Graph_Graph08->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetXaxis()->SetLabelSize(22);
   Graph_Graph08->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph08->GetXaxis()->SetTitleOffset(1);
   Graph_Graph08->GetXaxis()->SetTitleFont(42);
   Graph_Graph08->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph08->GetYaxis()->SetLabelFont(43);
   Graph_Graph08->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetYaxis()->SetLabelSize(22);
   Graph_Graph08->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph08->GetYaxis()->SetTickLength(0.02);
   Graph_Graph08->GetYaxis()->SetTitleOffset(1);
   Graph_Graph08->GetYaxis()->SetTitleFont(42);
   Graph_Graph08->GetZaxis()->SetLabelFont(42);
   Graph_Graph08->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph08->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph08->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph08->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph08->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph08);
   
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
   
   Double_t Graph1_fx9[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy9[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 65.59938, 100, 100, 19.77618, 86.52827, 26.2414, 54.31593, 5.960464e-06, 23.68468,
   50.35577, 75.47869, 59.10484, 16.90037, 22.87549, 27.38022, 31.62942, 34.29493, 5.292928, 9.466815, 1.659298, 2.762508, 4.894238, 5.908751, 8.553696, 1.878858,
   17.94295, 11.83543 };
   graph = new TGraph(35,Graph1_fx9,Graph1_fy9);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph19 = new TH1F("Graph_Graph19","Graph",100,0,3520);
   Graph_Graph19->SetMinimum(5.364418e-06);
   Graph_Graph19->SetMaximum(110);
   Graph_Graph19->SetDirectory(nullptr);
   Graph_Graph19->SetStats(0);
   Graph_Graph19->SetLineWidth(2);
   Graph_Graph19->SetMarkerStyle(20);
   Graph_Graph19->SetMarkerSize(0.9);
   Graph_Graph19->GetXaxis()->SetLabelFont(42);
   Graph_Graph19->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetXaxis()->SetTitleFont(42);
   Graph_Graph19->GetYaxis()->SetLabelFont(42);
   Graph_Graph19->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetYaxis()->SetTickLength(0.02);
   Graph_Graph19->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetYaxis()->SetTitleFont(42);
   Graph_Graph19->GetZaxis()->SetLabelFont(42);
   Graph_Graph19->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph19->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph19->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph19->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph19->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph19);
   
   graph->Draw("p");
   
   Double_t Graph2_fx10[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy10[35] = { 100, 100, 100, 100, 100, 100, 100, 27.03399, 35.32102, 100, 100, 20.02131, 23.19874, 24.19314, 17.18956, 13.84544, 14.02628,
   28.51107, 42.59614, 7.484424, 27.675, 25.01362, 6.087947, 11.92023, 7.235241, 15.11927, 4.807067, 5.541337, 2.302349, 0.7685065, 3.936511, 3.380811, 8.706474,
   8.956814, 12.79103 };
   graph = new TGraph(35,Graph2_fx10,Graph2_fy10);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph210 = new TH1F("Graph_Graph210","Graph",100,0,3520);
   Graph_Graph210->SetMinimum(0.6916559);
   Graph_Graph210->SetMaximum(109.9231);
   Graph_Graph210->SetDirectory(nullptr);
   Graph_Graph210->SetStats(0);
   Graph_Graph210->SetLineWidth(2);
   Graph_Graph210->SetMarkerStyle(20);
   Graph_Graph210->SetMarkerSize(0.9);
   Graph_Graph210->GetXaxis()->SetLabelFont(42);
   Graph_Graph210->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetXaxis()->SetTitleFont(42);
   Graph_Graph210->GetYaxis()->SetLabelFont(42);
   Graph_Graph210->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetYaxis()->SetTickLength(0.02);
   Graph_Graph210->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetYaxis()->SetTitleFont(42);
   Graph_Graph210->GetZaxis()->SetLabelFont(42);
   Graph_Graph210->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph210->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph210->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph210->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph210->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph210);
   
   graph->Draw("p");
   
   Double_t Graph3_fx11[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy11[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 3.548044, 0, 0, 29.39451, 8.027184, 7.600498, 6.04952, 2.364099, 3.87392, 3.722763, 2.539349, 2.817357, 2.348292, 0,
   6.519186, 2.358592 };
   graph = new TGraph(35,Graph3_fx11,Graph3_fy11);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph311 = new TH1F("Graph_Graph311","Graph",100,0,3520);
   Graph_Graph311->SetMinimum(99.9);
   Graph_Graph311->SetMaximum(101.1);
   Graph_Graph311->SetDirectory(nullptr);
   Graph_Graph311->SetStats(0);
   Graph_Graph311->SetLineWidth(2);
   Graph_Graph311->SetMarkerStyle(20);
   Graph_Graph311->SetMarkerSize(0.9);
   Graph_Graph311->GetXaxis()->SetLabelFont(42);
   Graph_Graph311->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetXaxis()->SetTitleFont(42);
   Graph_Graph311->GetYaxis()->SetLabelFont(42);
   Graph_Graph311->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetYaxis()->SetTickLength(0.02);
   Graph_Graph311->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetYaxis()->SetTitleFont(42);
   Graph_Graph311->GetZaxis()->SetLabelFont(42);
   Graph_Graph311->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph311->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph311->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph311->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph311->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph311);
   
   graph->Draw("p");
   
   Double_t Graph4_fx12[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy12[35] = { 100, 100, 100, 100, 100, 100, 100, 3.647757, 2.645302, 100, 100, 1.338339, 1.857972, 2.296519, 3.274274, 1.881361, 3.112805,
   2.632725, 1.184368, 2.755058, 2.003062, 2.500677, 3.323549, 3.076136, 3.003383, 2.83289, 3.054988, 3.209341, 3.385258, 3.300858, 3.297341, 3.044975, 3.656363,
   3.854537, 3.315735 };
   graph = new TGraph(35,Graph4_fx12,Graph4_fy12);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph412 = new TH1F("Graph_Graph412","Graph",100,0,3520);
   Graph_Graph412->SetMinimum(1.065931);
   Graph_Graph412->SetMaximum(109.8816);
   Graph_Graph412->SetDirectory(nullptr);
   Graph_Graph412->SetStats(0);
   Graph_Graph412->SetLineWidth(2);
   Graph_Graph412->SetMarkerStyle(20);
   Graph_Graph412->SetMarkerSize(0.9);
   Graph_Graph412->GetXaxis()->SetLabelFont(42);
   Graph_Graph412->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetXaxis()->SetTitleFont(42);
   Graph_Graph412->GetYaxis()->SetLabelFont(42);
   Graph_Graph412->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetYaxis()->SetTickLength(0.02);
   Graph_Graph412->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetYaxis()->SetTitleFont(42);
   Graph_Graph412->GetZaxis()->SetLabelFont(42);
   Graph_Graph412->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph412->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph412->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph412->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph412->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph412);
   
   graph->Draw("p");
   
   Double_t Graph5_fx13[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy13[35] = { 100, 100, 100, 100, 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 18.78004, 0,
   0, 0, 0, 44.98097, 0, 0, 5.458736, 0, 1.00944, 0.62657, 0.7269025, 0.5967021, 0.2509236, 0.5709946, 0, 0,
   0, 5.480587 };
   graph = new TGraph(35,Graph5_fx13,Graph5_fy13);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph513 = new TH1F("Graph_Graph513","Graph",100,0,3520);
   Graph_Graph513->SetMinimum(99.9);
   Graph_Graph513->SetMaximum(101.1);
   Graph_Graph513->SetDirectory(nullptr);
   Graph_Graph513->SetStats(0);
   Graph_Graph513->SetLineWidth(2);
   Graph_Graph513->SetMarkerStyle(20);
   Graph_Graph513->SetMarkerSize(0.9);
   Graph_Graph513->GetXaxis()->SetLabelFont(42);
   Graph_Graph513->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetXaxis()->SetTitleFont(42);
   Graph_Graph513->GetYaxis()->SetLabelFont(42);
   Graph_Graph513->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetYaxis()->SetTickLength(0.02);
   Graph_Graph513->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetYaxis()->SetTitleFont(42);
   Graph_Graph513->GetZaxis()->SetLabelFont(42);
   Graph_Graph513->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph513->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph513->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph513->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph513->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph513);
   
   graph->Draw("p");
   
   Double_t Graph6_fx14[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy14[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 103.654, 74.55099, 223.6068, 223.6068, 28.17341, 89.60343, 35.76581, 57.06507, 13.97268, 27.70181,
   57.92681, 86.6768, 67.15957, 61.65793, 33.98856, 42.0448, 34.94133, 36.04335, 17.68151, 11.72155, 7.911233, 6.586652, 6.627826, 8.510923, 11.26443, 9.739021,
   26.22753, 20.62761 };
   graph = new TGraph(35,Graph6_fx14,Graph6_fy14);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph614 = new TH1F("Graph_Graph614","Graph",100,0,3520);
   Graph_Graph614->SetMinimum(5.927987);
   Graph_Graph614->SetMaximum(245.3088);
   Graph_Graph614->SetDirectory(nullptr);
   Graph_Graph614->SetStats(0);
   Graph_Graph614->SetLineWidth(2);
   Graph_Graph614->SetMarkerStyle(20);
   Graph_Graph614->SetMarkerSize(0.9);
   Graph_Graph614->GetXaxis()->SetLabelFont(42);
   Graph_Graph614->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetXaxis()->SetTitleFont(42);
   Graph_Graph614->GetYaxis()->SetLabelFont(42);
   Graph_Graph614->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetYaxis()->SetTickLength(0.02);
   Graph_Graph614->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetYaxis()->SetTitleFont(42);
   Graph_Graph614->GetZaxis()->SetLabelFont(42);
   Graph_Graph614->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph614->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph614->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph614->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph614->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph614);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1200 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
