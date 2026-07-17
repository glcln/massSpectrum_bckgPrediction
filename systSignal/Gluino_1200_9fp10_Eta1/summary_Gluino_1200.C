#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1200()
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
   
   Double_t Graph0_fx8[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy8[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 100, 16.00609, 0, 0, 0, 100, 0, 15.38267, 100,
   0, 29.52772, 11.41744, 5.960464e-06, 6.050629, 1.128554, 8.604413, 1.421249, 5.563247, 6.31395, 2.181888, 1.488197, 3.654504, 3.679556, 2.101392, 0.239706,
   4.133624, 0 };
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
   Double_t Graph1_fy9[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 100, 16.00609, 0, 25.59274, 4.0205, 100, 0, 0, 0,
   0, 59.25233, 27.64891, 32.58535, 20.53813, 30.23638, 1.247412, 5.754137, 4.797208, 5.926359, 2.325177, 1.628882, 4.362965, 6.939131, 4.178786, 4.152232,
   4.133624, 0 };
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
   Graph_Graph19->SetMinimum(99.9);
   Graph_Graph19->SetMaximum(101.1);
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
   Double_t Graph2_fy10[35] = { 100, 100, 100, 100, 20.74046, 100, 20.79192, 22.16222, 17.02496, 28.08389, 34.30193, 43.68462, 74.09564, 100, 32.36191, 22.29773, 27.74897,
   17.30343, 35.53002, 13.61142, 14.51739, 8.320803, 25.96807, 5.35351, 5.255604, 4.684424, 0.9977698, 1.555419, 1.740324, 3.403485, 3.393281, 8.964193, 15.01254,
   12.50402, 7.068288 };
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
   Graph_Graph210->SetMinimum(0.8979928);
   Graph_Graph210->SetMaximum(109.9002);
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
   Double_t Graph3_fy11[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 3.37882, 2.795154, 2.285743, 1.603848, 1.868331, 2.173519, 2.903914, 2.488983, 4.310787,
   3.37106, 1.015484 };
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
   Double_t Graph4_fy12[35] = { 100, 100, 100, 100, 0.9231091, 100, 2.684402, 2.027321, 1.60315, 2.311772, 3.169811, 1.771796, 2.976382, 100, 1.073658, 1.739252, 3.647768,
   2.113044, 2.490997, 2.820539, 3.160357, 3.594947, 3.414226, 2.962112, 2.872229, 3.00318, 3.192484, 3.360581, 3.375757, 3.349328, 3.432393, 3.876638, 3.588569,
   3.477967, 3.792882 };
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
   Graph_Graph412->SetMinimum(0.8307981);
   Graph_Graph412->SetMaximum(109.9077);
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
   Double_t Graph5_fy13[35] = { 100, 100, 100, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 1.224458, 1.180363, 0.509268, 0.4441142, 0.3318012, 0.5279243, 0.7641375, 0, 0,
   0, 0 };
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
   Double_t Graph6_fy14[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 102.1324, 223.6068, 20.9645, 22.25475, 142.4515, 36.14469, 34.44807, 50.66038, 74.26431, 223.6068, 32.37971, 27.14481, 103.8427,
   17.43197, 75.17522, 32.98555, 35.81268, 23.25047, 40.01887, 10.63138, 9.07841, 9.630095, 9.560325, 5.143175, 4.773166, 7.740591, 9.665354, 11.11118, 16.55704,
   14.62821, 8.085658 };
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
   Graph_Graph614->SetMinimum(4.295849);
   Graph_Graph614->SetMaximum(245.4902);
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
