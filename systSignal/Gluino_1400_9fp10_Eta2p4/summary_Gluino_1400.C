#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1400()
{
//=========Macro generated from canvas: c2/c2
//=========  (Thu Jul 16 15:10:02 2026) by ROOT version 6.32.13
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
   
   Double_t Graph0_fx22[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy22[35] = { 100, 100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 22.63669, 20.05178, 31.69736, 44.78931, 8.791912, 6.74628,
   7.428622, 15.40499, 30.41431, 6.291127, 18.29327, 6.15583, 9.42849, 6.819701, 2.002287, 4.18607, 3.129047, 1.409566, 1.745915, 3.374195, 2.667844, 0.7467866,
   3.586113, 0 };
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
   Graph_Graph022->GetXaxis()->SetLabelSize(22);
   Graph_Graph022->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph022->GetXaxis()->SetTitleOffset(1);
   Graph_Graph022->GetXaxis()->SetTitleFont(42);
   Graph_Graph022->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph022->GetYaxis()->SetLabelFont(43);
   Graph_Graph022->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph022->GetYaxis()->SetLabelSize(22);
   Graph_Graph022->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph022->GetYaxis()->SetTickLength(0.02);
   Graph_Graph022->GetYaxis()->SetTitleOffset(1);
   Graph_Graph022->GetYaxis()->SetTitleFont(42);
   Graph_Graph022->GetZaxis()->SetLabelFont(42);
   Graph_Graph022->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph022->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph022->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph022->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph022->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph022);
   
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
   
   Double_t Graph1_fx23[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy23[35] = { 100, 100, 67.8756, 73.9323, 40.52021, 100, 132.9389, 75.42682, 84.73432, 34.52751, 86.28674, 1.192093e-05, 0, 76.0998, 62.96245, 2.331793, 0.116992,
   0.8072376, 9.856569, 51.20547, 20.02657, 45.23406, 3.029704, 20.47508, 9.81288, 4.818344, 7.446384, 4.535145, 1.398546, 1.770806, 3.678602, 2.901483, 3.252316,
   3.083444, 1.57547 };
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
   Graph_Graph123->SetMinimum(1.072884e-05);
   Graph_Graph123->SetMaximum(146.2328);
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
   Double_t Graph2_fy24[35] = { 100, 100, 43.9416, 24.34452, 43.26846, 100, 42.82398, 7.311451, 20.84246, 7.062716, 14.54738, 25.82039, 10.45403, 23.27374, 29.92421, 10.4367, 21.33442,
   12.75587, 13.0273, 12.5986, 13.26053, 15.91367, 13.27342, 11.68591, 0.6696463, 2.269441, 5.837476, 2.924526, 3.409493, 1.552749, 4.730797, 4.418206, 3.859425,
   3.938991, 4.616988 };
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
   Graph_Graph224->SetMinimum(0.6026816);
   Graph_Graph224->SetMaximum(109.933);
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
   Double_t Graph3_fy25[35] = { 100, 100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 7.48446, 0, 0, 0, 0, 0, 2.870059, 1.096314, 1.590955, 1.860666, 1.979339, 1.535833, 2.371156, 1.795328, 3.179049,
   2.615309, 4.163229 };
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
   Double_t Graph4_fy26[35] = { 100, 100, 2.507889, 4.396004, 3.647757, 100, 3.188342, 2.243125, 2.946782, 3.050661, 2.769864, 3.339934, 2.10036, 2.324045, 2.736497, 3.068829, 1.512444,
   1.175761, 3.151894, 2.10588, 3.139853, 2.808595, 3.266883, 2.885878, 3.573966, 2.764344, 2.977002, 3.216648, 3.347266, 3.412497, 3.388548, 3.496385, 3.421187,
   3.443944, 3.477418 };
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
   Graph_Graph426->SetMinimum(1.058185);
   Graph_Graph426->SetMaximum(109.8824);
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
   Double_t Graph5_fy27[35] = { 100, 100, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   38.80337, 0, 38.74492, 6.330031, 0, 0, 0, 0, 4.834855, 0.4640043, 0.3373921, 0.6105304, 0.4301906, 0.2192438, 0.3253639, 0,
   0.9370208, 0.8708954 };
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
   Double_t Graph6_fy28[35] = { 223.6068, 223.6068, 80.89655, 77.96131, 59.39153, 223.6068, 139.7026, 75.81355, 87.30978, 35.37425, 87.54827, 34.50025, 22.71061, 85.69113, 82.90536, 14.18017, 22.42702,
   14.83006, 23.87735, 60.91132, 25.02682, 51.39939, 15.29476, 25.55415, 12.81639, 6.420238, 10.88309, 7.26087, 5.539804, 4.753832, 8.025109, 7.106648, 6.916739,
   7.522239, 7.295451 };
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
   Graph_Graph628->SetMinimum(4.278449);
   Graph_Graph628->SetMaximum(245.4921);
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
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1400 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
