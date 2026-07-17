#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1300()
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
   
   Double_t Graph0_fx15[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy15[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 8.192503, 18.20425,
   5.276948, 11.36873, 0, 40.28271, 4.870582, 8.121133, 7.298446, 4.695714, 7.466316, 2.491522, 3.917253, 1.313657, 2.246249, 3.066552, 2.88592, 5.688977,
   0, 0 };
   TGraph *graph = new TGraph(35,Graph0_fx15,Graph0_fy15);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph015 = new TH1F("Graph_Graph015","Graph",100,0,3520);
   Graph_Graph015->SetMinimum(1);
   Graph_Graph015->SetMaximum(2000);
   Graph_Graph015->SetDirectory(nullptr);
   Graph_Graph015->SetStats(0);
   Graph_Graph015->SetLineWidth(2);
   Graph_Graph015->SetMarkerStyle(20);
   Graph_Graph015->SetMarkerSize(0.9);
   Graph_Graph015->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph015->GetXaxis()->SetRange(1,101);
   Graph_Graph015->GetXaxis()->SetLabelFont(43);
   Graph_Graph015->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetXaxis()->SetLabelSize(22);
   Graph_Graph015->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph015->GetXaxis()->SetTitleOffset(1);
   Graph_Graph015->GetXaxis()->SetTitleFont(42);
   Graph_Graph015->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph015->GetYaxis()->SetLabelFont(43);
   Graph_Graph015->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetYaxis()->SetLabelSize(22);
   Graph_Graph015->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph015->GetYaxis()->SetTickLength(0.02);
   Graph_Graph015->GetYaxis()->SetTitleOffset(1);
   Graph_Graph015->GetYaxis()->SetTitleFont(42);
   Graph_Graph015->GetZaxis()->SetLabelFont(42);
   Graph_Graph015->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph015->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph015->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph015->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph015->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph015);
   
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
   
   Double_t Graph1_fx16[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy16[35] = { 100, 100, 100, 100, 100, 100, 100, 100, 0, 0, 0, 100, 31.34479, 58.33184, 131.4963, 48.41817, 9.548545,
   28.6845, 67.00005, 0, 170.1121, 11.26991, 19.99598, 8.035469, 13.39455, 6.738889, 5.777061, 6.39863, 0.6334543, 5.66029, 4.047841, 2.379346, 5.937881,
   5.20668, 0 };
   graph = new TGraph(35,Graph1_fx16,Graph1_fy16);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph116 = new TH1F("Graph_Graph116","Graph",100,0,3520);
   Graph_Graph116->SetMinimum(92.98879);
   Graph_Graph116->SetMaximum(177.1233);
   Graph_Graph116->SetDirectory(nullptr);
   Graph_Graph116->SetStats(0);
   Graph_Graph116->SetLineWidth(2);
   Graph_Graph116->SetMarkerStyle(20);
   Graph_Graph116->SetMarkerSize(0.9);
   Graph_Graph116->GetXaxis()->SetLabelFont(42);
   Graph_Graph116->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetXaxis()->SetTitleFont(42);
   Graph_Graph116->GetYaxis()->SetLabelFont(42);
   Graph_Graph116->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetYaxis()->SetTickLength(0.02);
   Graph_Graph116->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetYaxis()->SetTitleFont(42);
   Graph_Graph116->GetZaxis()->SetLabelFont(42);
   Graph_Graph116->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph116->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph116->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph116->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph116->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph116);
   
   graph->Draw("p");
   
   Double_t Graph2_fx17[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy17[35] = { 100, 100, 52.13418, 100, 33.25831, 100, 100, 31.28117, 51.08392, 19.77257, 27.03399, 100, 26.36795, 20.74047, 35.1756, 9.365487, 38.03018,
   26.64644, 22.20877, 30.83687, 29.52177, 23.04423, 14.43214, 7.524741, 8.734221, 8.113867, 6.00642, 0.9299636, 4.560661, 3.840905, 1.953101, 9.555298, 20.80399,
   39.10294, 9.205246 };
   graph = new TGraph(35,Graph2_fx17,Graph2_fy17);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph217 = new TH1F("Graph_Graph217","Graph",100,0,3520);
   Graph_Graph217->SetMinimum(0.8369672);
   Graph_Graph217->SetMaximum(109.907);
   Graph_Graph217->SetDirectory(nullptr);
   Graph_Graph217->SetStats(0);
   Graph_Graph217->SetLineWidth(2);
   Graph_Graph217->SetMarkerStyle(20);
   Graph_Graph217->SetMarkerSize(0.9);
   Graph_Graph217->GetXaxis()->SetLabelFont(42);
   Graph_Graph217->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetXaxis()->SetTitleFont(42);
   Graph_Graph217->GetYaxis()->SetLabelFont(42);
   Graph_Graph217->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetYaxis()->SetTickLength(0.02);
   Graph_Graph217->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetYaxis()->SetTitleFont(42);
   Graph_Graph217->GetZaxis()->SetLabelFont(42);
   Graph_Graph217->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph217->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph217->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph217->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph217->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph217);
   
   graph->Draw("p");
   
   Double_t Graph3_fx18[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy18[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0,
   0, 0, 6.648266, 0, 0, 0, 0, 4.695714, 5.076331, 1.752335, 3.370118, 2.5998, 3.717053, 4.694486, 3.37671, 4.185593,
   3.797567, 8.153391 };
   graph = new TGraph(35,Graph3_fx18,Graph3_fy18);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph318 = new TH1F("Graph_Graph318","Graph",100,0,3520);
   Graph_Graph318->SetMinimum(99.9);
   Graph_Graph318->SetMaximum(101.1);
   Graph_Graph318->SetDirectory(nullptr);
   Graph_Graph318->SetStats(0);
   Graph_Graph318->SetLineWidth(2);
   Graph_Graph318->SetMarkerStyle(20);
   Graph_Graph318->SetMarkerSize(0.9);
   Graph_Graph318->GetXaxis()->SetLabelFont(42);
   Graph_Graph318->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetXaxis()->SetTitleFont(42);
   Graph_Graph318->GetYaxis()->SetLabelFont(42);
   Graph_Graph318->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetYaxis()->SetTickLength(0.02);
   Graph_Graph318->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetYaxis()->SetTitleFont(42);
   Graph_Graph318->GetZaxis()->SetLabelFont(42);
   Graph_Graph318->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph318->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph318->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph318->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph318->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph318);
   
   graph->Draw("p");
   
   Double_t Graph4_fx19[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy19[35] = { 100, 100, 3.722632, 100, 3.647768, 100, 100, 1.60315, 0.9231031, 5.533356, 1.234829, 100, 3.28567, 1.234829, 3.682387, 3.781223, 2.359426,
   2.123082, 1.919317, 2.427483, 2.505225, 2.285254, 3.17111, 2.00724, 2.971768, 2.865875, 2.882123, 3.280199, 3.274012, 3.334939, 3.295135, 3.360319, 3.071773,
   3.133512, 3.005075 };
   graph = new TGraph(35,Graph4_fx19,Graph4_fy19);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph419 = new TH1F("Graph_Graph419","Graph",100,0,3520);
   Graph_Graph419->SetMinimum(0.8307928);
   Graph_Graph419->SetMaximum(109.9077);
   Graph_Graph419->SetDirectory(nullptr);
   Graph_Graph419->SetStats(0);
   Graph_Graph419->SetLineWidth(2);
   Graph_Graph419->SetMarkerStyle(20);
   Graph_Graph419->SetMarkerSize(0.9);
   Graph_Graph419->GetXaxis()->SetLabelFont(42);
   Graph_Graph419->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetXaxis()->SetTitleFont(42);
   Graph_Graph419->GetYaxis()->SetLabelFont(42);
   Graph_Graph419->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetYaxis()->SetTickLength(0.02);
   Graph_Graph419->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetYaxis()->SetTitleFont(42);
   Graph_Graph419->GetZaxis()->SetLabelFont(42);
   Graph_Graph419->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph419->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph419->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph419->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph419->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph419);
   
   graph->Draw("p");
   
   Double_t Graph5_fx20[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy20[35] = { 100, 100, 0, 100, 0, 100, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 16.64824,
   0, 0, 0, 0, 24.50599, 0, 0, 0, 0, 0, 0.4009306, 0.6878138, 0.6267786, 0.3701031, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph5_fx20,Graph5_fy20);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph520 = new TH1F("Graph_Graph520","Graph",100,0,3520);
   Graph_Graph520->SetMinimum(99.9);
   Graph_Graph520->SetMaximum(101.1);
   Graph_Graph520->SetDirectory(nullptr);
   Graph_Graph520->SetStats(0);
   Graph_Graph520->SetLineWidth(2);
   Graph_Graph520->SetMarkerStyle(20);
   Graph_Graph520->SetMarkerSize(0.9);
   Graph_Graph520->GetXaxis()->SetLabelFont(42);
   Graph_Graph520->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetXaxis()->SetTitleFont(42);
   Graph_Graph520->GetYaxis()->SetLabelFont(42);
   Graph_Graph520->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetYaxis()->SetTickLength(0.02);
   Graph_Graph520->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetYaxis()->SetTitleFont(42);
   Graph_Graph520->GetZaxis()->SetLabelFont(42);
   Graph_Graph520->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph520->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph520->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph520->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph520->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph520);
   
   graph->Draw("p");
   
   Double_t Graph6_fx21[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy21[35] = { 223.6068, 223.6068, 112.8354, 223.6068, 105.4487, 223.6068, 223.6068, 104.7907, 51.09226, 20.53223, 27.06218, 223.6068, 41.09209, 61.92168, 136.1696, 50.13428, 43.2947,
   39.56245, 71.52041, 31.63866, 177.3094, 26.21054, 26.15596, 13.3599, 17.56791, 14.17659, 9.329335, 8.903346, 6.356469, 8.762156, 7.905661, 11.31315, 22.96484,
   39.75411, 12.65878 };
   graph = new TGraph(35,Graph6_fx21,Graph6_fy21);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph621 = new TH1F("Graph_Graph621","Graph",100,0,3520);
   Graph_Graph621->SetMinimum(5.720822);
   Graph_Graph621->SetMaximum(245.3318);
   Graph_Graph621->SetDirectory(nullptr);
   Graph_Graph621->SetStats(0);
   Graph_Graph621->SetLineWidth(2);
   Graph_Graph621->SetMarkerStyle(20);
   Graph_Graph621->SetMarkerSize(0.9);
   Graph_Graph621->GetXaxis()->SetLabelFont(42);
   Graph_Graph621->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetXaxis()->SetTitleFont(42);
   Graph_Graph621->GetYaxis()->SetLabelFont(42);
   Graph_Graph621->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetYaxis()->SetTickLength(0.02);
   Graph_Graph621->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetYaxis()->SetTitleFont(42);
   Graph_Graph621->GetZaxis()->SetLabelFont(42);
   Graph_Graph621->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph621->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph621->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph621->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph621->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph621);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1300 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
