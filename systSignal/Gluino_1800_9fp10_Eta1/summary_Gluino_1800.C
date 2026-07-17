#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1800()
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
   
   Double_t Graph0_fx36[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy36[35] = { 100, 100, 100, 100, 0, 0, 91.52025, 100, 0, 0, 0, 0, 0, 0, 100, 9.674478, 1.192093e-05,
   0, 0, 0, 29.16968, 16.1471, 27.79392, 10.27325, 0, 6.846368, 2.65128, 6.093991, 5.473095, 1.931608, 1.040268, 1.6366, 2.790272,
   3.432441, 1.81511 };
   TGraph *graph = new TGraph(35,Graph0_fx36,Graph0_fy36);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph036 = new TH1F("Graph_Graph036","Graph",100,0,3520);
   Graph_Graph036->SetMinimum(1);
   Graph_Graph036->SetMaximum(2000);
   Graph_Graph036->SetDirectory(nullptr);
   Graph_Graph036->SetStats(0);
   Graph_Graph036->SetLineWidth(2);
   Graph_Graph036->SetMarkerStyle(20);
   Graph_Graph036->SetMarkerSize(0.9);
   Graph_Graph036->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph036->GetXaxis()->SetRange(1,101);
   Graph_Graph036->GetXaxis()->SetLabelFont(43);
   Graph_Graph036->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetXaxis()->SetLabelSize(22);
   Graph_Graph036->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph036->GetXaxis()->SetTitleOffset(1);
   Graph_Graph036->GetXaxis()->SetTitleFont(42);
   Graph_Graph036->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph036->GetYaxis()->SetLabelFont(43);
   Graph_Graph036->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetYaxis()->SetLabelSize(22);
   Graph_Graph036->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph036->GetYaxis()->SetTickLength(0.02);
   Graph_Graph036->GetYaxis()->SetTitleOffset(1);
   Graph_Graph036->GetYaxis()->SetTitleFont(42);
   Graph_Graph036->GetZaxis()->SetLabelFont(42);
   Graph_Graph036->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph036->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph036->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph036->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph036->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph036);
   
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
   
   Double_t Graph1_fx37[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy37[35] = { 100, 100, 100, 100, 0, 0, 100, 100, 0, 29.3586, 0, 0, 0, 0, 100, 9.674478, 1.192093e-05,
   0, 0, 0, 15.46319, 8.55977, 37.05569, 9.340686, 9.084207, 7.384884, 2.024436, 4.589808, 3.954387, 1.798296, 0.6044865, 1.974136, 2.300024,
   5.037272, 0.6084979 };
   graph = new TGraph(35,Graph1_fx37,Graph1_fy37);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph137 = new TH1F("Graph_Graph137","Graph",100,0,3520);
   Graph_Graph137->SetMinimum(99.9);
   Graph_Graph137->SetMaximum(101.1);
   Graph_Graph137->SetDirectory(nullptr);
   Graph_Graph137->SetStats(0);
   Graph_Graph137->SetLineWidth(2);
   Graph_Graph137->SetMarkerStyle(20);
   Graph_Graph137->SetMarkerSize(0.9);
   Graph_Graph137->GetXaxis()->SetLabelFont(42);
   Graph_Graph137->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetXaxis()->SetTitleFont(42);
   Graph_Graph137->GetYaxis()->SetLabelFont(42);
   Graph_Graph137->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetYaxis()->SetTickLength(0.02);
   Graph_Graph137->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetYaxis()->SetTitleFont(42);
   Graph_Graph137->GetZaxis()->SetLabelFont(42);
   Graph_Graph137->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph137->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph137->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph137->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph137->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph137);
   
   graph->Draw("p");
   
   Double_t Graph2_fx38[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy38[35] = { 100, 100, 100, 100, 45.27677, 31.28117, 31.28117, 34.71781, 39.68093, 17.98391, 19.5901, 15.70612, 39.81429, 12.63813, 52.13418, 67.17352, 8.016128,
   38.16611, 10.38067, 17.02496, 108.8573, 17.48975, 28.73767, 7.428283, 7.974732, 10.72679, 6.540179, 29.09931, 1.368845, 1.358956, 0.3907979, 1.667345, 1.686966,
   9.609842, 6.116491 };
   graph = new TGraph(35,Graph2_fx38,Graph2_fy38);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph238 = new TH1F("Graph_Graph238","Graph",100,0,3520);
   Graph_Graph238->SetMinimum(0.3517181);
   Graph_Graph238->SetMaximum(119.7039);
   Graph_Graph238->SetDirectory(nullptr);
   Graph_Graph238->SetStats(0);
   Graph_Graph238->SetLineWidth(2);
   Graph_Graph238->SetMarkerStyle(20);
   Graph_Graph238->SetMarkerSize(0.9);
   Graph_Graph238->GetXaxis()->SetLabelFont(42);
   Graph_Graph238->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetXaxis()->SetTitleFont(42);
   Graph_Graph238->GetYaxis()->SetLabelFont(42);
   Graph_Graph238->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetYaxis()->SetTickLength(0.02);
   Graph_Graph238->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetYaxis()->SetTitleFont(42);
   Graph_Graph238->GetZaxis()->SetLabelFont(42);
   Graph_Graph238->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph238->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph238->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph238->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph238->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph238);
   
   graph->Draw("p");
   
   Double_t Graph3_fx39[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy39[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.05499, 1.689011, 0.8310199, 0.8578897, 1.463223, 2.308369,
   3.783691, 2.923024 };
   graph = new TGraph(35,Graph3_fx39,Graph3_fy39);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph339 = new TH1F("Graph_Graph339","Graph",100,0,3520);
   Graph_Graph339->SetMinimum(99.9);
   Graph_Graph339->SetMaximum(101.1);
   Graph_Graph339->SetDirectory(nullptr);
   Graph_Graph339->SetStats(0);
   Graph_Graph339->SetLineWidth(2);
   Graph_Graph339->SetMarkerStyle(20);
   Graph_Graph339->SetMarkerSize(0.9);
   Graph_Graph339->GetXaxis()->SetLabelFont(42);
   Graph_Graph339->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetXaxis()->SetTitleFont(42);
   Graph_Graph339->GetYaxis()->SetLabelFont(42);
   Graph_Graph339->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetYaxis()->SetTickLength(0.02);
   Graph_Graph339->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetYaxis()->SetTitleFont(42);
   Graph_Graph339->GetZaxis()->SetLabelFont(42);
   Graph_Graph339->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph339->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph339->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph339->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph339->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph339);
   
   graph->Draw("p");
   
   Double_t Graph4_fx40[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy40[35] = { 100, 100, 100, 100, 1.60315, 6.518424, 3.647768, 3.647768, 2.532613, 4.56531, 2.148604, 2.846038, 4.08327, 3.538215, 2.042818, 3.572798, 2.174318,
   2.978152, 2.70223, 2.042806, 1.339865, 3.034556, 3.072476, 2.976775, 2.556551, 3.166258, 2.577269, 3.291154, 3.277087, 3.348184, 3.449142, 3.406608, 3.505552,
   3.541446, 3.525245 };
   graph = new TGraph(35,Graph4_fx40,Graph4_fy40);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph440 = new TH1F("Graph_Graph440","Graph",100,0,3520);
   Graph_Graph440->SetMinimum(1.205878);
   Graph_Graph440->SetMaximum(109.866);
   Graph_Graph440->SetDirectory(nullptr);
   Graph_Graph440->SetStats(0);
   Graph_Graph440->SetLineWidth(2);
   Graph_Graph440->SetMarkerStyle(20);
   Graph_Graph440->SetMarkerSize(0.9);
   Graph_Graph440->GetXaxis()->SetLabelFont(42);
   Graph_Graph440->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetXaxis()->SetTitleFont(42);
   Graph_Graph440->GetYaxis()->SetLabelFont(42);
   Graph_Graph440->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetYaxis()->SetTickLength(0.02);
   Graph_Graph440->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetYaxis()->SetTitleFont(42);
   Graph_Graph440->GetZaxis()->SetLabelFont(42);
   Graph_Graph440->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph440->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph440->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph440->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph440->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph440);
   
   graph->Draw("p");
   
   Double_t Graph5_fx41[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy41[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 5.10549, 1.42405, 0.454855, 0.6579161, 0.2761424, 0.4017472, 0.2304077,
   0, 0 };
   graph = new TGraph(35,Graph5_fx41,Graph5_fy41);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph541 = new TH1F("Graph_Graph541","Graph",100,0,3520);
   Graph_Graph541->SetMinimum(99.9);
   Graph_Graph541->SetMaximum(101.1);
   Graph_Graph541->SetDirectory(nullptr);
   Graph_Graph541->SetStats(0);
   Graph_Graph541->SetLineWidth(2);
   Graph_Graph541->SetMarkerStyle(20);
   Graph_Graph541->SetMarkerSize(0.9);
   Graph_Graph541->GetXaxis()->SetLabelFont(42);
   Graph_Graph541->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetXaxis()->SetTitleFont(42);
   Graph_Graph541->GetYaxis()->SetLabelFont(42);
   Graph_Graph541->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetYaxis()->SetTickLength(0.02);
   Graph_Graph541->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetYaxis()->SetTitleFont(42);
   Graph_Graph541->GetZaxis()->SetLabelFont(42);
   Graph_Graph541->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph541->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph541->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph541->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph541->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph541);
   
   graph->Draw("p");
   
   Double_t Graph6_fx42[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy42[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 45.30515, 31.95312, 139.1682, 145.6662, 39.76167, 34.73026, 19.70758, 15.9619, 40.02313, 13.12408, 150.7387, 68.64574, 8.305779,
   38.28212, 10.72662, 17.14708, 113.7615, 25.47741, 54.59781, 16.02588, 12.35537, 15.04986, 7.780995, 30.28065, 7.813945, 4.551117, 3.772644, 4.806424, 5.79129,
   12.50452, 7.877041 };
   graph = new TGraph(35,Graph6_fx42,Graph6_fy42);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph642 = new TH1F("Graph_Graph642","Graph",100,0,3520);
   Graph_Graph642->SetMinimum(3.395379);
   Graph_Graph642->SetMaximum(245.5902);
   Graph_Graph642->SetDirectory(nullptr);
   Graph_Graph642->SetStats(0);
   Graph_Graph642->SetLineWidth(2);
   Graph_Graph642->SetMarkerStyle(20);
   Graph_Graph642->SetMarkerSize(0.9);
   Graph_Graph642->GetXaxis()->SetLabelFont(42);
   Graph_Graph642->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetXaxis()->SetTitleFont(42);
   Graph_Graph642->GetYaxis()->SetLabelFont(42);
   Graph_Graph642->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetYaxis()->SetTickLength(0.02);
   Graph_Graph642->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetYaxis()->SetTitleFont(42);
   Graph_Graph642->GetZaxis()->SetLabelFont(42);
   Graph_Graph642->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph642->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph642->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph642->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph642->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph642);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=1800 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
