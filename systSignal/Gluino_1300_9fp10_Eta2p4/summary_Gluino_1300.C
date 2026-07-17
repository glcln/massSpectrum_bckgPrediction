#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1300()
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
   
   Double_t Graph0_fx15[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy15[35] = { 100, 100, 0, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 12.29782, 14.7956, 6.722445, 7.258713,
   4.254085, 5.052829, 13.58404, 39.35975, 2.41785, 6.683284, 3.048384, 9.89772, 3.688931, 3.107214, 2.639627, 0.2644122, 2.854633, 2.489156, 1.720977, 3.906941,
   3.009009, 2.331638 };
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
   Double_t Graph1_fy16[35] = { 100, 100, 100, 100, 35.99813, 100, 51.57832, 117.2107, 0, 7.340801, 0, 29.05661, 27.9246, 28.93497, 77.65378, 39.73004, 13.59254,
   23.12442, 46.42174, 5.960464e-06, 98.65285, 13.38307, 7.781667, 8.862609, 13.24906, 4.72213, 5.088079, 3.859144, 0.5461335, 3.645301, 4.125542, 6.070304, 3.250933,
   6.327063, 2.331638 };
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
   Graph_Graph116->SetMinimum(27.87687);
   Graph_Graph116->SetMaximum(125.332);
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
   Double_t Graph2_fy17[35] = { 100, 100, 52.13418, 100, 28.49717, 100, 7.604379, 22.47303, 31.44739, 29.03594, 27.03399, 34.71781, 27.78487, 25.24585, 33.05026, 10.07112, 20.67612,
   12.96394, 24.88776, 31.14432, 33.41825, 24.62463, 15.05164, 10.81446, 11.50323, 13.64291, 1.593649, 3.848904, 2.035761, 4.684204, 7.253504, 2.901322, 15.03831,
   6.236303, 7.714725 };
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
   Graph_Graph217->SetMinimum(1.434284);
   Graph_Graph217->SetMaximum(109.8406);
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
   Double_t Graph3_fy18[35] = { 100, 100, 0, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 6.291044, 0, 0, 0, 1.302797, 2.692533, 3.467488, 2.270913, 1.977909, 1.548666, 2.111286, 2.767706, 1.346737, 1.73952,
   4.513651, 2.968144 };
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
   Double_t Graph4_fy19[35] = { 100, 100, 3.722632, 100, 2.103448, 100, 5.128348, 1.326168, 1.737231, 4.534692, 1.234829, 1.234829, 3.325188, 2.875352, 3.898442, 3.43678, 2.443528,
   2.188873, 2.440548, 2.196968, 2.793515, 2.578235, 3.095949, 2.346778, 3.22274, 2.770662, 3.098047, 3.325248, 3.367162, 3.408682, 3.45211, 3.397489, 3.132796,
   3.502989, 3.445864 };
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
   Graph_Graph419->SetMinimum(1.111346);
   Graph_Graph419->SetMaximum(109.8765);
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
   Double_t Graph5_fy20[35] = { 100, 100, 0, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 12.98427, 1.903564, 0, 0, 0.3402531, 0.5976617, 0.791961, 0.6626606, 0.3778219, 0.3578663, 0, 0,
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
   Double_t Graph6_fy21[35] = { 223.6068, 223.6068, 112.8354, 223.6068, 45.96062, 223.6068, 52.3875, 119.353, 31.49534, 30.29086, 27.06218, 45.28949, 39.53276, 40.42387, 85.77026, 41.67621, 25.9021,
   26.93865, 52.97043, 34.62512, 111.3829, 28.2484, 18.47586, 14.56007, 20.57821, 15.54785, 7.268974, 7.186382, 4.271843, 7.710743, 9.767629, 7.847589, 16.27341,
   10.98277, 9.543265 };
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
   Graph_Graph621->SetMinimum(3.844659);
   Graph_Graph621->SetMaximum(245.5403);
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
