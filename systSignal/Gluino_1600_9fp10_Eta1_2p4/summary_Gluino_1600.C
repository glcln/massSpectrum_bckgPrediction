#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy29[35] = { 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 0, 0, 100, 0, 100,
   21.25779, 33.85832, 20.15773, 38.32642, 0, 0, 13.68453, 6.382871, 9.692204, 2.722096, 3.535199, 3.923368, 0.5260885, 2.124631, 1.942325, 2.618206,
   2.307177, 2.589381 };
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
   Double_t Graph1_fy30[35] = { 100, 100, 100, 100, 100, 0, 100, 100, 100, 0, 100, 19.68638, 58.95676, 71.66038, 100, 22.42421, 100,
   79.17995, 126.1138, 20.15773, 67.42649, 22.29033, 0, 52.86162, 8.765483, 5.27339, 14.00719, 10.39905, 5.053163, 0.9661257, 3.382772, 3.22994, 8.675539,
   4.576719, 2.589381 };
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
   Graph_Graph130->SetMinimum(97.38862);
   Graph_Graph130->SetMaximum(128.7251);
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
   Double_t Graph2_fy31[35] = { 100, 54.12139, 100, 100, 100, 29.92773, 100, 100, 100, 36.49092, 100, 27.91953, 28.89853, 36.98486, 100, 33.17332, 100,
   5.485952, 22.81462, 33.45159, 26.13817, 19.46707, 8.023024, 36.23886, 25.87402, 18.59411, 8.022177, 1.014197, 1.295435, 5.918407, 3.393245, 4.631507, 13.42086,
   4.154486, 1.63489 };
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
   Graph_Graph231->SetMinimum(0.9127772);
   Graph_Graph231->SetMaximum(109.8986);
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
   Double_t Graph3_fy32[35] = { 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 0, 0, 100, 18.57879, 100,
   0, 0, 0, 0, 0, 0, 0, 0, 12.59104, 0.9933174, 2.246296, 2.136141, 1.975983, 2.51922, 3.809834, 4.16891,
   5.198484, 1.383936 };
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
   Graph_Graph332->SetMinimum(0.11);
   Graph_Graph332->SetMaximum(110);
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
   Double_t Graph4_fy33[35] = { 100, 2.507889, 100, 100, 100, 1.234829, 100, 100, 100, 3.396153, 100, 4.081607, 3.647768, 3.495669, 100, 1.16874, 100,
   2.720618, 2.386212, 2.82352, 1.366568, 1.535571, 2.71616, 2.03321, 2.403927, 2.577233, 2.597797, 2.852166, 3.174901, 3.351688, 3.328001, 3.35412, 3.279746,
   3.373206, 3.692949 };
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
   Graph_Graph433->SetMinimum(1.051866);
   Graph_Graph433->SetMaximum(109.8831);
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
   Double_t Graph5_fy34[35] = { 100, 0, 100, 100, 100, 0, 100, 100, 100, 0, 100, 0, 0, 0, 100, 0, 100,
   0, 32.02876, 0, 0, 0, 0, 0, 0, 0, 2.313244, 1.063222, 0.5443454, 0.406605, 0.7764697, 0.8736014, 0.8978605,
   0, 0 };
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
   Graph_Graph534->SetMinimum(0.11);
   Graph_Graph534->SetMaximum(110);
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
   Double_t Graph6_fy35[35] = { 223.6068, 113.734, 223.6068, 223.6068, 223.6068, 29.95319, 223.6068, 223.6068, 223.6068, 36.64862, 223.6068, 34.40514, 65.75964, 80.71747, 223.6068, 44.15713, 223.6068,
   82.21226, 132.5793, 44.04145, 81.85547, 29.63416, 8.470326, 65.5668, 28.15703, 25.15283, 16.60427, 11.61237, 7.566275, 7.167708, 6.700252, 7.837152, 17.04042,
   9.051625, 5.624571 };
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
   Graph_Graph635->SetMinimum(5.062114);
   Graph_Graph635->SetMaximum(245.405);
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
