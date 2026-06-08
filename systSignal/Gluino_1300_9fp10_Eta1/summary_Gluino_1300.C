#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1300()
{
//=========Macro generated from canvas: c1/c1
//=========  (Fri Jun  5 11:29:21 2026) by ROOT version 6.32.13
   TCanvas *c1 = new TCanvas("c1", "c1",0,0,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->SetHighLightColor(2);
   c1->Range(-720.0405,-0.668563,3780.213,3.509956);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLeftMargin(0.16);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   
   Double_t Graph0_fx13[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy13[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 1.192093e-05, 100, 78.26144, 68.87587, 39.28905, 21.32539, 113.5812, 13.57679,
   52.04154, 7.693768, 1.551235, 50.58646, 26.53029, 19.23171, 9.360516, 32.24822, 6.398159, 13.46195, 7.53088, 0.5851686, 6.598461, 9.070963, 8.848, 6.786048,
   3.801906, 4.139167 };
   TGraph *graph = new TGraph(35,Graph0_fx13,Graph0_fy13);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph013 = new TH1F("Graph_Graph013","Graph",100,0,3520);
   Graph_Graph013->SetMinimum(1);
   Graph_Graph013->SetMaximum(2000);
   Graph_Graph013->SetDirectory(nullptr);
   Graph_Graph013->SetStats(0);
   Graph_Graph013->SetLineWidth(2);
   Graph_Graph013->SetMarkerStyle(20);
   Graph_Graph013->SetMarkerSize(0.9);
   Graph_Graph013->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph013->GetXaxis()->SetRange(1,101);
   Graph_Graph013->GetXaxis()->SetLabelFont(43);
   Graph_Graph013->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetXaxis()->SetLabelSize(16);
   Graph_Graph013->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph013->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetXaxis()->SetTitleFont(42);
   Graph_Graph013->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph013->GetYaxis()->SetLabelFont(43);
   Graph_Graph013->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetYaxis()->SetLabelSize(16);
   Graph_Graph013->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph013->GetYaxis()->SetTickLength(0.02);
   Graph_Graph013->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetYaxis()->SetTitleFont(42);
   Graph_Graph013->GetZaxis()->SetLabelFont(42);
   Graph_Graph013->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph013->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph013->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph013->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph013->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph013);
   
   graph->Draw("ap");
   
   TLegend *leg = new TLegend(0.27,0.7,0.5,0.93,NULL,"brNDC");
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
   leg->Draw();
   
   Double_t Graph1_fx14[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy14[35] = { 100, 100, 100, 100, 0, 100, 48.46771, 0, 100, 0, 100, 171.414, 54.34362, 20.72161, 18.77546, 100, 0,
   52.04154, 26.36508, 17.47038, 41.53788, 5.960464e-06, 8.855724, 9.570312, 5.407357, 1.864707, 3.383791, 1.084882, 0.5438089, 1.643491, 2.670056, 1.061928, 4.313159,
   3.801906, 0 };
   graph = new TGraph(35,Graph1_fx14,Graph1_fy14);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph114 = new TH1F("Graph_Graph114","Graph",100,0,3520);
   Graph_Graph114->SetMinimum(92.8586);
   Graph_Graph114->SetMaximum(178.5555);
   Graph_Graph114->SetDirectory(nullptr);
   Graph_Graph114->SetStats(0);
   Graph_Graph114->SetLineWidth(2);
   Graph_Graph114->SetMarkerStyle(20);
   Graph_Graph114->SetMarkerSize(0.9);
   Graph_Graph114->GetXaxis()->SetLabelFont(42);
   Graph_Graph114->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetXaxis()->SetTitleFont(42);
   Graph_Graph114->GetYaxis()->SetLabelFont(42);
   Graph_Graph114->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetYaxis()->SetTickLength(0.02);
   Graph_Graph114->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetYaxis()->SetTitleFont(42);
   Graph_Graph114->GetZaxis()->SetLabelFont(42);
   Graph_Graph114->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph114->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph114->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph114->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph114->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph114);
   
   graph->Draw("p");
   
   Double_t Graph2_fx15[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy15[35] = { 100, 100, 100, 27.44938, 23.04946, 100, 30.59461, 28.24067, 23.59124, 136.2849, 100, 30.59461, 237.6894, 84.20563, 158.6791, 31.63878, 98.50137,
   20.16787, 32.88016, 28.95981, 30.64523, 220.202, 19.94328, 29.77039, 18.19582, 54.71046, 30.98527, 37.96197, 41.70592, 33.96362, 36.93254, 34.44015, 16.52872,
   13.64875, 15.28776 };
   graph = new TGraph(35,Graph2_fx15,Graph2_fy15);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph215 = new TH1F("Graph_Graph215","Graph",100,0,3520);
   Graph_Graph215->SetMinimum(12.28387);
   Graph_Graph215->SetMaximum(260.0934);
   Graph_Graph215->SetDirectory(nullptr);
   Graph_Graph215->SetStats(0);
   Graph_Graph215->SetLineWidth(2);
   Graph_Graph215->SetMarkerStyle(20);
   Graph_Graph215->SetMarkerSize(0.9);
   Graph_Graph215->GetXaxis()->SetLabelFont(42);
   Graph_Graph215->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetXaxis()->SetTitleFont(42);
   Graph_Graph215->GetYaxis()->SetLabelFont(42);
   Graph_Graph215->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetYaxis()->SetTickLength(0.02);
   Graph_Graph215->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetYaxis()->SetTitleFont(42);
   Graph_Graph215->GetZaxis()->SetLabelFont(42);
   Graph_Graph215->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph215->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph215->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph215->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph215->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph215);
   
   graph->Draw("p");
   
   Double_t Graph3_fx16[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy16[35] = { 100, 100, 100, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx16,Graph3_fy16);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph316 = new TH1F("Graph_Graph316","Graph",100,0,3520);
   Graph_Graph316->SetMinimum(99.9);
   Graph_Graph316->SetMaximum(101.1);
   Graph_Graph316->SetDirectory(nullptr);
   Graph_Graph316->SetStats(0);
   Graph_Graph316->SetLineWidth(2);
   Graph_Graph316->SetMarkerStyle(20);
   Graph_Graph316->SetMarkerSize(0.9);
   Graph_Graph316->GetXaxis()->SetLabelFont(42);
   Graph_Graph316->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetXaxis()->SetTitleFont(42);
   Graph_Graph316->GetYaxis()->SetLabelFont(42);
   Graph_Graph316->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetYaxis()->SetTickLength(0.02);
   Graph_Graph316->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetYaxis()->SetTitleFont(42);
   Graph_Graph316->GetZaxis()->SetLabelFont(42);
   Graph_Graph316->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph316->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph316->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph316->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph316->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph316);
   
   graph->Draw("p");
   
   Double_t Graph4_fx17[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy17[35] = { 100, 100, 100, 2.803653, 0.6861925, 100, 2.298272, 2.825499, 2.803665, 2.803659, 100, 0.9149969, 1.081431, 2.694821, 2.20083, 0.7658124, 2.436876,
   0.930953, 1.958025, 2.404273, 2.287471, 2.227938, 2.029544, 2.155697, 1.801032, 1.940709, 1.837093, 1.95666, 1.964509, 1.943159, 1.936275, 1.994902, 2.439207,
   1.915306, 1.656818 };
   graph = new TGraph(35,Graph4_fx17,Graph4_fy17);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph417 = new TH1F("Graph_Graph417","Graph",100,0,3520);
   Graph_Graph417->SetMinimum(0.6175733);
   Graph_Graph417->SetMaximum(109.9314);
   Graph_Graph417->SetDirectory(nullptr);
   Graph_Graph417->SetStats(0);
   Graph_Graph417->SetLineWidth(2);
   Graph_Graph417->SetMarkerStyle(20);
   Graph_Graph417->SetMarkerSize(0.9);
   Graph_Graph417->GetXaxis()->SetLabelFont(42);
   Graph_Graph417->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetXaxis()->SetTitleFont(42);
   Graph_Graph417->GetYaxis()->SetLabelFont(42);
   Graph_Graph417->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetYaxis()->SetTickLength(0.02);
   Graph_Graph417->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetYaxis()->SetTitleFont(42);
   Graph_Graph417->GetZaxis()->SetLabelFont(42);
   Graph_Graph417->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph417->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph417->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph417->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph417->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph417);
   
   graph->Draw("p");
   
   Double_t Graph5_fx18[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy18[35] = { 223.6068, 223.6068, 223.6068, 103.7368, 23.05967, 223.6068, 57.36228, 28.38166, 102.7833, 136.3138, 223.6068, 190.9044, 253.3664, 95.24109, 161.2178, 154.6037, 99.46249,
   76.3168, 42.88648, 33.94218, 72.31008, 221.8056, 29.15709, 32.71289, 37.46357, 55.14902, 34.002, 38.76636, 41.7598, 34.69213, 38.17294, 35.6303, 18.54191,
   14.79411, 15.92461 };
   graph = new TGraph(35,Graph5_fx18,Graph5_fy18);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph518 = new TH1F("Graph_Graph518","Graph",100,0,3520);
   Graph_Graph518->SetMinimum(13.3147);
   Graph_Graph518->SetMaximum(277.2236);
   Graph_Graph518->SetDirectory(nullptr);
   Graph_Graph518->SetStats(0);
   Graph_Graph518->SetLineWidth(2);
   Graph_Graph518->SetMarkerStyle(20);
   Graph_Graph518->SetMarkerSize(0.9);
   Graph_Graph518->GetXaxis()->SetLabelFont(42);
   Graph_Graph518->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetXaxis()->SetTitleFont(42);
   Graph_Graph518->GetYaxis()->SetLabelFont(42);
   Graph_Graph518->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetYaxis()->SetTickLength(0.02);
   Graph_Graph518->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetYaxis()->SetTitleFont(42);
   Graph_Graph518->GetZaxis()->SetLabelFont(42);
   Graph_Graph518->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph518->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph518->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph518->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph518->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph518);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.16,0.96,"#scale[1.3]{#bf{CMS}}#it{Simulation Work in progress}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c1->Modified();
   c1->SetSelected(c1);
}
