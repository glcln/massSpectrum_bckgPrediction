#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx25[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy25[35] = { 100, 100, 100, 100, 0, 0, 0, 39.96457, 221.241, 100, 101.8407, 0.896436, 9.280962, 100, 0, 29.41645, 100,
   0, 100, 22.52771, 0, 29.40916, 9.629226, 7.497609, 6.623155, 3.935075, 9.711123, 7.028615, 10.15071, 1.802343, 5.097365, 5.447745, 6.63718,
   6.344771, 5.20153 };
   TGraph *graph = new TGraph(35,Graph0_fx25,Graph0_fy25);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph025 = new TH1F("Graph_Graph025","Graph",100,0,3520);
   Graph_Graph025->SetMinimum(1);
   Graph_Graph025->SetMaximum(2000);
   Graph_Graph025->SetDirectory(nullptr);
   Graph_Graph025->SetStats(0);
   Graph_Graph025->SetLineWidth(2);
   Graph_Graph025->SetMarkerStyle(20);
   Graph_Graph025->SetMarkerSize(0.9);
   Graph_Graph025->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph025->GetXaxis()->SetRange(1,101);
   Graph_Graph025->GetXaxis()->SetLabelFont(43);
   Graph_Graph025->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetXaxis()->SetLabelSize(16);
   Graph_Graph025->GetXaxis()->SetTitleSize(0.04);
   Graph_Graph025->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetXaxis()->SetTitleFont(42);
   Graph_Graph025->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph025->GetYaxis()->SetLabelFont(43);
   Graph_Graph025->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetYaxis()->SetLabelSize(16);
   Graph_Graph025->GetYaxis()->SetTitleSize(0.04);
   Graph_Graph025->GetYaxis()->SetTickLength(0.02);
   Graph_Graph025->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetYaxis()->SetTitleFont(42);
   Graph_Graph025->GetZaxis()->SetLabelFont(42);
   Graph_Graph025->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph025->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph025->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph025->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph025->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph025);
   
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
   
   Double_t Graph1_fx26[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy26[35] = { 100, 100, 100, 100, 0, 46.98205, 0, 39.96457, 121.52, 100, 101.8407, 53.53054, 27.11737, 6.771088, 0, 30.15056, 102.4956,
   0, 100, 0, 1.192093e-05, 0, 0, 0, 0, 4.689562, 2.919817, 2.764809, 1.102912, 0.4244924, 0.8594632, 1.259816, 1.770782,
   1.039815, 0 };
   graph = new TGraph(35,Graph1_fx26,Graph1_fy26);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph126 = new TH1F("Graph_Graph126","Graph",100,0,3520);
   Graph_Graph126->SetMinimum(97.848);
   Graph_Graph126->SetMaximum(123.672);
   Graph_Graph126->SetDirectory(nullptr);
   Graph_Graph126->SetStats(0);
   Graph_Graph126->SetLineWidth(2);
   Graph_Graph126->SetMarkerStyle(20);
   Graph_Graph126->SetMarkerSize(0.9);
   Graph_Graph126->GetXaxis()->SetLabelFont(42);
   Graph_Graph126->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetXaxis()->SetTitleFont(42);
   Graph_Graph126->GetYaxis()->SetLabelFont(42);
   Graph_Graph126->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetYaxis()->SetTickLength(0.02);
   Graph_Graph126->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetYaxis()->SetTitleFont(42);
   Graph_Graph126->GetZaxis()->SetLabelFont(42);
   Graph_Graph126->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph126->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph126->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph126->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph126->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph126);
   
   graph->Draw("p");
   
   Double_t Graph2_fx27[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy27[35] = { 100, 100, 100, 100, 19.08832, 19.72026, 36.23116, 22.52525, 22.36394, 100, 28.98337, 22.55288, 112.8161, 16.41644, 22.36394, 32.17925, 30.59461,
   24.09071, 100, 25.3762, 15.99126, 106.1284, 41.22542, 43.71197, 85.51138, 35.77089, 55.02653, 26.53872, 35.07348, 39.3963, 37.80286, 34.30271, 53.0834,
   32.02887, 32.37941 };
   graph = new TGraph(35,Graph2_fx27,Graph2_fy27);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph227 = new TH1F("Graph_Graph227","Graph",100,0,3520);
   Graph_Graph227->SetMinimum(6.30877);
   Graph_Graph227->SetMaximum(122.4986);
   Graph_Graph227->SetDirectory(nullptr);
   Graph_Graph227->SetStats(0);
   Graph_Graph227->SetLineWidth(2);
   Graph_Graph227->SetMarkerStyle(20);
   Graph_Graph227->SetMarkerSize(0.9);
   Graph_Graph227->GetXaxis()->SetLabelFont(42);
   Graph_Graph227->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetXaxis()->SetTitleFont(42);
   Graph_Graph227->GetYaxis()->SetLabelFont(42);
   Graph_Graph227->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetYaxis()->SetTickLength(0.02);
   Graph_Graph227->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetYaxis()->SetTitleFont(42);
   Graph_Graph227->GetZaxis()->SetLabelFont(42);
   Graph_Graph227->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph227->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph227->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph227->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph227->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph227);
   
   graph->Draw("p");
   
   Double_t Graph3_fx28[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy28[35] = { 100, 100, 100, 100, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0,
   0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0 };
   graph = new TGraph(35,Graph3_fx28,Graph3_fy28);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph328 = new TH1F("Graph_Graph328","Graph",100,0,3520);
   Graph_Graph328->SetMinimum(99.9);
   Graph_Graph328->SetMaximum(101.1);
   Graph_Graph328->SetDirectory(nullptr);
   Graph_Graph328->SetStats(0);
   Graph_Graph328->SetLineWidth(2);
   Graph_Graph328->SetMarkerStyle(20);
   Graph_Graph328->SetMarkerSize(0.9);
   Graph_Graph328->GetXaxis()->SetLabelFont(42);
   Graph_Graph328->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetXaxis()->SetTitleFont(42);
   Graph_Graph328->GetYaxis()->SetLabelFont(42);
   Graph_Graph328->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetYaxis()->SetTickLength(0.02);
   Graph_Graph328->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetYaxis()->SetTitleFont(42);
   Graph_Graph328->GetZaxis()->SetLabelFont(42);
   Graph_Graph328->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph328->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph328->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph328->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph328->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph328);
   
   graph->Draw("p");
   
   Double_t Graph4_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy29[35] = { 100, 100, 100, 100, 2.803653, 1.60706, 0.9150028, 1.691902, 1.127869, 100, 2.868772, 2.089661, 2.284694, 0.7658124, 0.6861925, 1.499599, 0.6861925,
   1.756173, 100, 2.130258, 1.532018, 2.243418, 1.934737, 1.63635, 1.882058, 2.054149, 1.87645, 2.133435, 1.913714, 1.973605, 1.987058, 1.941288, 1.847667,
   2.023137, 1.853448 };
   graph = new TGraph(35,Graph4_fx29,Graph4_fy29);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph429 = new TH1F("Graph_Graph429","Graph",100,0,3520);
   Graph_Graph429->SetMinimum(0.6175733);
   Graph_Graph429->SetMaximum(109.9314);
   Graph_Graph429->SetDirectory(nullptr);
   Graph_Graph429->SetStats(0);
   Graph_Graph429->SetLineWidth(2);
   Graph_Graph429->SetMarkerStyle(20);
   Graph_Graph429->SetMarkerSize(0.9);
   Graph_Graph429->GetXaxis()->SetLabelFont(42);
   Graph_Graph429->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetXaxis()->SetTitleFont(42);
   Graph_Graph429->GetYaxis()->SetLabelFont(42);
   Graph_Graph429->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetYaxis()->SetTickLength(0.02);
   Graph_Graph429->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetYaxis()->SetTitleFont(42);
   Graph_Graph429->GetZaxis()->SetLabelFont(42);
   Graph_Graph429->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph429->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph429->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph429->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph429->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph429);
   
   graph->Draw("p");
   
   Double_t Graph5_fx30[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy30[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 19.29312, 50.97827, 36.24271, 60.86529, 253.409, 223.6068, 146.9398, 58.13194, 116.4225, 101.5674, 22.37446, 53.02958, 146.4303,
   24.15463, 223.6068, 33.99981, 16.06448, 110.1507, 42.37925, 44.38049, 85.78813, 36.34904, 55.98457, 27.67491, 36.57957, 39.48914, 38.20637, 34.80962, 53.5579,
   32.73039, 32.84688 };
   graph = new TGraph(35,Graph5_fx30,Graph5_fy30);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph530 = new TH1F("Graph_Graph530","Graph",100,0,3520);
   Graph_Graph530->SetMinimum(14.45803);
   Graph_Graph530->SetMaximum(277.1435);
   Graph_Graph530->SetDirectory(nullptr);
   Graph_Graph530->SetStats(0);
   Graph_Graph530->SetLineWidth(2);
   Graph_Graph530->SetMarkerStyle(20);
   Graph_Graph530->SetMarkerSize(0.9);
   Graph_Graph530->GetXaxis()->SetLabelFont(42);
   Graph_Graph530->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetXaxis()->SetTitleFont(42);
   Graph_Graph530->GetYaxis()->SetLabelFont(42);
   Graph_Graph530->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetYaxis()->SetTickLength(0.02);
   Graph_Graph530->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetYaxis()->SetTitleFont(42);
   Graph_Graph530->GetZaxis()->SetLabelFont(42);
   Graph_Graph530->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph530->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph530->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph530->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph530->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph530);
   
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
