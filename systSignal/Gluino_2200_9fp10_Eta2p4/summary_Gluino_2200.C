#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_2200()
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
   
   Double_t Graph0_fx50[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy50[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   11.86163, 24.73135, 0, 4.908586, 3.11197, 0, 19.61773, 7.986128, 3.614593, 2.383673, 2.875686, 5.551827, 2.779913, 2.055252, 0.5236864, 1.900566,
   2.252376, 1.521707 };
   TGraph *graph = new TGraph(35,Graph0_fx50,Graph0_fy50);
   graph->SetName("Graph0");
   graph->SetTitle("Graph");
   graph->SetFillColor(30);
   graph->SetFillStyle(1000);
   graph->SetLineColor(30);
   graph->SetMarkerColor(30);
   graph->SetMarkerStyle(21);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph050 = new TH1F("Graph_Graph050","Graph",100,0,3520);
   Graph_Graph050->SetMinimum(1);
   Graph_Graph050->SetMaximum(2000);
   Graph_Graph050->SetDirectory(nullptr);
   Graph_Graph050->SetStats(0);
   Graph_Graph050->SetLineWidth(2);
   Graph_Graph050->SetMarkerStyle(20);
   Graph_Graph050->SetMarkerSize(0.9);
   Graph_Graph050->GetXaxis()->SetTitle("Mass bin");
   Graph_Graph050->GetXaxis()->SetRange(1,101);
   Graph_Graph050->GetXaxis()->SetLabelFont(43);
   Graph_Graph050->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetXaxis()->SetLabelSize(22);
   Graph_Graph050->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph050->GetXaxis()->SetTitleOffset(1);
   Graph_Graph050->GetXaxis()->SetTitleFont(42);
   Graph_Graph050->GetYaxis()->SetTitle("Systematic Uncertainty [%]");
   Graph_Graph050->GetYaxis()->SetLabelFont(43);
   Graph_Graph050->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetYaxis()->SetLabelSize(22);
   Graph_Graph050->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph050->GetYaxis()->SetTickLength(0.02);
   Graph_Graph050->GetYaxis()->SetTitleOffset(1);
   Graph_Graph050->GetYaxis()->SetTitleFont(42);
   Graph_Graph050->GetZaxis()->SetLabelFont(42);
   Graph_Graph050->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph050->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph050->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph050->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph050->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph050);
   
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
   
   Double_t Graph1_fx51[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph1_fy51[35] = { 100, 100, 100, 100, 100, 100, 0, 30.09958, 412.6682, 5.960464e-06, 56.58397, 71.47784, 0, 37.13564, 41.2521, 15.36809, 31.01386,
   0.9622276, 24.73135, 0, 0, 17.15922, 1.823699, 17.58004, 4.903936, 3.614599, 3.410959, 0.1382351, 5.081391, 2.547789, 2.12177, 0.5658865, 1.560289,
   2.591372, 1.164913 };
   graph = new TGraph(35,Graph1_fx51,Graph1_fy51);
   graph->SetName("Graph1");
   graph->SetTitle("Graph");
   graph->SetFillColor(38);
   graph->SetFillStyle(1000);
   graph->SetLineColor(38);
   graph->SetMarkerColor(38);
   graph->SetMarkerStyle(22);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph151 = new TH1F("Graph_Graph151","Graph",100,0,3520);
   Graph_Graph151->SetMinimum(68.73318);
   Graph_Graph151->SetMaximum(443.935);
   Graph_Graph151->SetDirectory(nullptr);
   Graph_Graph151->SetStats(0);
   Graph_Graph151->SetLineWidth(2);
   Graph_Graph151->SetMarkerStyle(20);
   Graph_Graph151->SetMarkerSize(0.9);
   Graph_Graph151->GetXaxis()->SetLabelFont(42);
   Graph_Graph151->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetXaxis()->SetTitleFont(42);
   Graph_Graph151->GetYaxis()->SetLabelFont(42);
   Graph_Graph151->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetYaxis()->SetTickLength(0.02);
   Graph_Graph151->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetYaxis()->SetTitleFont(42);
   Graph_Graph151->GetZaxis()->SetLabelFont(42);
   Graph_Graph151->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph151->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph151->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph151->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph151->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph151);
   
   graph->Draw("p");
   
   Double_t Graph2_fx52[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph2_fy52[35] = { 100, 100, 100, 100, 100, 100, 91.52737, 27.84832, 26.5033, 22.02033, 55.43044, 30.00287, 24.57522, 29.71905, 27.03399, 15.08158, 34.31872,
   26.77629, 33.68576, 24.32543, 15.9551, 2.497375, 1.48201, 23.34278, 10.44331, 18.61055, 6.99054, 10.67054, 3.199589, 2.680469, 1.787376, 4.171968, 2.326119,
   1.321363, 5.11939 };
   graph = new TGraph(35,Graph2_fx52,Graph2_fy52);
   graph->SetName("Graph2");
   graph->SetTitle("Graph");
   graph->SetFillColor(46);
   graph->SetFillStyle(1000);
   graph->SetLineColor(46);
   graph->SetMarkerColor(46);
   graph->SetMarkerStyle(23);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph252 = new TH1F("Graph_Graph252","Graph",100,0,3520);
   Graph_Graph252->SetMinimum(1.189227);
   Graph_Graph252->SetMaximum(109.8679);
   Graph_Graph252->SetDirectory(nullptr);
   Graph_Graph252->SetStats(0);
   Graph_Graph252->SetLineWidth(2);
   Graph_Graph252->SetMarkerStyle(20);
   Graph_Graph252->SetMarkerSize(0.9);
   Graph_Graph252->GetXaxis()->SetLabelFont(42);
   Graph_Graph252->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetXaxis()->SetTitleFont(42);
   Graph_Graph252->GetYaxis()->SetLabelFont(42);
   Graph_Graph252->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetYaxis()->SetTickLength(0.02);
   Graph_Graph252->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetYaxis()->SetTitleFont(42);
   Graph_Graph252->GetZaxis()->SetLabelFont(42);
   Graph_Graph252->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph252->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph252->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph252->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph252->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph252);
   
   graph->Draw("p");
   
   Double_t Graph3_fx53[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph3_fy53[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 6.787884, 0, 0, 0, 2.274776, 0, 1.754886, 1.228577, 1.608866, 0.9610534, 0.9507835, 1.233935,
   1.900721, 1.558769 };
   graph = new TGraph(35,Graph3_fx53,Graph3_fy53);
   graph->SetName("Graph3");
   graph->SetTitle("Graph");
   graph->SetFillColor(43);
   graph->SetFillStyle(1000);
   graph->SetLineColor(43);
   graph->SetMarkerColor(43);
   graph->SetMarkerStyle(43);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph353 = new TH1F("Graph_Graph353","Graph",100,0,3520);
   Graph_Graph353->SetMinimum(99.9);
   Graph_Graph353->SetMaximum(101.1);
   Graph_Graph353->SetDirectory(nullptr);
   Graph_Graph353->SetStats(0);
   Graph_Graph353->SetLineWidth(2);
   Graph_Graph353->SetMarkerStyle(20);
   Graph_Graph353->SetMarkerSize(0.9);
   Graph_Graph353->GetXaxis()->SetLabelFont(42);
   Graph_Graph353->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetXaxis()->SetTitleFont(42);
   Graph_Graph353->GetYaxis()->SetLabelFont(42);
   Graph_Graph353->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetYaxis()->SetTickLength(0.02);
   Graph_Graph353->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetYaxis()->SetTitleFont(42);
   Graph_Graph353->GetZaxis()->SetLabelFont(42);
   Graph_Graph353->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph353->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph353->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph353->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph353->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph353);
   
   graph->Draw("p");
   
   Double_t Graph4_fx54[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph4_fy54[35] = { 100, 100, 100, 100, 100, 100, 3.647757, 2.761465, 3.111076, 3.477997, 1.549745, 2.424598, 2.881527, 3.230691, 5.113363, 1.700366, 3.359795,
   3.075111, 1.182175, 3.357863, 4.805887, 3.251612, 3.195584, 3.147769, 2.647221, 2.694595, 2.881873, 2.99319, 3.106201, 3.283167, 3.38701, 3.398645, 3.47935,
   3.55798, 3.497601 };
   graph = new TGraph(35,Graph4_fx54,Graph4_fy54);
   graph->SetName("Graph4");
   graph->SetTitle("Graph");
   graph->SetFillColor(40);
   graph->SetFillStyle(1000);
   graph->SetLineColor(40);
   graph->SetMarkerColor(40);
   graph->SetMarkerStyle(39);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph454 = new TH1F("Graph_Graph454","Graph",100,0,3520);
   Graph_Graph454->SetMinimum(1.063957);
   Graph_Graph454->SetMaximum(109.8818);
   Graph_Graph454->SetDirectory(nullptr);
   Graph_Graph454->SetStats(0);
   Graph_Graph454->SetLineWidth(2);
   Graph_Graph454->SetMarkerStyle(20);
   Graph_Graph454->SetMarkerSize(0.9);
   Graph_Graph454->GetXaxis()->SetLabelFont(42);
   Graph_Graph454->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetXaxis()->SetTitleFont(42);
   Graph_Graph454->GetYaxis()->SetLabelFont(42);
   Graph_Graph454->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetYaxis()->SetTickLength(0.02);
   Graph_Graph454->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetYaxis()->SetTitleFont(42);
   Graph_Graph454->GetZaxis()->SetLabelFont(42);
   Graph_Graph454->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph454->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph454->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph454->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph454->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph454);
   
   graph->Draw("p");
   
   Double_t Graph5_fx55[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph5_fy55[35] = { 100, 100, 100, 100, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 3.077102, 0, 0, 0, 0, 0, 0, 1.362079, 0.918138, 0.472337, 0.4318118, 0.3732145, 0.270617, 0.1830161,
   0.2362311, 0 };
   graph = new TGraph(35,Graph5_fx55,Graph5_fy55);
   graph->SetName("Graph5");
   graph->SetTitle("Graph");
   graph->SetFillColor(41);
   graph->SetFillStyle(1000);
   graph->SetLineColor(41);
   graph->SetMarkerColor(41);
   graph->SetMarkerStyle(42);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph555 = new TH1F("Graph_Graph555","Graph",100,0,3520);
   Graph_Graph555->SetMinimum(99.9);
   Graph_Graph555->SetMaximum(101.1);
   Graph_Graph555->SetDirectory(nullptr);
   Graph_Graph555->SetStats(0);
   Graph_Graph555->SetLineWidth(2);
   Graph_Graph555->SetMarkerStyle(20);
   Graph_Graph555->SetMarkerSize(0.9);
   Graph_Graph555->GetXaxis()->SetLabelFont(42);
   Graph_Graph555->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetXaxis()->SetTitleFont(42);
   Graph_Graph555->GetYaxis()->SetLabelFont(42);
   Graph_Graph555->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetYaxis()->SetTickLength(0.02);
   Graph_Graph555->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetYaxis()->SetTitleFont(42);
   Graph_Graph555->GetZaxis()->SetLabelFont(42);
   Graph_Graph555->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph555->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph555->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph555->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph555->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph555);
   
   graph->Draw("p");
   
   Double_t Graph6_fx56[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph6_fy56[35] = { 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 223.6068, 91.60003, 41.09914, 413.5301, 22.2933, 79.2255, 77.55728, 24.74357, 47.673, 49.58547, 21.59916, 46.37804,
   29.46269, 48.57373, 24.5561, 17.37113, 19.15747, 3.966609, 35.33705, 14.27926, 19.61934, 8.630723, 11.58395, 8.833944, 5.896752, 4.93114, 5.518566, 5.008588,
   5.459521, 6.674109 };
   graph = new TGraph(35,Graph6_fx56,Graph6_fy56);
   graph->SetName("Graph6");
   graph->SetTitle("Graph");
   graph->SetFillColor(28);
   graph->SetFillStyle(1000);
   graph->SetLineColor(28);
   graph->SetMarkerColor(28);
   graph->SetMarkerStyle(34);
   graph->SetMarkerSize(0.9);
   
   TH1F *Graph_Graph656 = new TH1F("Graph_Graph656","Graph",100,0,3520);
   Graph_Graph656->SetMinimum(3.569948);
   Graph_Graph656->SetMaximum(454.4864);
   Graph_Graph656->SetDirectory(nullptr);
   Graph_Graph656->SetStats(0);
   Graph_Graph656->SetLineWidth(2);
   Graph_Graph656->SetMarkerStyle(20);
   Graph_Graph656->SetMarkerSize(0.9);
   Graph_Graph656->GetXaxis()->SetLabelFont(42);
   Graph_Graph656->GetXaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetXaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetXaxis()->SetTitleFont(42);
   Graph_Graph656->GetYaxis()->SetLabelFont(42);
   Graph_Graph656->GetYaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetYaxis()->SetTickLength(0.02);
   Graph_Graph656->GetYaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetYaxis()->SetTitleFont(42);
   Graph_Graph656->GetZaxis()->SetLabelFont(42);
   Graph_Graph656->GetZaxis()->SetLabelOffset(0.015);
   Graph_Graph656->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph656->GetZaxis()->SetTitleSize(0.065);
   Graph_Graph656->GetZaxis()->SetTitleOffset(1.1);
   Graph_Graph656->GetZaxis()->SetTitleFont(42);
   graph->SetHistogram(Graph_Graph656);
   
   graph->Draw("p");
   TLatex *   tex = new TLatex(0.1,0.96,"#scale[1.3]{#it{Private work (CMS simulation)}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.73,0.96,"#scale[1.3]{#bf{m_{#tilde{g}}=2200 GeV}}");
   tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   c2->Modified();
   c2->SetSelected(c2);
}
