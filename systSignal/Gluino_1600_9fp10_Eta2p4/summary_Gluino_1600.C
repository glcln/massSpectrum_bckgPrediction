#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1600()
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
   
   Double_t Graph0_fx29[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy29[35] = { 100, 0, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 13.79598, 19.36631,
   0.03204346, 19.07293, 16.03257, 31.39418, 7.803178, 0, 12.54374, 3.00231, 5.944431, 3.825843, 3.213549, 3.242135, 0.116539, 2.008843, 2.572143, 1.667833,
   2.842033, 1.001132 };
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
   Double_t Graph1_fy30[35] = { 100, 100, 100, 100, 0, 68.16498, 320.2316, 36.78538, 0, 29.03384, 0, 8.964276, 0, 71.66038, 100, 15.24549, 29.03975,
   10.71893, 51.744, 8.08329, 43.3125, 24.0157, 6.62579, 27.04868, 12.62243, 2.833223, 7.849193, 5.019498, 2.944577, 0.3974438, 2.129245, 3.792763, 1.942098,
   3.878689, 1.001137 };
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
   Graph_Graph130->SetMinimum(77.97684);
   Graph_Graph130->SetMaximum(342.2547);
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
   Double_t Graph2_fy31[35] = { 100, 54.12139, 100, 100, 30.75011, 31.88645, 43.9416, 24.15559, 24.34925, 36.49092, 47.97527, 16.53525, 15.38796, 36.98486, 100, 26.96863, 43.31357,
   18.203, 23.39746, 18.24149, 31.72739, 8.316982, 11.65035, 5.920267, 15.66396, 5.83396, 2.654892, 6.701517, 2.747375, 1.661551, 1.754546, 2.771199, 4.918361,
   5.210644, 6.081176 };
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
   Graph_Graph231->SetMinimum(1.495396);
   Graph_Graph231->SetMaximum(109.8338);
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
   Double_t Graph3_fy32[35] = { 100, 0, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 9.263444, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 2.931708, 1.365638, 1.335174, 1.399118, 1.558816, 2.026236, 3.283024,
   4.302454, 2.851826 };
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
   Double_t Graph4_fy33[35] = { 100, 2.507889, 100, 100, 4.904199, 2.653927, 1.234829, 2.92716, 1.300681, 3.396153, 3.118026, 3.190589, 2.89042, 3.495669, 100, 1.772869, 2.800655,
   2.827299, 2.833545, 2.633405, 1.308167, 2.047372, 2.662873, 2.766407, 2.395475, 3.157651, 2.99201, 3.12233, 3.30286, 3.416574, 3.436971, 3.420269, 3.451252,
   3.415215, 3.520358 };
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
   Graph_Graph433->SetMinimum(1.111346);
   Graph_Graph433->SetMaximum(109.8765);
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
   Double_t Graph5_fy34[35] = { 100, 0, 100, 100, 0, 0, 0, 0, 28.58557, 0, 0, 0, 13.19151, 0, 100, 0, 0,
   0, 18.04231, 0, 0, 0, 12.25576, 0, 0, 7.962781, 1.013613, 0.8849204, 0.1628458, 0.2885282, 0.2795756, 0.5437374, 0.5302131,
   0, 0.3386676 };
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
   Double_t Graph6_fy35[35] = { 223.6068, 113.734, 223.6068, 223.6068, 31.13873, 75.30109, 323.2346, 44.1047, 24.38396, 46.75559, 48.07648, 19.07754, 15.65707, 80.71747, 223.6068, 35.19965, 55.69802,
   21.31289, 59.97239, 25.73071, 62.20857, 26.66472, 13.66465, 30.52342, 20.48017, 9.347153, 10.04203, 9.594087, 6.278664, 4.069742, 5.088246, 6.669662, 7.309793,
   8.969225, 7.714342 };
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
   Graph_Graph635->SetMinimum(3.662768);
   Graph_Graph635->SetMaximum(355.1511);
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
