#ifdef __CLING__
#pragma cling optimize(0)
#endif
void summary_Gluino_1400()
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
   
   Double_t Graph0_fx22[35] = { 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320,
   340, 360, 380, 410, 440, 480, 530, 590, 660, 760, 880, 1030, 1210, 1440, 1730, 2000,
   2500, 3200 };
   Double_t Graph0_fy22[35] = { 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 31.22379, 36.13875, 57.01214, 100.4947, 17.02057, 28.73288,
   24.36884, 31.34621, 86.86649, 0, 16.54417, 4.150498, 31.41952, 11.99589, 1.632869, 3.546971, 2.906984, 1.773393, 1.861167, 2.828103, 3.109241, 0.8568287,
   5.519354, 0 };
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
   Double_t Graph1_fy23[35] = { 100, 100, 100, 73.9323, 40.52021, 100, 100, 47.99643, 84.73432, 69.50414, 0, 0, 0, 57.01214, 100, 22.02649, 0.4982352,
   24.36884, 20.05627, 86.64072, 16.30896, 20.10764, 4.150498, 11.31714, 4.722619, 4.621661, 5.318594, 3.992152, 1.364917, 1.572597, 3.195977, 3.091872, 3.47954,
   5.519354, 1.408398 };
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
   Graph_Graph123->SetMinimum(34.57224);
   Graph_Graph123->SetMaximum(105.948);
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
   Double_t Graph2_fy24[35] = { 100, 100, 100, 24.34452, 43.26846, 100, 100, 51.63009, 20.84246, 12.35838, 27.03399, 32.44827, 12.48886, 34.66217, 34.71781, 6.717807, 47.71345,
   10.21199, 21.51796, 21.97553, 11.15347, 21.52226, 10.37251, 23.90988, 22.34188, 4.871833, 4.725826, 1.971, 2.341354, 1.223862, 3.375459, 2.678001, 5.402124,
   3.680336, 4.788125 };
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
   Graph_Graph224->SetMinimum(1.101476);
   Graph_Graph224->SetMaximum(109.8776);
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
   Double_t Graph3_fy25[35] = { 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 5.048454, 1.632869, 1.845622, 1.348686, 1.308179, 1.439571, 2.401006, 1.851726, 2.03892,
   1.184464, 5.203247 };
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
   Double_t Graph4_fy26[35] = { 100, 100, 100, 4.396004, 3.647757, 100, 100, 3.368664, 2.946782, 2.963948, 2.507889, 3.832, 2.794743, 1.695478, 1.603138, 2.510232, 1.77393,
   1.165712, 3.297019, 2.259111, 3.249466, 2.891445, 2.898729, 2.980924, 3.401387, 2.765501, 2.980638, 3.224528, 3.339577, 3.414237, 3.344262, 3.51088, 3.447092,
   3.488696, 3.575957 };
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
   Graph_Graph426->SetMinimum(1.049141);
   Graph_Graph426->SetMaximum(109.8834);
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
   Double_t Graph5_fy27[35] = { 100, 100, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   50.96358, 0, 0, 0, 0, 0, 0, 0, 7.201111, 0.6596684, 0.3250062, 0.6934941, 0.4032433, 0.2904832, 0.4435301, 0,
   0, 1.32134 };
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
   Double_t Graph6_fy28[35] = { 223.6068, 223.6068, 223.6068, 77.96131, 59.39153, 223.6068, 223.6068, 70.57388, 87.30978, 70.65649, 27.15007, 45.19402, 38.33786, 87.77876, 145.9695, 28.74537, 55.72743,
   35.96281, 43.11301, 124.6613, 20.02351, 33.90568, 12.2656, 41.18048, 26.50322, 7.620697, 8.688636, 6.363113, 4.832594, 4.600459, 6.823128, 6.492604, 7.619969,
   9.383249, 8.048041 };
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
   Graph_Graph628->SetMinimum(4.140413);
   Graph_Graph628->SetMaximum(245.5074);
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
