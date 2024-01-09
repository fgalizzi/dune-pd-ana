#include "flc.hpp"

#include "/Users/federico/G_Class/G_Func.hpp"
#include "/Users/federico/G_Class/G_Read.hpp"
#include "/Users/federico/G_Class/G_WF.hpp"
#include "/Users/federico/G_Class/G_Utility.hpp"



void SPE() {
  double y0, y1;
  vector<double> x, y, y2, y3, avg_wf, int_wf;
    
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i*TICK_LEN);
  
  //CompleteWFfloat_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  CAEN_WF_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  SubBaseline(y, MEMORYDEPTH, PREPULSE_TICKS);
  int nnn =  NonSat_WF(y, y2, MEMORYDEPTH, SAT_LOW, SAT_UP);
  TH1D* hI = BuildRawChargeHisto(y2, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  Avg_Sel_WF (y2, y3, avg_wf, int_wf, MU0_LOW, MU0_UP);
 
  Print_RiseFallTime(avg_wf, TICK_LEN);
  MovingAverageWF(avg_wf, avg_wf, 100);

  TCanvas* cTime = new TCanvas("wavedec","wavedec");
 
  TGraph *g1 = new TGraph(avg_wf.size(), &x[0], &avg_wf[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  g1->Draw();
  c2->Modified();
  c2->Update();
  gPad->Modified();
  gPad->Update();
}
