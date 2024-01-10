#include "flc.hpp"


void AverageWF() {
  double y0, y1;
  vector<double> x , y , y2 , avg_wf, avg_wf2;
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i*TICK_LEN);
  
  CAEN_WF_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  //CompleteWFfloat_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  //SubBaseline_Invert(y, MEMORYDEPTH, PREPULSE_TICKS);
  SubBaseline(y, MEMORYDEPTH, PREPULSE_TICKS);
  //SubBaseline_Range(y, MEMORYDEPTH, 2200, 3000);
  //int nnn =  NonSat_WF(y, y2, MEMORYDEPTH, SAT_LOW, SAT_UP);
  //avgWF(y2 , avg_wf, MEMORYDEPTH);
  avgWF(y , avg_wf, MEMORYDEPTH);
  MovingAverageWF(avg_wf, avg_wf, 10);
  std::cout << "size " << y.size() << " #WF " << y.size()/MEMORYDEPTH << std::endl;  
  std::cout << "avg size " << avg_wf.size() << std::endl;
  Print_RiseFallTime(avg_wf, TICK_LEN);
  TCanvas* cTime = new TCanvas("wavedec","wavedec");
 
  DisplayWFs (y, MEMORYDEPTH, TICK_LEN, 10);
 
  y0 = *min_element(std::begin(y), std::end(y));
  y1 = *max_element(std::begin(y), std::end(y));
  //y0 = *min_element(std::begin(y2), std::end(y2));
  //y1 = *max_element(std::begin(y2), std::end(y2));

  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), MEMORYDEPTH/2, 0., MEMORYDEPTH, 120, y0, y1);
  //for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y[i*MEMORYDEPTH+j]);
  //for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y2[i*MEMORYDEPTH+j]);
  //y.erase(y.begin(), y.end());
  //h2->Draw("COLZ");
  //return;
  TGraph *g1 = new TGraph(avg_wf.size(), &x[0], &avg_wf[0]);
  g1->GetXaxis()->SetTitle("Time [#mus]");
  g1->GetYaxis()->SetTitle("ADC counts");
  g1->SetLineColor(2);
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,800);
  c2->cd();
  g1->Draw();
  c2->Modified();
  c2->Update();
  gPad->Modified();
  gPad->Update();
//  VecDouble_in_Binary("Template50L.dat",avg_wf);
}
