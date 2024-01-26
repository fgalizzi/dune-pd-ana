#include "flc.hpp"


void Persistence() {
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);

  double y0, y1;
  vector<double> x , y , y2 , avg_wf, avg_wf2;
  
  for (size_t i = 0; i < MEMORYDEPTH; i++) x.push_back( (double) i);
  
  CAEN_WF_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  //CompleteWFfloat_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  SubBaseline(y, MEMORYDEPTH, PREPULSE_TICKS);

  std::cout << "size " << y.size() << " #WF " << y.size()/MEMORYDEPTH << std::endl;  
  std::cout << "avg size " << avg_wf.size() << std::endl;
  
  TCanvas* cTime = new TCanvas("wavedec","wavedec");
  
  //MovingAverageWF(y, y2, 100);
  DisplayWFs (y, MEMORYDEPTH, TICK_LEN, 15);


  //return;
  y0 = *min_element(std::begin(y), std::end(y));
  y1 = *max_element(std::begin(y), std::end(y));
  std::cout << "Max - Min " << y1-y0 << "\n";
  //y0 = *min_element(std::begin(y2), std::end(y2));
  //y1 = *max_element(std::begin(y2), std::end(y2));

  //TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), MEMORYDEPTH/2, 0., MEMORYDEPTH, 200, y0, y1);
  TH2D* h2 = new TH2D("h2", Form("%s;%s;%s", "DecoWF", "Ticks", "ADC Counts"), MEMORYDEPTH/2, 0., MEMORYDEPTH, 200, -30, 300);
  for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y[i*MEMORYDEPTH+j]);
  //for (int i=0; i<N_WF; i++) for (int j=0; j<MEMORYDEPTH; j=j+2) h2->Fill(j, y2[i*MEMORYDEPTH+j]);
  y.erase(y.begin(), y.end());
  h2->Draw("COLZ");

  return;
}
