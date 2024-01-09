
#include "flc.hpp"

#include "/Users/federico/G_Class/G_Func.hpp"
#include "/Users/federico/G_Class/G_Read.hpp"
#include "/Users/federico/G_Class/G_WF.hpp"
#include "/Users/federico/G_Class/G_Utility.hpp"

void Noise_PSD(){
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  gStyle->SetPalette(kSunset);
  
  double t;
  std::vector<Double_t> avg_wf, noise2, noise, y3, int_wf;
  
  //Read, subtract the baseline to the wfs and build the charge histogram
  CAEN_WF_Binary(WF_FILE, noise, N_WF, MEMORYDEPTH);
  //CompleteWFfloat_Binary(WF_FILE, noise, N_WF, MEMORYDEPTH);
  SubBaseline(noise, MEMORYDEPTH, PREPULSE_TICKS);
  int nnn =  NonSat_WF(noise, noise2, MEMORYDEPTH, SAT_LOW, SAT_UP);
  TH1D* hI = BuildRawChargeHisto(noise2, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  Avg_Sel_WF (noise2, y3, avg_wf, int_wf, MU0_LOW, MU0_UP);

  std::cout << "N wf for FFT : " << nnn << std::endl;

  TGraph* gNoise_spectral_density = build_avg_spectral_density(MEMORYDEPTH, TICK_LEN*MEMORYDEPTH, TICK_LEN, y3, RES);
  return;  
  std::ofstream OutFile ("FFT_short.dat", ios::binary);
  for(int i = 0; i < gNoise_spectral_density->GetN(); i++){
    t = gNoise_spectral_density->GetPointY(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
    t = gNoise_spectral_density->GetPointX(i);
    OutFile.write(reinterpret_cast<char*> (&t), sizeof(t));
}
  
  std::cout << "Vector saved in ---> " << std::endl;
  OutFile.close();
}
 
