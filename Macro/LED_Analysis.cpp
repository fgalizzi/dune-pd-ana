#include "flc.hpp"


//*** MAIN ************************************
void LED_Analysis(){
//*********************************************
  gStyle->SetOptFit(1111); gStyle->SetOptTitle(0);
  gStyle->SetStatX(0.9); gStyle->SetStatY(0.9);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit"); 

  std::vector<double> y, sel_wf, int_wf;

  //CompleteWFfloat_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  CAEN_WF_Binary(WF_FILE, y, N_WF, MEMORYDEPTH);
  //SubBaseline_Invert(all_wf, MEMORYDEPTH, PREPULSE_TICKS);
  SubBaseline(y, MEMORYDEPTH, PREPULSE_TICKS);
  SelCalib_WF(y, sel_wf, MEMORYDEPTH, PREPULSE_TICKS, SAT_LOW, SAT_UP, BSL);  
  
  TH1D* hI = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP);
  //TH1D* hI = BuildRawChargeHisto(sel_wf, int_wf, MEMORYDEPTH, INT_LOW, INT_UP, HMIN, HMAX);

  double par[NMAXPEAKS+4] = {0}; //Norm_Consts, sigma_0 , sigma_cell, G

  hI->Rebin(4);
  hI->Draw();gPad->Update();//gPad->WaitPrimitive();

  //go for peaks: create an instance of TSpectrum
  TSpectrum *s = new TSpectrum(NMAXPEAKS);
  int npeaks = s->Search(hI,1,"goff",0.05)+1;
  par[0] = 0;     //peakx[0];
  par[1] = (SPE_LOW+SPE_UP)/2;   //peakx[1]-peakx[0];
  par[2] = (S0_LOW+S0_UP)/2;   //sigma_0
  par[3] = (SC_LOW+SC_UP)/2;    //sigma_cell
  for(int i = 0 ; i < npeaks ; i++){
    par[i+4] = 140;//peaky[i];
  }

  //Gaussian sumatory
  TF1* fgaus = new TF1("fgaus", fRandomName, FIT_LOW, FIT_UP, npeaks+3);
  fRandomName_set(fgaus);
  fgaus->SetParameters(par);
  fgaus->SetParLimits(0, MU0_LOW , MU0_UP);
  //fgaus->FixParameter(0,0);
  fgaus->SetParLimits(1,SPE_LOW,SPE_UP);
  fgaus->SetParLimits(2,S0_LOW,S0_UP);
  fgaus->SetParLimits(3,SC_LOW,SC_UP);
  for(int i = 0 ; i < npeaks ; i++) fgaus->SetParLimits(i+4,0,2700);

  TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,900);
  c1->cd();

  //fit histogram
  hI->Draw();
  cout << "Fit ... " << endl;
  hI->Fit("fgaus", "R");
  cout << "... end fit. " << endl;
  
  c1->Modified();
  c1->Update();
  
  //fit CX
  auto g_CX = Build_CX_Graph(fgaus, hI);
  TF1* f_CX = new TF1("f_CX", fCX, -0.5, 5.5, 2);
  fCX_set(f_CX);
  
  g_CX->Fit("f_CX", "R");
  
  TCanvas *c2 = new TCanvas("c2","c2",20,20,1000,900);
  c2->cd();
  g_CX->Draw();
  c2->Modified();
  c2->Update();
  
}
