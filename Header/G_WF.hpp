//
//  G_WF.hpp
//
//  Created by Federico Galizzi on 02/10/23
//

#ifndef G_WF_hpp
#define G_WF_hpp

#include <stdio.h>
#include "G_Read.hpp"
#include "G_Utility.hpp"


// Subtract baseline to all the len-long waveform in all_wf. Baseline
// computed in pre-trigger
//*********************************************
void SubBaseline(std::vector<double>& all_wf, int len, int pre){
//*********************************************
  if (pre<20) {
    std::cout << "Too few ticks to subtract the baseline \t \t " ;
    exit(0);
  }
  
  if (all_wf.size()%len != 0) {
    std::cout << "Non-integer number of WFs \t \t " ;
    exit(0);
  }
  
  int WFs = (int) all_wf.size()/len;
  double baseline = 0.;
  
  for (int i=0; i<WFs; i++) {
    baseline = 0.;
    for (int j = 0; j<pre-10; j++) baseline += all_wf[i*len+j];

    baseline /= (double) pre-10;
    for (int j=0; j<len; j++) all_wf[i*len+j] = all_wf[i*len+j]-baseline;
  }
}

// Fixing a value
//*********************************************
void SubBaseline(std::vector<double>& all_wf, int len, int pre, double bsl){
//*********************************************
  if (pre<20) {
    std::cout << "Too few ticks to subtract the baseline \t \t " ;
    exit(0);
  }
  
  if (all_wf.size()%len != 0) {
    std::cout << "Non-integer number of WFs \t \t " ;
    exit(0);
  }
  
  int WFs = (int) all_wf.size()/len;
  double baseline = bsl;
  
  for (int i=0; i<WFs; i++) {
    for (int j=0; j<len; j++) all_wf[i*len+j] = all_wf[i*len+j]-baseline;
  }
}

//*********************************************
void SubBaseline_Range(std::vector<double>& all_wf, int len, int low, int up){
//*********************************************
  if (up-low<20) {
    std::cout << "Too few ticks to subtract the baseline \t \t " ;
    exit(0);
  }
  
  if (all_wf.size()%len != 0) {
    std::cout << "Non-integer number of WFs \t \t " ;
    exit(0);
  }
  
  int WFs = (int) all_wf.size()/len;
  double baseline = 0.;
  
  for (int i=0; i<WFs; i++) {
    baseline = 0.;
    for (int j=low; j<up; j++) baseline += all_wf[i*len+j];
    baseline /= double(up-low);
    for (int j=0; j<len; j++) all_wf[i*len+j] = all_wf[i*len+j]-baseline;
  }
}

//*********************************************
void SubBaseline_Invert(std::vector<double>& all_wf, int len, int pre){
//*********************************************
  if (pre<20) {
    std::cout << "Too few ticks to subtract the baseline \t \t " ;
    exit(0);
  }
  
  if (all_wf.size()%len != 0) {
    std::cout << "Non-integer number of WFs \t \t " ;
    exit(0);
  }
  
  int WFs = (int) all_wf.size()/len;
  double baseline = 0.;
  
  for (int i=0; i<WFs; i++) {
    baseline = 0.;
    for (int j = 0; j<pre-10; j++) baseline += all_wf[i*len+j];
    
    baseline /= (double) pre-10;
    for (int j=0; j<len; j++) all_wf[i*len+j] = -all_wf[i*len+j]+baseline;
    
  }
}

// Compute the integral of a WF in [I_low ; I_up] range
//*********************************************
template <typename T, typename U, typename H>
void ComputeIntegral(H* h, std::vector<T>& int_wf, U len, U I_low, U I_up){
//*********************************************
  T I = 0;
  for(U i = I_low; i < I_up; i++) I += h->GetBinContent(i+1);

  int_wf.push_back(I);
}
//*********************************************
template <typename T, typename U, typename H>
double ComputeFprompt(H* h, std::vector<T>& f_wf, U len, U I_low, U I_up, U I_pr){
//*********************************************
  T I = 0;
  T P = 0;
  for(int ibin = I_low; ibin < I_up; ibin++){
    I += h->GetBinContent(ibin+1);
    if(ibin < I_pr) P += h->GetBinContent(ibin+1);
  }

  return P/I;
}


// With the entire set of WFs (all_wf) it  build the calibration histogram integrating [I_low;I_up]
//*********************************************
TH1D* BuildRawChargeHisto(std::vector<Double_t>& all_wf , std::vector<double>& int_wf, int len, int I_low, int I_up){
//*********************************************

  double hmin, hmax;
  int ibin = 1, j=0;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
    
  while( j < all_wf.size()){
    hwf->SetBinContent(ibin,all_wf[j]);
    ibin++;
    j++;
    if(ibin%(len+1)==0){
      ComputeIntegral(hwf, int_wf, len, I_low, I_up);
      hwf->Reset();
      ibin=1;
    }
  }
  
  hmax = *max_element(std::begin(int_wf), std::end(int_wf));
  hmin = *min_element(std::begin(int_wf), std::end(int_wf));
  TH1D* hI  = new TH1D("hI" ,"hI", 2000, hmin, hmax);
  
  for (int i=0; i<int_wf.size(); i++) hI->Fill(int_wf[i]);
  
  return hI;
}

// With the entire set of WFs (all_wf) it  build the calibration histogram integrating [I_low;I_up]
//*********************************************
TH1D* BuildRawChargeHisto(std::vector<Double_t>& all_wf , std::vector<double>& int_wf, int len, int I_low, int I_up, double hmin, double hmax){
//*********************************************
  
  int ibin = 1, j=0;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
    
  while( j < all_wf.size()){
    hwf->SetBinContent(ibin,all_wf[j]);
    ibin++;
    j++;
    if(ibin%(len+1)==0){
      ComputeIntegral(hwf, int_wf, len, I_low, I_up);
      hwf->Reset();
      ibin=1;
    }
  }
  
  TH1D* hI  = new TH1D("hI" ,"hI", 2000, hmin, hmax);
  
  for (int i=0; i<int_wf.size(); i++) hI->Fill(int_wf[i]);
  
  return hI;
}


// With the non-saturating WFs (ns_wf) it  build the F_prompt histogram integrating [I_low;I_up] and [I_low;I_prompt]
//*********************************************
TH1D* BuildFpromptHisto(std::vector<double>& ns_wf, std::vector<double>& mu_wf, int len, int I_low, int I_up, int I_pr, double f_pr){
//*********************************************
  double t;
  int ibin = 1, j=0;
  TH1D* hwf = new TH1D("hwf","hwf",len,0,len);
  TH1D* hI  = new TH1D("hI" ,"hI", 2000, 0., 1.);
    
  while( j < ns_wf.size()){
    hwf->SetBinContent(ibin,ns_wf[j]);
    ibin++; j++;
   
    if(ibin%(len+1)==0){
      t = ComputeFprompt(hwf, mu_wf, len, I_low, I_up, I_pr);
      hI->Fill(t);
      if (t < f_pr) for (int i=0; i<len; i++) mu_wf.push_back(ns_wf[j+i*len]);
      hwf->Reset();
      ibin=1;
    }
  }
  
  return hI;
}

// Select WFs on their integral basis, store them in sel_wf and compute their average spe_wf
//*********************************************
void Avg_Sel_WF (std::vector<Double_t>& all_wf, std::vector<Double_t>& sel_wf, std::vector<Double_t>& spe_wf, const std::vector<double>& int_wf, double I_low, double I_up){
//*********************************************
  int nspe_wf=0;
  int len = all_wf.size()/int_wf.size();
  
  if (spe_wf.size()<4) for (int i = 0; i < len; i++) spe_wf.push_back(0.);
 
  for (int i = 0; i < int_wf.size(); i++) {
    if (int_wf[i] > I_low && int_wf[i] < I_up) {
      nspe_wf += 1;
      for (int j = 0; j < len; j++) sel_wf.push_back(all_wf[i*len+j]);
    }
  }
  
  //Averaging the waveforms
  for (int i = 0; i < nspe_wf; i++) {
    for (unsigned int j = 0; j < len; j++) spe_wf[j] += sel_wf[i*len+j];
  }
  
  for (int i = 0; i < len; i++) spe_wf[i] /= (double) nspe_wf;
  std::cout << "N_sel " << nspe_wf << std::endl;
}

// Build a TH2D with all the PSD of the waveform and compute the average
//*********************************************
TGraph* build_avg_spectral_density(int nsample, double t1, double t0, std::vector<Double_t>& wf) {
//*********************************************
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;
  double nwindow = wf.size()/nsample;
  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;
  

  int nsample_fft = 0.5*nsample+1;
  TGraph* g_avg_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++)
    g_avg_spectral_density->SetPoint(j, j/t1, 0.);
  double ymin = 1e-6;
  double ymax = 9e+1;
  int    nbinsy = 100;
  
  std::vector<double> ybin_exponent = linspace(TMath::Log10(ymin), TMath::Log10(ymax), nbinsy);
  std::vector<double> ybins(100, 0);
  for (int iy=0; iy<100; iy++) ybins[iy] = TMath::Power(10., ybin_exponent[iy]);

  TH2D* h2_spectral_density = new TH2D("h2_spectral_density",
      Form("%s;%s;%s", "Noise spectral density", "Frequency [MHz]", "Amplitude [A.U.]"),
      nsample_fft,
      0.,
      (double) 0.5*nsample/t1,
      nbinsy-1, &ybins.at(0));
  h2_spectral_density->GetXaxis()->CenterTitle();
  h2_spectral_density->GetYaxis()->CenterTitle();
  
  //  FFT
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");

  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) xn[ip] = (double) wf[iw*nsample+ip];
    
    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j])*c_scale;
      h2_spectral_density->Fill(j/t1, xN[j].Rho2());
      g_avg_spectral_density->GetY()[j] += (xN[j].Rho2()*scale);
    }
  }

  TCanvas* cNoise = new TCanvas("cNoise", "Noise spectral density", 100, 100, 800, 600);
  cNoise->SetLogy(1);
  cNoise->SetLogx(1);
  cNoise->SetTicks(1, 1);
  cNoise->SetGrid(1, 1);
  h2_spectral_density->Draw("colz");
  g_avg_spectral_density->SetLineColor(kGray+2);
  g_avg_spectral_density->SetLineWidth(2);
  g_avg_spectral_density->Draw("l");
  
  return g_avg_spectral_density;
}

// Build a TH2D with all the PSD of the waveform and compute the average
//*********************************************
TGraph* build_avg_spectral_density(int nsample, double t1, double t0, std::vector<double>& wf, double res) {
//*********************************************
  double dt = (t1-t0)/nsample;
  const int nsample_ = nsample;
  double nwindow = wf.size()/nsample;
  double    xn[nsample_];
  TComplex  xN[nsample_];
  double    xN_re[nsample_];
  double    xN_im[nsample_];
  double scale = 1./nwindow;
  double c_scale = 1./nsample;
  double t;
  

  int nsample_fft = 0.5*nsample+1;
  TGraph* g_avg_spectral_density = new TGraph(nsample_fft);
  for (int j=0; j<nsample_fft; j++)
    g_avg_spectral_density->SetPoint(j, j/t1, 0.);
  double ymin = -100;
  double ymax = -20;
  int    nbinsy = 100;

  TH2D* h2_spectral_density = new TH2D("h2_spectral_density",
      Form("%s;%s;%s", "Noise spectral density", "Frequency [MHz]", "Power Spectral Density [db]"),
      nsample_fft,
      0.,
      (double) 0.5*nsample/t1,
      nbinsy-1, ymin, ymax);
  
  //  FFT
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &nsample, "M R2C");

  for (int iw=0; iw<nwindow; iw++) {
    for (int ip=0; ip<nsample_; ip++) xn[ip] = (double) wf[iw*nsample+ip];
    
    fft->SetPoints(xn);
    fft->Transform();
    fft->GetPointsComplex(xN_re, xN_im);

    for (int j=0; j<nsample_fft; j++) {
      xN[j] = TComplex(xN_re[j], xN_im[j]);
      t = 10*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res*2)));
      //t = 20*TMath::Log10(xN[j].Rho2()*c_scale/(pow(2,res)));
      h2_spectral_density->Fill(j/t1, t);
      g_avg_spectral_density->GetY()[j] += (t*scale);
    }
  }

  TCanvas* cNoise = new TCanvas("cNoise", "Power Spectral Density", 100, 100, 800, 600);
  cNoise->SetLogz(1);
  cNoise->SetLogx(1);
  cNoise->SetTicks(1, 1);
  cNoise->SetGrid(1, 1);
  h2_spectral_density->Draw("colz");
  g_avg_spectral_density->SetLineColor(kGray+2);
  g_avg_spectral_density->SetLineWidth(2);
  g_avg_spectral_density->Draw("l");
  
  return g_avg_spectral_density;
}


//*********************************************
template <typename T>
void SubVec_to_WFs(vector<T>& y, vector<T>& sub, int len){
//*********************************************
  for (int i = 0; i*len < y.size(); i++) {
    for (unsigned int j = 0; j < len; j++) y[i*len+j] -= sub[j];
  }
}

//*********************************************
template <typename T>
void avgWF (const vector<T>& y, vector<T>& avg_wf, int len){
//*********************************************
  T wfs = (T) y.size()/len;
  
  if (avg_wf.size() < 5) for (int i=0; i<len; i++) avg_wf.push_back(0.);
  
  for (int i = 0; i < wfs; i++) {
    for (unsigned int j = 0; j < len; j++) avg_wf[j] += y[i*len+j];
  }

  for (int i = 0; i < len; i++) avg_wf[i] /= wfs;
}

// Average the non-saturating waveforms (specify the saturating value with sat_low and up)
//*********************************************
template <typename T>
void Avg_NonSat_WF (vector<T>& y, vector<T>& avg_wf, int len, T sat_low, T sat_up){
//*********************************************
  T max_el, min_el;
  int norm=0;
  
  if (avg_wf.size()<10){
    for (size_t i = 0; i < len; i++) avg_wf.push_back(0.);
  }
  
  for (int i = 0; i*len < y.size(); i++) {
    max_el = *max_element( &y[i*len], &y[(i+1)*len-1]);
    min_el = *min_element( &y[i*len], &y[(i+1)*len-1]);
    if (max_el<sat_up && min_el > sat_low) {
     for (unsigned int j = 0; j < len; j++) avg_wf[j] += y[i*len+j];
     norm +=1;
    }
  }
  
  for (int i = 0; i < len; i++) avg_wf[i] /= (T) norm;
  std::cout << "N avg WF " << norm << std::endl;
}

// Pick up the first saturating WF 
//*********************************************
template <typename T>
void Sat_WF(vector<T>& y, vector<T>& y2, int len, T sat_up){
//*********************************************
  T max_el, min_el;
    
  for (int i = 0; i*len < y.size(); i++) {
    max_el = *max_element( &y[i*len], &y[(i+1)*len-1]);
    
    if (max_el>sat_up) {
      for (unsigned int j = 0; j < len; j++) y2.push_back(y[i*len+j]);
      return;
    }
  }
    
  std::cout << "\n \n !!! \n No saturating WF \n Threshold was set to " << sat_up << "\n \n" ;
}

// Return the waveform where the prepulse-trigger ticks are within the
// sat_low-sat_up range
//*********************************************
template <typename T>
void SelCalib_WF(vector<T>& y, vector<T>& y2, int len, int pre, T sat_low, T sat_up, T bsl){
//*********************************************
  T max_el, min_el;
    
  for (int i = 0; i*len < y.size(); i++) {
    max_el = *max_element( &y[i*len], &y[i*len+pre]);
    min_el = *min_element( &y[i*len], &y[i*len+pre]);
    
    if (max_el<bsl && min_el > -bsl) {
      max_el = *max_element(&y[i*len+pre], &y[(i+1)*len-1]);
      min_el = *min_element(&y[i*len+pre], &y[(i+1)*len-1]);
      
      if (max_el<sat_up && min_el > sat_low) {
      for (unsigned int j = 0; j < len; j++) y2.push_back(y[i*len+j]);
      }
    }
  }  
  return;
}

// Return the waveform where the prepulse-trigger ticks are within the
// sat_low-sat_up range
//*********************************************
template <typename T>
void GoodBaseline_WF(vector<T>& y, vector<T>& y2, int len, int pre, T sat_low, T sat_up){
//*********************************************
  T max_el, min_el;
    
  for (int i = 0; i*len < y.size(); i++) {
    max_el = *max_element( &y[i*len], &y[i*len+pre]);
    min_el = *min_element( &y[i*len], &y[i*len+pre]);
    
    if (max_el<sat_up && min_el > sat_low) {
      for (unsigned int j = 0; j < len; j++) y2.push_back(y[i*len+j]);
    }
  }
    
  return;
}


// Averages just the non saturating waveform and gives you the number of those
//*********************************************
template <typename T>
int NonSat_WF(vector<T>& y, vector<T>& y2, int len, T sat_low, T sat_up){
//*********************************************
  T max_el, min_el;
  int N=0;
    
  for (int i = 0; i*len < y.size(); i++) {
    max_el = *max_element( &y[i*len], &y[(i+1)*len-1]);
    min_el = *min_element( &y[i*len], &y[(i+1)*len-1]);
    
    if (max_el<sat_up && min_el > sat_low) {
      for (unsigned int j = 0; j < len; j++) y2.push_back(y[i*len+j]);
      N +=1;
    }
  }
    
  return N;
}

// Display num waveforms belonging to a single vector
//*********************************************
template <typename T>
void DisplayWFs (const std::vector<T>& y, int len, T tt, int num){
//*********************************************
  for (int i = 0; i<num; i++) {
    TH1D *waa = new TH1D("Waveform", "Waveform",len,0,len*tt);
    for(int bin = 0; bin < len; bin++){
      waa->SetBinContent(bin+1, y[i*len+bin]);}
    waa->Draw();gPad->Update();gPad->WaitPrimitive();
    delete waa;
  }
}

// Display num waveforms belonging to two vectors
//*********************************************
template <typename T>
void DisplayWFs (const std::vector<T>& y, const std::vector<T>& y2, int len, T tt, int num){
//*********************************************
  for (int i = 0; i<num; i++) {
    TH1D *h1 = new TH1D("h1", "h1",len,0,len*tt);
    TH1D *h2 = new TH1D("h2", "h2",len,0,len*tt);
    h1->SetLineColor(2);
    h2->SetLineColor(4);
    for(int bin = 0; bin < len; bin++){
      h1->SetBinContent(bin+1, y[i*len+bin]);
      h2->SetBinContent(bin+1, y2[i*len+bin]);
    }
    h1->Draw();h2->Draw("same");gPad->Update();gPad->WaitPrimitive();
    delete h1; delete h2;
  }
}

// Filter all the wfs according to the G filter
//*********************************************
void FilterAllWF(const std::vector<double>& all_wf, std::vector<double>& filt_wf, TComplex G[], int len){
//*********************************************
  double xv[len];
  double* xy;
  double c_scale = 1./len;
  TComplex xV[len]; double xV_re[len]; double xV_im[len];
  TComplex xY[len]; double xY_re[len]; double xY_im[len];
  TVirtualFFT* fft = TVirtualFFT::FFT(1, &len, "M R2C");
  
  for (int i=0; i<all_wf.size()/len; i++) {
    for (int j=0; j<len; j++) xv[j]=all_wf[i*len+j];
    
    fft = TVirtualFFT::FFT(1, &len, "M R2C");
    fft->SetPoints(xv);
    fft->Transform();
    fft->GetPointsComplex(xV_re, xV_im);
    
    for (int j=0; j<len*0.5+1; j++) {
      xV[j] = TComplex(xV_re[j], xV_im[j]);
      xY[j] = G[j]*xV[j]; //G[j]*P[j]*xV[j]
      xY_re[j] = xY[j].Re(); xY_im[j] = xY[j].Im();
    }
    
    fft = TVirtualFFT::FFT(1, &len, "M C2R");
    fft->SetPointsComplex(xY_re, xY_im);
    fft->Transform();
    xy = fft->GetPointsReal();
    
    for (int j=0; j<len; j++) filt_wf.push_back(xy[j]*c_scale);
     
  }
}


 
template <typename T>
void MovingAverageWF (std::vector<T> in, std::vector<T>& out, int w){
  T sum = 0.;
 
  if (out.size()>0) out.erase(out.begin(), out.end());
  if (w <= 0) std::cout << "Change window \n" ;  
  
  for (int i=0; i<in.size(); i++) {
    sum += in[i];

    if (i >= w){
      sum -= in[i-w];
      out.push_back(sum/w);
    } else out.push_back( sum/((T)i+1.) );
 
  } 
}



/*
template <typename T>
void MovingAverageWF (const std::vector<T>& in, std::vector<T>& out, int w){
  T sum = 0.;
 
  if (out.size()>0) out.erase(out.begin(), out.end());
  if (w <= 0) std::cout << "Change window \n" ;  
  
  for (int i=0; i<w; i++){
    sum += in[i];
    out.push_back(sum/((T)i+1.));
  } 

  for (int i=w; i<in.size(); i++) {
    sum += in[i+w]-in[i-w];
    out.push_back(sum/w);
  }

}
*/



#endif /* G_WF_hpp */
