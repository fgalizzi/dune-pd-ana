//
//  G_Read.hpp
//
//  Created by Federico Galizzi on 28/09/23
//

#ifndef G_Read_hpp
#define G_Read_hpp

#include <stdio.h>

// Read Binary file with #WF=WFs len-tick long
//**********************************************************
void CompleteWF_Binary(std::string fileName, vector<double>& y, int WFs, int len){
//**********************************************************
  double t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < len; i++) {
      file.read( reinterpret_cast<char*>(&t), sizeof(t) );
      y.push_back((double)t);
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

// Read Binary file with #WF=WFs len-tick long
//**********************************************************
void CompleteWFfloat_Binary(std::string fileName, vector<double>& y, int WFs, int len){
//**********************************************************
  float t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < len; i++) {
      file.read( reinterpret_cast<char*>(&t), sizeof(t) );
      y.push_back((double)t);
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

// Read CAEN-Wavedump Binary file with #WF=WFs len-tick long
//**********************************************************
void CAEN_WF_Binary(std::string fileName, vector<double>& y, int WFs, int len){
//**********************************************************
  int skip;
  uint16_t t;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  
  for (int n_wf=0; n_wf < WFs; n_wf++) {
    for (int i=0; i < 6; i++) file.read( reinterpret_cast<char*>(&skip), sizeof(skip) );
    for (int i=0; i < len; i++) {
      file.read( reinterpret_cast<char*>(&t), sizeof(t) );
      y.push_back((double)t);
    }
  }
  std::cout << "The file has been correctly read \t \t" << std::endl;
}

//*********************************************
void CompleteWF_Binary_Swap(std::string fileName, vector<double>& y, int WFs, int len){
//*********************************************
  int skip;
  uint16_t t;
  uint8_t t_8_msb, t_8_lsb;
  std::ifstream file;

  file.open( fileName, std::ios::binary );
  if (!file.is_open()) {
    std::cout << "Error opening file " << std::endl;
    return;
  }
  std::cout << "Running on " << fileName << std::endl;
  
  for (int n_wf=0; n_wf<WFs; n_wf++) {
    for (int i=0; i<len; i++) {
      file.read( reinterpret_cast<char*>( &t ), sizeof( t ) );
      t = OSSwapBigToHostInt16(t);
      y.push_back((double) t);
    }
  }
  
  std::cout << "The file has been correctly read \t \t" << std::endl;
  
}

// Read a .txt and put in a vector
//**********************************************************
void CompleteWF(std::string fileName, vector<double>& y){
//**********************************************************
  double a;
  std::fstream file (fileName.c_str());
  if (file.is_open()) std::cout << "Running on " << fileName << std::endl;
  else {
    std::cout << "Error opening file" << std::endl;
    return;
  }
  while (true) {
    file >> a;
    if (file.eof() == true ){
      cout << "The file has been correctly read \t \t" <<std::endl;
      return;}
    y.push_back(a);
  }
  return;
}



#endif /* G_Read_hpp */
