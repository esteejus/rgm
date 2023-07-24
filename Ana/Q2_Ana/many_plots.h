#ifndef MANYPLOT_HH
#define MANYPLOT_HH
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"

class many_plots{
 public:
  
  many_plots();
  many_plots(std::string temp_name, std::string temp_title, double xmin, double xmax);
  ~many_plots();

  void Fill_hist_set(bool is_epp, double Q2, double x);
  void Write_hist_set(TFile *f, char fileName[100], TCanvas * myCanvas);
  int binQ2(double q2);
  
 private:
  
  TH1D * h_ep;
  TH1D * h_epp;
  TH1D * h_ep_Q2bin[10];
  TH1D * h_epp_Q2bin[10];

};

#endif

