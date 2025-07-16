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
  many_plots(std::string temp_name, std::string temp_title, TFile * inFile);
  ~many_plots();

  void Fill_hist_set(bool is_epp, double Q2, double x, double wep, double wepp);
  void Write_hist_set(TFile *f, char fileName[100], TCanvas * myCanvas);
  void Write_hist_set_epp(TFile *f, char fileName[100], TCanvas * myCanvas);
  void Write_ratio_set(TFile *f, char fileName[100], TCanvas * myCanvas);
  void Draw_hist_set(char fileName[100], TCanvas * myCanvas);
  void Draw_hist_set_same(char fileName[100], TCanvas * myCanvas, many_plots scaler);
  void Draw_hist_set_epp(char fileName[100], TCanvas * myCanvas);
  void Draw_hist_set_epp_same(char fileName[100], TCanvas * myCanvas, many_plots scaler);
  double getScale_ep();
  double getScale_epp();
  double getScale_epQ2(int i);
  double getScale_eppQ2(int i);
  int binQ2(double q2);
  
 private:
  
  TH1D * h_ep;
  TH1D * h_epp;
  TH1D * h_ep_Q2bin[10];
  TH1D * h_epp_Q2bin[10];

};

#endif

