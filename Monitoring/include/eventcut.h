#ifndef EVENTCUT_HH
#define EVENTCUT_HH
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "HipoChain.h"

#include "TF1.h"
#include "TCanvas.h"

enum cutName{e_cuts,e_nphe,e_calv,e_calw,e_SF,e_mom,e_vtze,l_cuts,l_pid,l_scint,l_theta,l_thetalq,l_chipid,l_vtzdiff,l_phidiff,lsrc_cuts,lsrc_Q2,lsrc_xB,lsrc_pmiss,lsrc_mmiss,lsrc_loq,rsrc_cuts,rsrc_pid,rsrc_mom,rsrc_chipid,fake};

struct cutInfo{
  bool docut;
  double min;
  double max;
  int count;
  std::string label;
  };

class eventcut{
 public:
  
  eventcut(double E, char * filename);
  eventcut(double E);
  ~eventcut();

  bool getDoCut(cutName thisCut);
  double getCutMin(cutName thisCut);
  double getCutMax(cutName thisCut);
  int getCutCount(cutName thisCut);
  std::string getCutLabel(cutName thisCut);
  std::string getCutName(cutName thisCut);
  std::string getCutInformation(cutName thisCut);

  void set_cuts(char * filename);
  void print_cuts();
  void print_cut_loop(cutName startCut, cutName endCut);
  void print_cut_onPDF(TLatex &myText, cutName thisCut, double &line);
  
  bool electroncut(const std::unique_ptr<clas12::clas12reader>& c12);
  int leadnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool leadSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L);
  int recoilSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int index_L);

  
 private:
  
  cutName hashit(std::string cut_name);

  //Electron Cuts
  bool e_nphecut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_calvcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_calwcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_SFcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_momcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_vtzecut(const std::unique_ptr<clas12::clas12reader>& c12);


  //Lead Nucleon Cuts
  bool l_pidcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool l_scintcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool l_thetacut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool l_thetalqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool l_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool l_vtzdiffcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool l_phidiffcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);

  //SRC (e,e'N) Cuts
  bool lsrc_Q2cut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool lsrc_xBcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool lsrc_pmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool lsrc_mmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  bool lsrc_loqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i);
  
  //SRC (e,e'NN) Cuts
  bool rsrc_momcut(const std::unique_ptr<clas12::clas12reader>& c12,int j);
  bool rsrc_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12,int j);

  //General Cut
  bool inRange(double x, cutName thisCut);

  const double mN = 0.939;
  const double mD = 1.8756;
  const int proton_number = 2212;
  const int neutron_number = 2112;

  char name;
  double Ebeam;
  TVector3 vbeam;
  std::map<cutName,cutInfo> cutmap;

};

#endif

