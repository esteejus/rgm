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
#include "eventcut.h"

#include "TF1.h"
#include "TCanvas.h"

enum ecuts {e_nphe,e_calv,e_calw,e_SF,e_mom,l_pid,l_scint,l_theta,l_thetalq,l_chipid,lsrc_Q2,lsrc_xB,lsrc_pmiss,lsrc_mmiss,lsrc_loq};

struct cutInfo{
  bool docut;
  double min;
  double max;
  int count;
  };

class eventcut{
 public:
  
  eventcut(double E);
  ~eventcut();
  bool electroncut(const std::unique_ptr<clas12::clas12reader>& c12);
  int leadnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool leadSRCnucleoncut(const std::unique_ptr<clas12::clas12reader>& c12, int i);

  void setl_scint(int i);

 private:

  bool e_nphecut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_calvcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_calwcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_SFcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool e_momcut(const std::unique_ptr<clas12::clas12reader>& c12);


  bool l_pidcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool l_scintcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  bool l_thetacut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  bool l_thetalqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  bool l_chipidcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);


  bool lsrc_Q2cut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool lsrc_xBcut(const std::unique_ptr<clas12::clas12reader>& c12);
  bool lsrc_pmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  bool lsrc_mmisscut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  bool lsrc_loqcut(const std::unique_ptr<clas12::clas12reader>& c12, int i, int pid);
  

  const double mN = 0.939;
  const double mD = 1.8756;
  char name;
  double Ebeam;
  TVector3 vbeam;
  std::map<ecuts,cutInfo> cutmap;
};

#endif

