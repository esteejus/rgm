#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH
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

const double mN = 0.939;
const double mD = 1.8756;

double get_mmiss(TVector3 vbeam, TVector3 ve, TVector3 vp){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vp.Mag2());

  TVector3 vmiss = vbeam - ve - vp;
  double emiss = Ebeam + mD - Ee - Ep;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mpi(TVector3 vbeam, TVector3 ve, TVector3 vp){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vp.Mag2());

  TVector3 vpi = ve + vp - vbeam;

  double epi = Ep + Ee - Ebeam - mN;
  double mpi = sqrt((epi * epi) - vpi.Mag2());

  return mpi;
}

double get_phi_diff(TVector3 p_e, TVector3 p_p){

  double e_phi = p_e.Phi()*180/M_PI;
  double p_phi = p_p.Phi()*180/M_PI;

  if(e_phi>p_phi){
    if((e_phi-p_phi)<=180){
      return (e_phi-p_phi);
    }
    else{
      return 360 - (e_phi-p_phi);
    }
  }
  else{
    if((p_phi-e_phi)<=180){
      return (p_phi-e_phi);
    }
    else{
      return 360 - (p_phi-e_phi);
    }
  }
}


#endif

