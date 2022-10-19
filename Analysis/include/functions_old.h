#ifndef FUNCTIONS_old_HH
#define FUNCTIONS_old_HH
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

const double me = 0.000511;
const double mU = 0.9314941024;
const double mPion = 0.139570;
const double mN = 0.939;
const double mD = 1.8756;
const double mT = 3*mN - 0.007;
const double m_4He = 4.00260325415 * mU - 2*me;
const double c = 29.9792458;

double sq(double x){
  return x*x;
}

double isq(double x){
  return 1/sq(x);
}

double guass(double x, double mu, double sigma){
  return (1/(sigma*sqrt(2*M_PI))) * exp(-sq((x-mu)/sigma)/2);
}

double getMin(double mom){
  if(mom<0.4){return -0.5;}
  if(mom<0.5){return -0.2;}
  if(mom<0.6){return -0.1;}
  if(mom<0.7){return  0.0;}
  if(mom<0.8){return -0.4;}
  return -0.7;
}

double get_mmiss(TVector3 vbeam, TVector3 ve, TVector3 vp){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vp.Mag2());

  TVector3 vmiss = vbeam - ve - vp;
  double emiss = Ebeam + mD - Ee - Ep;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mmiss_initial3(TVector3 vbeam, TVector3 ve, TVector3 vL, TVector3 vR1, TVector3 vR2){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double EL  = sqrt((mN * mN) + vL.Mag2());
  double ER1 = sqrt((mN * mN) + vR1.Mag2());
  double ER2 = sqrt((mN * mN) + vR2.Mag2());

  TVector3 vmiss = ve + vL + vR1 + vR2 - vbeam;
  double emiss = Ee + EL + ER1 + ER2 - Ebeam;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mmiss_T(TVector3 vbeam, TVector3 ve, TVector3 vR1, TVector3 vR2){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double ER1 = sqrt((mN * mN) + vR1.Mag2());
  double ER2 = sqrt((mN * mN) + vR2.Mag2());

  TVector3 vmiss = vbeam - ve - vR1 - vR2;
  double emiss = Ebeam + mT - Ee - ER1 - ER2;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mmiss_alpha(TVector3 vbeam, TVector3 ve, TVector3 v1, TVector3 v2, TVector3 v3){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double E1 = sqrt((mN * mN) + v1.Mag2());
  double E2 = sqrt((mN * mN) + v2.Mag2());
  double E3 = sqrt((mN * mN) + v3.Mag2());

  TVector3 vmiss = vbeam - ve - v1 - v2 - v3;
  double emiss = Ebeam + (m_4He) - Ee - E1 - E2 - E3;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_m2(TVector3 v1, TVector3 v2){
  
  double E1 = sqrt((mN * mN) + v1.Mag2());
  double E2 = sqrt((mN * mN) + v2.Mag2());

  TVector3 vmiss = v1 + v2;
  double emiss = E1 + E2;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mass_second_recoil(TVector3 vbeam, TVector3 ve, TVector3 v1, TVector3 v2, TVector3 v3){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double E1 = sqrt((mN * mN) + v1.Mag2());
  double E2 = sqrt((mN * mN) + v2.Mag2());
  double E3 = sqrt((mN * mN) + v3.Mag2());

  TVector3 vmiss = - vbeam + ve + v1 + v2 + v3;
  double emiss = - Ebeam - mD + Ee + E1 + E2 + E3;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_mpi(TVector3 vbeam, TVector3 ve, TVector3 vp){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vp.Mag2());

  TVector3 vpi = vbeam - ve - vp;

  double epi = Ebeam + mN - Ep - Ee;
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

bool LeadCDProton_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index){
  TVector3 p_b(0,0,Ebeam);

  auto electrons=c12->getByID(11);
  TVector3 p_e;
  p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  TVector3 p_q = p_b - p_e;
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();
  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);

  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  TVector3 p_1 = p_p - p_q;
  TVector3 p_miss = p_1;
  double E_p = sqrt(mN*mN + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = (allParticles[index]->getTime()-c12->event()->getStartTime()) - time_frommom_p;
  //double td = time_frombeta_p - time_frommom_p;
  
  //Angle for finding lead
  double theta_p = p_p.Theta() * 180 / M_PI;
  double theta_pq = p_p.Angle(p_q) * 180 / M_PI;
  
  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
  
  if(allParticles[index]->par()->getCharge()<1){return false;}
  if(!CTOF){return false;}
  
  if(theta_p < 25){return false;}
  if(theta_p > 125){return false;}
  if(theta_pq > 35){return false;}  
  
  //if(vtzdiff < -3){return false;}
  //if(vtzdiff > 3){return false;}
  
  if(td < -0.25){return false;}
  if(td > 0.4){return false;}
  
  if(p_p.Mag() < 0.7){return false;}
  if(beta_p > 1.0){return false;}
  
  return true;
}

bool LeadFDProton_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index){
  TVector3 p_b(0,0,Ebeam);

  auto electrons=c12->getByID(11);
  TVector3 p_e;
  p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
  TVector3 p_q = p_b - p_e;
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();
  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);

  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  TVector3 p_1 = p_p - p_q;
  TVector3 p_miss = p_1;
  double E_p = sqrt(mN*mN + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = time_frombeta_p - time_frommom_p;

  //Angle for finding lead
  double theta_p = p_p.Theta() * 180 / M_PI;
  double theta_pq = p_p.Angle(p_q) * 180 / M_PI;

  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
	
  if(allParticles[index]->par()->getCharge()<1){return false;}
  if(!(FTOF1A || FTOF1B || FTOF2)){return false;}
  if(theta_p < 0){return false;}
  if(theta_p > 50){return false;}
  if(theta_pq > 35){return false;}
  if(td < -1.0){return false;}
  if(td > 1.0){return false;}
  if(vtzdiff < -3){return false;}
  if(vtzdiff > 3){return false;}
  if(p_p.Mag() < 0.7){return false;}
  if(beta_p > 1.0){return false;}

  return true;
}

bool AllProton_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index)
{
  TVector3 p_b(0,0,Ebeam);
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();

  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);
  
  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  double E_p = sqrt(mN*mN + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = time_frombeta_p - time_frommom_p;
    
  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
    
  if(allParticles[index]->par()->getCharge()<1){return false;}
  
  if(p_p.Mag()<0.3){return false;}
  if(beta_p<0.1){return false;}
  if(beta_p>0.98){return false;}
  if(vtzdiff < -3){return false;}
  if(vtzdiff > 3){return false;}
  	
  if(FTOF1A || FTOF1B || FTOF2){
    if(td<-5){return false;}
    if(td>5){return false;}
  }
  else if(CTOF){
    if(td<getMin(p_p.Mag())){return false;}
    if(td>0.75){return false;}
  }
  else{
    return false;
  }

  return true;
}

bool AllProton_LooseCut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index)
{
  TVector3 p_b(0,0,Ebeam);
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();

  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);
  
  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  double E_p = sqrt(mN*mN + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = time_frombeta_p - time_frommom_p;
    
  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
    
  if(allParticles[index]->par()->getCharge()<1){return false;}
  if(p_p.Mag()<0.3){return false;}
  if(beta_p<0.1){return false;}
  if(beta_p>0.98){return false;}
  //if(vtzdiff < -3){return false;}
  //if(vtzdiff > 3){return false;}
	
  if(FTOF1A || FTOF1B || FTOF2){
    if(td<-5){return false;}
    if(td>5){return false;}
  }
  else if(CTOF){
    if(td<getMin(p_p.Mag())){return false;}
    if(td>0.75){return false;}
  }
  else{
    return false;
  }
  return true;
}

bool edepProtonCD_Cut(double edep, double beta){
  if((edep>isq(beta)*2.5) && (edep<(100.0/6.0)) && (beta<1.0) && (beta>0.6)){
    return true;
  }
  if((edep>isq(beta)*2.5) && (edep<isq(beta)*7.0) && (beta<0.7) && (beta>0.425)){
    return true;
  }
  if((edep>isq(beta)*2.5) && (edep<550*(beta-0.24)) && (beta>0.35) && (beta<0.425)){
    return true;    
  }
  if((edep>400*(beta-0.275)) && (edep<550*(beta-0.24)) && (beta>0.0) && (beta<0.35)){
    return true;
  }
  return false;
}

TVector3 getNeutronMom(const std::unique_ptr<clas12::clas12reader>& c12, int index)
{

  auto allParticles = c12->getDetParticles();

  //Scintillator hit location
  bool C1 = (allParticles[index]->sci(clas12::CND1)->getDetector() == 3);
  bool C2 = (allParticles[index]->sci(clas12::CND2)->getDetector() == 3);
  bool C3 = (allParticles[index]->sci(clas12::CND3)->getDetector() == 3);
  
  //Neutrals with beta>0.8 don't get directions assigned
  //We need to get them ourselves
  double nvtx_x = allParticles[index]->par()->getVx();
  double nvtx_y = allParticles[index]->par()->getVy();
  double nvtx_z = allParticles[index]->par()->getVz();
  TVector3 v_nvtx(nvtx_x,nvtx_y,nvtx_z);
  TVector3 v_hit;
  if(C1){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND1)->getX(),allParticles[index]->sci(clas12::CND1)->getY(),allParticles[index]->sci(clas12::CND1)->getZ());}
  else if(C2){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND2)->getX(),allParticles[index]->sci(clas12::CND2)->getY(),allParticles[index]->sci(clas12::CND2)->getZ());}
  else if(C3){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND3)->getX(),allParticles[index]->sci(clas12::CND3)->getY(),allParticles[index]->sci(clas12::CND3)->getZ());}
  TVector3 v_path = v_hit - v_nvtx;

  //Now you can define the momentum
  double beta = allParticles[index]->par()->getBeta();
  double gamma = 1/sqrt(1-(beta*beta));
  double mom = gamma * beta * mN;
  TVector3 v_n;
  v_n.SetMagThetaPhi(mom,v_path.Theta(),v_path.Phi());

  return v_n;
}

double getNeutronTP(const std::unique_ptr<clas12::clas12reader>& c12, int index)
{

  auto allParticles = c12->getDetParticles();

  //Scintillator hit location
  bool C1 = (allParticles[index]->sci(clas12::CND1)->getDetector() == 3);
  bool C2 = (allParticles[index]->sci(clas12::CND2)->getDetector() == 3);
  bool C3 = (allParticles[index]->sci(clas12::CND3)->getDetector() == 3);
  
  //Neutrals with beta>0.8 don't get directions assigned
  //We need to get them ourselves
  double nvtx_x = allParticles[index]->par()->getVx();
  double nvtx_y = allParticles[index]->par()->getVy();
  double nvtx_z = allParticles[index]->par()->getVz();
  TVector3 v_nvtx(nvtx_x,nvtx_y,nvtx_z);
  TVector3 v_hit;
  if(C1){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND1)->getX(),allParticles[index]->sci(clas12::CND1)->getY(),allParticles[index]->sci(clas12::CND1)->getZ());}
  else if(C2){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND2)->getX(),allParticles[index]->sci(clas12::CND2)->getY(),allParticles[index]->sci(clas12::CND2)->getZ());}
  else if(C3){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND3)->getX(),allParticles[index]->sci(clas12::CND3)->getY(),allParticles[index]->sci(clas12::CND3)->getZ());}
  TVector3 v_path = v_hit - v_nvtx;

  double ToF = allParticles[index]->getTime() - c12->event()->getStartTime();

  return (100*ToF)/(v_path.Mag());
}

bool AllNeutron_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index)
{
  TVector3 p_b(0,0,Ebeam);
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();
  if(allParticles[index]->par()->getCharge()!=0){return false;}

  //Scintillator hit location
  bool C1 = (allParticles[index]->sci(clas12::CND1)->getDetector() == 3);
  bool C2 = (allParticles[index]->sci(clas12::CND2)->getDetector() == 3);
  bool C3 = (allParticles[index]->sci(clas12::CND3)->getDetector() == 3);
  
  //Getting the Edep is straight forward
  double edep = allParticles[index]->sci(clas12::CND1)->getEnergy() + allParticles[index]->sci(clas12::CND2)->getEnergy() + allParticles[index]->sci(clas12::CND3)->getEnergy();
  

  //Neutrals with beta>0.8 don't get directions assigned
  //We need to get them ourselves
  double nvtx_x = allParticles[index]->par()->getVx();
  double nvtx_y = allParticles[index]->par()->getVy();
  double nvtx_z = allParticles[index]->par()->getVz();
  TVector3 v_nvtx(nvtx_x,nvtx_y,nvtx_z);
  TVector3 v_hit;
  if(C1){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND1)->getX(),allParticles[index]->sci(clas12::CND1)->getY(),allParticles[index]->sci(clas12::CND1)->getZ());}
  else if(C2){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND2)->getX(),allParticles[index]->sci(clas12::CND2)->getY(),allParticles[index]->sci(clas12::CND2)->getZ());}
  else if(C3){v_hit.SetXYZ(allParticles[index]->sci(clas12::CND3)->getX(),allParticles[index]->sci(clas12::CND3)->getY(),allParticles[index]->sci(clas12::CND3)->getZ());}
  TVector3 v_path = v_hit - v_nvtx;

  //The time of flight is not always in line with beta
  //We want to recalculate beta
  double ToF = allParticles[index]->getTime() - c12->event()->getStartTime();
  //double beta_TOF = (v_path.Mag())/(ToF*c);
  double beta_TOF = (v_path.Mag())/(ToF*c);

  //Now you can define the momentum
  double beta = allParticles[index]->par()->getBeta();
  double gamma = 1/sqrt(1-(beta*beta));
  double mom = gamma * beta * mN;
  TVector3 v_n;
  v_n.SetMagThetaPhi(mom,v_path.Theta(),v_path.Phi());

  if(!(C1||C2||C3)){return false;}
  if((beta-beta_TOF)> 0.01){return false;}
  if((beta-beta_TOF)<-0.01){return false;}
  if(v_hit.Z()>40){return false;}	
  if(v_hit.Z()<-30){return false;}

  //Exclude the forward part of the detector
  if(C3 && (v_hit.Z()>25)){return false;}
  else if(C2 && (v_hit.Z()>20)){return false;}
  else if(C1 && (v_hit.Z()>10)){return false;}

  if(edep<5){return false;}
  if(ToF<0){return false;}
  if(ToF>7){return false;}
  
  return true;
  /*
  if(!(C1||C2||C3)){return false;}
  if(beta_p<0.1){return false;}
  if(beta_p>0.8){return false;}
  if(allParticles[index]->getP()<0.1){return false;}
  if(allParticles[index]->getTheta()*180/M_PI > 160){return false;}
  */
}

bool AllPionPlus_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index)
{
  TVector3 p_b(0,0,Ebeam);
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();
  if(allParticles[index]->par()->getCharge()<1){return false;}

  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);
  
  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  double E_p = sqrt(mPion*mPion + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = time_frombeta_p - time_frommom_p;
    
  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
    
  if(beta_p<0.1){return false;}
  if(beta_p>1.2){return false;}
  if(vtzdiff < -3){return false;}
  if(vtzdiff > 3){return false;}
	
  if(FTOF1A || FTOF1B || FTOF2){
    if(p_p.Mag()<0.1){return false;}
    if(p_p.Mag()>2.0){return false;}
    if(td<-1){return false;}
    if(td>1){return false;}
  }
  if(CTOF){
    if(p_p.Mag()<0.1){return false;}
    if(p_p.Mag()>0.7){return false;}
    if(td<-0.5){return false;}
    if(td>0.75){return false;}
  }
  return true;
}

bool AllPionMinus_Cut(const std::unique_ptr<clas12::clas12reader>& c12, double Ebeam, int index)
{
  TVector3 p_b(0,0,Ebeam);
  auto electrons=c12->getByID(11);
  double vtze = electrons[0]->par()->getVz();  

  auto allParticles = c12->getDetParticles();
  if(allParticles[index]->par()->getCharge()>-1){return false;}

  //Scintillator hit location
  bool FTOF1A = (allParticles[index]->sci(clas12::FTOF1A)->getDetector() == 12);
  bool FTOF1B = (allParticles[index]->sci(clas12::FTOF1B)->getDetector() == 12);
  bool FTOF2 = (allParticles[index]->sci(clas12::FTOF2)->getDetector() == 12);
  bool CTOF = (allParticles[index]->sci(clas12::CTOF)->getDetector() == 4);
  
  //Calculate the time difference
  double path_p = allParticles[index]->getPath();
  TVector3 p_p;
  p_p.SetMagThetaPhi(allParticles[index]->getP(),allParticles[index]->getTheta(),allParticles[index]->getPhi());
  double E_p = sqrt(mPion*mPion + p_p.Mag2());
  double beta_p = allParticles[index]->par()->getBeta();
  double beta_frommom_p = p_p.Mag()/E_p;
  double time_frommom_p = path_p / (c*beta_frommom_p);
  double time_frombeta_p = path_p / (c*beta_p);
  double td = time_frombeta_p - time_frommom_p;
    
  //Vertex information
  double vtzl = allParticles[index]->par()->getVz();
  double vtzdiff = vtze - vtzl;
    
  if(beta_p<0.1){return false;}
  if(beta_p>1.2){return false;}
  if(vtzdiff < -3){return false;}
  if(vtzdiff > 3){return false;}
  if(p_p.Mag()<0.1){return false;}
  if(p_p.Mag()>4.0){return false;}
  if(allParticles[index]->getPid()==11){return false;}
	
  if(FTOF1A || FTOF1B || FTOF2){
    if(td<-1){return false;}
    if(td>1){return false;}
  }
  if(CTOF){
    if(td<-1.5){return false;}
    if(td>1.5){return false;}
  }
  return true;
}


#endif

