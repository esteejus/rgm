#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TLine.h>
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
#include "clas12ana.h"
#include "many_plots.h"
#include "reweighter.h"
#include "TGraphErrors.h"
#include "Corrections.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
/*
TRandom3 * thisRand = new TRandom3(0);;

double SmearFD[6][6]={{1.34173,-0.110811,0.00388463,0.28124,-0.00377903,0.00100301},
		      {0.756057,-0.0425847,0.00210147,0.534038,-0.0298776,0.00157995},
		      {1.14129,-0.0889379,0.00342942,-0.382759,0.065036,-0.000616358},
		      {1.0249,-0.075999,0.00293772,0.418168,-0.0225906,0.00151777},
		      {1.29883,-0.111417,0.0038588,-0.328985,0.0601144,-0.000506229},
		      {0.484246,-0.0154072,0.00152646,-0.162174,0.0382665,0.000149443}};

double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}

double DoubQuad(double x, double A, double B, double C, double D, double E, double F){
  double X = Quad(x,A,B,C)*Quad(x,A,B,C) - Quad(x,D,E,F)*Quad(x,D,E,F);
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

void SetLorentzVector_MomentumSimulationSmear(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double smear = 0.01*DoubQuad(theta,SmearFD[sector-1][0],SmearFD[sector-1][1],SmearFD[sector-1][2],SmearFD[sector-1][3],SmearFD[sector-1][4],SmearFD[sector-1][5]);
    mom*=thisRand->Gaus(1.0,smear);
  }
  else if(p->getRegion()==CD){
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}
*/

auto db=TDatabasePDG::Instance();
double mass_n = db->GetParticle(2112)->Mass();
double mass_p = db->GetParticle(2212)->Mass();
double mass_pi = db->GetParticle(-211)->Mass();
double mD = 1.8756;

double beam_E = 5.98636;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

vector<double> bE_xB = {1.1,1.2,1.3,1.4,1.5,1.6,1.7};

int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

void Usage()
{
  std::cerr << "Usage: ./code A outputfile.root outputfile.pdf DataFile.hipo SimFile.hipo \n\n\n";
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

int main(int argc, char ** argv)
{

  if(argc < 6)
    {
      Usage();
      return -1;
    }



  int nucleus_A = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();

  ////////////////////////////
  clas12root::HipoChain chain;
  cout<<"Input file "<<argv[4]<<endl;
  chain.Add(argv[4]);
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();
  ////////////////////////////
  clas12root::HipoChain chain_Sim;
  cout<<"Input file "<<argv[5]<<endl;
  chain_Sim.Add(argv[5]);
  chain_Sim.SetReaderTags({0});
  chain_Sim.db()->turnOffQADB();
  auto config_c12_Sim=chain_Sim.GetC12Reader();
  auto &c12_Sim=chain_Sim.C12ref();

  ////////////////////////////  
  //For both sets
  ////////////////////////////  
  char temp_name[100];
  char temp_title[100];

  int counter = 0;
  int cutcounter = 0;
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  int Z=nucleus_A/2;
  int N=nucleus_A/2;
  //////////////////////////////////

  vector<TH1*> hist_list_Data;
  TH1D * h_xB_Data_ep = new TH1D("xB","xB (e,e'p);xB;Counts",100,0.9,2.5);
  hist_list_Data.push_back(h_xB_Data_ep);
  TH1D * h_Q2_Data_ep = new TH1D("Q2","Q2 (e,e'p);Q2;Counts",100,1.0,5);
  hist_list_Data.push_back(h_Q2_Data_ep);
  TH1D * h_pLead_Data_ep = new TH1D("pLead_Data_ep","p_{Lead} (e,e'p);p_{Lead};Counts",100,0,3);
  hist_list_Data.push_back(h_pLead_Data_ep);
  TH1D * h_pMiss_Data_ep = new TH1D("pMiss_Data_ep","p_{Miss} (e,e'p);p_{Miss};Counts",50,0.2,1.2);
  hist_list_Data.push_back(h_pMiss_Data_ep);
  TH2D * h_xB_mMiss_Data_ep = new TH2D("xB_mMiss_Data_ep","x_{B} vs. m_{Miss} (e,e'p);x_{B};m_{Miss}",50,1.0,2.0,50,0.2,1.8);
  hist_list_Data.push_back(h_xB_mMiss_Data_ep);
  TH1D * h_mMiss_Data_ep[6];
  for(int i=0; i<6; i++){
    double min = bE_xB[i];
    double max = bE_xB[i+1];
    std::string te = "Counts vs. m_{Miss} ("+std::to_string(min).substr(0,4)+"< x_{B} < "+std::to_string(max).substr(0,4)+");m_{Miss} [GeV];Counts";
    sprintf(temp_name,"mMiss_Data_ep_%d",i+1);
    h_mMiss_Data_ep[i] = new TH1D(temp_name,te.c_str(),50,0.2,1.8);    
    hist_list_Data.push_back(h_mMiss_Data_ep[i]);    
  }
  TH1D * h_thetalead_Data_ep = new TH1D("thetalead_Data_ep","Counts vs. #theta_{p} (e,e'p);#theta_{p};Counts",50,5,45);
  hist_list_Data.push_back(h_thetalead_Data_ep);
  TH1D * h_thetalead_Data_epp = new TH1D("thetalead_Data_epp","Counts vs. #theta_{p} (e,e'pp);#theta_{p};Counts",50,5,45);
  hist_list_Data.push_back(h_thetalead_Data_epp);

  TH1D * h_thetaqlead_Data_ep = new TH1D("thetaqlead_Data_ep","Counts vs. #theta_{q,lead} (e,e'p);#theta_{q,lead};Counts",50,0,40);
  hist_list_Data.push_back(h_thetaqlead_Data_ep);
  TH1D * h_thetaqmiss_Data_ep = new TH1D("thetaqmiss_Data_ep","Counts vs. #theta_{q,miss} (e,e'p);#theta_{q,miss};Counts",50,100,170);
  hist_list_Data.push_back(h_thetaqmiss_Data_ep);
  TH1D * h_poq_Data_ep = new TH1D("poq_Data_ep","Counts vs. p_{lead}/q (e,e'p);p_{lead}/q;Counts",50,0.4,1.0);
  hist_list_Data.push_back(h_poq_Data_ep);

  TH1D * h_precoil_Data_epp = new TH1D("precoil_Data_epp","Counts vs. p_{recoil} (e,e'pp);p_{recoil};Counts",50,0.15,1.0);
  hist_list_Data.push_back(h_precoil_Data_epp);
  TH1D * h_prelative_Data_epp = new TH1D("prelative_Data_epp","Counts vs. p_{relative} (e,e'pp);p_{relative};Counts",50,0.15,1.0);
  hist_list_Data.push_back(h_prelative_Data_epp);
  TH1D * h_thetamissrecoil_Data_epp = new TH1D("thetamissrecoil_Data_epp","Counts vs. #theta_{miss,recoil} (e,e'pp);#theta_{miss,recoil};Counts",50,0,180);
  hist_list_Data.push_back(h_thetamissrecoil_Data_epp);
  TH1D * h_thetacmrelative_Data_epp = new TH1D("thetacmrelative_Data_epp","Counts vs. #theta_{cm,relative} (e,e'pp);#theta_{cm,relative};Counts",50,0,180);
  hist_list_Data.push_back(h_thetacmrelative_Data_epp);
  TH1D * h_E2miss_Data_epp = new TH1D("E2miss_Data_epp","Counts vs. E_{2miss} (e,e'pp);E_{2miss};Counts",50,-0.2,0.4);
  hist_list_Data.push_back(h_E2miss_Data_epp);

  
  for(int i=0; i<hist_list_Data.size(); i++){
    hist_list_Data[i]->Sumw2();
    hist_list_Data[i]->GetXaxis()->CenterTitle();
    hist_list_Data[i]->GetYaxis()->CenterTitle();
  }

  reweighter newWeight_kelly(beam_E,Z,N,kelly,"AV18");
  vector<TH1*> hist_list_SimKelly;
  TH1D * h_xB_SimKelly_ep = new TH1D("xB","xB (e,e'p);xB;Counts",100,0.9,2.5);
  hist_list_SimKelly.push_back(h_xB_SimKelly_ep);
  TH1D * h_Q2_SimKelly_ep = new TH1D("Q2","Q2 (e,e'p);Q2;Counts",100,1.0,5);
  hist_list_SimKelly.push_back(h_Q2_SimKelly_ep);
  TH1D * h_pLead_SimKelly_ep = new TH1D("pLead_SimKelly_ep","p_{Lead} (e,e'p);p_{Lead};Counts",100,0,3);
  hist_list_SimKelly.push_back(h_pLead_SimKelly_ep);
  TH1D * h_pMiss_SimKelly_ep = new TH1D("pMiss_SimKelly_ep","p_{Miss} (e,e'p);p_{Miss};Counts",50,0.2,1.2);
  hist_list_SimKelly.push_back(h_pMiss_SimKelly_ep);
  TH2D * h_xB_mMiss_SimKelly_ep = new TH2D("xB_mMiss_SimKelly_ep","x_{B} vs. m_{Miss} (e,e'p);x_{B};m_{Miss}",50,1.0,2.0,50,0.2,1.8);
  hist_list_SimKelly.push_back(h_xB_mMiss_SimKelly_ep);
  TH1D * h_mMiss_SimKelly_ep[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mMiss_SimKelly_ep_%d",i+1);
    sprintf(temp_title,"mMiss_SimKelly_ep_%d",i+1);
    h_mMiss_SimKelly_ep[i] = new TH1D(temp_name,temp_title,50,0.2,1.8);    
    hist_list_SimKelly.push_back(h_mMiss_SimKelly_ep[i]);    
  }
  TH1D * h_thetalead_SimKelly_ep = new TH1D("thetalead_SimKelly_ep","Counts vs. #theta_{p} (e,e'p);#theta_{p};Counts",50,5,45);
  hist_list_SimKelly.push_back(h_thetalead_SimKelly_ep);
  TH1D * h_thetalead_SimKelly_epp = new TH1D("thetalead_SimKelly_epp","Counts vs. #theta_{p} (e,e'p)p;#theta_{p};Counts",50,5,45);
  hist_list_SimKelly.push_back(h_thetalead_SimKelly_epp);

  TH1D * h_thetaqlead_SimKelly_ep = new TH1D("thetaqlead_SimKelly_ep","Counts vs. #theta_{q,lead} (e,e'p);#theta_{q,lead};Counts",50,0,40);
  hist_list_SimKelly.push_back(h_thetaqlead_SimKelly_ep);
  TH1D * h_thetaqmiss_SimKelly_ep = new TH1D("thetaqmiss_SimKelly_ep","Counts vs. #theta_{q,miss} (e,e'p);#theta_{q,miss};Counts",50,100,170);
  hist_list_SimKelly.push_back(h_thetaqmiss_SimKelly_ep);
  TH1D * h_poq_SimKelly_ep = new TH1D("poq_SimKelly_ep","Counts vs. p_{lead}/q (e,e'p);p_{lead}/q;Counts",50,0.4,1.0);
  hist_list_SimKelly.push_back(h_poq_SimKelly_ep);

  TH1D * h_precoil_SimKelly_epp = new TH1D("precoil_SimKelly_epp","Counts vs. p_{recoil} (e,e'pp);p_{recoil};Counts",50,0.15,1.0);
  hist_list_SimKelly.push_back(h_precoil_SimKelly_epp);
  TH1D * h_prelative_SimKelly_epp = new TH1D("prelative_SimKelly_epp","Counts vs. p_{relative} (e,e'pp);p_{relative};Counts",50,0.15,1.0);
  hist_list_SimKelly.push_back(h_prelative_SimKelly_epp);
  TH1D * h_thetamissrecoil_SimKelly_epp = new TH1D("thetamissrecoil_SimKelly_epp","Counts vs. #theta_{miss,recoil} (e,e'pp);#theta_{miss,recoil};Counts",50,0,180);
  hist_list_SimKelly.push_back(h_thetamissrecoil_SimKelly_epp);
  TH1D * h_thetacmrelative_SimKelly_epp = new TH1D("thetacmrelative_SimKelly_epp","Counts vs. #theta_{cm,relative} (e,e'pp);#theta_{cm,relative};Counts",50,0,180);
  hist_list_SimKelly.push_back(h_thetacmrelative_SimKelly_epp);
  TH1D * h_E2miss_SimKelly_epp = new TH1D("E2miss_SimKelly_epp","Counts vs. E_{2miss} (e,e'pp);E_{2miss};Counts",50,-0.2,0.4);
  hist_list_SimKelly.push_back(h_E2miss_SimKelly_epp);

  
  for(int i=0; i<hist_list_SimKelly.size(); i++){
    hist_list_SimKelly[i]->SetLineColor(2);
    hist_list_SimKelly[i]->Sumw2();
    hist_list_SimKelly[i]->GetXaxis()->CenterTitle();
    hist_list_SimKelly[i]->GetYaxis()->CenterTitle();
  }
  
  reweighter newWeight_dipole(beam_E,Z,N,dipole,"AV18");
  vector<TH1*> hist_list_SimDipole;
  TH1D * h_xB_SimDipole_ep = new TH1D("xB","xB ep;xB;Counts",100,0.9,2.5);
  hist_list_SimDipole.push_back(h_xB_SimDipole_ep);
  TH1D * h_Q2_SimDipole_ep = new TH1D("Q2","Q2 ep;Q2;Counts",100,1.0,5);
  hist_list_SimDipole.push_back(h_Q2_SimDipole_ep);
  TH1D * h_pLead_SimDipole_ep = new TH1D("pLead_SimDipole_ep","p_{Lead} ep;p_{Lead};Counts",100,0,3);
  hist_list_SimDipole.push_back(h_pLead_SimDipole_ep);
  TH1D * h_pMiss_SimDipole_ep = new TH1D("pMiss_SimDipole_ep","p_{Miss} ep;p_{Miss};Counts",50,0.2,1.2);
  hist_list_SimDipole.push_back(h_pMiss_SimDipole_ep);
  TH1D * h_mMiss_SimDipole_ep[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mMiss_SimDipole_ep_%d",i+1);
    sprintf(temp_title,"mMiss_SimDipole_ep_%d",i+1);
    h_mMiss_SimDipole_ep[i] = new TH1D(temp_name,temp_title,50,0.7,1.2);    
    hist_list_SimDipole.push_back(h_mMiss_SimDipole_ep[i]);    
  }

  for(int i=0; i<hist_list_SimDipole.size(); i++){
    hist_list_SimDipole[i]->SetLineColor(3);
    hist_list_SimDipole[i]->Sumw2();
    hist_list_SimDipole[i]->GetXaxis()->CenterTitle();
    hist_list_SimDipole[i]->GetYaxis()->CenterTitle();
  }  

  reweighter newWeight_n2lo(beam_E,Z,N,kelly,"N2LO");
  vector<TH1*> hist_list_SimN2lo;
  TH1D * h_xB_SimN2lo_ep = new TH1D("xB","xB ep;xB;Counts",100,0.9,2.5);
  hist_list_SimN2lo.push_back(h_xB_SimN2lo_ep);
  TH1D * h_Q2_SimN2lo_ep = new TH1D("Q2","Q2 ep;Q2;Counts",100,1.0,5);
  hist_list_SimN2lo.push_back(h_Q2_SimN2lo_ep);
  TH1D * h_pLead_SimN2lo_ep = new TH1D("pLead_SimN2lo_ep","p_{Lead} ep;p_{Lead};Counts",100,0,3);
  hist_list_SimN2lo.push_back(h_pLead_SimN2lo_ep);
  TH1D * h_pMiss_SimN2lo_ep = new TH1D("pMiss_SimN2lo_ep","p_{Miss} ep;p_{Miss};Counts",50,0.2,1.2);
  hist_list_SimN2lo.push_back(h_pMiss_SimN2lo_ep);
  TH1D * h_mMiss_SimN2lo_ep[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mMiss_SimN2lo_ep_%d",i+1);
    sprintf(temp_title,"mMiss_SimN2lo_ep_%d",i+1);
    h_mMiss_SimN2lo_ep[i] = new TH1D(temp_name,temp_title,50,0.7,1.2);    
    hist_list_SimN2lo.push_back(h_mMiss_SimN2lo_ep[i]);    
  }

  for(int i=0; i<hist_list_SimN2lo.size(); i++){
    hist_list_SimN2lo[i]->SetLineColor(4);
    hist_list_SimN2lo[i]->Sumw2();
    hist_list_SimN2lo[i]->GetXaxis()->CenterTitle();
    hist_list_SimN2lo[i]->GetYaxis()->CenterTitle();
  }  

  
  while(chain.Next() && counter <100000000000)
    {

      double wep = 1;

      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pims = clasAna.getByPid(-211);
      auto pips = clasAna.getByPid(211);
      if(electrons.size() != 1){continue;}

      GetLorentzVector_ReconVector(el,electrons[0]);
      SetLorentzVector_ThetaCorrection(el,electrons[0]);
      SetLorentzVector_MomentumCorrection(el,electrons[0]);
      TLorentzVector q = beam - el;
      double Q2        = -q.M2();
      double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );

      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      if(lead.size()!=1){continue;}
      if((lead[0]->getRegion()!=FD)){continue;}
      GetLorentzVector_ReconVector(lead_ptr,lead[0]);
      SetLorentzVector_ThetaCorrection(lead_ptr,lead[0]);
      SetLorentzVector_EnergyLossCorrection(lead_ptr,lead[0]);
      SetLorentzVector_MomentumCorrection(lead_ptr,lead[0]);
      
      TLorentzVector miss = q + deut_ptr - lead_ptr;
      TVector3 neg_miss = -miss.Vect();
      double mmiss = miss.M();
      double mom_miss = miss.P();
      double mom_lead = lead_ptr.P();
      double theta_lead = lead_ptr.Theta()*180/M_PI;
      TLorentzVector miss_LC = lead_ptr - q;
      TVector3 u_ZQ = q.Vect().Unit();
      double pmm_ZQ = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ = miss_LC.Vect().Perp(u_ZQ);
      double kmiss_ZQ = sqrt(mN*mN*((pmp_ZQ*pmp_ZQ+mN*mN)/(pmm_ZQ*(2*mN-pmm_ZQ))) - mN*mN);
      
      if(kmiss_ZQ<0.3){continue;}            
      if(Q2<1.5){continue;}
      if(Q2>5){continue;}

      if(binX(bE_xB,xB)!=-1){
	h_mMiss_Data_ep[binX(bE_xB,xB)]->Fill(mmiss,wep);
      }
      h_xB_mMiss_Data_ep->Fill(xB,mmiss,wep);
      
      if(xB<1.2){continue;}
      if(mmiss<0.65){continue;}
      if(mmiss>1.1){continue;}

      h_pLead_Data_ep->Fill(mom_lead,wep);	  
      h_pMiss_Data_ep->Fill(mom_miss,wep);	  
      h_thetalead_Data_ep->Fill(lead_ptr.Theta()*180/M_PI,wep);

      if(mom_lead<1.0){continue;}
      if(theta_lead>37){continue;}
      h_xB_Data_ep->Fill(xB,wep);
      h_Q2_Data_ep->Fill(Q2,wep);
      h_thetaqlead_Data_ep->Fill(lead_ptr.Vect().Angle(q.Vect())*180/M_PI,wep);
      h_thetaqmiss_Data_ep->Fill(neg_miss.Angle(q.Vect())*180/M_PI,wep);
      h_poq_Data_ep->Fill(lead_ptr.P()/q.P(),wep);
      if(recoil.size()!=1){continue;}
      h_thetalead_Data_epp->Fill(lead_ptr.Theta()*180/M_PI,wep);


      GetLorentzVector_ReconVector(recoil_ptr,recoil[0]);
      SetLorentzVector_ThetaCorrection(recoil_ptr,recoil[0]);
      SetLorentzVector_EnergyLossCorrection(recoil_ptr,recoil[0]);
      SetLorentzVector_MomentumCorrection(recoil_ptr,recoil[0]);
      double TP1 = lead_ptr.E() - lead_ptr.M();
      double TP2 = recoil_ptr.E() - recoil_ptr.M();
      TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
      double TB2 = miss_Am2.E() - miss_Am2.M();
      double E2miss = q.E() - TP1 - TP2 - TB2;
      
      TVector3 v_rec  = recoil_ptr.Vect();
      TVector3 v_rel  = (neg_miss - v_rec) * 0.5;
      TVector3 v_cm   = neg_miss + v_rec;
      double vtz_recoil = recoil[0]->par()->getVz();
      
      double precoil = v_rec.Mag();
      double prelative = v_rel.Mag();
      double thetamissrecoil = neg_miss.Angle(v_rec)*180/M_PI;
      double thetacmrelative = v_cm.Angle(v_rel)*180/M_PI;

      h_precoil_Data_epp->Fill(precoil,wep);
      if(precoil<0.3){continue;}
      h_prelative_Data_epp->Fill(prelative,wep);
      h_thetamissrecoil_Data_epp->Fill(thetamissrecoil,wep);
      h_thetacmrelative_Data_epp->Fill(thetacmrelative,wep);
      h_E2miss_Data_epp->Fill(E2miss,wep);
      
    }



  while(chain_Sim.Next() && counter <100000000)
    {

      double wep_SimKelly = c12_Sim->mcevent()->getWeight() * newWeight_kelly.get_weight_ep(c12_Sim->mcparts());
      double wepp_SimKelly = c12_Sim->mcevent()->getWeight() * newWeight_kelly.get_weight_epp(c12_Sim->mcparts());
      double wep_SimDipole = c12_Sim->mcevent()->getWeight() * newWeight_dipole.get_weight_ep(c12_Sim->mcparts());
      double wep_SimN2lo = c12_Sim->mcevent()->getWeight() * newWeight_n2lo.get_weight_ep(c12_Sim->mcparts());

      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12_Sim);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pims = clasAna.getByPid(-211);
      auto pips = clasAna.getByPid(211);
      if(electrons.size() != 1){continue;}

      GetLorentzVector_ReconVector(el,electrons[0]);
      SetLorentzVector_MomentumSimulationSmear(el,electrons[0]);
      TLorentzVector q = beam - el;
      double Q2        = -q.M2();
      double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );

      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      if(lead.size()!=1){continue;}
      if((lead[0]->getRegion()!=FD)){continue;}
      GetLorentzVector_ReconVector(lead_ptr,lead[0]);
      SetLorentzVector_EnergyLossCorrection(lead_ptr,lead[0]);
      SetLorentzVector_MomentumSimulationSmear(lead_ptr,lead[0]);    
      TLorentzVector miss = q + deut_ptr - lead_ptr;
      TVector3 neg_miss = -miss.Vect();
      double mmiss = miss.M();
      double mom_miss = miss.P();
      double mom_lead = lead_ptr.P();
      double theta_lead = lead_ptr.Theta()*180/M_PI;
      TLorentzVector miss_LC = lead_ptr - q;
      TVector3 u_ZQ = q.Vect().Unit();
      double pmm_ZQ = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ = miss_LC.Vect().Perp(u_ZQ);
      double kmiss_ZQ = sqrt(mN*mN*((pmp_ZQ*pmp_ZQ+mN*mN)/(pmm_ZQ*(2*mN-pmm_ZQ))) - mN*mN);
      
      if(kmiss_ZQ<0.3){continue;}            
      if(Q2<1.5){continue;}
      if(Q2>5){continue;}

      if(binX(bE_xB,xB)!=-1){
	h_mMiss_SimKelly_ep[binX(bE_xB,xB)]->Fill(mmiss,wep_SimKelly);
      }
      h_xB_mMiss_SimKelly_ep->Fill(xB,mmiss,wep_SimKelly);

      if(xB<1.2){continue;}
      if(mmiss<0.65){continue;}
      if(mmiss>1.1){continue;}
      
      h_pLead_SimKelly_ep->Fill(mom_lead,wep_SimKelly);	  
      h_pMiss_SimKelly_ep->Fill(mom_miss,wep_SimKelly);	  
      h_thetalead_SimKelly_ep->Fill(lead_ptr.Theta()*180/M_PI,wep_SimKelly);

      if(mom_lead<1.0){continue;}
      if(theta_lead>37){continue;}
      h_xB_SimKelly_ep->Fill(xB,wep_SimKelly);
      h_Q2_SimKelly_ep->Fill(Q2,wep_SimKelly);
      h_thetaqlead_SimKelly_ep->Fill(lead_ptr.Vect().Angle(q.Vect())*180/M_PI,wep_SimKelly);
      h_thetaqmiss_SimKelly_ep->Fill(neg_miss.Angle(q.Vect())*180/M_PI,wep_SimKelly);
      h_poq_SimKelly_ep->Fill(lead_ptr.P()/q.P(),wep_SimKelly);
      if(recoil.size()!=1){continue;}
      h_thetalead_SimKelly_epp->Fill(lead_ptr.Theta()*180/M_PI,wep_SimKelly);


      GetLorentzVector_ReconVector(recoil_ptr,recoil[0]);
      SetLorentzVector_EnergyLossCorrection(recoil_ptr,recoil[0]);
      SetLorentzVector_MomentumSimulationSmear(recoil_ptr,recoil[0]);    
      double TP1 = lead_ptr.E() - lead_ptr.M();
      double TP2 = recoil_ptr.E() - recoil_ptr.M();
      TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
      double TB2 = miss_Am2.E() - miss_Am2.M();
      double E2miss = q.E() - TP1 - TP2 - TB2;
      
      TVector3 v_rec  = recoil_ptr.Vect();
      TVector3 v_rel  = (neg_miss - v_rec) * 0.5;
      TVector3 v_cm   = neg_miss + v_rec;
      double vtz_recoil = recoil[0]->par()->getVz();
      
      double precoil = v_rec.Mag();
      double prelative = v_rel.Mag();
      double thetamissrecoil = neg_miss.Angle(v_rec)*180/M_PI;
      double thetacmrelative = v_cm.Angle(v_rel)*180/M_PI;

      h_precoil_SimKelly_epp->Fill(precoil,wepp_SimKelly);
      if(precoil<0.3){continue;}
      h_prelative_SimKelly_epp->Fill(prelative,wepp_SimKelly);
      h_thetamissrecoil_SimKelly_epp->Fill(thetamissrecoil,wepp_SimKelly);
      h_thetacmrelative_SimKelly_epp->Fill(thetacmrelative,wepp_SimKelly);
      h_E2miss_SimKelly_epp->Fill(E2miss,wepp_SimKelly);

      /*
      h_xB_SimDipole_ep->Fill(xB,wep_SimDipole);
      h_Q2_SimDipole_ep->Fill(Q2,wep_SimDipole);
      h_pLead_SimDipole_ep->Fill(mom_lead,wep_SimDipole);
      h_pMiss_SimDipole_ep->Fill(mom_miss,wep_SimDipole);  

      h_xB_SimN2lo_ep->Fill(xB,wep_SimN2lo);
      h_Q2_SimN2lo_ep->Fill(Q2,wep_SimN2lo);
      h_pLead_SimN2lo_ep->Fill(mom_lead,wep_SimN2lo);	  
      h_pMiss_SimN2lo_ep->Fill(mom_miss,wep_SimN2lo);	  
      */
    }

  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list_Data.size(); i++){
    hist_list_Data[i]->Write();
  }
  for(int i=0; i<hist_list_SimKelly.size(); i++){
    hist_list_SimKelly[i]->Write();
  }

  for(int i=0; i<hist_list_SimKelly.size(); i++){
    hist_list_SimKelly[i]->Scale(hist_list_Data[i]->GetMaximum()/hist_list_SimKelly[i]->GetMaximum());
  }
  /*
  for(int i=0; i<hist_list_SimDipole.size(); i++){
    hist_list_SimDipole[i]->Scale(hist_list_Data[i]->GetMaximum()/hist_list_SimDipole[i]->GetMaximum());
  }
  for(int i=0; i<hist_list_SimN2lo.size(); i++){
    hist_list_SimN2lo[i]->Scale(hist_list_Data[i]->GetMaximum()/hist_list_SimN2lo[i]->GetMaximum());
  }
  */
  
  //Plot on pdf
  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");
  // from ROOT plain style
  myStyle->SetPalette(1,0);
  myStyle->SetOptStat(0);
  // myStyle->SetOptTitle(0);
  myStyle->SetOptDate(0);
  myStyle->SetLabelSize(0.5, "xyz");
  myStyle->SetTitleSize(1.07, "xyz");
  //myStyle->SetTitleOffset(1.0, "xyz");
  //myStyle->SetTitleOffset(4.0, "y");
  myStyle->SetTitleFont(132, "xyz");
  myStyle->SetLabelFont(82, "xyz"); // size of axis values
  myStyle->SetNdivisions(4, "xyz");
  // Some canvas borders and stuff
  myStyle->SetCanvasColor(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasBorderSize(0);
  myStyle->SetPadBottomMargin(0.1);
  myStyle->SetPadTopMargin(0.1);
  myStyle->SetPadLeftMargin(0.1);
  myStyle->SetPadRightMargin(0.1);
  myStyle->SetPadGridX(0);
  myStyle->SetPadGridY(0);
  //myStyle->SetPadTickX(1);
  //myStyle->SetPadTickY(1);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetPaperSize(20, 24);
  myStyle->SetPadBorderMode(0);
  // Title styles and histogram styles
  //myStyle->SetTitleStyle(0000);
  myStyle->SetHistFillStyle(3001);
  myStyle->SetHistLineColor(kBlack);
  myStyle->SetHistLineWidth(2); //Style option to make the plots look a certain way
  myStyle->SetHistFillColor(kYellow);
  
  //myStyle->SetTitleSize(0.0002, "t");
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleFillColor(0);
  //myStyle->SetTitleTextColor(0);
  myStyle->SetTitleAlign(18);
  
  myStyle->SetPalette(kBird);
  //myStyle->SetLabelSize(3.5);
  myStyle->cd();

  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);

  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);
  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_Data_ep->Draw();
  h_xB_SimKelly_ep->Draw("SAME");
  h_xB_SimDipole_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_Q2_Data_ep->Draw();
  h_Q2_SimKelly_ep->Draw("SAME");
  h_Q2_SimDipole_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pLead_Data_ep->Draw();
  h_pLead_SimKelly_ep->Draw("SAME");
  h_pLead_SimN2lo_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMiss_Data_ep->Draw();
  h_pMiss_SimKelly_ep->Draw("SAME");
  h_pMiss_SimN2lo_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_mMiss_Data_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_mMiss_SimKelly_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  
  myCanvas->Divide(3,2);
  for(int i = 0; i<6; i++){
    myCanvas->cd(i+1);    
    h_mMiss_Data_ep[i]->Draw();
    TLine *line = new TLine(1.1,0,1.1,h_mMiss_Data_ep[i]->GetMaximum());
    line->SetLineColor(3);
    line->Draw("SAME");
    TLine *lin2 = new TLine(0.65,0,0.65,h_mMiss_Data_ep[i]->GetMaximum());
    lin2->SetLineColor(3);
    lin2->Draw("SAME");
    h_mMiss_SimKelly_ep[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetalead_Data_ep->Draw();
  h_thetalead_SimKelly_ep->Draw("SAME");
  TLine *line1 = new TLine(37,0,37,h_thetalead_Data_ep->GetMaximum());
  line1->SetLineColor(3);
  line1->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetalead_Data_epp->Draw();
  h_thetalead_SimKelly_epp->Draw("SAME");
  TLine *line2 = new TLine(37,0,37,h_thetalead_Data_epp->GetMaximum());
  line2->SetLineColor(3);
  line2->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetalead_Data_ep->Draw();
  h_thetalead_Data_epp->SetLineColor(1);
  h_thetalead_Data_ep->Scale(h_thetalead_Data_epp->Integral()/h_thetalead_Data_ep->Integral());
  h_thetalead_Data_epp->Draw("SAME");
  line2->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaqlead_Data_ep->Draw();
  h_thetaqlead_SimKelly_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaqmiss_Data_ep->Draw();
  h_thetaqmiss_SimKelly_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_poq_Data_ep->Draw();
  h_poq_SimKelly_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_Data_ep->Draw();
  h_xB_SimKelly_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_Q2_Data_ep->Draw();
  h_Q2_SimKelly_ep->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  /////////////////

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_precoil_Data_epp->Draw();
  h_precoil_SimKelly_epp->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_prelative_Data_epp->Draw();
  h_prelative_SimKelly_epp->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetamissrecoil_Data_epp->Draw();
  h_thetamissrecoil_SimKelly_epp->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetacmrelative_Data_epp->Draw();
  h_thetacmrelative_SimKelly_epp->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_E2miss_Data_epp->Draw();
  h_E2miss_SimKelly_epp->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");
  
  f->Close();

  return 0;
}

  /*
  TH1D * h_mMiss_Data_ep_SomeCuts[10];
  for(int i=0; i<10; i++){
    sprintf(temp_name,"mMiss_Data_ep_SomeCuts_%d",i+1);
    sprintf(temp_title,"mMiss__Dataep_SomeCuts_%d",i+1);
    h_mMiss_Data_ep_SomeCuts[i] = new TH1D(temp_name,temp_title,100,0,2);    
  }
  */
