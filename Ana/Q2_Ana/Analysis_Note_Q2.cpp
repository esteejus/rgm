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
#include "clas12ana.h"
#include "many_plots.h"
#include "reweighter.h"
#include "TGraphErrors.h"
#include "Corrections.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
double mass_n = db->GetParticle(2112)->Mass();
double mass_p = db->GetParticle(2212)->Mass();
double mass_pi = db->GetParticle(-211)->Mass();
double mD = 1.8756;

double beam_E = 5.98636;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

double binEdges_Q2[] = {1.5,1.65,1.80,1.95,2.10,2.25,2.40,2.70,3.00,3.50,5.0};
int binEdgeslength_Q2 = sizeof(binEdges_Q2)/sizeof(binEdges_Q2[0]) -1;

int binQ2(double q2){
  if(q2<1.65){return 0;}
  else if(q2<1.80){return 1;}
  else if(q2<1.95){return 2;}
  else if(q2<2.10){return 3;}
  else if(q2<2.25){return 4;}
  else if(q2<2.40){return 5;}
  else if(q2<2.70){return 6;}
  else if(q2<3.00){return 7;}
  else if(q2<3.50){return 8;}
  else{return 9;}
}

double binEdges_pmiss[] = {0.4,0.5,0.6,0.75,1.0};
int binEdgeslength_pmiss = sizeof(binEdges_pmiss)/sizeof(binEdges_pmiss[0]) -1;
int binPMISS(double pmiss){
  for(int i = 0; i <= binEdgeslength_pmiss; i++){
    if(pmiss<binEdges_pmiss[i]){
      return i-1;
    }
  }
  return -1;
}



void Usage()
{
  std::cerr << "Usage: ./code isMC A outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

int main(int argc, char ** argv)
{

  if(argc < 5)
    {
      Usage();
      return -1;
    }



  int isMC = atoi(argv[1]);
  int nucleus_A = atoi(argv[2]);
  TString outFile = argv[3];
  char * pdfFile = argv[4];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  int Z=2;
  int N=2;
  if(isMC){
    Z=nucleus_A/2;
    N=nucleus_A/2;
  }
  reweighter newWeight(beam_E,Z,N,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;

  TH1D * h_multlead_ep = new TH1D("mult","mult (e,e'p);mult;Counts",5,-0.5,4.5);
  hist_list.push_back(h_multlead_ep);
  TH1D * h_xB_ep = new TH1D("xB","x_{B} (e,e'p);x_{B};Counts",50,1.2,2.0);
  hist_list.push_back(h_xB_ep);
  TH1D * h_Q2_ep = new TH1D("Q2","Q^{2} (e,e'p);Q^{2} [GeV];Counts",50,1.5,5);
  hist_list.push_back(h_Q2_ep);
  TH1D * h_pLead_ep = new TH1D("pLead_ep","p_{Lead} (e,e'p);p_{Lead} [GeV];Counts",50,1.0,3);
  hist_list.push_back(h_pLead_ep);
  TH1D * h_kMiss_ep = new TH1D("kMiss_ep","k_{Miss} (e,e'p);k_{Miss} [GeV];Counts",50,0.3,1.0);
  hist_list.push_back(h_kMiss_ep);
  TH1D * h_pMiss_ep = new TH1D("pMiss_ep","p_{Miss} (e,e'p);p_{Miss} [GeV];Counts",50,0.2,1.0);
  hist_list.push_back(h_pMiss_ep);
  TH1D * h_mMiss_ep = new TH1D("mMiss_ep","M_{Miss} (e,e'p);M_{Miss} [GeV];Counts",50,0,2);
  hist_list.push_back(h_mMiss_ep);
  TH1D * h_thetamissq_ep = new TH1D("thetamissq_ep","#theta_{miss,q} (e,e'p);#theta_{miss,q} [degrees];Counts",100,100,180);
  hist_list.push_back(h_thetamissq_ep);
  TH1D * h_EMiss_ep = new TH1D("EMiss_ep","E_{Miss} (e,e'p);E_{Miss} [GeV];Counts",50,-0.1,0.5);
  hist_list.push_back(h_EMiss_ep);
  TH1D * h_vtz_e_ep = new TH1D("vtz_e_ep","Z Vertex e (e,e'p);V_{Z}^{e} [cm];Counts",50,-7,1);
  hist_list.push_back(h_vtz_e_ep);
  TH1D * h_vtz_lead_ep = new TH1D("vtz_lead_ep","Z Vertex p (e,e'p);V_{Z}^{p} [cm];Counts",50,-7,1);
  hist_list.push_back(h_vtz_lead_ep);
  TH1D * h_vtz_diff_ep = new TH1D("vtz_diff_ep","Z Vertex Difference (e,e'p);(V_{Z}^{e}-V_{Z}^{p})/#sqrt{#sigma_{Z}^{e}-#sigma_{Z}^{p}};Counts",50,-5,5);
  hist_list.push_back(h_vtz_diff_ep);
  TH2D * h_mom_SF_ep = new TH2D("mom_SF_ep","Sampling Fraction vs. p for Electrons (e,e'p);p [GeV];Sampling Fraction;Counts",50,3.5,6.0,50,0.17,0.33);
  TH2D * h_mom_beta_ep = new TH2D("mom_beta_ep","#beta vs. p for Protons (e,e'p);p [GeV];#beta;Counts",50,1.0,3.0,50,0.6,1.0);

  TH2D * h_vtz_FD = new TH2D("vtz_FD","vtz_FD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_FD);
  TH2D * h_vtz_CD = new TH2D("vtz_CD","vtz_CD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_CD);
  TH2D * h_vtz_epp_FD_FD = new TH2D("vtz_epp_FD_FD","vtz_epp_FD_FD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_epp_FD_FD);
  TH2D * h_vtz_epp_FD_CD = new TH2D("vtz_epp_FD_CD","vtz_epp_FD_CD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_epp_FD_CD);
  TH2D * h_vtz_epp_CD_FD = new TH2D("vtz_epp_CD_FD","vtz_epp_CD_FD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_epp_CD_FD);
  TH2D * h_vtz_epp_CD_CD = new TH2D("vtz_epp_CD_CD","vtz_epp_CD_CD",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_epp_CD_CD);
  TH1D * h_thetaleadrec_epp = new TH1D("thetaleadrec_epp","thetaleadrec_epp",100,0,180);
  hist_list.push_back(h_thetaleadrec_epp);
  TH1D * h_thetaqrec_epp = new TH1D("thetaqrec_epp","thetaqrec_epp",100,0,180);
  hist_list.push_back(h_thetaqrec_epp);
  TH1D * h_thetaleadq_epp = new TH1D("thetaleadq_epp","thetaleadq_epp",100,0,180);
  hist_list.push_back(h_thetaleadq_epp);
  
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","mmiss_pmiss",100,0.0,2,100,0.3,1.0);
  hist_list.push_back(h_mmiss_pmiss);

  TH2D * h_mmiss_kmiss = new TH2D("mmiss_kmiss","mmiss_kmiss",100,0.0,2,100,0.0,1.0);
  hist_list.push_back(h_mmiss_kmiss);

  
  TH1D * h_mMiss_ep_SomeCuts[10];
  for(int i=0; i<10; i++){
    sprintf(temp_name,"mMiss_ep_SomeCuts_%d",i+1);
    sprintf(temp_title,"mMiss_ep_SomeCuts_%d",i+1);
    h_mMiss_ep_SomeCuts[i] = new TH1D(temp_name,temp_title,100,0,2);    
  }

  
  TH1D * h_xB_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"xB_ep_SRC_%d",i+1);
    sprintf(temp_title,"xB_ep_SRC_%d",i+1);
    h_xB_ep_SRC[i] = new TH1D(temp_name,temp_title,100,1.0,2.0);    
  }
  TH1D * h_thetamissq_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"thetamissq_ep_SRC_%d",i+1);
    sprintf(temp_title,"thetamissq_ep_SRC_%d",i+1);
    h_thetamissq_ep_SRC[i] = new TH1D(temp_name,temp_title,100,100,180);    
  }
  TH1D * h_thetapq_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"thetapq_ep_SRC_%d",i+1);
    sprintf(temp_title,"thetapq_ep_SRC_%d",i+1);
    h_thetapq_ep_SRC[i] = new TH1D(temp_name,temp_title,100,0,50);    
  }
  TH1D * h_poq_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"poq_ep_SRC_%d",i+1);
    sprintf(temp_title,"poq_ep_SRC_%d",i+1);
    h_poq_ep_SRC[i] = new TH1D(temp_name,temp_title,100,0.4,1.0);
  }
  TH1D * h_pmiss_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"pmiss_ep_SRC_%d",i+1);
    sprintf(temp_title,"pmiss_ep_SRC_%d",i+1);
    h_pmiss_ep_SRC[i] = new TH1D(temp_name,temp_title,100,0.4,1.0);
  }
  TH1D * h_emiss_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"emiss_ep_SRC_%d",i+1);
    sprintf(temp_title,"emiss_ep_SRC_%d",i+1);
    h_emiss_ep_SRC[i] = new TH1D(temp_name,temp_title,100,-0.1,0.5);
  }
  TH1D * h_mmiss_ep_SRC[9];
  for(int i=0; i<9; i++){
    sprintf(temp_name,"mmiss_ep_SRC_%d",i+1);
    sprintf(temp_title,"mmiss_ep_SRC_%d",i+1);
    h_mmiss_ep_SRC[i] = new TH1D(temp_name,temp_title,100,0.5,1.6);
  }



  vector<many_plots> hist_list_ep;

  many_plots h_xB("xB","x_{B}",1.1,2);
  hist_list_ep.push_back(h_xB);
  many_plots h_Q2("Q2","Q^{2}",1.5,5.0);
  hist_list_ep.push_back(h_Q2);
  many_plots h_omega("omega","#omega",0.5,2.5);
  hist_list_ep.push_back(h_omega);
  many_plots h_thetae("thetae","#theta_{e}",10,35);
  hist_list_ep.push_back(h_thetae);
  many_plots h_phie("phie","#phi_{e}",-180,180);
  hist_list_ep.push_back(h_phie);

  many_plots h_plead("plead","p_{Lead}",0.9,3.5);
  hist_list_ep.push_back(h_plead);
  many_plots h_thetalead("thetalead","#theta_{Lead}",0.0,50);
  hist_list_ep.push_back(h_thetalead);
  many_plots h_thetalead_FD("thetalead_FD","FD #theta_{Lead}",0.0,90);
  hist_list_ep.push_back(h_thetalead_FD);
  many_plots h_thetalead_CD("thetalead_CD","CD #theta_{Lead}",0.0,90);
  hist_list_ep.push_back(h_thetalead_CD);
  many_plots h_philead("philead","#phi_{Lead}",-180,180);
  hist_list_ep.push_back(h_philead);
  many_plots h_pmiss("pmiss","p_{miss}",0.2,1.0);
  hist_list_ep.push_back(h_pmiss);
  many_plots h_mmiss("mmiss","m_{miss}",0.5,1.6);
  hist_list_ep.push_back(h_mmiss);
  many_plots h_emiss("emiss","E_{miss}",-0.1,0.6);
  hist_list_ep.push_back(h_emiss);
  many_plots h_thetapq("thetapq","#theta_{Lead,q}",0,40);
  hist_list_ep.push_back(h_thetapq);
  many_plots h_thetamissq("thetamissq","#theta_{miss,q}",100,170);
  hist_list_ep.push_back(h_thetamissq);
  many_plots h_poq("poq","p/q",0.5,1.0);
  hist_list_ep.push_back(h_poq);
  many_plots h_kmissZQ("kmissZQ","k_{missZQ}",0.0,1.0);
  hist_list_ep.push_back(h_kmissZQ);
  many_plots h_kmissZB("kmissZB","k_{missZB}",0.0,1.0);
  hist_list_ep.push_back(h_kmissZB);

  vector<many_plots> hist_list_epp;
  many_plots h_precoil("precoil","p_{recoil}",0.15,1.0);
  hist_list_epp.push_back(h_precoil);
  many_plots h_prel("prel","p_{rel}",0.15,1.0);
  hist_list_epp.push_back(h_prel);
  many_plots h_thetamissrecoil("thetamissrecoil","#theta_{miss,recoil}",0,180);
  hist_list_epp.push_back(h_thetamissrecoil);
  many_plots h_thetacmrel("thetacmrel","#theta_{cm,rel}",0,180);
  hist_list_epp.push_back(h_thetacmrel);
  many_plots h_pcm("pcm","p_{cm}",0.0,1.0);
  hist_list_epp.push_back(h_pcm);
  many_plots h_pcmx("pcmx","p_{X,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmx);
  many_plots h_pcmy("pcmy","p_{||,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmy);
  many_plots h_pcmz("pcmz","p_{miss,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmz);
  many_plots h_E2miss("E2miss","E_{2,miss}",-0.2,0.4);
  hist_list_epp.push_back(h_E2miss);


  TH1D * h_Q2_ep_SRC_pmiss[4];
  TH1D * h_Q2_epp_SRC_pmiss[4];

  TH1D * h_E1miss_ep_SRC_pmiss[4];
  TH1D * h_E1miss_epp_SRC_pmiss[4];
  TH1D * h_E2miss_epp_SRC_pmiss[4];
  for(int i=0; i<4; i++){
    sprintf(temp_title,"%f - %f",binEdges_pmiss[i],binEdges_pmiss[i+1]);

    sprintf(temp_name,"h_Q2_ep_SRC_pmiss_%d",i+1);
    h_Q2_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,binEdgeslength_Q2,binEdges_Q2);
    sprintf(temp_name,"h_Q2_epp_SRC_pmiss_%d",i+1);
    h_Q2_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,binEdgeslength_Q2,binEdges_Q2);

    sprintf(temp_name,"h_E1miss_ep_SRC_pmiss_%d",i+1);
    h_E1miss_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,50,-0.2,0.4);
    sprintf(temp_name,"h_E1miss_epp_SRC_pmiss_%d",i+1);
    h_E1miss_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,50,-0.2,0.4);
    sprintf(temp_name,"h_E2miss_epp_SRC_pmiss_%d",i+1);
    h_E2miss_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,50,-0.2,0.4);
  }


  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next() && counter <100000000000)
    {

      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight;// * newWeight.get_weight_ep(c12->mcparts());
	wepp = original_weight;// * newWeight.get_weight_epp(c12->mcparts());
      }

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
      if(!isMC){
	SetLorentzVector_ThetaCorrection(el,electrons[0]);
	SetLorentzVector_MomentumCorrection(el,electrons[0]);
      }
      if(isMC){SetLorentzVector_MomentumSimulationSmear(el,electrons[0]);}
      
      TLorentzVector q = beam - el;
      double Q2        = -q.M2();
      double omega = q.E();
      double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
      double WSq = (mN*mN) - Q2 + (2*omega*mN);
      double W = sqrt(WSq);
      double vtz_e = electrons[0]->par()->getVz();
	  
      int sector_e = electrons[0]->getSector()-1;
      double mom_e = el.P();
      double theta_e = el.Theta() * 180 / M_PI;
      double phi_e = el.Phi() * 180 / M_PI;
      double shift_e = 7.5;
      shift_e += (sector_e==0)?0:(sector_e==1)?60:(sector_e==2)?120:(sector_e==3 && phi_e>0)?180:(sector_e==3 && phi_e<0)?-180:(sector_e==4)?-120:(sector_e==5)?-60:0;
      phi_e -= shift_e;

      double SF = (electrons[0]->cal(clas12::PCAL)->getEnergy() +  electrons[0]->cal(clas12::ECIN)->getEnergy() +  electrons[0]->cal(clas12::ECOUT)->getEnergy()) / electrons[0]->par()->getP();
      //if(vtz_e<-5.5){continue;}
      //if(vtz_e>0){continue;}
      
      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      
      h_multlead_ep->Fill(lead.size(),wep);
      //if(lead.size()==2){cout<<"here\n\n\n\n";}
      if(lead.size()!=1){continue;}

      GetLorentzVector_ReconVector(lead_ptr,lead[0]);
      if(!isMC){SetLorentzVector_ThetaCorrection(lead_ptr,lead[0]);}
      SetLorentzVector_EnergyLossCorrection(lead_ptr,lead[0]);
      if(!isMC){SetLorentzVector_MomentumCorrection(lead_ptr,lead[0]);}
      if(isMC){SetLorentzVector_MomentumSimulationSmear(lead_ptr,protons[0]);}
      
      //SetLorentzVector(lead_ptr,lead[0]);

      TLorentzVector miss = q + deut_ptr - lead_ptr;
      double mmiss2 = miss.M2();
      double mmiss= sqrt(mmiss2);
      double alphamiss = (miss.E() - miss.Vect().Dot(q.Vect().Unit()))/mN;
      TVector3 miss_neg = -miss.Vect();
      double mom_miss = miss.P();
      
      
      double beta_lead = lead[0]->par()->getBeta();
      double mom_lead = lead_ptr.P();
      double momT_lead = lead_ptr.Vect().Perp();
      double theta_lead = lead_ptr.Theta() * 180 / M_PI;
      double phi_lead = lead_ptr.Phi() * 180 / M_PI;
      double vtz_lead = lead[0]->par()->getVz();
      double EP = lead_ptr.E();
      double EB = omega + nucleus_ptr.M() - EP;
      double TB = EB - sqrt(EB*EB - mom_miss*mom_miss);
      double TP = EP - sqrt(EP*EP - mom_lead*mom_lead);
      double Emiss = omega - TP - TB;
      double thetamissq = miss_neg.Angle(q.Vect())*180/M_PI;
      double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
      double poq = mom_lead/q.P();

      TLorentzVector miss_LC = lead_ptr - q;

      TVector3 u_ZQ = q.Vect().Unit();
      double pmm_ZQ = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
      double pmp_ZQ = miss_LC.Vect().Perp(u_ZQ);
      double kmiss_ZQ = sqrt(mN*mN*((pmp_ZQ*pmp_ZQ+mN*mN)/(pmm_ZQ*(2*mN-pmm_ZQ))) - mN*mN);
      
      TVector3 u_ZB(0,0,1.0);
      double pmm_ZB = miss_LC.E() - miss_LC.Vect().Dot(u_ZB);
      double pmp_ZB = miss_LC.Vect().Perp(u_ZB);
      double kmiss_ZB = sqrt(mN*mN*((pmp_ZB*pmp_ZB+mN*mN)/(pmm_ZB*(2*mN-pmm_ZB))) - mN*mN);
      
      if(xB<1.2){continue;}
      if(Q2<1.5){continue;}
      if(Q2>5){continue;}

      if(mom_lead<1.0){continue;}
      if(kmiss_ZQ<0.3){continue;}      
      if(mmiss<0.65){continue;}
      if(mmiss>1.1){continue;}
      if(theta_lead>37){continue;}
      
      //if(fabs(vtz_e-vtz_lead)>1.5){continue;}
      if((lead[0]->getRegion()==FD)){
      }
      else if((lead[0]->getRegion()==CD)){
	h_vtz_CD->Fill(vtz_e,vtz_lead,wep);
	continue;
      }
      h_vtz_FD->Fill(vtz_e,vtz_lead,wep);
      h_vtz_diff_ep->Fill((vtz_e-vtz_lead)/sqrt(0.75*0.75 + 1.4*1.4),wep);
      h_vtz_e_ep->Fill(vtz_e,wep);
      h_vtz_lead_ep->Fill(vtz_lead,wep);
      h_mom_SF_ep->Fill(el.P(),SF,wep);
      h_mom_beta_ep->Fill(lead_ptr.P(),beta_lead,wep);
      
      h_xB_ep->Fill(xB,wep);
      h_Q2_ep->Fill(Q2,wep);
      h_pLead_ep->Fill(mom_lead,wep);	  
      h_kMiss_ep->Fill(kmiss_ZQ,wep);	  
      h_pMiss_ep->Fill(mom_miss,wep);	  
      h_mMiss_ep->Fill(mmiss,wep);	        
      h_thetamissq_ep->Fill(thetamissq,wep);	        
      h_EMiss_ep->Fill(Emiss,wep);	        
      h_mmiss_pmiss->Fill(mmiss,mom_miss,wep);
      h_mmiss_kmiss->Fill(mmiss,kmiss_ZQ,wep);
      
      int fbin = (xB-1.0)*10;
      if((fbin>=0) && (fbin<10)){
	h_mMiss_ep_SomeCuts[fbin]->Fill(mmiss,wep);
      }
      
      int sbin = fbin-1;
      /*
      if((sbin>=1) && (fbin<10)){
	h_xB_ep_SRC[0]->Fill(xB,wep);
	h_xB_ep_SRC[sbin]->Fill(xB,wep);

	h_thetamissq_ep_SRC[0]->Fill(thetamissq,wep);
	h_thetamissq_ep_SRC[sbin]->Fill(thetamissq,wep);

      	h_thetapq_ep_SRC[0]->Fill(thetapq,wep);
	h_thetapq_ep_SRC[sbin]->Fill(thetapq,wep);

	h_poq_ep_SRC[0]->Fill(poq,wep);
	h_poq_ep_SRC[sbin]->Fill(poq,wep);

      	h_pmiss_ep_SRC[0]->Fill(mom_miss,wep);
	h_pmiss_ep_SRC[sbin]->Fill(mom_miss,wep);

	h_emiss_ep_SRC[0]->Fill(Emiss,wep);
	h_emiss_ep_SRC[sbin]->Fill(Emiss,wep);

	h_mmiss_ep_SRC[0]->Fill(mmiss,wep);
	h_mmiss_ep_SRC[sbin]->Fill(mmiss,wep);

      }      
      */
      bool rec = false;
      if(recoil.size()==1){rec = true;}

      /*
      if(binPMISS(mom_miss)!=-1){
	h_Q2_ep_SRC_pmiss[binPMISS(mom_miss)]->Fill(Q2,wep);
	h_E1miss_ep_SRC_pmiss[binPMISS(mom_miss)]->Fill(Emiss,wep);
	if(rec){
	  h_Q2_epp_SRC_pmiss[binPMISS(mom_miss)]->Fill(Q2,wepp);
	  h_E1miss_epp_SRC_pmiss[binPMISS(mom_miss)]->Fill(Emiss,wep);
	}
      }
      */
      h_xB.Fill_hist_set(rec,Q2,xB,wep,wepp);
      h_Q2.Fill_hist_set(rec,Q2,Q2,wep,wepp);
      h_omega.Fill_hist_set(rec,Q2,omega,wep,wepp);
      h_thetae.Fill_hist_set(rec,Q2,el.Theta()*180/M_PI,wep,wepp);
      h_phie.Fill_hist_set(rec,Q2,el.Phi()*180/M_PI,wep,wepp);
      h_plead.Fill_hist_set(rec,Q2,lead_ptr.P(),wep,wepp);
      h_thetalead.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI,wep,wepp);
      if(lead[0]->getRegion()==FD){
	h_thetalead_FD.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI,wep,wepp);
      }
      else if(lead[0]->getRegion()==CD){
	h_thetalead_CD.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI,wep,wepp);	      
      }
      h_philead.Fill_hist_set(rec,Q2,lead_ptr.Phi()*180/M_PI,wep,wepp);
      h_pmiss.Fill_hist_set(rec,Q2,miss.P(),wep,wepp);
      h_mmiss.Fill_hist_set(rec,Q2,miss.M(),wep,wepp);
      h_emiss.Fill_hist_set(rec,Q2,Emiss,wep,wepp);	      	 
      h_thetapq.Fill_hist_set(rec,Q2,lead_ptr.Angle(q.Vect())*180/M_PI,wep,wepp);
      h_thetamissq.Fill_hist_set(rec,Q2,thetamissq,wep,wepp);
      h_poq.Fill_hist_set(rec,Q2,lead_ptr.P()/q.P(),wep,wepp);
      h_kmissZQ.Fill_hist_set(rec,Q2,kmiss_ZQ,wep,wepp);
      h_kmissZB.Fill_hist_set(rec,Q2,kmiss_ZB,wep,wepp);

      if(!rec){continue;}

      GetLorentzVector_ReconVector(recoil_ptr,recoil[0]);
      SetLorentzVector_ThetaCorrection(recoil_ptr,recoil[0]);
      if(!isMC){SetLorentzVector_EnergyLossCorrection(recoil_ptr,recoil[0]);}
      SetLorentzVector_MomentumCorrection(recoil_ptr,recoil[0]);
      //SetLorentzVector(recoil_ptr,recoil[0]);
      double TP2 = recoil_ptr.E() - recoil_ptr.M();
      TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
      double TB2 = miss_Am2.E() - miss_Am2.M();
      double E2miss = q.E() - TP - TP2 - TB2;
      
      TVector3 v_rec  = recoil_ptr.Vect();
      TVector3 v_rel  = (miss_neg - v_rec) * 0.5;
      TVector3 v_cm   = miss_neg + v_rec;
      
      TVector3 vz = miss_neg.Unit();
      TVector3 vy = miss_neg.Cross(q.Vect()).Unit();
      TVector3 vx = vz.Cross(vy).Unit();
      
      double vtz_recoil = recoil[0]->par()->getVz();

      /*
      if((lead[0]->getRegion()==FD) && (recoil[0]->getRegion()==FD)){
	h_vtz_epp_FD_FD->Fill(vtz_lead,vtz_recoil,wep);
      }
      else if((lead[0]->getRegion()==FD) && (recoil[0]->getRegion()==CD)){
	h_vtz_epp_FD_CD->Fill(vtz_lead,vtz_recoil,wep);
      }
      else if((lead[0]->getRegion()==CD) && (recoil[0]->getRegion()==FD)){
	h_vtz_epp_CD_FD->Fill(vtz_lead,vtz_recoil,wep);
      }
      else if((lead[0]->getRegion()==CD) && (recoil[0]->getRegion()==CD)){
	h_vtz_epp_CD_CD->Fill(vtz_lead,vtz_recoil,wep);
      }
      h_thetaleadrec_epp->Fill(lead_ptr.Vect().Angle(recoil_ptr.Vect())*180/M_PI,wep);
      h_thetaqrec_epp->Fill(q.Vect().Angle(recoil_ptr.Vect())*180/M_PI,wep);
      h_thetaleadq_epp->Fill(lead_ptr.Vect().Angle(q.Vect())*180/M_PI,wep);

      
      if(binPMISS(mom_miss)!=-1){
	h_E2miss_epp_SRC_pmiss[binPMISS(mom_miss)]->Fill(E2miss,wep);
      }
      */
      
      h_precoil.Fill_hist_set(rec,Q2,recoil_ptr.P(),wep,wepp);
      h_prel.Fill_hist_set(rec,Q2,v_rel.Mag(),wep,wepp);
      h_thetamissrecoil.Fill_hist_set(rec,Q2,miss_neg.Angle(v_rec)*180/M_PI,wep,wepp);
      h_thetacmrel.Fill_hist_set(rec,Q2,v_cm.Angle(v_rel)*180/M_PI,wep,wepp);
      h_pcm.Fill_hist_set(rec,Q2,v_cm.Mag(),wep,wepp);
      h_pcmx.Fill_hist_set(rec,Q2,v_cm.Dot(vx),wep,wepp);
      h_pcmy.Fill_hist_set(rec,Q2,v_cm.Dot(vy),wep,wepp);
      h_pcmz.Fill_hist_set(rec,Q2,v_cm.Dot(vz),wep,wepp);
      h_E2miss.Fill_hist_set(rec,Q2,E2miss,wep,wepp);

      
    }

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }

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
  h_multlead_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_Q2_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pLead_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_kMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_EMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_vtz_diff_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_vtz_e_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_vtz_lead_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mom_SF_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mom_beta_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetamissq_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */

  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_vtz_FD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_vtz_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */

  /*
  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  h_vtz_epp_FD_FD->Draw("colz");
  myCanvas->cd(2);    
  h_vtz_epp_FD_CD->Draw("colz");
  myCanvas->cd(3);    
  h_vtz_epp_CD_FD->Draw("colz");
  myCanvas->cd(4);    
  h_vtz_epp_CD_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
    
  myCanvas->Divide(2,2);
  myCanvas->cd(1);    
  h_thetaleadrec_epp->Draw("colz");
  myCanvas->cd(2);    
  h_thetaqrec_epp->Draw("colz");
  myCanvas->cd(3);    
  h_thetaleadq_epp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mmiss_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mmiss_kmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */

  /*
  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_mMiss_ep_SomeCuts[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  
  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_xB_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_thetamissq_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_thetapq_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_poq_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_pmiss_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_emiss_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i =0; i<9; i++){
    myCanvas->cd(i+1);    
    h_mmiss_ep_SRC[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  



  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_Q2_ep_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_Q2_epp_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_Q2_epp_SRC_pmiss[i]->Divide(h_Q2_ep_SRC_pmiss[i]);
    h_Q2_epp_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_E1miss_ep_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_E1miss_epp_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  for(int i =0; i<4; i++){
    myCanvas->cd(i+1);    
    h_E2miss_epp_SRC_pmiss[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  
  //h_pmiss.Write_ratio_set(f,fileName,myCanvas);
    
  for(int i=0; i<hist_list_ep.size(); i++){
    hist_list_ep[i].Write_hist_set(f,fileName,myCanvas);
  } 

  for(int i=0; i<hist_list_epp.size(); i++){
    hist_list_epp[i].Write_hist_set_epp(f,fileName,myCanvas);
  } 
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
}

/*
	  TLorentzVector balance_ptr = beam + target_ptr - proton_ptr_p1 - proton_ptr_p2;
	  double theta_be = balance_ptr.Vect().Angle(el.Vect());
	  double Eb = balance_ptr.E();
	  double Pb = balance_ptr.P();
	  double K = (pim_ptr.M2() - balance_ptr.M2() - el.M2()) / 2;
	  double a = Pb*Pb*cos(theta_be)*cos(theta_be) - Eb*Eb;
	  double b = -2 * K * Pb * cos(theta_be);
	  double c = K*K - Eb*Eb*el.M2();
	  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
	  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
*/
