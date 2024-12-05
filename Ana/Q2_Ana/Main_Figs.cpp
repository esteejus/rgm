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

vector<double> bE_Q2 = {1.5,1.65,1.80,1.95,2.10,2.25,2.40,2.70,3.00,3.50,5.0}; 
vector<double> bE_pmiss = {0.4,0.5,0.6,0.75,1.0};
vector<double> bE_kmiss = {0.2,0.35,0.45,0.6,0.85};

int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

void getGraph(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH2D * h_myhist, double min, double max);


void Usage()
{
  std::cerr << "Usage: ./code isMC A outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";
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
  TLorentzVector elcpy(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptrcpy(0,0,0,db->GetParticle(2212)->Mass());
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

  TH1D * h_Q2_ep = new TH1D("Q2","Q2 ep;Q2;Counts",500,1.5,5);
  hist_list.push_back(h_Q2_ep);
  TH1D * h_Q2_epp = new TH1D("Q2","Q2 epp;Q2;Counts",500,1.5,5);
  hist_list.push_back(h_Q2_epp);

  TH2D * h_Q2_pcmx_epp = new TH2D("Q2_pcmx_epp","Q2_pcmx_epp",bE_Q2.size()-1,&bE_Q2[0],50,-0.75,0.75);
  hist_list.push_back(h_Q2_pcmx_epp);
  TH2D * h_Q2_pcmy_epp = new TH2D("Q2_pcmy_epp","Q2_pcmy_epp",bE_Q2.size()-1,&bE_Q2[0],50,-0.75,0.75);
  hist_list.push_back(h_Q2_pcmy_epp);
  TH2D * h_Q2_pcmz_epp = new TH2D("Q2_pcmz_epp","Q2_pcmz_epp",bE_Q2.size()-1,&bE_Q2[0],50,-0.75,0.75);
  hist_list.push_back(h_Q2_pcmz_epp);
  
  TH1D * h_pMiss_ep = new TH1D("pMiss_ep","p_{Miss} ep;p_{Miss};Counts",50,0.0,1.0);
  hist_list.push_back(h_pMiss_ep);
  TH1D * h_pMiss_epp = new TH1D("pMiss_epp","p_{Miss} epp;p_{Miss};Counts",50,0.0,1.0);
  hist_list.push_back(h_pMiss_epp);

  TH1D * h_kMiss_ep = new TH1D("kMiss_ep","k_{Miss} ep;k_{Miss};Counts",50,bE_kmiss.front(),bE_kmiss.back());
  hist_list.push_back(h_kMiss_ep);
  TH1D * h_kMiss_epp = new TH1D("kMiss_epp","k_{Miss} epp;k_{Miss};Counts",50,bE_kmiss.front(),bE_kmiss.back());
  hist_list.push_back(h_kMiss_epp);

  TH1D * h_phidiff_ep = new TH1D("phidiff_ep","phidiff ep;phidiff;Counts",100,-45,45);
  hist_list.push_back(h_phidiff_ep);
  TH2D * h_phidiff_plead_ep = new TH2D("phidiff_plead_ep","phidiff plead ep;phidiff;plead;Counts",100,-45,45,100,0.5,2.5);
  hist_list.push_back(h_phidiff_plead_ep);
  TH2D * h_phidiff_thetalead_ep = new TH2D("phidiff_thetalead_ep","phidiff thetalead ep;phidiff;thetalead;Counts",100,-45,45,100,30,100);
  hist_list.push_back(h_phidiff_thetalead_ep);
  TH2D * h_philead_thetalead_ep = new TH2D("philead_thetalead_ep","philead thetalead ep;philead;thetalead;Counts",100,-180,180,100,30,100);
  hist_list.push_back(h_philead_thetalead_ep);
  TH1D * h_thetalead_ep = new TH1D("thetalead_ep","thetalead ep;thetalead;Counts",100,0,100);
  hist_list.push_back(h_thetalead_ep);
  TH1D * h_philead_ep = new TH1D("philead_ep","philead ep;philead;Counts",100,-180,180);
  hist_list.push_back(h_philead_ep);
  TH1D * h_thetalead_epp = new TH1D("thetalead_epp","thetalead epp;thetalead;Counts",100,0,100);
  hist_list.push_back(h_thetalead_epp);
  TH1D * h_philead_epp = new TH1D("philead_epp","philead epp;philead;Counts",100,-180,180);
  hist_list.push_back(h_philead_epp);
  TH1D * h_precoil_epp = new TH1D("precoil_epp","precoil epp;precoil;Counts",50,0.3,1.0);
  hist_list.push_back(h_precoil_epp);
  TH1D * h_thetaleadrecoil_epp = new TH1D("thetaleadrecoil_epp","thetaleadrecoil epp;thetaleadrecoil;Counts",50,0,180);
  hist_list.push_back(h_thetaleadrecoil_epp);
  TH1D * h_thetaqrecoil_epp = new TH1D("thetaqrecoil_epp","thetaqrecoil epp;thetaqrecoil;Counts",50,0,180);
  hist_list.push_back(h_thetaqrecoil_epp);
  TH1D * h_thetarecoil_epp = new TH1D("thetarecoil_epp","thetarecoil epp;thetarecoil;Counts",100,0,180);
  hist_list.push_back(h_thetarecoil_epp);
  TH1D * h_phirecoil_epp = new TH1D("phirecoil_epp","phirecoil epp;phirecoil;Counts",100,-180,180);
  hist_list.push_back(h_phirecoil_epp);
  TH1D * h_thetapp_epp = new TH1D("thetapp_epp","thetapp epp;thetapp;Counts",100,0,180);
  hist_list.push_back(h_thetapp_epp);
  TH1D * h_phidiffpp_epp = new TH1D("phidiffpp_epp","phidiffpp epp;phidiffpp;Counts",100,-180,180);
  hist_list.push_back(h_phidiffpp_epp);
 

  TH1D * h_Q2_ep_SRC_pmiss[4];
  TH1D * h_Q2_epp_SRC_pmiss[4];
  TH1D * h_Q2_ep_SRC_kmiss[4];
  TH1D * h_Q2_epp_SRC_kmiss[4];

  for(int i=0; i<4; i++){
    sprintf(temp_title,"%f - %f",bE_pmiss[i],bE_pmiss[i+1]);
    sprintf(temp_name,"h_Q2_ep_SRC_pmiss_%d",i+1);
    h_Q2_ep_SRC_pmiss[i] = new TH1D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0]);
    hist_list.push_back(h_Q2_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_Q2_epp_SRC_pmiss_%d",i+1);
    h_Q2_epp_SRC_pmiss[i] = new TH1D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0]);
    hist_list.push_back(h_Q2_epp_SRC_pmiss[i]);

    sprintf(temp_title,"%f - %f",bE_kmiss[i],bE_kmiss[i+1]);
    sprintf(temp_name,"h_Q2_ep_SRC_kmiss_%d",i+1);
    h_Q2_ep_SRC_kmiss[i] = new TH1D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0]);
    hist_list.push_back(h_Q2_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_Q2_epp_SRC_kmiss_%d",i+1);
    h_Q2_epp_SRC_kmiss[i] = new TH1D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0]);
    hist_list.push_back(h_Q2_epp_SRC_kmiss[i]);
  }

  
  TH2D * h_Q2_E0miss_ep_SRC_pmiss[4];
  TH2D * h_Q2_E1miss_ep_SRC_pmiss[4];
  TH2D * h_Q2_E1miss_epp_SRC_pmiss[4];
  TH2D * h_Q2_E2miss_epp_SRC_pmiss[4];

  TH2D * h_Q2_E0miss_ep_SRC_kmiss[4];
  TH2D * h_Q2_E1miss_ep_SRC_kmiss[4];
  TH2D * h_Q2_E1miss_epp_SRC_kmiss[4];
  TH2D * h_Q2_E2miss_epp_SRC_kmiss[4];
  for(int i=0; i<4; i++){
    sprintf(temp_title,"%f - %f",bE_pmiss[i],bE_pmiss[i+1]);

    sprintf(temp_name,"h_Q2_E0miss_ep_SRC_pmiss_%d",i+1);
    h_Q2_E0miss_ep_SRC_pmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],100,0.05,0.45);
    hist_list.push_back(h_Q2_E0miss_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_Q2_E1miss_ep_SRC_pmiss_%d",i+1);
    h_Q2_E1miss_ep_SRC_pmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.1,0.4);
    hist_list.push_back(h_Q2_E1miss_ep_SRC_pmiss[i]);
    sprintf(temp_name,"h_Q2_E1miss_epp_SRC_pmiss_%d",i+1);
    h_Q2_E1miss_epp_SRC_pmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.2,0.4);
    hist_list.push_back(h_Q2_E1miss_epp_SRC_pmiss[i]);
    sprintf(temp_name,"h_Q2_E2miss_epp_SRC_pmiss_%d",i+1);
    h_Q2_E2miss_epp_SRC_pmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.2,0.4);
    hist_list.push_back(h_Q2_E2miss_epp_SRC_pmiss[i]);


    sprintf(temp_title,"%f - %f",bE_kmiss[i],bE_kmiss[i+1]);

    sprintf(temp_name,"h_Q2_E0miss_ep_SRC_kmiss_%d",i+1);
    h_Q2_E0miss_ep_SRC_kmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],100,0.0,0.55);
    hist_list.push_back(h_Q2_E0miss_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_Q2_E1miss_ep_SRC_kmiss_%d",i+1);
    h_Q2_E1miss_ep_SRC_kmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.1,0.4);
    hist_list.push_back(h_Q2_E1miss_ep_SRC_kmiss[i]);
    sprintf(temp_name,"h_Q2_E1miss_epp_SRC_kmiss_%d",i+1);
    h_Q2_E1miss_epp_SRC_kmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.2,0.4);
    hist_list.push_back(h_Q2_E1miss_epp_SRC_kmiss[i]);
    sprintf(temp_name,"h_Q2_E2miss_epp_SRC_kmiss_%d",i+1);
    h_Q2_E2miss_epp_SRC_kmiss[i] = new TH2D(temp_name,temp_title,bE_Q2.size()-1,&bE_Q2[0],50,-0.2,0.4);
    hist_list.push_back(h_Q2_E2miss_epp_SRC_kmiss[i]);
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
      if(isMC){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
	wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
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

      GetLorentzVector_Corrected(el,electrons[0],isMC);
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
      double phi_e = el.Vect().Phi() * 180 / M_PI;
      double shift_e = 7.5;

      //shift_e += (sector_e==0)?0:(sector_e==1)?60:(sector_e==2)?120:(sector_e==3 && phi_e>0)?180:(sector_e==3 && phi_e<0)?-180:(sector_e==4)?-120:(sector_e==5)?-60:0;
      //phi_e -= shift_e;

      
      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      
      if(lead.size()!=1){continue;}
      GetLorentzVector_Corrected(lead_ptr,lead[0],isMC);

      TLorentzVector miss = q + deut_ptr - lead_ptr;
      double mmiss2 = miss.M2();
      double mmiss= sqrt(mmiss2);
      double alphamiss = (miss.E() - miss.Vect().Dot(q.Vect().Unit()))/mN;
      TVector3 miss_neg = -miss.Vect();
      double mom_miss = miss.P();
      
      double E0miss = sqrt(mom_miss*mom_miss + mN*mN)-mN;
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
      double E1miss = omega - TP - TB;
      double thetamissq = miss_neg.Angle(q.Vect())*180/M_PI;
      double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
      double poq = mom_lead/q.P();

      TLorentzVector miss_LC = lead_ptr - q;

      TVector3 u = q.Vect().Unit();
      double pmm = miss_LC.E() - miss_LC.Vect().Dot(u);
      double pmp = miss_LC.Vect().Perp(u);
      double pmiss = miss.P();
      double kmiss = sqrt(mN*mN*((pmp*pmp+mN*mN)/(pmm*(2*mN-pmm))) - mN*mN);

      double phidiff = (q.Vect().Phi()*180/M_PI)-phi_lead;
      if(phidiff<-180){phidiff+=360;}
      if(phidiff>180){phidiff-=360;}
            
      bool rec = false;
      if(recoil.size()==1){rec = true;}

      if(xB<1.2){continue;}
      if(Q2<1.5){continue;}
      if(Q2>5){continue;}
      if(mmiss<0.65){continue;}
      if(mmiss>1.10){continue;}
      if(lead[0]->getRegion()!=FD){continue;}

      if(kmiss<0.3){continue;}            
      if(mom_lead<1.0){continue;}
      if(theta_lead>37){continue;}

      h_Q2_ep->Fill(Q2,wep);
      h_pMiss_ep->Fill(pmiss,wep);
      h_kMiss_ep->Fill(kmiss,wep);
      if(binX(bE_pmiss,pmiss)!=-1){
	h_Q2_ep_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,wep);
	h_Q2_E0miss_ep_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,E0miss,wep);
	h_Q2_E1miss_ep_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,E1miss,wep);
      }
      if(binX(bE_kmiss,kmiss)!=-1){
	h_Q2_ep_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,wep);
	h_Q2_E0miss_ep_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,E0miss,wep);
	h_Q2_E1miss_ep_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,E1miss,wep);
      }

      //checks
      h_phidiff_ep->Fill(phidiff,wep);
      h_phidiff_plead_ep->Fill(phidiff,lead_ptr.P(),wep);
      h_phidiff_thetalead_ep->Fill(phidiff,lead_ptr.Vect().Theta()*180/M_PI,wep);
      h_philead_thetalead_ep->Fill(phi_lead,lead_ptr.Vect().Theta()*180/M_PI,wep);
      h_thetalead_ep->Fill(lead_ptr.Vect().Theta()*180/M_PI,wep);
      h_philead_ep->Fill(lead_ptr.Vect().Phi()*180/M_PI,wep);
      //

      
      if(!rec){continue;}
      GetLorentzVector_Corrected(recoil_ptr,recoil[0],isMC);

      double vtz_recoil = recoil[0]->par()->getVz();
      double phi_rec = recoil_ptr.Phi() * 180 / M_PI;
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

      double phidiffpp = phi_rec-phi_lead;
      if(phidiffpp<-180){phidiffpp+=360;}
      if(phidiffpp>180){phidiffpp-=360;}
      //if(phidiffpp>0){continue;}
      h_phidiffpp_epp->Fill(phidiffpp,wepp);


      
      h_Q2_epp->Fill(Q2,wepp);
      h_Q2_pcmx_epp->Fill(Q2,v_cm.Dot(vx),wepp);
      h_Q2_pcmy_epp->Fill(Q2,v_cm.Dot(vy),wepp);
      h_Q2_pcmz_epp->Fill(Q2,v_cm.Dot(vz),wepp);
      h_pMiss_epp->Fill(pmiss,wepp);
      h_kMiss_epp->Fill(kmiss,wepp);

      //checks
      h_thetalead_epp->Fill(lead_ptr.Vect().Theta()*180/M_PI,wepp);
      h_philead_epp->Fill(lead_ptr.Vect().Phi()*180/M_PI,wepp);
      h_precoil_epp->Fill(recoil_ptr.Vect().Mag(),wepp);
      h_thetaleadrecoil_epp->Fill(lead_ptr.Vect().Angle(recoil_ptr.Vect())*180/M_PI,wepp);
      h_thetaqrecoil_epp->Fill(q.Vect().Angle(recoil_ptr.Vect())*180/M_PI,wepp);
      h_thetarecoil_epp->Fill(recoil_ptr.Vect().Theta()*180/M_PI,wepp);
      h_phirecoil_epp->Fill(recoil_ptr.Vect().Phi()*180/M_PI,wepp);
      h_thetapp_epp->Fill(lead_ptr.Vect().Angle(recoil_ptr.Vect())*180/M_PI,wepp);
      //

      
      if(binX(bE_pmiss,pmiss)!=-1){
	h_Q2_epp_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,wepp);
	h_Q2_E1miss_epp_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,E1miss,wepp);
	h_Q2_E2miss_epp_SRC_pmiss[binX(bE_pmiss,pmiss)]->Fill(Q2,E2miss,wepp);
      }
      if(binX(bE_kmiss,kmiss)!=-1){
	h_Q2_epp_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,wepp);
	h_Q2_E1miss_epp_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,E1miss,wepp);
	h_Q2_E2miss_epp_SRC_kmiss[binX(bE_kmiss,kmiss)]->Fill(Q2,E2miss,wepp);
      }

      
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

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  TH1D * h_pMiss_epp_clone = (TH1D*) h_pMiss_epp->Clone();
  h_pMiss_epp_clone->Divide(h_pMiss_ep);
  h_pMiss_epp_clone->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  TH1D * h_kMiss_epp_clone = (TH1D*) h_kMiss_epp->Clone();
  h_kMiss_epp_clone->Divide(h_kMiss_ep);
  h_kMiss_epp_clone->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  getGraph(f,myCanvas,fileName,"Q2_pcmx",h_Q2_pcmy_epp,-0.3,0.3);
  getGraph(f,myCanvas,fileName,"Q2_pcmy",h_Q2_pcmy_epp,-0.3,0.3);
  getGraph(f,myCanvas,fileName,"Q2_pcmz",h_Q2_pcmy_epp,-0.3,0.3);

  for(int i = 0; i < 4; i++){
    getGraph(f,myCanvas,fileName,"Q2_E0miss_binPmiss_"+to_string(i)+"_ep",h_Q2_E0miss_ep_SRC_pmiss[i],0.05,0.4);
  }

  for(int i = 0; i < 4; i++){
    getGraph(f,myCanvas,fileName,"Q2_E1miss_binPmiss_"+to_string(i)+"_ep",h_Q2_E1miss_ep_SRC_pmiss[i],0.0,0.5);
  }

  for(int i = 0; i < 4; i++){
    getGraph(f,myCanvas,fileName,"Q2_E0miss_binKmiss_"+to_string(i)+"_ep",h_Q2_E0miss_ep_SRC_kmiss[i],0.0,0.45);
  }

  for(int i = 0; i < 4; i++){
    getGraph(f,myCanvas,fileName,"Q2_E1miss_binKmiss_"+to_string(i)+"_ep",h_Q2_E1miss_ep_SRC_kmiss[i],0.0,0.5);
  }


  for(int i=0; i<hist_list.size(); i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    hist_list[i]->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
}

void getGraph(TFile *f, TCanvas * myCanvas, char fileName[100], string objectName, TH2D * h_myhist, double min, double max){

  TGraphErrors * g_mu = new TGraphErrors();
  g_mu->SetName(("g_mu_"+objectName).c_str());
  TGraphErrors * g_sigma = new TGraphErrors();
  g_sigma->SetName(("g_sigma_"+objectName).c_str());

  int ctr = 0;
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    TH1D * proj = h_myhist->ProjectionY(("h_proj_"+objectName+"_"+to_string(ctr)).c_str(),j+1,j+1);

    //Now preform a guassian fit
    if(proj->GetEntries()<15){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },min,max,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,(max+min)/2);
    gFit->SetParLimits(1,min,max);
    gFit->SetParameter(2,(max-min)/4);
    gFit->SetParLimits(2,0.00,max-min);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",min,max);
    if(gPoint == 0){
      g_mu->SetPoint(g_mu->GetN(),x,gPoint->Parameter(1));
      g_mu->SetPointError(g_mu->GetN()-1,0,gPoint->ParError(1));

      g_sigma->SetPoint(g_sigma->GetN(),x,gPoint->Parameter(2));
      g_sigma->SetPointError(g_sigma->GetN()-1,0,gPoint->ParError(2));
    }
    proj->Write();

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    proj->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_mu->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  g_sigma->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  g_mu->Write();
  g_sigma->Write();
  
}
