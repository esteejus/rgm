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
#include "TRandom3.h"


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

void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  int isMC = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
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
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //some MC particles
  TLorentzVector MCbeam(0,0,beam_E,beam_E);
  TLorentzVector MCnucleus_ptr(0,0,0,m_4He);
  TLorentzVector MCdeut_ptr(0,0,0,mD);
  TLorentzVector MCel(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector MClead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector MCrecoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector MCntr(0,0,0,db->GetParticle(2112)->Mass());

  
  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;
  /*
  TH1D * h_xB_ep = new TH1D("xB","xB ep;xB;Counts",100,0.9,2.5);
  hist_list.push_back(h_xB_ep);
  TH1D * h_Q2_ep = new TH1D("Q2","Q2 ep;Q2;Counts",100,1.0,5);
  hist_list.push_back(h_Q2_ep);
  TH1D * h_pLead_ep = new TH1D("pLead_ep","p_{Lead} ep;p_{Lead};Counts",100,0,3);
  hist_list.push_back(h_pLead_ep);
  TH1D * h_mMiss_ep = new TH1D("mMiss_ep","m_{Miss} ep;m_{Miss};Counts",100,0,2);
  hist_list.push_back(h_mMiss_ep);
  TH1D * h_EMiss_ep = new TH1D("EMiss_ep","E_{Miss} ep;E_{Miss};Counts",100,-0.1,0.5);
  hist_list.push_back(h_EMiss_ep);
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","mmiss_pmiss",100,0.0,2,100,0.0,1.0);
  hist_list.push_back(h_mmiss_pmiss);
  TH2D * h_mmiss_kmiss = new TH2D("mmiss_kmiss","mmiss_kmiss",100,0.0,2,100,0.0,1.0);
  hist_list.push_back(h_mmiss_kmiss);
  TH1D * h_pMiss_ep = new TH1D("pMiss_ep","p_{Miss} ep;p_{Miss};Counts",100,0,1.0);
  hist_list.push_back(h_pMiss_ep);
  TH1D * h_kMissZQ_ep = new TH1D("kMissZQ_ep","k_{Miss,ZQ} ep;k_{Miss,ZQ};Counts",100,0.0,1.0);
  hist_list.push_back(h_kMissZQ_ep);
  TH1D * h_kMissZB_ep = new TH1D("kMissZB_ep","k_{Miss,ZB} ep;k_{Miss,ZB};Counts",100,0.0,1.0);
  hist_list.push_back(h_kMissZB_ep);
  */
  TH1D * h_kMissZQ_ep = new TH1D("kMissZQ_ep","(e,e'p) k_{Miss};k_{Miss};Counts",100,0.0,1.0);
  hist_list.push_back(h_kMissZQ_ep);

  TH1D * h_DpMiss_ep = new TH1D("DpMiss_ep","#Delta p_{Miss} ep;#Delta p_{Miss};Counts",100,-0.15,0.15);
  hist_list.push_back(h_DpMiss_ep);
  TH1D * h_DkMissZQ_ep = new TH1D("DkMissZQ_ep","#Delta k_{Miss} ep;#Delta k_{Miss};Counts",100,-0.15,0.15);
  hist_list.push_back(h_DkMissZQ_ep);
  //TH1D * h_DkMissZB_ep = new TH1D("DkMissZB_ep","Dk_{Miss,ZB} ep;Dk_{Miss,ZB};Counts",100,-0.15,0.15);
  //hist_list.push_back(h_DkMissZB_ep);

  TH2D * h_pMiss_kMissZQ_ep = new TH2D("pMiss_kMissZQ_ep","Reconstructed Values;p_{Miss} [GeV]; k_{Miss} [GeV];Counts",100,0.2,1.0,100,0.2,1.0);
  hist_list.push_back(h_pMiss_kMissZQ_ep);

  TH2D * h_pMiss_kMissZQ_gen_ep = new TH2D("pMiss_kMissZQ_gen_ep","Generated Values;p_{Miss} [GeV]; k_{Miss} [GeV];Counts",100,0.2,1.0,100,0.2,1.0);
  hist_list.push_back(h_pMiss_kMissZQ_gen_ep);
  
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next() && counter <10000000000)
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

      //if(vtz_e<-5.5){continue;}
      //if(vtz_e>0){continue;}
      
      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
      
      if(lead.size()!=1){continue;}
      GetLorentzVector_ReconVector(lead_ptr,lead[0]);
      if(!isMC){SetLorentzVector_ThetaCorrection(lead_ptr,lead[0]);}
      SetLorentzVector_EnergyLossCorrection(lead_ptr,lead[0]);
      if(!isMC){SetLorentzVector_MomentumCorrection(lead_ptr,lead[0]);}
      //if(lead[0]->getRegion()!=FD){continue;}
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
      
      
      //if(kmiss_ZQ<0.3){continue;}
      //if(mmiss<0.85){continue;}
      //if(mmiss>1.05){continue;}
      if(xB<1.2){continue;}
      if(Q2<1.5){continue;}
      if(Q2>5){continue;}
      if(mom_lead<1.0){continue;}
      
      /*
      h_xB_ep->Fill(xB,wep);
      h_Q2_ep->Fill(Q2,wep);
      h_pLead_ep->Fill(mom_lead,wep);	  
      h_mMiss_ep->Fill(mmiss,wep);	        
      h_EMiss_ep->Fill(Emiss,wep);	        
      h_mmiss_pmiss->Fill(mmiss,mom_miss,wep);
      h_mmiss_kmiss->Fill(mmiss,kmiss_ZQ,wep);
      h_pMiss_ep->Fill(miss.P(),wep);	  
      */
      h_kMissZQ_ep->Fill(kmiss_ZQ,wep);
      if(mmiss<1.2){
      }
      if(kmiss_ZQ>0.3){
	//h_mMiss_ep->Fill(mmiss,wep);	        
      }
      //h_kMissZB_ep->Fill(kmiss_ZB,wep);	  

      h_pMiss_kMissZQ_ep->Fill(miss.P(),kmiss_ZQ,wep);
      
      if(isMC)
	{
	  MCel.SetXYZM(c12->mcparts()->getPx(0),c12->mcparts()->getPy(0),c12->mcparts()->getPz(0),MCel.M());
	  MClead_ptr.SetXYZM(c12->mcparts()->getPx(1),c12->mcparts()->getPy(1),c12->mcparts()->getPz(1),MClead_ptr.M());
	  TLorentzVector MCq = beam - MCel;
	  TLorentzVector MCmiss = MCq + deut_ptr - MClead_ptr;
	  TLorentzVector MCmiss_LC = MClead_ptr - MCq;      
	  
	  TVector3 MCu_ZQ = MCq.Vect().Unit();
	  double MCpmm_ZQ = MCmiss_LC.E() - MCmiss_LC.Vect().Dot(MCu_ZQ);
	  double MCpmp_ZQ = MCmiss_LC.Vect().Perp(MCu_ZQ);
	  double MCkmiss_ZQ = sqrt(mN*mN*((MCpmp_ZQ*MCpmp_ZQ+mN*mN)/(MCpmm_ZQ*(2*mN-MCpmm_ZQ))) - mN*mN);
	  
	  TVector3 MCu_ZB(0,0,1.0);
	  double MCpmm_ZB = MCmiss_LC.E() - MCmiss_LC.Vect().Dot(MCu_ZB);
	  double MCpmp_ZB = MCmiss_LC.Vect().Perp(MCu_ZB);
	  double MCkmiss_ZB = sqrt(mN*mN*((MCpmp_ZB*MCpmp_ZB+mN*mN)/(MCpmm_ZB*(2*mN-MCpmm_ZB))) - mN*mN);

	  h_DpMiss_ep->Fill(miss.P()-MCmiss.P(),wep);	  
	  h_DkMissZQ_ep->Fill(kmiss_ZQ-MCkmiss_ZQ,wep);	  
	  //h_DkMissZB_ep->Fill(kmiss_ZB-MCkmiss_ZB,wep);	  
	  h_pMiss_kMissZQ_gen_ep->Fill(MCmiss.P(),MCkmiss_ZQ,wep);
	 
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
  /*
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
  h_mmiss_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mmiss_kmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_kMissZQ_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
 
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_kMissZB_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_DpMiss_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_DkMissZQ_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_DkMissZB_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMiss_kMissZQ_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pMiss_kMissZQ_gen_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

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
