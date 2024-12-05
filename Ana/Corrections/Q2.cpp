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
  TLorentzVector MCel(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector MClead(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector MCrec(0,0,0,db->GetParticle(2212)->Mass());

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

  TH2D * h_thetap_Q2 = new TH2D("thetap_Q2","thetap_Q2;thetap;Q2",100,5,45,100,0.0,5.0);
  TH2D * h_xB_Q2 = new TH2D("xB_Q2","xB_Q2;xB;Q2",100,0.8,2.0,100,0.0,5.0);
  TH2D * h_mom_theta = new TH2D("mom_theta","mom_theta;mom;theta",100,0.0,6.0,100,0.0,45.0);

  TH1D * h_Delta_Eprime = new TH1D("Delta_Eprime","#Delta p_{e} (e,e'p);#Delta p_{e} [GeV];Counts",100,-0.2,0.2);
  
  TH1D * h_thetaeFD_e = new TH1D("thetaeFD_e","#theta_{e} (e,e');#theta_{e};Counts",100,5,55);
  hist_list.push_back(h_thetaeFD_e);
  TH1D * h_thetaeFD_eMC = new TH1D("thetaeFD_eMC","#theta_{e} (e,e') MC;#theta_{e};Counts",100,5,55);
  hist_list.push_back(h_thetaeFD_eMC);
  TH1D * h_thetaeFD_ep = new TH1D("thetaFD_ep","#theta_{e} (e,e'p);#theta_{e};Counts",100,5,55);
  hist_list.push_back(h_thetaeFD_ep);

  TH1D * h_thetapFD_ep = new TH1D("thetapFD_ep","thetapFD_ep",100,5,155);
  hist_list.push_back(h_thetapFD_ep);
  
  TH1D * h_pcmyMC = new TH1D("pcmyMC","pcmy",100,-0.5,0.5);
  hist_list.push_back(h_pcmyMC);
  
  TH1D * h_Q2FD_e = new TH1D("Q2FD_e","Q^{2} (e,e');Q^{2};Counts",100,0,5);
  hist_list.push_back(h_Q2FD_e);
  TH1D * h_Q2FD_ep = new TH1D("Q2FD_ep","Q^{2} (e,e'p);Q^{2};Counts",100,0,5);
  hist_list.push_back(h_Q2FD_ep);

  TH1D * h_thetaFD_e = new TH1D("thetaFD_e","#theta_{q} (e,e');#theta_{q};Counts",100,5,155);
  hist_list.push_back(h_thetaFD_e);
  TH1D * h_thetaFD_ep = new TH1D("thetaFD_ep","#theta_{q} (e,e'p);#theta_{q};Counts",100,5,155);
  hist_list.push_back(h_thetaFD_ep);

  TH1D * h_phidiffFD_ep = new TH1D("phidiffFD_ep","#Delta #phi (e,e'p);#Delta #phi;Counts",100,-180,180);
  hist_list.push_back(h_phidiffFD_ep);
  TH1D * h_angleFD_ep = new TH1D("angleFD_ep","#theta (e,e'p);#theta;Counts",100,0,50);
  hist_list.push_back(h_angleFD_ep);
  
  
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
	wep = original_weight* newWeight.get_weight_ep(c12->mcparts());
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
      double theta_q = q.Theta()*180/M_PI;
      double phi_e = el.Phi() * 180 / M_PI;
      double shift_e = 7.5;
      shift_e += (sector_e==0)?0:(sector_e==1)?60:(sector_e==2)?120:(sector_e==3 && phi_e>0)?180:(sector_e==3 && phi_e<0)?-180:(sector_e==4)?-120:(sector_e==5)?-60:0;
      phi_e -= shift_e;

      double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
      if(Delta_Eprime>0.1){continue;}
      //if(xB<1.2){continue;}
      //if(Q2<1.5){continue;}
      //if(Q2>5){continue;}

      if(isMC){
	MCel.SetXYZM(c12->mcparts()->getPx(0),c12->mcparts()->getPy(0),c12->mcparts()->getPz(0),MCel.M());
	MClead.SetXYZM(c12->mcparts()->getPx(1),c12->mcparts()->getPy(1),c12->mcparts()->getPz(1),MClead.M());
	MCrec.SetXYZM(c12->mcparts()->getPx(2),c12->mcparts()->getPy(2),c12->mcparts()->getPz(2),MCrec.M());
	TLorentzVector MCq = beam - MCel;
	double MCQ2        = -MCq.M2();
	double MCtheta_e     = MCel.Theta() * 180 / M_PI;
	TVector3 qMC = MCq.Vect();
	TVector3 leadMC = MClead.Vect();
	TVector3 recMC = MCrec.Vect();
	TVector3 missMC = leadMC-qMC; 
	TVector3 cmMC = missMC+recMC;

	TVector3 vzMC = missMC.Unit();
	TVector3 vyMC = missMC.Cross(qMC).Unit();
	TVector3 vxMC = vzMC.Cross(vyMC).Unit();

	h_pcmyMC->Fill(cmMC.Dot(vyMC),wep);
	h_thetaeFD_eMC->Fill(MCtheta_e,wep);
      }
      
      h_Delta_Eprime->Fill(Delta_Eprime,wep);
      h_thetaeFD_e->Fill(theta_e,wep);
      h_thetaFD_e->Fill(theta_q,wep);
      h_Q2FD_e->Fill(Q2,wep);

      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      auto recoil  = clasAna.getRecoilSRC();
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
     
      double phi_p = lead_ptr.Phi()*180/M_PI;
      double phi_q = q.Phi()*180/M_PI;
      double Delta_phi = (phi_p-phi_q);
      if(Delta_phi<-180){Delta_phi+=360;}
      else if(Delta_phi>180){Delta_phi-=360;}

      h_thetapFD_ep->Fill(theta_lead,wep);
      
      if(mom_lead<1.0){continue;}
      if(kmiss_ZQ<0.3){continue;}      
      if(mmiss<0.7){continue;}
      if(mmiss>1.2){continue;}
      if(theta_lead<37){continue;}
      
      //if(fabs(vtz_e-vtz_lead)>1.5){continue;}
      if((lead[0]->getRegion()!=FD)){ continue; }

      h_thetap_Q2->Fill(theta_lead,Q2,wep);
      h_xB_Q2->Fill(xB,Q2,wep);
      h_mom_theta->Fill(mom_e,theta_e,wep);

      
      h_phidiffFD_ep->Fill(Delta_phi,wep);
      h_angleFD_ep->Fill(q.Vect().Angle(lead_ptr.Vect())*180/M_PI,wep);
      //if(fabs(Delta_phi)>3){continue;}
      //if((q.Vect().Angle(lead_ptr.Vect())*180/M_PI)>5){continue;}
      h_thetaeFD_ep->Fill(theta_e,wep);
      h_Q2FD_ep->Fill(Q2,wep);	  
      h_thetaFD_ep->Fill(theta_q,wep);

      
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
  h_Delta_Eprime->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_xB_Q2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_mom_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_Q2FD_e->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaFD_e->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaeFD_e->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  if(isMC){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_thetaeFD_eMC->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    h_pcmyMC->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }
  
  //

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetapFD_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_phidiffFD_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_angleFD_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  //
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_Q2FD_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaFD_ep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_thetaeFD_ep->Draw("colz");
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
