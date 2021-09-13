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

using namespace std;
using namespace clas12;

const double mN = 0.939;
const double mD = 1.8756;

void printProgress(double percentage);

double get_mmiss(TVector3 vbeam, TVector3 ve, TVector3 vp){
  
  double Ebeam = vbeam.Mag();
  double Ee = ve.Mag();
  double Ep = sqrt((mN * mN) + vp.Mag2());

  TVector3 vmiss = vbeam - ve - vp;
  double emiss = Ebeam + mD - Ee - Ep;
  double mmiss = sqrt((emiss * emiss) - vmiss.Mag2());

  return mmiss;
}

double get_phi_diff(double e_phi, double p_phi){

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

bool lowThetaCut(double theta, double chi2PID, double vtzDiff){
  
  if(theta > (50 * M_PI / 180)){
    return false;
  }
  if(fabs(chi2PID-0.459179)>(3*1.2085)){
    return false;
  }
  if(fabs(vtzDiff-0.484268)>(3*1.30286)){
    return false;
  }
  
  return true;
}


void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> [scintillator number (4 or 12)] <path/to/input.hipo> \n";
}


int main(int argc, char ** argv)
{

  if(argc < 7)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  
  bool isMC = false;
  if(atoi(argv[1]) == 1){isMC=true;}

  double Ebeam = atof(argv[2]);
  
  TFile * outFile = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];
  int TOFID = atoi(argv[5]);
  clas12root::HipoChain chain;
  for(int k = 6; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  

  /////////////////////////////////////
  //Prepare histograms
  /////////////////////////////////////
  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(0.8);

  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////
  TH2D * h_Vcal_EoP = new TH2D("Vcal_EoP","ECAL V coordinate vs. Sampling Fraction ;ECAL V coordinate;Sampling Fraction",60,0,30,150,0.05,0.40);
  hist_list_2.push_back(h_Vcal_EoP);
  TH2D * h_Wcal_EoP = new TH2D("Wcal_EoP","ECAL W coordinate vs. Sampling Fraction ;ECAL W coordinate;Sampling Fraction",60,0,30,150,0.05,0.40);
  hist_list_2.push_back(h_Wcal_EoP);
  TH2D * h_phi_theta = new TH2D("phi_theta","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);
  hist_list_2.push_back(h_phi_theta);
  TH1D * h_sector = new TH1D("sector","ECAL Sector;Sector;Counts",6,1,7);
  hist_list_1.push_back(h_sector);

  /////////////////////////////////////
  //Electron Pid
  /////////////////////////////////////
  TH2D * h_P_EoP = new TH2D("P_EoP","p_{e} vs. Sampling Fraction ;p_{e};Sampling Faction",100,0,7,100,0.15,0.35);
  hist_list_2.push_back(h_P_EoP);
  TH1D * h_nphe = new TH1D("nphe","#Photo-electrons in HTCC;#Photo-electrons;Counts",100,0,50);
  hist_list_1.push_back(h_nphe);
  
  /////////////////////////////////////
  //Electron Kinematics  
  /////////////////////////////////////
  TH1D * h_xB = new TH1D("xB","x_{B};x_{B};Counts",100,0,2);
  hist_list_1.push_back(h_xB);
  TH1D * h_QSq = new TH1D("QSq","Q^{2};Q^{2};Counts",100,0,3);
  hist_list_1.push_back(h_QSq);
  TH1D * h_WSq = new TH1D("WSq","W^{2};W^{2}",100,0,7);
  hist_list_1.push_back(h_WSq);
  TH2D * h_xB_QSq = new TH2D("xB_QSq","x_{B} vs. Q^{2} ;x_{B};Q^{2}",100,0,2,100,0,3);
  hist_list_2.push_back(h_xB_QSq);
  TH2D * h_xB_WSq = new TH2D("xB_WSq","x_{B} vs. W^{2} ;x_{B};W^{2}",100,0,2,100,0,7);
  hist_list_2.push_back(h_xB_WSq);
  TH2D * h_QSq_WSq = new TH2D("QSq_WSq","Q^{2} vs. W^{2} ;Q^{2};W^{2}",100,0,3,100,0,7);
  hist_list_2.push_back(h_QSq_WSq);

  /////////////////////////////////////
  //All Proton Angles
  /////////////////////////////////////
  TH1D * h_theta_L = new TH1D("theta_L","#theta_{proton};#theta_{proton};Counts",180,0,180);
  hist_list_1.push_back(h_theta_L);
  TH1D * h_theta_Lq = new TH1D("theta_Lq","#theta_{pq};#theta_{pq};Counts",180,0,180);
  hist_list_1.push_back(h_theta_Lq);

  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
  TH1D * h_theta_L_FTOF = new TH1D("theta_L_FTOF","#theta_{proton} Lead;#theta_{proton};Counts",180,0,180);
  hist_list_1.push_back(h_theta_L_FTOF);
  TH1D * h_theta_Lq_FTOF = new TH1D("theta_Lq_FTOF","#theta_{pq} Lead;#theta_{pq};Counts",180,0,180);
  hist_list_1.push_back(h_theta_Lq_FTOF);

  TH1D * h_phi_e_L = new TH1D("phi_e_L","|#phi_{e} - #phi_{p}|;|#phi_{e} - #phi_{p}|,Counts",180,0,180);
  hist_list_1.push_back(h_phi_e_L);
  TH1D * h_mmiss_FTOF = new TH1D("mmiss_FTOF","m_{miss};m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_FTOF);
  TH2D * h_mmiss_phi_e_L = new TH2D("mmiss_phi_e_L","m_{miss} vs. |#phi_{e} - #phi_{p}|;m_{miss};|#phi_{e} - #phi_{p};Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_phi_e_L);
  TH2D * h_xB_mmiss = new TH2D("xB_mmiss","x_{B} vs. m_{miss};x_{B};m_{miss};Counts",100,0,2,100,0.4,1.4);
  hist_list_2.push_back(h_xB_mmiss);
  TH2D * h_pmiss_mmiss = new TH2D("pmiss_mmiss","p_{miss} vs. m_{miss};p_{miss};m_{miss};Counts",100,0,1.5,100,0.4,1.4);
  hist_list_2.push_back(h_pmiss_mmiss);
  TH2D * h_xB_theta_1q = new TH2D("xB_theta_1q","x_{B} vs. #theta_{miss,q};x_{B};#theta_{miss,q};Counts",100,0,2,180,0,180);
  hist_list_2.push_back(h_xB_theta_1q);
  TH2D * h_Loq_theta_1q = new TH2D("Loq_theta_1q","|p|/|q| vs. #theta_{miss,q};|p|/|q|;#theta_{miss,q}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_Loq_theta_1q);
  TH2D * h_pmiss_theta_miss = new TH2D("pmiss_theta_miss","p_{miss} vs. #theta_{miss};p_{miss};#theta_{miss}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss);


  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  TH1D * h_pmiss = new TH1D("pmiss","p_{miss};p_{miss};Counts",100,0,1.5);
  hist_list_1.push_back(h_pmiss);
  TH1D * h_mmiss = new TH1D("mmiss","m_{miss};m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss);
  TH2D * h_pmiss_theta_miss_SRC = new TH2D("pmiss_theta_miss_SRC","p_{miss} vs. #theta_{miss};p_{miss};theta_{miss};Counts",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_xB_Loq_SRC = new TH2D("xB_Loq","x_{B} vs |p|/|q|;x_{B};|p|/|q|",100,0,2,100,0,1.5);
  hist_list_2.push_back(h_xB_Loq_SRC);

  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Sumw2();
    hist_list_1[i]->GetXaxis()->CenterTitle();
    hist_list_1[i]->GetYaxis()->CenterTitle();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Sumw2();
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->CenterTitle();
  }


  int counter = 0;

  //Define cut class
  eventcut myCut(Ebeam);
  myCut.setl_scint(TOFID);
  while(chain.Next()==true){
      //Display completed  
      counter++;
      if((counter%100000) == 0){
	//	printProgress(.5);
	cerr << counter <<" completed \n";
      }    
      // get particles by type
      auto electrons=c12->getByID(11);
      auto protons=c12->getByID(2212);
      auto neutrons=c12->getByID(2112);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);


      if(electrons.size()!=1){ continue;}
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double EoP_e =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
      int nphe = electrons[0]->che(HTCC)->getNphe();

  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////      
      h_Vcal_EoP->Fill(electrons[0]->cal(PCAL)->getLv(),EoP_e,weight);
      h_Wcal_EoP->Fill(electrons[0]->cal(PCAL)->getLw(),EoP_e,weight);
      h_phi_theta->Fill(p_e.Phi()*180/M_PI,p_e.Theta()*180/M_PI,weight);
      h_sector->Fill(electrons[0]->getSector(),weight);
	  
  /////////////////////////////////////
  //Electron Pid
  /////////////////////////////////////
      h_P_EoP->Fill(p_e.Mag(),EoP_e,weight);
      h_nphe->Fill(nphe,weight);


      if(!myCut.electroncut(c12)){continue;}      
      if(!isMC){
	if(electrons[0]->par()->getVz() < -5){continue;}
	if(electrons[0]->par()->getVz() > -1){continue;}
      }
  /////////////////////////////////////
  //Electron Kinematics  
  /////////////////////////////////////
      TVector3	p_q = p_b - p_e;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;
      //if(WSq>1.25){continue;}
      
      h_xB->Fill(xB,weight);
      h_QSq->Fill(QSq,weight);
      h_WSq->Fill(WSq,weight);
      h_xB_QSq->Fill(xB,QSq,weight);
      h_xB_WSq->Fill(xB,WSq,weight);
      h_QSq_WSq->Fill(QSq,WSq,weight);

  /////////////////////////////////////
  //All Proton Angles
  /////////////////////////////////////
      for(int j = 0; j < protons.size(); j++){
	TVector3 p_L;
	p_L.SetMagThetaPhi(protons[j]->getP(),protons[j]->getTheta(),protons[j]->getPhi());
	double theta_L = p_L.Theta() * 180 / M_PI;
	double theta_Lq = p_L.Angle(p_q) * 180 / M_PI;
	h_theta_L->Fill(theta_L,weight);
	h_theta_Lq->Fill(theta_Lq,weight);
      }

  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
      int index_L = myCut.leadnucleoncut(c12);
      if(index_L < 0){ continue; }
      TVector3 p_L;
      p_L.SetMagThetaPhi(protons[index_L]->getP(),protons[index_L]->getTheta(),protons[index_L]->getPhi());
      TVector3 p_1 = p_L - p_q;
      TVector3 p_miss = -p_1;
      double mmiss = get_mmiss(p_b,p_e,p_L);
      double phi_diff = get_phi_diff(p_e.Phi()*180/M_PI,p_L.Phi()*180/M_PI);
      double theta_L = p_L.Theta() * 180 / M_PI;
      double theta_miss = p_miss.Theta() * 180 / M_PI;
      double theta_Lq = p_L.Angle(p_q) * 180 / M_PI;
      double Loq = p_L.Mag() / p_q.Mag();
      double theta_1q = p_1.Angle(p_q) * 180 / M_PI;

      h_theta_L_FTOF->Fill(theta_L,weight);
      h_theta_Lq_FTOF->Fill(theta_Lq,weight);
      h_phi_e_L->Fill(phi_diff,weight);
      h_mmiss_FTOF->Fill(mmiss,weight);
      h_mmiss_phi_e_L->Fill(mmiss,phi_diff,weight);
      h_xB_mmiss->Fill(xB,mmiss,weight);
      h_pmiss_mmiss->Fill(p_miss.Mag(),mmiss,weight);
      h_xB_theta_1q->Fill(xB,theta_1q,weight);
      h_Loq_theta_1q->Fill(Loq,theta_1q,weight);
      h_pmiss_theta_miss->Fill(p_miss.Mag(),theta_miss,weight);

      
  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
      if(!myCut.leadSRCnucleoncut(c12,index_L)){continue;}      
      h_pmiss->Fill(p_miss.Mag(),weight);
      h_mmiss->Fill(mmiss,weight);
      h_pmiss_theta_miss_SRC->Fill(p_miss.Mag(),theta_miss,weight);
      h_xB_Loq_SRC->Fill(xB,Loq,weight);

  }
  cout<<counter<<endl;

  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }


  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
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

  /////////////////////////////////////
  //Electron fiducials and Pid
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Candidates:");
  text.DrawLatex(0.2,0.8,"No Cuts");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_Vcal_EoP->Draw("colz");
  myCanvas->cd(2);
  h_Wcal_EoP->Draw("colz");
  myCanvas->cd(3);
  h_phi_theta->Draw("colz");
  myCanvas->cd(4);
  h_sector->Draw("colz");
  myCanvas->cd(5);
  h_P_EoP->Draw("colz");
  myCanvas->cd(6);
  h_nphe->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  /////////////////////////////////////
  //Electron Kinematics  
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts:");
  text.DrawLatex(0.2,0.8,"V_{cal} and W_{cal} > 14 [cm]");
  text.DrawLatex(0.2,0.7,"0.18 < SF < 0.28");
  text.DrawLatex(0.2,0.6,"1 [GeV] < p_{e} < E_{beam}");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_xB->Draw();
  myCanvas->cd(2);
  h_QSq->Draw();
  myCanvas->cd(3);
  h_WSq->Draw();
  myCanvas->cd(4);
  h_xB_QSq->Draw("colz");
  myCanvas->cd(5);
  h_xB_WSq->Draw("colz");
  myCanvas->cd(6);
  h_QSq_WSq->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  /////////////////////////////////////
  //All Proton Angles
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Proton Detected");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_L->Draw("colz");
  myCanvas->cd(2);
  h_theta_Lq->Draw("colz");
  myCanvas->cd(3);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  char temp[100];
  sprintf(temp,"Scintillator = %d",TOFID);
  text.DrawLatex(0.2,0.7,temp);
  text.DrawLatex(0.2,0.6,"#theta_{p,q}<25^{o}");
  text.DrawLatex(0.2,0.5,"-3 < #chi^{2} PID<3 ");
  //text.DrawLatex(0.2,0.4,"#theta_{p} <50^{o}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_L_FTOF->Draw();
  myCanvas->cd(2);
  h_theta_Lq_FTOF->Draw();
  myCanvas->cd(3);
  h_phi_e_L->Draw();
  myCanvas->cd(4);
  h_mmiss_FTOF->Draw();
  myCanvas->cd(5);
  h_mmiss_phi_e_L->Draw("colz");
  myCanvas->cd(6);
  h_xB_mmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_mmiss->Draw("colz");
  myCanvas->cd(2);
  h_xB_theta_1q->Draw("colz");
  myCanvas->cd(3);
  h_Loq_theta_1q->Draw("colz");
  myCanvas->cd(4);
  h_pmiss_theta_miss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD}) Cuts");
  text.DrawLatex(0.2,0.7,"1.5 < Q^{2} [GeV]");
  text.DrawLatex(0.2,0.6,"0.3 [GeV] < p_{miss}");
  text.DrawLatex(0.2,0.5,"0.84 [GeV] < m_{mmiss} < 1.04 [GeV]");
  text.DrawLatex(0.2,0.4,"0.62 < |p|/|q| < 0.96");
  text.DrawLatex(0.2,0.3,"1.2 < x_{B}");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss->Draw();
  myCanvas->cd(2);
  h_mmiss->Draw();
  myCanvas->cd(3);
  h_pmiss_theta_miss_SRC->Draw("colz");
  myCanvas->cd(4);
  h_xB_Loq_SRC->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
}


void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

