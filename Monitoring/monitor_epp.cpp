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
#include "functions.h"

using namespace std;
using namespace clas12;


void printProgress(double percentage);

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/cutfile.txt> <path/to/input.hipo> \n";
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
  eventcut myCut(Ebeam,argv[5]);
  myCut.print_cuts();
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

  char temp_name[100];
  char temp_title[100];

  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////
  TH2D * h_phi_theta = new TH2D("phi_theta","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);
  hist_list_2.push_back(h_phi_theta);

  TH1D * h_sector = new TH1D("sector","ECAL Sector;Sector;Counts",6,1,7);
  hist_list_1.push_back(h_sector);

  TH2D * h_Vcal_EoP[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"Vcal_EoP_%d",i+1);
    sprintf(temp_title,"ECAL V coordinate vs. Sampling Fraction Sector=%d;ECAL V coordinate;Sampling Fraction",i+1);
    h_Vcal_EoP[i] = new TH2D(temp_name,temp_title,60,0,30,100,0.1,0.35);
    hist_list_2.push_back(h_Vcal_EoP[i]);
  }

  TH2D * h_Wcal_EoP[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"Wcal_EoP_%d",i+1);
    sprintf(temp_title,"ECAL W coordinate vs. Sampling Fraction Sector=%d;ECAL W coordinate;Sampling Fraction",i+1);
    h_Wcal_EoP[i] = new TH2D(temp_name,temp_title,60,0,30,100,0.1,0.35);
    hist_list_2.push_back(h_Wcal_EoP[i]);
  }

  /////////////////////////////////////
  //Electron Pid and Vertex
  /////////////////////////////////////
  TH1D * h_nphe = new TH1D("nphe","#Photo-electrons in HTCC;#Photo-electrons;Counts",100,0,50);
  hist_list_1.push_back(h_nphe);

  TH1D * h_vtz_e = new TH1D("vtz_e","Electron Z Vertex;vertex;Counts",100,-15,15);
  hist_list_1.push_back(h_vtz_e);

  TH2D * h_mom_EoP[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_EoP_%d",i+1);
    sprintf(temp_title,"p_{e} vs. Sampling Faction Sector=%d;Momentum;Theta",i+1);
    h_mom_EoP[i] = new TH2D(temp_name,temp_title,100,0,7,100,0.1,0.35);
    hist_list_2.push_back(h_mom_EoP[i]);
  }

  
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

  TH2D * h_mom_theta[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_theta_%d",i+1);
    sprintf(temp_title,"p_{e} vs. #theta_{e} Sector=%d;Momentum;Theta",i+1);
    h_mom_theta[i] = new TH2D(temp_name,temp_title,100,0,7,100,5,40);
    hist_list_2.push_back(h_mom_theta[i]);
  }

  /////////////////////////////////////
  //FTOF Proton Information
  /////////////////////////////////////
  TH1D * h_chiSq_p_FTOF = new TH1D("chiSq_p_FTOF","#chi^{2}_{p};#chi^{2}_{p};Counts",100,-5,5);
  hist_list_1.push_back(h_chiSq_p_FTOF);

  TH1D * h_vtz_p_FTOF = new TH1D("vtz_p_FTOF","Proton Z Vertex;vertex;Counts",100,-15,15);
  hist_list_1.push_back(h_vtz_p_FTOF);
  TH1D * h_vtz_ep_delta_FTOF = new TH1D("vtz_ep_delta_FTOF","#Delta Vertex;#Delta Vertex;Counts",100,-5,5);
  hist_list_1.push_back(h_vtz_ep_delta_FTOF);
  TH2D * h_vtz_e_vtz_p_FTOF = new TH2D("vtz_e_vtz_p_FTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-15,15,100,-15,15);
  hist_list_2.push_back(h_vtz_e_vtz_p_FTOF);


  TH1D * h_theta_p_FTOF = new TH1D("theta_p_FTOF","#theta_{proton};#theta_{proton};Counts",180,0,180);
  hist_list_1.push_back(h_theta_p_FTOF);
  TH1D * h_theta_pq_FTOF = new TH1D("theta_pq_FTOF","#theta_{pq};#theta_{pq};Counts",180,0,180);
  hist_list_1.push_back(h_theta_pq_FTOF);
  TH2D * h_phi_theta_p_FTOF = new TH2D("phi_theta_p_FTOF","#phi_{p} vs. #theta_{p} ;#phi_{p};#theta_{p}",100,-180,180,100,5,45);
  hist_list_2.push_back(h_phi_theta_p_FTOF);
  TH2D * h_mom_beta_p_FTOF = new TH2D("mom_beta_p_FTOF","p_{p} vs. #beta_{p} ;p_{p};#beta_{p}",100,0,4,100,0.7,1);
  hist_list_2.push_back(h_mom_beta_p_FTOF);


  TH2D * h_mom_theta_p_FTOF[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_theta_p_%d",i+1);
    sprintf(temp_title,"p_{p} vs. #theta_{p} Sector=%d;Momentum;Theta",i+1);
    h_mom_theta_p_FTOF[i] = new TH2D(temp_name,temp_title,100,0,4,100,5,50);
    hist_list_2.push_back(h_mom_theta_p_FTOF[i]);
  }

  /////////////////////////////////////
  //CTOF Proton Information
  /////////////////////////////////////
  TH1D * h_chiSq_p_CTOF = new TH1D("chiSq_p_CTOF","#chi^{2}_{p};#chi^{2}_{p};Counts",100,-5,5);
  hist_list_1.push_back(h_chiSq_p_CTOF);

  TH1D * h_vtz_p_CTOF = new TH1D("vtz_p_CTOF","Proton Z Vertex;vertex;Counts",100,-15,15);
  hist_list_1.push_back(h_vtz_p_CTOF);
  TH1D * h_vtz_ep_delta_CTOF = new TH1D("vtz_ep_delta_CTOF","#Delta Vertex;#Delta Vertex;Counts",100,-5,5);
  hist_list_1.push_back(h_vtz_ep_delta_CTOF);
  TH2D * h_vtz_e_vtz_p_CTOF = new TH2D("vtz_e_vtz_p_CTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-15,15,100,-15,15);
  hist_list_2.push_back(h_vtz_e_vtz_p_CTOF);


  TH1D * h_theta_p_CTOF = new TH1D("theta_p_CTOF","#theta_{proton};#theta_{proton};Counts",180,0,180);
  hist_list_1.push_back(h_theta_p_CTOF);
  TH1D * h_theta_pq_CTOF = new TH1D("theta_pq_CTOF","#theta_{pq};#theta_{pq};Counts",180,0,180);
  hist_list_1.push_back(h_theta_pq_CTOF);
  TH2D * h_phi_theta_p_CTOF = new TH2D("phi_theta_p_CTOF","#phi_{p} vs. #theta_{p} ;#phi_{p};#theta_{p}",100,-180,180,100,5,45);
  hist_list_2.push_back(h_phi_theta_p_CTOF);
  TH2D * h_mom_beta_p_CTOF = new TH2D("mom_beta_p_CTOF","p_{p} vs. #beta_{p} ;p_{p};#beta_{p}",100,0,4,100,0.7,1);
  hist_list_2.push_back(h_mom_beta_p_CTOF);
  TH2D * h_mom_theta_p_CTOF = new TH2D("mom_theta_p_CTOF","#p_{p} vs. #theta_{p} ;#p_{p};#theta_{p}",100,0,4,100,30,135);
  hist_list_2.push_back(h_mom_theta_p_CTOF);


  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
  TH1D * h_theta_p_Lead = new TH1D("theta_p_Lead","#theta_{p,Lead};#theta_{p,Lead};Counts",180,0,180);
  hist_list_1.push_back(h_theta_p_Lead);
  TH1D * h_theta_pq_Lead = new TH1D("theta_pq_Lead","#theta_{pq} Lead;#theta_{pq};Counts",180,0,90);
  hist_list_1.push_back(h_theta_pq_Lead);
  TH2D * h_mom_theta_p_Lead = new TH2D("mom_theta_p_Lead","#p_{p,Lead} vs. #theta_{p,Lead} ;#p_{p,Lead};#theta_{p,Lead}",100,0,4,100,0,135);
  hist_list_2.push_back(h_mom_theta_p_Lead);
  TH1D * h_phi_e_p_Lead = new TH1D("phi_e_p_Lead","|#phi_{e} - #phi_{p,Lead}|;|#phi_{e} - #phi_{p,Lead}|,Counts",100,120,180);
  hist_list_1.push_back(h_phi_e_p_Lead);
  TH1D * h_xB_Lead = new TH1D("xB_Lead","x_{B} Lead;x_{B};Counts",100,0.0,2.0);
  hist_list_1.push_back(h_xB_Lead);
  TH2D * h_vtz_e_vtz_p_Lead = new TH2D("vtz_e_vtz_p_Lead","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-15,15,100,-15,15);
  hist_list_2.push_back(h_vtz_e_vtz_p_Lead);


  TH1D * h_pmiss_Lead = new TH1D("pmiss_Lead","p_{miss} Lead;p_{miss};Counts",100,0,1.5);
  hist_list_1.push_back(h_pmiss_Lead);
  TH2D * h_pmiss_thetamiss_Lead = new TH2D("pmiss_thetamiss_Lead","p_{miss} vs. #theta_{miss} Lead;p_{miss};#theta_{miss}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_thetamiss_Lead);
  TH2D * h_xB_theta_1q_Lead = new TH2D("xB_theta_1q_Lead","x_{B} vs. #theta_{miss,q} Lead;x_{B};#theta_{miss,q};Counts",100,0,2,180,0,180);
  hist_list_2.push_back(h_xB_theta_1q_Lead);
  TH2D * h_Loq_theta_1q_Lead = new TH2D("Loq_theta_1q_Lead","|p|/|q| vs. #theta_{miss,q} Lead;|p|/|q|;#theta_{miss,q}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_Loq_theta_1q_Lead);


  TH1D * h_mmiss_Lead = new TH1D("mmiss_Lead","m_{miss} Lead;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_Lead);
  TH2D * h_mmiss_phi_e_p_Lead = new TH2D("mmiss_phi_e_p_Lead","m_{miss} vs. |#phi_{e} - #phi_{p}| Lead;m_{miss};|#phi_{e} - #phi_{p};Counts",100,0.4,1.4,100,120,180);
  hist_list_2.push_back(h_mmiss_phi_e_p_Lead);
  TH2D * h_mmiss_xB_Lead = new TH2D("mmiss_xB_Lead","m_{miss} vs. x_{B} Lead;m_{miss};x_{B};Counts",100,0.4,1.4,100,0.0,2.0);
  hist_list_2.push_back(h_mmiss_xB_Lead);
  TH2D * h_mmiss_pmiss_Lead = new TH2D("mmiss_pmiss_Lead","m_{miss} vs. p_{miss} Lead;m_{miss};p_{miss};Counts",100,0.4,1.4,100,0.0,1.5);
  hist_list_2.push_back(h_mmiss_pmiss_Lead);
  TH2D * h_mmiss_theta_1q_Lead = new TH2D("mmiss_theta_1q_Lead","m_{miss} vs. #theta_{miss,q} Lead;m_{miss};#theta_{miss,q};Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_theta_1q_Lead);
  TH2D * h_mmiss_theta_p_Lead = new TH2D("mmiss_theta_p_Lead","m_{miss} vs. #theta_{p,Lead} Lead;m_{miss};#theta_{p,Lead};Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_theta_p_Lead);
  TH2D * h_mmiss_mom_p_Lead = new TH2D("mmiss_mom_p_Lead","m_{miss} vs. p_{p,Lead} Lead;m_{miss};p_{p,Lead};Counts",100,0.4,1.4,100,0,4);
  hist_list_2.push_back(h_mmiss_mom_p_Lead);
  TH2D * h_mmiss_momT_p_Lead = new TH2D("mmiss_momT_p_Lead","m_{miss} vs. p_{p,T,Lead} Lead;m_{miss};p_{p,T,Lead};Counts",100,0.4,1.4,100,0,2.5);
  hist_list_2.push_back(h_mmiss_momT_p_Lead);


  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  TH1D * h_xB_SRC = new TH1D("xB_SRC","x_{B} SRC;x_{B};Counts",100,1.0,2.0);
  hist_list_1.push_back(h_xB_SRC);
  TH1D * h_pmiss_SRC = new TH1D("pmiss_SRC","p_{miss} SRC;p_{miss};Counts",100,0,1.5);
  hist_list_1.push_back(h_pmiss_SRC);
  TH1D * h_mmiss_SRC = new TH1D("mmiss_SRC","m_{miss} SRC;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_SRC);
  TH2D * h_pmiss_theta_miss_SRC = new TH2D("pmiss_theta_miss_SRC","p_{miss} vs. #theta_{miss};p_{miss};#theta_{miss};Counts",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_pmiss_theta_L_SRC = new TH2D("pmiss_theta_L_SRC","p_{miss} vs. #theta_{L};p_{miss};#theta_{L};Counts",100,0,1.5,100,5,45);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_xB_Loq_SRC = new TH2D("xB_Loq","x_{B} vs |p|/|q|;x_{B};|p|/|q|",100,0,2,100,0,1.5);
  hist_list_2.push_back(h_xB_Loq_SRC);

  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
  TH1D * h_p_2_AllRec = new TH1D("p_2_AllRec","p All Recoils;p_2",100,0,1.5);
  hist_list_1.push_back(h_p_2_AllRec);
  TH1D * h_chiSq_rec_AllRec = new TH1D("chiSq_rec_AllRec","#chi^{2}_{rec} All Recoils;#chi^{2}_{rec}",100,-5,5);
  hist_list_1.push_back(h_chiSq_rec_AllRec);
  TH2D * h_mom_beta_rec_AllRec = new TH2D("mom_beta_rec_AllRec","p_{rec} vs. #beta_{rec} ;p_{rec};#beta_{rec}",100,0,4,100,0.7,1);
  hist_list_2.push_back(h_mom_beta_rec_AllRec);
  TH1D * h_count_AllRec = new TH1D("count_AllRec","Number of Recoils;Multiplicity",5,0,5);
  hist_list_1.push_back(h_count_AllRec);

  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////
  TH1D * h_p_2_Rec = new TH1D("p_2_Rec","p_{rec};p_{rec};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_2_Rec);
  TH1D * h_p_rel_Rec = new TH1D("p_rel_Rec","p_{rel};p_{rel};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_rel_Rec);
  TH1D * h_p_cm_Rec = new TH1D("p_cm_Rec","p_{C.M.};p_{C.M.};Counts",100,0,0.5);
  hist_list_1.push_back(h_p_cm_Rec);
  TH1D * h_p_t_cm_Rec = new TH1D("p_t_cm_Rec","p_{t,C.M.};p_{t,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_t_cm_Rec);
  TH1D * h_p_y_cm_Rec = new TH1D("p_y_cm_Rec","p_{y,C.M.};p_{y,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_y_cm_Rec);
  TH1D * h_p_x_cm_Rec = new TH1D("p_x_cm_Rec","p_{x,C.M.};p_{x,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_x_cm_Rec);
  TH1D * h_theta_rel_Rec = new TH1D("theta_rel_Rec","#theta_{rel};#theta_{rel};Counts",180,0,180);
  hist_list_1.push_back(h_theta_rel_Rec);
  TH2D * h_p_cm_theta_rel_Rec = new TH2D("p_cm_theta_rel_Rec","p_{C.M.} vs. #theta_{rel};p_{C.M.};#theta_{rel}",100,0,0.5,180,100,180);
  hist_list_2.push_back(h_p_cm_theta_rel_Rec);


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
  while(chain.Next()==true){
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
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
      int esector = electrons[0]->getSector();
      if(isMC){esector = (p_e.Phi()*3/M_PI)+4;}
      double vtz_e = electrons[0]->par()->getVz();

  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////      
      h_phi_theta->Fill(p_e.Phi()*180/M_PI,p_e.Theta()*180/M_PI,weight);
      h_sector->Fill(esector,weight);

      h_Vcal_EoP[esector-1]->Fill(electrons[0]->cal(PCAL)->getLv(),EoP_e,weight);
      h_Wcal_EoP[esector-1]->Fill(electrons[0]->cal(PCAL)->getLw(),EoP_e,weight);


  /////////////////////////////////////
  //Electron Pid
  /////////////////////////////////////
      h_mom_EoP[esector-1]->Fill(p_e.Mag(),EoP_e,weight);
      h_nphe->Fill(nphe,weight);
      h_vtz_e->Fill(vtz_e,weight);


      if(!myCut.electroncut(c12)){continue;}      
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
      h_mom_theta[esector-1]->Fill(p_e.Mag(),p_e.Theta()*180/M_PI,weight);
  /////////////////////////////////////
  //FTOF and CTOF Proton Information
  /////////////////////////////////////
      for(int j = 0; j < protons.size(); j++){
	TVector3 p_p;
	p_p.SetMagThetaPhi(protons[j]->getP(),protons[j]->getTheta(),protons[j]->getPhi());
	double theta_p = p_p.Theta() * 180 / M_PI;
	double phi_p = p_p.Phi() * 180 / M_PI;
	double theta_pq = p_p.Angle(p_q) * 180 / M_PI;
	double beta_p = protons[j]->par()->getBeta();
	double Chi2Pid_p = protons[j]->par()->getChi2Pid();
	double vtz_p = protons[j]->par()->getVz();
	double vtz_ep_delta = vtz_e - vtz_p;

	bool FTOF1A = (protons[j]->sci(clas12::FTOF1A)->getDetector() == 12);
	bool FTOF1B = (protons[j]->sci(clas12::FTOF1B)->getDetector() == 12);
	bool FTOF2 = (protons[j]->sci(clas12::FTOF2)->getDetector() == 12);
	bool CTOF = (protons[j]->sci(clas12::CTOF)->getDetector() == 4);

	if(FTOF1A || FTOF1B || FTOF2){
	  h_chiSq_p_FTOF->Fill(Chi2Pid_p,weight);
	  h_vtz_p_FTOF->Fill(vtz_p,weight);
	  h_vtz_ep_delta_FTOF->Fill(vtz_ep_delta,weight);
	  h_vtz_e_vtz_p_FTOF->Fill(vtz_e,vtz_p,weight);

	  h_theta_p_FTOF->Fill(theta_p,weight);
	  h_theta_pq_FTOF->Fill(theta_pq,weight);	
	  h_phi_theta_p_FTOF->Fill(phi_p,theta_p,weight);
	  h_mom_beta_p_FTOF->Fill(p_p.Mag(),beta_p,weight);
	  
	  h_mom_theta_p_FTOF[protons[j]->getSector()-1]->Fill(p_p.Mag(),theta_p,weight);
	}

	if(CTOF){
	  h_chiSq_p_CTOF->Fill(Chi2Pid_p,weight);
	  h_vtz_p_CTOF->Fill(vtz_p,weight);
	  h_vtz_ep_delta_CTOF->Fill(vtz_ep_delta,weight);
	  h_vtz_e_vtz_p_CTOF->Fill(vtz_e,vtz_p,weight);

	  h_theta_p_CTOF->Fill(theta_p,weight);
	  h_theta_pq_CTOF->Fill(theta_pq,weight);	
	  h_phi_theta_p_CTOF->Fill(phi_p,theta_p,weight);
	  h_mom_beta_p_CTOF->Fill(p_p.Mag(),beta_p,weight);
	  
	  h_mom_theta_p_CTOF->Fill(p_p.Mag(),theta_p,weight);
	}
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
      double phi_diff = get_phi_diff(p_e,p_L);
      double theta_L = p_L.Theta() * 180 / M_PI;
      double phi_L = p_L.Phi() * 180 / M_PI;
      double theta_miss = p_miss.Theta() * 180 / M_PI;
      double theta_Lq = p_L.Angle(p_q) * 180 / M_PI;
      double Loq = p_L.Mag() / p_q.Mag();
      double theta_1q = p_1.Angle(p_q) * 180 / M_PI;
      double vtz_p = protons[index_L]->par()->getVz();

      h_theta_p_Lead->Fill(theta_L,weight);
      h_theta_pq_Lead->Fill(theta_Lq,weight);
      h_mom_theta_p_Lead->Fill(p_L.Mag(),theta_L,weight);
      h_phi_e_p_Lead->Fill(phi_diff,weight);
      h_xB_Lead->Fill(xB,weight);
      h_vtz_e_vtz_p_Lead->Fill(vtz_e,vtz_p,weight);
      
      h_pmiss_Lead->Fill(p_miss.Mag(),weight);
      h_pmiss_thetamiss_Lead->Fill(p_miss.Mag(),theta_miss,weight);
      h_xB_theta_1q_Lead->Fill(xB,theta_1q,weight);
      h_Loq_theta_1q_Lead->Fill(Loq,theta_1q,weight);
      
      h_mmiss_Lead->Fill(mmiss,weight);
      h_mmiss_phi_e_p_Lead->Fill(mmiss,phi_diff,weight);
      h_mmiss_xB_Lead->Fill(mmiss,xB,weight);
      h_mmiss_pmiss_Lead->Fill(mmiss,p_miss.Mag(),weight);
      h_mmiss_theta_1q_Lead->Fill(mmiss,theta_1q,weight);
      h_mmiss_theta_p_Lead->Fill(mmiss,theta_L,weight);
      h_mmiss_mom_p_Lead->Fill(mmiss,p_L.Mag(),weight);
      h_mmiss_momT_p_Lead->Fill(mmiss,p_L.Perp(),weight);
        
  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
      if(!myCut.leadSRCnucleoncut(c12,index_L)){continue;}

      h_xB_SRC->Fill(xB,weight);
      h_pmiss_SRC->Fill(p_miss.Mag(),weight);
      h_mmiss_SRC->Fill(mmiss,weight);
      h_pmiss_theta_miss_SRC->Fill(p_miss.Mag(),theta_miss,weight);
      h_pmiss_theta_L_SRC->Fill(p_miss.Mag(),theta_L,weight);
      h_xB_Loq_SRC->Fill(xB,Loq,weight);


  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
      for(int j = 0; j < protons.size(); j++){
	if(j==index_L){continue;}
	h_p_2_AllRec->Fill(protons[j]->getP(),weight);
	h_chiSq_rec_AllRec->Fill(protons[j]->par()->getChi2Pid(),weight);
	h_mom_beta_rec_AllRec->Fill(protons[j]->getP(),protons[j]->par()->getBeta(),weight);
      }
      h_count_AllRec->Fill(protons.size()-1,weight);

  /////////////////////////////////////
  //Recoil SRC Proton Checks
  /////////////////////////////////////
      int index_R = myCut.recoilSRCnucleoncut(c12,index_L);
      if(index_R < 0){ continue; }
      TVector3 p_2;
      p_2.SetMagThetaPhi(protons[index_R]->getP(),protons[index_R]->getTheta(),protons[index_R]->getPhi());
      TVector3 p_rel = p_1-p_2;
      p_rel.SetMag(p_rel.Mag()/2);
      TVector3 p_cm = p_1+p_2;
      double theta_rel = p_1.Angle(p_2) * 180 / M_PI;
      
      //Create new reference frame
      TVector3 vt = p_2.Unit();
      TVector3 vy = p_2.Cross(p_q).Unit();
      TVector3 vx = vt.Cross(vy);
      
      
      h_p_2_Rec->Fill(protons[index_R]->getP(),weight);
      h_p_rel_Rec->Fill(p_rel.Mag(),weight);
      h_p_cm_Rec->Fill(p_cm.Mag(),weight);
      h_p_t_cm_Rec->Fill(p_cm.Dot(vt),weight);
      h_p_y_cm_Rec->Fill(p_cm.Dot(vy),weight);
      h_p_x_cm_Rec->Fill(p_cm.Dot(vx),weight);
      h_theta_rel_Rec->Fill(theta_rel,weight);
      h_p_cm_theta_rel_Rec->Fill(p_cm.Mag(),theta_rel,weight);
      
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
  h_phi_theta->Draw("colz");
  myCanvas->cd(2);
  h_sector->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_Vcal_EoP[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_Wcal_EoP[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_nphe->Draw();
  myCanvas->cd(2);
  h_vtz_e->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_EoP[i]->Draw("colz");  
  }
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
  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  /////////////////////////////////////
  //FTOF Proton Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Proton Detected in FTOF");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_chiSq_p_FTOF->Draw();
  myCanvas->cd(2);
  h_vtz_p_FTOF->Draw();
  myCanvas->cd(3);
  h_vtz_ep_delta_FTOF->Draw();
  myCanvas->cd(4);
  h_vtz_e_vtz_p_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_FTOF->Draw();
  myCanvas->cd(2);
  h_theta_pq_FTOF->Draw();
  myCanvas->cd(3);
  h_phi_theta_p_FTOF->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_p_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_p_FTOF[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //CTOF Proton Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Proton Detected in CTOF");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_chiSq_p_CTOF->Draw();
  myCanvas->cd(2);
  h_vtz_p_CTOF->Draw();
  myCanvas->cd(3);
  h_vtz_ep_delta_CTOF->Draw();
  myCanvas->cd(4);
  h_vtz_e_vtz_p_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_CTOF->Draw();
  myCanvas->cd(2);
  h_theta_pq_CTOF->Draw();
  myCanvas->cd(3);
  h_phi_theta_p_CTOF->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_p_CTOF->Draw("colz");
  myCanvas->cd(5);
  h_mom_theta_p_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  char temp[100];
  sprintf(temp,"Scintillator = %d",0);
  text.DrawLatex(0.2,0.7,temp);
  text.DrawLatex(0.2,0.6,"#theta_{p,q}<25^{o}");
  text.DrawLatex(0.2,0.5,"-3 < #chi^{2} PID<3 ");
  //text.DrawLatex(0.2,0.4,"#theta_{p} <50^{o}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_Lead->Draw();
  myCanvas->cd(2);
  h_theta_pq_Lead->Draw();
  myCanvas->cd(3);
  h_mom_theta_p_Lead->Draw("colz");
  myCanvas->cd(4);
  h_phi_e_p_Lead->Draw();
  myCanvas->cd(5);
  h_xB_Lead->Draw();
  myCanvas->cd(6);
  h_vtz_e_vtz_p_Lead->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_Lead->Draw();
  myCanvas->cd(2);
  h_pmiss_thetamiss_Lead->Draw("colz");
  myCanvas->cd(3);
  h_xB_theta_1q_Lead->Draw("colz");
  myCanvas->cd(4);
  h_Loq_theta_1q_Lead->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_Lead->Draw();
  myCanvas->cd(2);
  h_mmiss_phi_e_p_Lead->Draw("colz");
  myCanvas->cd(3);
  h_mmiss_xB_Lead->Draw("colz");
  myCanvas->cd(4);
  h_mmiss_pmiss_Lead->Draw("colz");
  myCanvas->cd(5);
  h_mmiss_theta_1q_Lead->Draw("colz");
  myCanvas->cd(6);
  h_mmiss_theta_p_Lead->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_mom_p_Lead->Draw("colz");
  myCanvas->cd(2);
  h_mmiss_momT_p_Lead->Draw("colz");
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
  
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_xB_SRC->Draw();
  myCanvas->cd(2);
  h_pmiss_SRC->Draw();
  myCanvas->cd(3);
  h_mmiss_SRC->Draw();
  myCanvas->cd(4);
  h_pmiss_theta_miss_SRC->Draw("colz");
  myCanvas->cd(5);
  h_pmiss_theta_L_SRC->Draw("colz");
  myCanvas->cd(6);
  h_xB_Loq_SRC->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////

  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC}) Cuts");
  text.DrawLatex(0.2,0.7,"Second Proton Detected");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_p_2_AllRec->Draw();
  myCanvas->cd(2);
  h_chiSq_rec_AllRec->Draw();
  myCanvas->cd(3);
  h_mom_beta_rec_AllRec->Draw("colz");
  myCanvas->cd(4);
  h_count_AllRec->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////

  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC},p_{Rec}) Cuts");
  text.DrawLatex(0.2,0.7,"0.35 [GeV] < p_{Rec}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_p_2_Rec->Draw();
  myCanvas->cd(2);
  h_p_rel_Rec->Draw();
  myCanvas->cd(3);
  h_p_cm_Rec->Draw();
  myCanvas->cd(4);
  h_p_t_cm_Rec->Draw();
  myCanvas->cd(5);
  h_p_y_cm_Rec->Draw();
  myCanvas->cd(6);
  h_p_x_cm_Rec->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_rel_Rec->Draw();
  myCanvas->cd(2);
  h_p_cm_theta_rel_Rec->Draw("colz");
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

