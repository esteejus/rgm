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

const double c = 29.9792458;

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
  chain.db()->turnOffQADB();
  

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
  TH1D * h_nphe = new TH1D("nphe","#Photo-electrons in HTCC;#Photo-electrons;Counts",40,0,40);
  hist_list_1.push_back(h_nphe);

  TH1D * h_vtz_e = new TH1D("vtz_e","Electron Z Vertex;vertex;Counts",100,-10,10);
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
  TH1D * h_W = new TH1D("W","W;W",100,0.5,2);
  hist_list_1.push_back(h_W);
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
  //Electron Kinematics with SRC kinematics
  /////////////////////////////////////
  TH2D * h_mom_theta_eSRC = new TH2D("mom_theta_eSRC","p_{e} vs. #theta_{e} ;p_{e};#theta_{e}",100,0,7,100,5,40);
  hist_list_2.push_back(h_mom_theta_eSRC);  
  TH2D * h_xB_theta_eSRC = new TH2D("xB_theta_eSRC","x_{B} vs. #theta_{e} ;x_{B};#theta_{e}",100,1,2,100,5,40);
  hist_list_2.push_back(h_xB_theta_eSRC);  
  TH2D * h_thetaq_theta_eSRC = new TH2D("thetaq_theta_eSRC","#theta_{q} vs. #theta_{e} ;#theta_{q};#theta_{e}",100,0,90,100,5,40);
  hist_list_2.push_back(h_thetaq_theta_eSRC);  
  TH2D * h_thetaq_xB_eSRC = new TH2D("thetaq_xB_eSRC","#theta_{q} vs. x_{B} ;#theta_{q};x_{B}",100,0,90,100,1,2);
  hist_list_2.push_back(h_thetaq_xB_eSRC);  
  TH1D * h_thetaq_eSRC = new TH1D("thetaq_eSRC","#theta_{q};#theta_{q}",100,0,90);
  hist_list_1.push_back(h_thetaq_eSRC);  
  

  /////////////////////////////////////
  //Hadron Information
  /////////////////////////////////////
  TH2D * h_mom_beta_hadplus_FTOF = new TH2D("mom_beta_hadplus_FTOF","p vs. #beta FTOF (All +hadrons);p;#beta",100,0,4,100,0.3,1.2);
  hist_list_2.push_back(h_mom_beta_hadplus_FTOF);
  TH2D * h_mom_beta_hadplus_CTOF = new TH2D("mom_beta_hadplus_CTOF","p vs. #beta CTOF (All +hadrons);p;#beta",100,0,4,100,0.0,5);
  hist_list_2.push_back(h_mom_beta_hadplus_CTOF);
  TH2D * h_mom_beta_hadplus_CTOF_zoom = new TH2D("mom_beta_hadplus_CTOF_zoom","p vs. #beta CTOF (All +hadrons);p;#beta",100,0,4,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_hadplus_CTOF_zoom);

  TH2D * h_mom_beta_proton_CTOF_zoom = new TH2D("mom_beta_proton_CTOF","p vs. #beta CTOF (Protons);p;#beta",100,0,4,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_proton_CTOF_zoom);
  TH2D * h_mom_beta_pion_CTOF_zoom = new TH2D("mom_beta_pion_CTOF","p vs. #beta CTOF (Pions);p;#beta",100,0,4,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_pion_CTOF_zoom);
  TH2D * h_mom_beta_kaon_CTOF_zoom = new TH2D("mom_beta_kaon_CTOF","p vs. #beta CTOF (Kaons);p;#beta",100,0,4,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_kaon_CTOF_zoom);
  TH2D * h_mom_beta_deuteron_CTOF_zoom = new TH2D("mom_beta_deuteron_CTOF","p vs. #beta CTOF (Deuterons);p;#beta",100,0,4,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_deuteron_CTOF_zoom);
  TH2D * h_mom_beta_0_CTOF_zoom = new TH2D("mom_beta_0_CTOF","p vs. #beta CTOF (PID=0);p;#beta",100,0,4,100,-1.0,1.2);
  hist_list_2.push_back(h_mom_beta_0_CTOF_zoom);

  /////////////////////////////////////
  //ECAL Neutron Information
  /////////////////////////////////////
  TH1D * h_chiSq_p_FTOF = new TH1D("chiSq_p_FTOF","#chi^{2}_{PID};#chi^{2}_{PID};Counts",100,-10,10);
  hist_list_1.push_back(h_chiSq_p_FTOF);

  TH1D * h_vtz_p_FTOF = new TH1D("vtz_p_FTOF","Proton Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_p_FTOF);
  TH1D * h_vtz_ep_delta_FTOF = new TH1D("vtz_ep_delta_FTOF","#Delta Vertex;#Delta Vertex;Counts",100,-5,5);
  hist_list_1.push_back(h_vtz_ep_delta_FTOF);
  TH2D * h_vtz_e_vtz_p_FTOF = new TH2D("vtz_e_vtz_p_FTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_vtz_e_vtz_p_FTOF);
  TH1D * h_phi_e_p_FTOF = new TH1D("phi_e_p_FTOF","|#phi_{e} - #phi_{p}|;|#phi_{e} - #phi_{p}|,Counts",100,100,180);
  hist_list_1.push_back(h_phi_e_p_FTOF);


  TH1D * h_theta_p_FTOF = new TH1D("theta_p_FTOF","#theta_{proton};#theta_{proton};Counts",100,0,60);
  hist_list_1.push_back(h_theta_p_FTOF);
  TH1D * h_theta_pq_FTOF = new TH1D("theta_pq_FTOF","#theta_{pq};#theta_{pq};Counts",100,0,80);
  hist_list_1.push_back(h_theta_pq_FTOF);
  TH2D * h_phi_theta_p_FTOF = new TH2D("phi_theta_p_FTOF","#phi_{p} vs. #theta_{p} ;#phi_{p};#theta_{p}",100,-180,180,100,5,45);
  hist_list_2.push_back(h_phi_theta_p_FTOF);
  TH2D * h_mom_beta_p_FTOF = new TH2D("mom_beta_p_FTOF","p_{p} vs. #beta_{p} ;p_{p};#beta_{p}",100,0,4,100,0.3,1);
  hist_list_2.push_back(h_mom_beta_p_FTOF);
  TH1D * h_timediff_p_FTOF = new TH1D("timediff_p_FTOF","ToF-ToF_{|p|} ;ToF-ToF_{|p|};Counts",100,-2,2);
  hist_list_1.push_back(h_timediff_p_FTOF);


  TH2D * h_mom_theta_p_FTOF[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_theta_p_%d",i+1);
    sprintf(temp_title,"p_{p} vs. #theta_{p} Sector=%d;Momentum;Theta",i+1);
    h_mom_theta_p_FTOF[i] = new TH2D(temp_name,temp_title,100,0,4,100,5,50);
    hist_list_2.push_back(h_mom_theta_p_FTOF[i]);
  }

  /////////////////////////////////////
  //CND Proton Information
  /////////////////////////////////////
  TH1D * h_chiSq_p_CTOF = new TH1D("chiSq_p_CTOF","#chi^{2}_{PID};#chi^{2}_{PID};Counts",100,-10,10);
  hist_list_1.push_back(h_chiSq_p_CTOF);

  TH1D * h_vtz_p_CTOF = new TH1D("vtz_p_CTOF","Proton Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_p_CTOF);
  TH1D * h_vtz_ep_delta_CTOF = new TH1D("vtz_ep_delta_CTOF","#Delta Vertex;#Delta Vertex;Counts",100,-5,5);
  hist_list_1.push_back(h_vtz_ep_delta_CTOF);
  TH2D * h_vtz_e_vtz_p_CTOF = new TH2D("vtz_e_vtz_p_CTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_vtz_e_vtz_p_CTOF);
  TH1D * h_phi_e_p_CTOF = new TH1D("phi_e_p_CTOF","|#phi_{e} - #phi_{p}|;|#phi_{e} - #phi_{p}|,Counts",100,100,180);
  hist_list_1.push_back(h_phi_e_p_CTOF);


  TH1D * h_theta_p_CTOF = new TH1D("theta_p_CTOF","#theta_{proton};#theta_{proton};Counts",100,20,140);
  hist_list_1.push_back(h_theta_p_CTOF);
  TH1D * h_theta_pq_CTOF = new TH1D("theta_pq_CTOF","#theta_{pq};#theta_{pq};Counts",100,0,80);
  hist_list_1.push_back(h_theta_pq_CTOF);
  TH2D * h_phi_theta_p_CTOF = new TH2D("phi_theta_p_CTOF","#phi_{p} vs. #theta_{p} ;#phi_{p};#theta_{p}",100,-180,180,100,20,140);
  hist_list_2.push_back(h_phi_theta_p_CTOF);
  TH2D * h_mom_beta_p_CTOF = new TH2D("mom_beta_p_CTOF","p_{p} vs. #beta_{p} ;p_{p};#beta_{p}",100,0,4,100,0.3,1);
  hist_list_2.push_back(h_mom_beta_p_CTOF);
  TH1D * h_timediff_p_CTOF = new TH1D("timediff_p_CTOF","ToF-ToF_{|p|} ;ToF-ToF_{|p|};Counts",100,-2,2);
  hist_list_1.push_back(h_timediff_p_CTOF);
  TH2D * h_mom_theta_p_CTOF = new TH2D("mom_theta_p_CTOF","p_{p} vs. #theta_{p} ;p_{p};#theta_{p}",100,0,4,100,30,135);
  hist_list_2.push_back(h_mom_theta_p_CTOF);

/*
  /////////////////////////////////////
  //Lead Proton FTOF  Checks
  /////////////////////////////////////
  TH1D * h_theta_p_Lead_FTOF = new TH1D("theta_p_Lead_FTOF","#theta_{p,Lead};#theta_{p,Lead};Counts",100,0,60);
  hist_list_1.push_back(h_theta_p_Lead_FTOF);
  TH1D * h_theta_pq_Lead_FTOF = new TH1D("theta_pq_Lead_FTOF","#theta_{pq};#theta_{pq};Counts",100,0,30);
  hist_list_1.push_back(h_theta_pq_Lead_FTOF);
  TH2D * h_mom_theta_p_Lead_FTOF = new TH2D("mom_theta_p_Lead_FTOF","p_{p,Lead} vs. #theta_{p,Lead} ;p_{p,Lead};#theta_{p,Lead}",100,0,4,100,0,60);
  hist_list_2.push_back(h_mom_theta_p_Lead_FTOF);
  TH1D * h_phi_e_p_Lead_FTOF = new TH1D("phi_e_p_Lead_FTOF","|#phi_{e} - #phi_{p,Lead}|;|#phi_{e} - #phi_{p,Lead}|,Counts",100,140,180);
  hist_list_1.push_back(h_phi_e_p_Lead_FTOF);
  TH1D * h_xB_Lead_FTOF = new TH1D("xB_Lead_FTOF","x_{B} Lead;x_{B};Counts",100,0.0,2.0);
  hist_list_1.push_back(h_xB_Lead_FTOF);
  TH2D * h_vtz_e_vtz_p_Lead_FTOF = new TH2D("vtz_e_vtz_p_Lead_FTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_vtz_e_vtz_p_Lead_FTOF);


  TH1D * h_pmiss_Lead_FTOF = new TH1D("pmiss_Lead_FTOF","p_{miss};p_{miss};Counts",100,0,2.0);
  hist_list_1.push_back(h_pmiss_Lead_FTOF);
  TH2D * h_pmiss_thetamiss_Lead_FTOF = new TH2D("pmiss_thetamiss_Lead_FTOF","p_{miss} vs. #theta_{miss};p_{miss};#theta_{miss}",100,0,2.0,180,0,180);
  hist_list_2.push_back(h_pmiss_thetamiss_Lead_FTOF);
  TH2D * h_xB_theta_1q_Lead_FTOF = new TH2D("xB_theta_1q_Lead_FTOF","x_{B} vs. #theta_{miss,q};x_{B};#theta_{miss,q};Counts",100,0,2,180,0,180);
  hist_list_2.push_back(h_xB_theta_1q_Lead_FTOF);
  TH2D * h_Loq_theta_1q_Lead_FTOF = new TH2D("Loq_theta_1q_Lead_FTOF","|p|/|q| vs. #theta_{miss,q};|p|/|q|;#theta_{miss,q}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_Loq_theta_1q_Lead_FTOF);


  TH1D * h_mmiss_Lead_FTOF = new TH1D("mmiss_Lead_FTOF","m_{miss} Lead;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_Lead_FTOF);
  TH2D * h_mmiss_phi_e_p_Lead_FTOF = new TH2D("mmiss_phi_e_p_Lead_FTOF","m_{miss} vs. |#phi_{e} - #phi_{p}|;m_{miss};|#phi_{e} - #phi_{p};Counts",100,0.4,1.4,100,140,180);
  hist_list_2.push_back(h_mmiss_phi_e_p_Lead_FTOF);
  TH2D * h_mmiss_xB_Lead_FTOF = new TH2D("mmiss_xB_Lead_FTOF","m_{miss} vs. x_{B};m_{miss};x_{B};Counts",100,0.4,1.4,100,0.0,2.0);
  hist_list_2.push_back(h_mmiss_xB_Lead_FTOF);
  TH2D * h_mmiss_pmiss_Lead_FTOF = new TH2D("mmiss_pmiss_Lead_FTOF","m_{miss} vs. p_{miss};m_{miss};p_{miss};Counts",100,0.4,1.4,100,0.0,2.0);
  hist_list_2.push_back(h_mmiss_pmiss_Lead_FTOF);
  TH2D * h_mmiss_theta_1q_Lead_FTOF = new TH2D("mmiss_theta_1q_Lead_FTOF","m_{miss} vs. #theta_{miss,q};m_{miss};#theta_{miss,q};Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_theta_1q_Lead_FTOF);
  TH2D * h_mmiss_theta_p_Lead_FTOF = new TH2D("mmiss_theta_p_Lead_FTOF","m_{miss} vs. #theta_{p,Lead};m_{miss};#theta_{p,Lead};Counts",100,0.4,1.4,100,0,60);
  hist_list_2.push_back(h_mmiss_theta_p_Lead_FTOF);
  TH2D * h_mmiss_mom_p_Lead_FTOF = new TH2D("mmiss_mom_p_Lead_FTOF","m_{miss} vs. p_{p,Lead};m_{miss};p_{p,Lead};Counts",100,0.4,1.4,100,0,4);
  hist_list_2.push_back(h_mmiss_mom_p_Lead_FTOF);
  TH2D * h_mmiss_momT_p_Lead_FTOF = new TH2D("mmiss_momT_p_Lead_FTOF","m_{miss} vs. p_{p,T,Lead};m_{miss};p_{p,T,Lead};Counts",100,0.4,1.4,100,0,2.5);
  hist_list_2.push_back(h_mmiss_momT_p_Lead_FTOF);

  /////////////////////////////////////
  //Lead Proton CTOF Checks
  /////////////////////////////////////
  TH1D * h_theta_p_Lead_CTOF = new TH1D("theta_p_Lead_CTOF","#theta_{p,Lead};#theta_{p,Lead};Counts",100,20,140);
  hist_list_1.push_back(h_theta_p_Lead_CTOF);
  TH1D * h_theta_pq_Lead_CTOF = new TH1D("theta_pq_Lead_CTOF","#theta_{pq};#theta_{pq};Counts",100,0,30);
  hist_list_1.push_back(h_theta_pq_Lead_CTOF);
  TH2D * h_mom_theta_p_Lead_CTOF = new TH2D("mom_theta_p_Lead_CTOF","p_{p,Lead} vs. #theta_{p,Lead} ;p_{p,Lead};#theta_{p,Lead}",100,0,4,100,20,140);
  hist_list_2.push_back(h_mom_theta_p_Lead_CTOF);
  TH1D * h_phi_e_p_Lead_CTOF = new TH1D("phi_e_p_Lead_CTOF","|#phi_{e} - #phi_{p,Lead}|;|#phi_{e} - #phi_{p,Lead}|,Counts",100,140,180);
  hist_list_1.push_back(h_phi_e_p_Lead_CTOF);
  TH1D * h_xB_Lead_CTOF = new TH1D("xB_Lead_CTOF","x_{B} Lead;x_{B};Counts",100,0.0,2.0);
  hist_list_1.push_back(h_xB_Lead_CTOF);
  TH2D * h_vtz_e_vtz_p_Lead_CTOF = new TH2D("vtz_e_vtz_p_Lead_CTOF","Electron Z Vertex vs. Proton Z Vertex;vertex e;vertex p",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_vtz_e_vtz_p_Lead_CTOF);


  TH1D * h_pmiss_Lead_CTOF = new TH1D("pmiss_Lead_CTOF","p_{miss};p_{miss};Counts",100,0,2.0);
  hist_list_1.push_back(h_pmiss_Lead_CTOF);
  TH2D * h_pmiss_thetamiss_Lead_CTOF = new TH2D("pmiss_thetamiss_Lead_CTOF","p_{miss} vs. #theta_{miss};p_{miss};#theta_{miss}",100,0,2.0,180,0,180);
  hist_list_2.push_back(h_pmiss_thetamiss_Lead_CTOF);
  TH2D * h_xB_theta_1q_Lead_CTOF = new TH2D("xB_theta_1q_Lead_CTOF","x_{B} vs. #theta_{miss,q};x_{B};#theta_{miss,q};Counts",100,0,2,180,0,180);
  hist_list_2.push_back(h_xB_theta_1q_Lead_CTOF);
  TH2D * h_Loq_theta_1q_Lead_CTOF = new TH2D("Loq_theta_1q_Lead_CTOF","|p|/|q| vs. #theta_{miss,q};|p|/|q|;#theta_{miss,q}",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_Loq_theta_1q_Lead_CTOF);


  TH1D * h_mmiss_Lead_CTOF = new TH1D("mmiss_Lead_CTOF","m_{miss} Lead;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_Lead_CTOF);
  TH2D * h_mmiss_phi_e_p_Lead_CTOF = new TH2D("mmiss_phi_e_p_Lead_CTOF","m_{miss} vs. |#phi_{e} - #phi_{p}|;m_{miss};|#phi_{e} - #phi_{p};Counts",100,0.4,1.4,100,140,180);
  hist_list_2.push_back(h_mmiss_phi_e_p_Lead_CTOF);
  TH2D * h_mmiss_xB_Lead_CTOF = new TH2D("mmiss_xB_Lead_CTOF","m_{miss} vs. x_{B};m_{miss};x_{B};Counts",100,0.4,1.4,100,0.0,2.0);
  hist_list_2.push_back(h_mmiss_xB_Lead_CTOF);
  TH2D * h_mmiss_pmiss_Lead_CTOF = new TH2D("mmiss_pmiss_Lead_CTOF","m_{miss} vs. p_{miss};m_{miss};p_{miss};Counts",100,0.4,1.4,100,0.0,2.0);
  hist_list_2.push_back(h_mmiss_pmiss_Lead_CTOF);
  TH2D * h_mmiss_theta_1q_Lead_CTOF = new TH2D("mmiss_theta_1q_Lead_CTOF","m_{miss} vs. #theta_{miss,q};m_{miss};#theta_{miss,q};Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_theta_1q_Lead_CTOF);
  TH2D * h_mmiss_theta_p_Lead_CTOF = new TH2D("mmiss_theta_p_Lead_CTOF","m_{miss} vs. #theta_{p,Lead};m_{miss};#theta_{p,Lead};Counts",100,0.4,1.4,100,20,140);
  hist_list_2.push_back(h_mmiss_theta_p_Lead_CTOF);
  TH2D * h_mmiss_mom_p_Lead_CTOF = new TH2D("mmiss_mom_p_Lead_CTOF","m_{miss} vs. p_{p,Lead};m_{miss};p_{p,Lead};Counts",100,0.4,1.4,100,0,4);
  hist_list_2.push_back(h_mmiss_mom_p_Lead_CTOF);
  TH2D * h_mmiss_momT_p_Lead_CTOF = new TH2D("mmiss_momT_p_Lead_CTOF","m_{miss} vs. p_{p,T,Lead};m_{miss};p_{p,T,Lead};Counts",100,0.4,1.4,100,0,2.5);
  hist_list_2.push_back(h_mmiss_momT_p_Lead_CTOF);
*/

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
      auto allParticles = c12->getDetParticles();
      auto electrons=c12->getByID(11);
      auto protons=c12->getByID(2212);
      auto pion0=c12->getByID(111);
      auto kaonlong=c12->getByID(130);
      auto kaonshort=c12->getByID(310);
      auto neutrons=c12->getByID(2112);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);
      if(electrons.size()!=1){ continue;}
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double EoP_e =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
      int nphe = electrons[0]->che(HTCC)->getNphe();
      double vtz_e = electrons[0]->par()->getVz();
  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////      
      h_phi_theta->Fill(p_e.Phi()*180/M_PI,p_e.Theta()*180/M_PI,weight);
      if(EoP_e<=0){ continue; }
      int esector = electrons[0]->getSector();
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
      double theta_q =  p_q.Theta() * 180 / M_PI;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;
      //if(WSq>1.25){continue;}

      h_xB->Fill(xB,weight);
      h_QSq->Fill(QSq,weight);
      h_WSq->Fill(WSq,weight);
      h_W->Fill(sqrt(WSq),weight);
      h_xB_QSq->Fill(xB,QSq,weight);
      h_xB_WSq->Fill(xB,WSq,weight);
      h_QSq_WSq->Fill(QSq,WSq,weight);
      h_mom_theta[esector-1]->Fill(p_e.Mag(),theta_e,weight);

  /////////////////////////////////////
  //Electron Kinematics with SRC kinematics
  /////////////////////////////////////
      if((xB>1) && (QSq>1.5)){
	h_mom_theta_eSRC->Fill(p_e.Mag(),theta_e,weight);
	h_xB_theta_eSRC->Fill(xB,theta_e,weight);
	h_thetaq_theta_eSRC->Fill(theta_q,theta_e,weight);
	h_thetaq_xB_eSRC->Fill(theta_q,xB,weight);
	h_thetaq_eSRC->Fill(theta_q,weight);
      }
      
  /////////////////////////////////////
  //Hadron Information
  /////////////////////////////////////
      for(int j = 0; j < allParticles.size(); j ++){
	if(allParticles[j]->par()->getCharge()!=0){continue;}
	//cout<<"PID="<<allParticles[j]->getPid()<<endl;
	bool PCAL = (allParticles[j]->cal(clas12::PCAL)->getDetector() == 7);
	bool ECin = (allParticles[j]->cal(clas12::ECIN)->getDetector() == 7);
        bool ECout = (allParticles[j]->cal(clas12::ECOUT)->getDetector() == 7);
        bool CND1 = (allParticles[j]->sci(clas12::CND1)->getDetector() == 3);
        bool CND2 = (allParticles[j]->sci(clas12::CND2)->getDetector() == 3);
        bool CND3 = (allParticles[j]->sci(clas12::CND3)->getDetector() == 3);
	//bool FTOF1A = (allParticles[j]->sci(clas12::FTOF1A)->getDetector() == 12);
	//bool FTOF1B = (allParticles[j]->sci(clas12::FTOF1B)->getDetector() == 12);
	//bool FTOF2 = (allParticles[j]->sci(clas12::FTOF2)->getDetector() == 12);
	//bool CTOF = (allParticles[j]->sci(clas12::CTOF)->getDetector() == 4);

	if(PCAL || ECin || ECout){
	  h_mom_beta_hadplus_FTOF->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
	}
	if(CND1 || CND2 || CND3){
	  h_mom_beta_hadplus_CTOF->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
	  h_mom_beta_hadplus_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
 	  
	  int PID = allParticles[j]->getPid();
	  if(PID==2112){
	    h_mom_beta_proton_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }
	  if(PID==111){
	    h_mom_beta_pion_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }
	  if(PID==130 || PID==310){
	    h_mom_beta_kaon_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }
	  if(PID==22){
	    h_mom_beta_deuteron_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
	  }
	  if(PID==0){
	    h_mom_beta_0_CTOF_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }

	}
      }
  /////////////////////////////////////
  //FTOF and CTOF Proton Information
  /////////////////////////////////////
      for(int j = 0; j < neutrons.size(); j++){
	TVector3 p_p;
	p_p.SetMagThetaPhi(neutrons[j]->getP(),neutrons[j]->getTheta(),neutrons[j]->getPhi());
	double E_p = sqrt(mN*mN + p_p.Mag2());
	double theta_p = p_p.Theta() * 180 / M_PI;
	double phi_p = p_p.Phi() * 180 / M_PI;
	double theta_pq = p_p.Angle(p_q) * 180 / M_PI;
	double beta_p = neutrons[j]->par()->getBeta();
	double Chi2Pid_p = neutrons[j]->par()->getChi2Pid();
	double phi_diff = get_phi_diff(p_e,p_p);
	double vtz_p = neutrons[j]->par()->getVz();
	double vtz_ep_delta = vtz_e - vtz_p;

	double path_p = neutrons[j]->getPath();
	double beta_frommom_p = p_p.Mag()/E_p;
	double time_frommom_p = path_p / (c*beta_frommom_p);
	double time_frombeta_p = path_p / (c*beta_p);

	bool PCAL = (neutrons[j]->cal(clas12::PCAL)->getDetector() == 7);
	bool ECin = (neutrons[j]->cal(clas12::ECIN)->getDetector() == 7);
        bool ECout = (neutrons[j]->cal(clas12::ECOUT)->getDetector() == 7);
        bool CND1 = (neutrons[j]->sci(clas12::CND1)->getDetector() == 3);
        bool CND2 = (neutrons[j]->sci(clas12::CND2)->getDetector() == 3);
        bool CND3 = (neutrons[j]->sci(clas12::CND3)->getDetector() == 3);

	//bool FTOF1A = (protons[j]->sci(clas12::FTOF1A)->getDetector() == 12);
	//bool FTOF1B = (protons[j]->sci(clas12::FTOF1B)->getDetector() == 12);
	//bool FTOF2 = (protons[j]->sci(clas12::FTOF2)->getDetector() == 12);
	//bool CTOF = (protons[j]->sci(clas12::CTOF)->getDetector() == 4);

	if(PCAL || ECin || ECout){
	  h_chiSq_p_FTOF->Fill(Chi2Pid_p,weight);
	  h_vtz_p_FTOF->Fill(vtz_p,weight);
	  h_vtz_ep_delta_FTOF->Fill(vtz_ep_delta,weight);
	  h_vtz_e_vtz_p_FTOF->Fill(vtz_e,vtz_p,weight);
	  h_phi_e_p_FTOF->Fill(phi_diff,weight);

	  h_theta_p_FTOF->Fill(theta_p,weight);
	  h_theta_pq_FTOF->Fill(theta_pq,weight);	
	  h_phi_theta_p_FTOF->Fill(phi_p,theta_p,weight);
	  h_mom_beta_p_FTOF->Fill(p_p.Mag(),beta_p,weight);
	  h_timediff_p_FTOF->Fill(time_frommom_p-time_frombeta_p,weight);

	  h_mom_theta_p_FTOF[neutrons[j]->getSector()-1]->Fill(p_p.Mag(),theta_p,weight);
	}

	if(CND){
	  h_chiSq_p_CTOF->Fill(Chi2Pid_p,weight);
	  h_vtz_p_CTOF->Fill(vtz_p,weight);
	  h_vtz_ep_delta_CTOF->Fill(vtz_ep_delta,weight);
	  h_vtz_e_vtz_p_CTOF->Fill(vtz_e,vtz_p,weight);
	  h_phi_e_p_CTOF->Fill(phi_diff,weight);

	  h_theta_p_CTOF->Fill(theta_p,weight);
	  h_theta_pq_CTOF->Fill(theta_pq,weight);	
	  h_phi_theta_p_CTOF->Fill(phi_p,theta_p,weight);
	  h_mom_beta_p_CTOF->Fill(p_p.Mag(),beta_p,weight);
	  h_timediff_p_CTOF->Fill(time_frommom_p-time_frombeta_p,weight);
	  
	  h_mom_theta_p_CTOF->Fill(p_p.Mag(),theta_p,weight);
	}
      }

/*  /////////////////////////////////////
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

      bool PCAL = (allParticles[j]->cal(clas12::PCAL)->getDetector() == 7);
      bool ECin = (allParticles[j]->cal(clas12::ECIN)->getDetector() == 7);
      bool ECout = (allParticles[j]->cal(clas12::ECOUT)->getDetector() == 7);
      bool CND1 = (allParticles[j]->sci(clas12::CND1)->getDetector() == 3);
      bool CND2 = (allParticles[j]->sci(clas12::CND2)->getDetector() == 3);
      bool CND3 = (allParticles[j]->sci(clas12::CND3)->getDetector() == 3);

      //bool FTOF1A = (protons[index_L]->sci(clas12::FTOF1A)->getDetector() == 12);
      //bool FTOF1B = (protons[index_L]->sci(clas12::FTOF1B)->getDetector() == 12);
      //bool FTOF2 =  (protons[index_L]->sci(clas12::FTOF2)->getDetector() == 12);
      //bool CTOF =   (protons[index_L]->sci(clas12::CTOF)->getDetector() == 4);

      if(PCAL || ECin || ECout){
	h_theta_p_Lead_FTOF->Fill(theta_L,weight);
	h_theta_pq_Lead_FTOF->Fill(theta_Lq,weight);
	h_mom_theta_p_Lead_FTOF->Fill(p_L.Mag(),theta_L,weight);
	h_phi_e_p_Lead_FTOF->Fill(phi_diff,weight);
	h_xB_Lead_FTOF->Fill(xB,weight);
	h_vtz_e_vtz_p_Lead_FTOF->Fill(vtz_e,vtz_p,weight);
	
	h_pmiss_Lead_FTOF->Fill(p_miss.Mag(),weight);
	h_pmiss_thetamiss_Lead_FTOF->Fill(p_miss.Mag(),theta_miss,weight);
	h_xB_theta_1q_Lead_FTOF->Fill(xB,theta_1q,weight);
	h_Loq_theta_1q_Lead_FTOF->Fill(Loq,theta_1q,weight);
	
	h_mmiss_Lead_FTOF->Fill(mmiss,weight);
	h_mmiss_phi_e_p_Lead_FTOF->Fill(mmiss,phi_diff,weight);
	h_mmiss_xB_Lead_FTOF->Fill(mmiss,xB,weight);
	h_mmiss_pmiss_Lead_FTOF->Fill(mmiss,p_miss.Mag(),weight);
	h_mmiss_theta_1q_Lead_FTOF->Fill(mmiss,theta_1q,weight);
	h_mmiss_theta_p_Lead_FTOF->Fill(mmiss,theta_L,weight);
	h_mmiss_mom_p_Lead_FTOF->Fill(mmiss,p_L.Mag(),weight);
	h_mmiss_momT_p_Lead_FTOF->Fill(mmiss,p_L.Perp(),weight);
      }
      if(CND){
	h_theta_p_Lead_CTOF->Fill(theta_L,weight);
	h_theta_pq_Lead_CTOF->Fill(theta_Lq,weight);
	h_mom_theta_p_Lead_CTOF->Fill(p_L.Mag(),theta_L,weight);
	h_phi_e_p_Lead_CTOF->Fill(phi_diff,weight);
	h_xB_Lead_CTOF->Fill(xB,weight);
	h_vtz_e_vtz_p_Lead_CTOF->Fill(vtz_e,vtz_p,weight);
	
	h_pmiss_Lead_CTOF->Fill(p_miss.Mag(),weight);
	h_pmiss_thetamiss_Lead_CTOF->Fill(p_miss.Mag(),theta_miss,weight);
	h_xB_theta_1q_Lead_CTOF->Fill(xB,theta_1q,weight);
	h_Loq_theta_1q_Lead_CTOF->Fill(Loq,theta_1q,weight);
	
	h_mmiss_Lead_CTOF->Fill(mmiss,weight);
	h_mmiss_phi_e_p_Lead_CTOF->Fill(mmiss,phi_diff,weight);
	h_mmiss_xB_Lead_CTOF->Fill(mmiss,xB,weight);
	h_mmiss_pmiss_Lead_CTOF->Fill(mmiss,p_miss.Mag(),weight);
	h_mmiss_theta_1q_Lead_CTOF->Fill(mmiss,theta_1q,weight);
	h_mmiss_theta_p_Lead_CTOF->Fill(mmiss,theta_L,weight);
	h_mmiss_mom_p_Lead_CTOF->Fill(mmiss,p_L.Mag(),weight);
	h_mmiss_momT_p_Lead_CTOF->Fill(mmiss,p_L.Perp(),weight);
      }
*/
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
  double line = 0.8;
  if(myCut.getDoCut(e_cuts)){
    myCut.print_cut_onPDF(text,e_nphe,line);
    myCut.print_cut_onPDF(text,e_calv,line);
    myCut.print_cut_onPDF(text,e_calw,line);
    myCut.print_cut_onPDF(text,e_SF,line);
    myCut.print_cut_onPDF(text,e_mom,line);
    myCut.print_cut_onPDF(text,e_vtze,line);
  }
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

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_W->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Electron Kinematics with SRC kinematics
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts");
  text.DrawLatex(0.2,0.8,"x_{B}: min=1. max=2.");
  text.DrawLatex(0.2,0.7,"Q^{2}: min=1.5 max=10.");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mom_theta_eSRC->Draw("colz");
  myCanvas->cd(2);
  h_xB_theta_eSRC->Draw("colz");
  myCanvas->cd(3);
  h_thetaq_theta_eSRC->Draw("colz");
  myCanvas->cd(4);
  h_thetaq_xB_eSRC->Draw("colz");
  myCanvas->cd(5);
  h_thetaq_eSRC->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Hadron Detected
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts with positive Hadron");
  myText->Print(fileName,"pdf");  
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mom_beta_hadplus_FTOF->Draw("colz");
  myCanvas->cd(2);
  h_mom_beta_hadplus_CTOF->Draw("colz");
  myCanvas->cd(3);
  h_mom_beta_hadplus_CTOF_zoom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mom_beta_proton_CTOF_zoom->Draw("colz");
  myCanvas->cd(2);
  h_mom_beta_pion_CTOF_zoom->Draw("colz");
  myCanvas->cd(3);
  h_mom_beta_kaon_CTOF_zoom->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_deuteron_CTOF_zoom->Draw("colz");
  myCanvas->cd(5);
  h_mom_beta_0_CTOF_zoom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //ECAL Neutron Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in ECAL");
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
  myCanvas->cd(5);
  h_phi_e_p_FTOF->Draw();
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
  myCanvas->cd(5);
  h_timediff_p_FTOF->Draw();
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
  //CND Neutron Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in CND");
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
  myCanvas->cd(5);
  h_phi_e_p_CTOF->Draw();
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
  h_timediff_p_CTOF->Draw();
  myCanvas->cd(6);
  h_mom_theta_p_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

/*  /////////////////////////////////////
  //Lead Proton FTOF Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  text.DrawLatex(0.2,0.7,"Lead Scintillator: FTOF");
  line = 0.6;
  if(myCut.getDoCut(l_cuts)){
    myCut.print_cut_onPDF(text,l_pid,line);
    myCut.print_cut_onPDF(text,l_scint,line);
    myCut.print_cut_onPDF(text,l_theta,line);
    myCut.print_cut_onPDF(text,l_thetalq,line);
    myCut.print_cut_onPDF(text,l_chipid,line);
    myCut.print_cut_onPDF(text,l_vtzdiff,line);
    myCut.print_cut_onPDF(text,l_phidiff,line);
  }
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_Lead_FTOF->Draw();
  myCanvas->cd(2);
  h_theta_pq_Lead_FTOF->Draw();
  myCanvas->cd(3);
  h_mom_theta_p_Lead_FTOF->Draw("colz");
  myCanvas->cd(4);
  h_phi_e_p_Lead_FTOF->Draw();
  myCanvas->cd(5);
  h_xB_Lead_FTOF->Draw();
  myCanvas->cd(6);
  h_vtz_e_vtz_p_Lead_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_Lead_FTOF->Draw();
  myCanvas->cd(2);
  h_pmiss_thetamiss_Lead_FTOF->Draw("colz");
  myCanvas->cd(3);
  h_xB_theta_1q_Lead_FTOF->Draw("colz");
  myCanvas->cd(4);
  h_Loq_theta_1q_Lead_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_Lead_FTOF->Draw();
  myCanvas->cd(2);
  h_mmiss_phi_e_p_Lead_FTOF->Draw("colz");
  myCanvas->cd(3);
  h_mmiss_xB_Lead_FTOF->Draw("colz");
  myCanvas->cd(4);
  h_mmiss_pmiss_Lead_FTOF->Draw("colz");
  myCanvas->cd(5);
  h_mmiss_theta_1q_Lead_FTOF->Draw("colz");
  myCanvas->cd(6);
  h_mmiss_theta_p_Lead_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_mom_p_Lead_FTOF->Draw("colz");
  myCanvas->cd(2);
  h_mmiss_momT_p_Lead_FTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Lead Proton CTOF Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  text.DrawLatex(0.2,0.7,"Lead Scintillator: CTOF");
  line = 0.6;
  if(myCut.getDoCut(l_cuts)){
    myCut.print_cut_onPDF(text,l_pid,line);
    myCut.print_cut_onPDF(text,l_scint,line);
    myCut.print_cut_onPDF(text,l_theta,line);
    myCut.print_cut_onPDF(text,l_thetalq,line);
    myCut.print_cut_onPDF(text,l_chipid,line);
    myCut.print_cut_onPDF(text,l_vtzdiff,line);
    myCut.print_cut_onPDF(text,l_phidiff,line);
  }
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_Lead_CTOF->Draw();
  myCanvas->cd(2);
  h_theta_pq_Lead_CTOF->Draw();
  myCanvas->cd(3);
  h_mom_theta_p_Lead_CTOF->Draw("colz");
  myCanvas->cd(4);
  h_phi_e_p_Lead_CTOF->Draw();
  myCanvas->cd(5);
  h_xB_Lead_CTOF->Draw();
  myCanvas->cd(6);
  h_vtz_e_vtz_p_Lead_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_Lead_CTOF->Draw();
  myCanvas->cd(2);
  h_pmiss_thetamiss_Lead_CTOF->Draw("colz");
  myCanvas->cd(3);
  h_xB_theta_1q_Lead_CTOF->Draw("colz");
  myCanvas->cd(4);
  h_Loq_theta_1q_Lead_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_Lead_CTOF->Draw();
  myCanvas->cd(2);
  h_mmiss_phi_e_p_Lead_CTOF->Draw("colz");
  myCanvas->cd(3);
  h_mmiss_xB_Lead_CTOF->Draw("colz");
  myCanvas->cd(4);
  h_mmiss_pmiss_Lead_CTOF->Draw("colz");
  myCanvas->cd(5);
  h_mmiss_theta_1q_Lead_CTOF->Draw("colz");
  myCanvas->cd(6);
  h_mmiss_theta_p_Lead_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_mmiss_mom_p_Lead_CTOF->Draw("colz");
  myCanvas->cd(2);
  h_mmiss_momT_p_Lead_CTOF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
*/
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

