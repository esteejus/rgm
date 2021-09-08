#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>

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
  std::cerr << "Usage: ./monitorPID <path/to/ouput.root> <path/to/ouput.pdf>  <path/to/input.hipo> \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }

  eventcut myCut();
  /////////////////////////////////////
  TFile * outFile = new TFile(argv[1],"RECREATE");
  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  bool isMC = false;
  bool input_energy = false; 
  double beam_energyIN = 0; //input beam energy which overides any info in RCDB or for simulation

  for(;;)
    {
      switch(getopt(argc-3, &argv[3], "SE:")) // note the colon (:) to indicate that 'b' has a parameter and is not a switch
	{
	case 'S':
	  printf("Simulation(MC) option 'S' specified\n");
	  isMC = true;
	  continue;
	case 'E':
	  printf("Beam Energy 'E' specified with the value %s\n (GeV)", optarg);
	  input_energy = true;
	  beam_energyIN = atof(optarg);
	  continue;

	case '?':
	case 'h':
	default :
	  printf("Help/Usage Example\n");
	  break;

	case -1:
	  break;
	}

      break;
    }

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
  //All Neutron Angles
  /////////////////////////////////////
  TH2D * h_ToM_eDep_nL = new TH2D("ToM_eDep_nL","ToM vs.eDep;ToM;eDep;events",100,0,15,100,0,50);
  hist_list_2.push_back(h_ToM_eDep_nL);
  TH1D * h_ToM_nL = new TH1D("ToM_nL","ToM;ToM;events",100,0,15);
  hist_list_1.push_back(h_ToM_nL);
  TH2D * h_mom_ToM_nL = new TH2D("mom_ToM_nL","mom vs. ToM;mom;ToM;events",100,0,2,100,0,15);
  hist_list_2.push_back(h_mom_ToM_nL);

  TH1D * h_mom_nL = new TH1D("mom_nL","mom;mom;events",100,0,2);
  hist_list_1.push_back(h_mom_nL);
  TH1D * h_theta_nL = new TH1D("theta_nL","theta_nL;theta_nL",180,0,180);
  hist_list_1.push_back(h_theta_nL);
  TH1D * h_theta_nLq = new TH1D("theta_nLq","theta_nLq;theta_nLq",180,0,180);
  hist_list_1.push_back(h_theta_nLq);
  TH1D * h_phi_e_nL = new TH1D("phi_e_nL","phi_e minus phi_nL;phi_e_nL",180,0,180);
  hist_list_1.push_back(h_phi_e_nL);


  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
  TH1D * h_theta_L_FTOF = new TH1D("theta_L_FTOF","#theta_{proton} Lead;theta_{proton}",180,0,180);
  hist_list_1.push_back(h_theta_L_FTOF);
  TH1D * h_theta_Lq_FTOF = new TH1D("theta_Lq_FTOF","#theta_{pq} Lead;#theta_{pq}",180,0,180);
  hist_list_1.push_back(h_theta_Lq_FTOF);
  TH1D * h_phi_e_L = new TH1D("phi_e_L","|#phi_{e} - #phi_{p}|;|#phi_{e} - #phi_{p}|,Counts",180,0,180);
  hist_list_1.push_back(h_phi_e_L);
  TH1D * h_mmiss_FTOF = new TH1D("mmiss_FTOF","m_{miss};m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_FTOF);
  TH2D * h_mmiss_phi_e_L = new TH2D("mmiss_phi_e_L","m_{miss} vs. |#phi_{e} - #phi_{p}|;m_{miss};|#phi_{e} - #phi_{p}",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_phi_e_L);
  TH2D * h_xB_mmiss = new TH2D("xB_mmiss","x_{B} vs. m_{miss};xB;mmiss",100,0,2,100,0.4,1.4);
  hist_list_2.push_back(h_xB_mmiss);
  TH2D * h_pmiss_mmiss = new TH2D("pmiss_mmiss","p_{miss} vs. m_{miss};p_{miss};m_{miss}",100,0,1.5,100,0.4,1.4);
  hist_list_2.push_back(h_pmiss_mmiss);
  TH2D * h_xB_theta_1q = new TH2D("xB_theta_1q","x_{B} vs. #theta_{p_{i},q};x_{B};#theta_{p_{i},q}",100,0,2,180,0,180);
  hist_list_2.push_back(h_xB_theta_1q);
  TH2D * h_Loq_theta_1q = new TH2D("Loq_theta_1q","|p|/|q| vs. #theta_{p_{i},q};Loq;theta_1q",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_Loq_theta_1q);
  TH2D * h_pmiss_theta_miss = new TH2D("pmiss_theta_miss","p_{miss} vs. #theta_{miss};pmiss;theta_miss",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss);


  h_theta_L_FTOF -> GetXaxis()->SetTitle("#theta_{proton} FTOF");
  h_theta_L_FTOF -> GetYaxis()->SetTitle("Counts");
  h_theta_L_FTOF -> GetXaxis()->CenterTitle();
  h_theta_L_FTOF -> GetYaxis()->CenterTitle();

  h_theta_Lq_FTOF -> GetXaxis()->SetTitle("#theta_{pq} FTOF");
  h_theta_Lq_FTOF -> GetYaxis()->SetTitle("Counts");
  h_theta_Lq_FTOF -> GetXaxis()->CenterTitle();
  h_theta_Lq_FTOF -> GetYaxis()->CenterTitle();

  h_phi_e_L -> GetXaxis()->SetTitle("#phi_{e^{-}} - #phi_{proton}");
  h_phi_e_L -> GetYaxis()->SetTitle("Counts");
  h_phi_e_L -> GetXaxis()->CenterTitle();
  h_phi_e_L -> GetYaxis()->CenterTitle();

  h_mmiss_phi_e_L -> GetXaxis()->SetTitle("Missing Mass (GeV)");
  h_mmiss_phi_e_L -> GetYaxis()->SetTitle("#phi_{e^{-}} - #phi_{proton}");
  h_mmiss_phi_e_L -> GetXaxis()->CenterTitle();
  h_mmiss_phi_e_L -> GetYaxis()->CenterTitle();

  h_xB_mmiss -> GetXaxis()->SetTitle("x_{B}");
  h_xB_mmiss -> GetYaxis()->SetTitle("Missing Mass (GeV)");
  h_xB_mmiss -> GetXaxis()->CenterTitle();
  h_xB_mmiss -> GetYaxis()->CenterTitle();

  h_pmiss_mmiss -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss_mmiss -> GetYaxis()->SetTitle("Missing Mass (GeV)");
  h_pmiss_mmiss -> GetXaxis()->CenterTitle();
  h_pmiss_mmiss -> GetYaxis()->CenterTitle();

  h_xB_theta_1q -> GetXaxis()->SetTitle("x_{B}");
  h_xB_theta_1q -> GetYaxis()->SetTitle("#theta_{Recoil,pq}");
  h_xB_theta_1q -> GetXaxis()->CenterTitle();
 h_xB_theta_1q -> GetYaxis()->CenterTitle();

  h_pmiss_theta_miss -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss_theta_miss -> GetYaxis()->SetTitle("#theta_{miss}");
  h_pmiss_theta_miss -> GetXaxis()->CenterTitle();
  h_pmiss_theta_miss -> GetYaxis()->CenterTitle();

  h_Loq_theta_1q -> GetXaxis()->SetTitle("|p|/|q|");
  h_Loq_theta_1q -> GetYaxis()->SetTitle("#theta_{p_{recoil}q}");
  h_Loq_theta_1q -> GetXaxis()->CenterTitle();
  h_Loq_theta_1q -> GetYaxis()->CenterTitle();


  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  TH1D * h_pmiss = new TH1D("pmiss","p_{miss};pmiss",100,0,1.5);
  hist_list_1.push_back(h_pmiss);
  TH1D * h_mmiss = new TH1D("mmiss","m_{miss};mmiss",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss);
  TH2D * h_pmiss_theta_miss_SRC = new TH2D("pmiss_theta_miss_SRC","p_{miss} vs. #theta_{miss};pmiss;theta_1",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_xB_Loq_SRC = new TH2D("xB_Loq","x_{B} vs |p|/|q|;xB;Loq",100,0,2,100,0,1.5);
  hist_list_2.push_back(h_xB_Loq_SRC);


  TH1D * h_pmiss_tight = new TH1D("pmiss_tight","p_{miss};pmiss",100,0,1.5);
  hist_list_1.push_back(h_pmiss_tight);
  TH1D * h_mmiss_tight = new TH1D("mmiss_tight","m_{miss};mmiss",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_tight);
  TH2D * h_pmiss_theta_miss_SRC_tight = new TH2D("pmiss_theta_miss_SRC_tight","p_{miss} vs. #theta_{miss};pmiss;theta_1",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC_tight);
  

  h_pmiss -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss -> GetYaxis()->SetTitle("Counts");
  h_pmiss -> GetXaxis()->CenterTitle();
  h_pmiss -> GetYaxis()->CenterTitle();

  h_mmiss -> GetXaxis()->SetTitle("Missing Mass (GeV)");
  h_mmiss -> GetYaxis()->SetTitle("Counts");
  h_mmiss -> GetXaxis()->CenterTitle();
  h_mmiss -> GetYaxis()->CenterTitle();

  h_pmiss_theta_miss_SRC -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss_theta_miss_SRC -> GetYaxis()->SetTitle("#theta_{miss}");
  h_pmiss_theta_miss_SRC -> GetXaxis()->CenterTitle();
  h_pmiss_theta_miss_SRC -> GetYaxis()->CenterTitle();

  h_xB_Loq_SRC -> GetXaxis()->SetTitle("x_{B}");
  h_xB_Loq_SRC -> GetYaxis()->SetTitle("#theta_{pq}");
  h_xB_Loq_SRC -> GetXaxis()->CenterTitle();
  h_xB_Loq_SRC -> GetYaxis()->CenterTitle();

  h_mmiss_tight -> GetXaxis()->SetTitle("Missing Mass (GeV)");
  h_mmiss_tight -> GetYaxis()->SetTitle("Counts");
  h_mmiss_tight -> GetXaxis()->CenterTitle();
  h_mmiss_tight -> GetYaxis()->CenterTitle();

  h_pmiss_tight -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss_tight -> GetYaxis()->SetTitle("Counts");
  h_pmiss_tight -> GetXaxis()->CenterTitle();
  h_pmiss_tight -> GetYaxis()->CenterTitle();

  h_pmiss_theta_miss_SRC_tight -> GetXaxis()->SetTitle("Missing Momentum (GeV)");
  h_pmiss_theta_miss_SRC_tight -> GetYaxis()->SetTitle("#theta_{miss}");
  h_pmiss_theta_miss_SRC_tight -> GetXaxis()->CenterTitle();
  h_pmiss_theta_miss_SRC_tight -> GetYaxis()->CenterTitle();



  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
  TH1D * h_num_R = new TH1D("num_R","num Recoil;num_R",5,0,5);
  hist_list_1.push_back(h_num_R);
  TH1D * h_p_2 = new TH1D("p_2","p Recoil;p_2",100,0,1.5);
  hist_list_1.push_back(h_p_2);

  TH1D * h_num_R_tight = new TH1D("num_R_tight","#Recoil Protons;#Recoil Protons;Counts",5,0,5);
  hist_list_1.push_back(h_num_R_tight);
  TH1D * h_p_2_tight = new TH1D("p_2_tight","p_{rec};p_{rec};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_2_tight);
  
  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////
  TH1D * h_p_2_high = new TH1D("p_2_high","p_{rec};p_{rec};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_2_high);
  TH1D * h_p_rel = new TH1D("p_rel","p_{rel};p_{rel};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_rel);
  TH1D * h_p_cm = new TH1D("p_cm","p_{C.M.};p_{C.M.};Counts",100,0,0.5);
  hist_list_1.push_back(h_p_cm);
  TH1D * h_p_t_cm = new TH1D("p_t_cm","p_{t,C.M.};p_{t,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_t_cm);
  TH1D * h_p_y_cm = new TH1D("p_y_cm","p_{y,C.M.};p_{y,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_y_cm);
  TH1D * h_p_x_cm = new TH1D("p_x_cm","p_{x,C.M.};p_{x,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_x_cm);
  TH1D * h_theta_rel = new TH1D("theta_rel","#theta_{rel};#theta_{rel};Counts",180,0,180);
  hist_list_1.push_back(h_theta_rel);
  TH2D * h_p_cm_theta_rel = new TH2D("p_cm_theta_rel","p_{C.M.} vs. #theta_{rel};p_{C.M.};#theta_{rel}",100,0,0.5,180,0,180);
  hist_list_2.push_back(h_p_cm_theta_rel);

  TH1D * h_p_2_high_tight = new TH1D("p_2_high_tight","p_{rec};p_{rec};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_2_high_tight);
  TH1D * h_p_rel_tight = new TH1D("p_rel_tight","p_{rel};p_{rel};Counts",100,0,1.5);
  hist_list_1.push_back(h_p_rel_tight);
  TH1D * h_p_cm_tight = new TH1D("p_cm_tight","p_{C.M.};p_{C.M.};Counts",100,0,0.5);
  hist_list_1.push_back(h_p_cm_tight);
  TH1D * h_p_t_cm_tight = new TH1D("p_t_cm_tight","p_{t,C.M.};p_{t,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_t_cm_tight);
  TH1D * h_p_y_cm_tight = new TH1D("p_y_cm_tight","p_{y,C.M.};p_{y,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_y_cm_tight);
  TH1D * h_p_x_cm_tight = new TH1D("p_x_cm_tight","p_{x,C.M.};p_{x,C.M.};Counts",100,-0.5,0.5);
  hist_list_1.push_back(h_p_x_cm_tight);
  TH1D * h_theta_rel_tight = new TH1D("theta_rel_tight","#theta_{rel};#theta_{rel};Counts",180,0,180);
  hist_list_1.push_back(h_theta_rel_tight);
  TH2D * h_p_cm_theta_rel_tight = new TH2D("p_cm_theta_rel_tight","p_{C.M.} vs. #theta_{rel};p_{C.M.};#theta_{rel}",100,0,0.5,180,0,180);
  hist_list_2.push_back(h_p_cm_theta_rel_tight);

  /////////////////////////////////////
  //Testers
  /////////////////////////////////////
  TH1D * h_mmiss_all = new TH1D("mmiss_all","m_{miss} All;Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_all);
  TH1D * h_mmiss_FTOF1 = new TH1D("mmiss_FTOF1","m_{miss} FTOF1;Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_FTOF1);
  TH1D * h_mmiss_FTOF2 = new TH1D("mmiss_FTOF2","m_{miss} FTOF2;Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_FTOF2);
  TH1D * h_mmiss_noFTOF = new TH1D("mmiss_noFTOF","m_{miss} noFTOF;Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_noFTOF);
  TH2D * h_mmiss_FTOF2_mag = new TH2D("mmiss_FTOF2_mag","m_{miss} FTOF2; Mag;Counts",100,0.4,1.4,100,0,4);
  hist_list_2.push_back(h_mmiss_FTOF2_mag);
  TH2D * h_mmiss_FTOF2_theta = new TH2D("mmiss_FTOF2_theta","m_{miss} FTOF2; Theta;Counts",100,0.4,1.4,100,30,50);
  hist_list_2.push_back(h_mmiss_FTOF2_theta);
  TH2D * h_mmiss_FTOF2_phi = new TH2D("mmiss_FTOF2_phi","m_{miss} FTOF2; Phi;Counts",100,0.4,1.4,360,-180,180);
  hist_list_2.push_back(h_mmiss_FTOF2_phi);
  TH2D * h_mmiss_FTOF2_chi = new TH2D("mmiss_FTOF2_chi","m_{miss} FTOF2; Chi;Counts",100,0.4,1.4,100,-10,10);
  hist_list_2.push_back(h_mmiss_FTOF2_chi);

  TH2D * h_mmiss_noFTOF_theta = new TH2D("mmiss_noFTOF_theta","m_{miss} noFTOF; Theta;Counts",100,0.4,1.4,180,0,180);
  hist_list_2.push_back(h_mmiss_noFTOF_theta);


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

  //  clas12databases::SetCCDBRemoteConnection(); //not recommended
  //  clas12databases::SetRCDBRemoteConnection(); //not recommended
  
  //create a local database with PrepareDatabases.C supplied in ../ or $CLAS12ROOT/RunRoot/
  clas12databases::SetCCDBLocalConnection("../ccdb.sqlite"); 
  clas12databases::SetRCDBRootConnection("../rcdb.root");

  clas12root::HipoChain chain;
  //  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[3]<<endl;
    chain.Add(argv[3]);
    //  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});
  auto& c12=chain.C12ref();

  auto& rcdbData= config_c12->rcdb()->current();//struct with all relevent rcdb values
  auto Ebeam = rcdbData.beam_energy/1000;

  if( input_energy)
    Ebeam = beam_energyIN;

  while(chain.Next()==true)
    {
      if(!input_energy && counter == 0)
	Ebeam = rcdbData.beam_energy/1000;
      if(counter == 0)
	cout<<"Beam Energy "<<Ebeam<<endl;

      //Display completed  
      counter++;
      if((counter%100000) == 0)
	cerr << counter <<" completed \n";

      // get particles by type
      auto electrons=c12->getByID(11);
      auto protons=c12->getByID(2212);
      auto neutrons=c12->getByID(2112);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);
      
  /////////////////////////////////////
  //Electron fiducials
  /////////////////////////////////////
      if(electrons.size()!=1){continue;}
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());

      double EoP_e =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
      
      h_Vcal_EoP->Fill(electrons[0]->cal(PCAL)->getLv(),EoP_e,weight);
      h_Wcal_EoP->Fill(electrons[0]->cal(PCAL)->getLw(),EoP_e,weight);
      h_phi_theta->Fill(p_e.Phi()*180/M_PI,p_e.Theta()*180/M_PI,weight);
      h_sector->Fill(electrons[0]->getSector(),weight);

	  
  /////////////////////////////////////
  //Electron Pid
  /////////////////////////////////////
      h_P_EoP->Fill(p_e.Mag(),EoP_e,weight);
      int nphe = electrons[0]->che(HTCC)->getNphe();
      if(nphe < 2){ continue; }
      h_nphe->Fill(nphe,weight);
      if(electrons[0]->cal(ECIN)->getLv() < 14){ continue; }
      if(electrons[0]->cal(ECIN)->getLw() < 14){ continue; }
      if(EoP_e < 0.18){ continue; }
      if(EoP_e > 0.28){ continue; }
      if(p_e.Mag() < 1){ continue; }
      if(p_e.Mag() > 6.6){ continue; }
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
      //if(WSq>1.25){continue;}
      
      h_xB->Fill(xB,weight);
      h_QSq->Fill(QSq,weight);
      h_WSq->Fill(WSq,weight);
      h_xB_QSq->Fill(xB,QSq,weight);
      h_xB_WSq->Fill(xB,WSq,weight);
      h_QSq_WSq->Fill(QSq,WSq,weight);

      /*
      if((electrons.size()>0) && (num_L==1)){
	TVector3 p_e_t;
	p_e_t.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
	TVector3 p_p_t;
	p_p_t.SetMagThetaPhi(protons[index_L]->getP(),protons[index_L]->getTheta(),protons[index_L]->getPhi());
	double mmiss_t = get_mmiss(p_b,p_e_t,p_p_t);
	h_mmiss_t->Fill(mmiss_t,weight);
	}*/
  /////////////////////////////////////
  //All Proton Angles
  /////////////////////////////////////
      int num_L = 0;
      int index_L = -1;      
      for(int j = 0; j < protons.size(); j++){
	TVector3 p_L;
	p_L.SetMagThetaPhi(protons[j]->getP(),protons[j]->getTheta(),protons[j]->getPhi());
	double theta_L = p_L.Theta() * 180 / M_PI;
	double phi_L = p_L.Phi() * 180 / M_PI;
	double theta_Lq = p_L.Angle(p_q) * 180 / M_PI;
	double vtz_L = protons[j]->par()->getVz();
	double Chi2Pid_L = protons[j]->par()->getChi2Pid();
	double vtz_diff = electrons[0]->par()->getVz()-vtz_L;
	if(isMC){vtz_diff=0.5;}
	h_theta_L->Fill(theta_L,weight);
	h_theta_Lq->Fill(theta_Lq,weight);

	double mmiss_t = get_mmiss(p_b,p_e,p_L);
	h_mmiss_all->Fill(mmiss_t,weight);
	if((protons[j]->sci(FTOF1A)->getDetector()==12) && (protons[j]->sci(FTOF1B)->getDetector()==12)){
	  h_mmiss_FTOF1->Fill(mmiss_t,weight);
	}
	else if(protons[j]->sci(FTOF2)->getDetector()==12){
	  h_mmiss_FTOF2->Fill(mmiss_t,weight);
	  h_mmiss_FTOF2_mag->Fill(mmiss_t,p_L.Mag(),weight);
	  h_mmiss_FTOF2_theta->Fill(mmiss_t,theta_L,weight);
	  h_mmiss_FTOF2_phi->Fill(mmiss_t,phi_L,weight);
	  h_mmiss_FTOF2_chi->Fill(mmiss_t,Chi2Pid_L,weight);

	  /*
	  h_mmiss_FTOF2_PCAL->Fill(mmiss_t,protons[j]->cal(PCAL)->getEnergy(),weight);
	  h_mmiss_FTOF2_ECIN->Fill(mmiss_t,protons[j]->cal(ECIN)->getEnergy(),weight);
	  h_mmiss_FTOF2_ECOUT->Fill(mmiss_t,protons[j]->cal(ECOUT)->getEnergy(),weight);
	  */
	}
	else{
	  h_mmiss_noFTOF->Fill(mmiss_t,weight);
	  h_mmiss_noFTOF_theta->Fill(mmiss_t,theta_L,weight);
	}

	if(((protons[j]->sci(FTOF1A)->getDetector()==12) || (protons[j]->sci(FTOF1B)->getDetector()==12)) || (protons[j]->sci(FTOF2)->getDetector()==12))
	//if(((protons[j]->sci(FTOF1A)->getDetector()==12) || (protons[j]->sci(FTOF1B)->getDetector()==12)))
	  {
	    if(theta_Lq<25){
	      if(lowThetaCut(p_L.Theta(),Chi2Pid_L,vtz_diff)){
		num_L++;
		index_L=j;
	      }	  
	    }
	  }
      }      

  /////////////////////////////////////
  //All Neutron Angles
  /////////////////////////////////////
      for(int j = 0; j < neutrons.size(); j++){
	TVector3 p_nL;
	p_nL.SetMagThetaPhi(neutrons[j]->getP(),neutrons[j]->getTheta(),neutrons[j]->getPhi());
	double theta_nL = p_nL.Theta() * 180 / M_PI;
	double theta_nLq = p_nL.Angle(p_q) * 180 / M_PI;
	double phi_diff_n = get_phi_diff(p_e.Phi()*180/M_PI,p_nL.Phi()*180/M_PI);
	if(p_nL.Mag()<0.0001){continue;}
	
	int DET = 0;
	if(neutrons[j]->sci(FTOF1A)->getDetector()==12){DET=FTOF1A;}
	if(neutrons[j]->sci(FTOF1B)->getDetector()==12){DET=FTOF1B;}
	if(neutrons[j]->sci(FTOF2)->getDetector()==12){DET=FTOF2;}
	
	if(DET!=0){
	  double t_nL = neutrons[j]->sci(DET)->getTime() - c12->event()->getStartTime();;
	  double l_nL = neutrons[j]->sci(DET)->getPath();
	  double eDep_nL = neutrons[j]->sci(DET)->getEnergy();
	  double ToM = t_nL*100/l_nL;
	  h_ToM_eDep_nL->Fill(ToM,eDep_nL,weight);
	  if(eDep_nL>5){
	    h_ToM_nL->Fill(ToM,weight);
	    h_mom_ToM_nL->Fill(p_nL.Mag(),ToM,weight);
	    if(ToM>3.5){
	      h_mom_nL->Fill(p_nL.Mag(),weight);
	      h_theta_nL->Fill(theta_nL,weight);
	      if(theta_nL<35){
		h_theta_nLq->Fill(theta_nLq,weight);
		if(theta_nLq<25){
		  h_phi_e_nL->Fill(phi_diff_n,weight);
		}
	      }
	    }
	  }
	}
      }      
      
  /////////////////////////////////////
  //Lead Proton Checks
  /////////////////////////////////////
      bool isSRC_loose = false;
      bool isSRC_tight = false;
      if(num_L==1){
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
	if(QSq>1.5){
	  if(p_miss.Mag()>0.25){
	    if((mmiss>0.84) && (mmiss<1.04)){
	      h_pmiss->Fill(p_miss.Mag(),weight);
	      h_mmiss->Fill(mmiss,weight);
	      h_pmiss_theta_miss_SRC->Fill(p_miss.Mag(),theta_miss,weight);
	      h_xB_Loq_SRC->Fill(xB,Loq,weight);
	      isSRC_loose=true;	      
	      //if((Loq > 0.62) && (Loq < 0.96)){
	      //if(xB > 1.2){
		  h_pmiss_tight->Fill(p_miss.Mag(),weight);
		  h_mmiss_tight->Fill(mmiss,weight);
		  h_pmiss_theta_miss_SRC_tight->Fill(p_miss.Mag(),theta_miss,weight);	       
		  isSRC_tight=true;
		  //}
		  //}
	    }
	  }
	}
      }

      
  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
      if(isSRC_loose){
	h_num_R->Fill(protons.size()-1,1);
	if(isSRC_tight){
	  h_num_R_tight->Fill(protons.size()-1,1);
	}
      }		
      
      double mom_2 = 0;
      double mom_2_tight = 0;
      int index_2 = -1;
      int index_2_tight = -1;
      for(int j = 0; j < protons.size(); j++){
	if(j==index_L){continue;}
	TVector3 p_R;
	p_R.SetMagThetaPhi(protons[j]->getP(),protons[j]->getTheta(),protons[j]->getPhi());
	if(isSRC_loose){
	  h_p_2->Fill(p_R.Mag(),weight);
	  if(p_R.Mag()>mom_2){
	    index_2=j;
	    mom_2=p_R.Mag();
	  }
	  if(isSRC_tight){
	    h_p_2_tight->Fill(p_R.Mag(),weight);
	    if(p_R.Mag()>mom_2_tight){
	      index_2_tight=j;
	      mom_2_tight=p_R.Mag();
	    }
	  }
	}	
      }
   
  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////
      if(index_2!=-1){
	//cout<< index_2 <<"\n";
	TVector3 p_L;
	p_L.SetMagThetaPhi(protons[index_L]->getP(),protons[index_L]->getTheta(),protons[index_L]->getPhi());
	TVector3 p_1 = p_L - p_q;
	TVector3 p_2;
	p_2.SetMagThetaPhi(protons[index_2]->getP(),protons[index_2]->getTheta(),protons[index_2]->getPhi());
	TVector3 p_rel = p_1-p_2;
	p_rel.SetMag(p_rel.Mag()/2);
	TVector3 p_cm = p_1+p_2;
	double theta_rel = p_1.Angle(p_2) * 180 / M_PI;
	
	//Create new reference frame
	TVector3 vt = p_2.Unit();
	TVector3 vy = p_2.Cross(p_q).Unit();
	TVector3 vx = vt.Cross(vy);

	h_p_2_high->Fill(p_2.Mag(),weight);
	if(p_2.Mag()>0.25){
	  h_p_rel->Fill(p_rel.Mag(),weight);
	  h_p_cm->Fill(p_cm.Mag(),weight);
	  h_p_t_cm->Fill(p_cm.Dot(vt),weight);
	  h_p_y_cm->Fill(p_cm.Dot(vy),weight);
	  h_p_x_cm->Fill(p_cm.Dot(vx),weight);
	  h_theta_rel->Fill(theta_rel,weight);
	  h_p_cm_theta_rel->Fill(p_cm.Mag(),theta_rel,weight);
	}
      }

      if(index_2_tight!=-1){
	TVector3 p_L;
	p_L.SetMagThetaPhi(protons[index_L]->getP(),protons[index_L]->getTheta(),protons[index_L]->getPhi());
	TVector3 p_1 = p_L - p_q;
	TVector3 p_2;
	p_2.SetMagThetaPhi(protons[index_2_tight]->getP(),protons[index_2_tight]->getTheta(),protons[index_2_tight]->getPhi());
	TVector3 p_rel = p_1-p_2;
	p_rel.SetMag(p_rel.Mag()/2);
	TVector3 p_cm = p_1+p_2;
	double theta_rel = p_1.Angle(p_2) * 180 / M_PI;	
	
	//Create new reference frame
	TVector3 vt = p_2.Unit();
	TVector3 vy = p_2.Cross(p_q).Unit();
	TVector3 vx = vt.Cross(vy);

	h_p_2_high_tight->Fill(p_2.Mag(),weight);
	if(p_2.Mag()>0.25){
	  h_p_rel_tight->Fill(p_rel.Mag(),weight);
	  h_p_cm_tight->Fill(p_cm.Mag(),weight);
	  h_p_t_cm_tight->Fill(p_cm.Dot(vt),weight);
	  h_p_y_cm_tight->Fill(p_cm.Dot(vy),weight);
	  h_p_x_cm_tight->Fill(p_cm.Dot(vx),weight);
	  h_theta_rel_tight->Fill(theta_rel,weight);
	  h_p_cm_theta_rel_tight->Fill(p_cm.Mag(),theta_rel,weight);
	}
      }
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

  TCanvas *cvs_electron = new TCanvas("cvs_electron","cvs_electron",pixelx,pixely);
  cvs_electron->Divide(2,3);
  cvs_electron->cd(1);
  h_Vcal_EoP->Draw("colz");
  cvs_electron->cd(2);
  h_Wcal_EoP->Draw("colz");
  cvs_electron->cd(3);
  h_phi_theta->Draw("colz");
  cvs_electron->cd(4);
  h_sector->Draw("colz");
  cvs_electron->cd(5);
  h_P_EoP->Draw("colz");
  cvs_electron->cd(6);
  h_nphe->Draw();

  //  cvs_electron->SaveAs("electron_pid.pdf");


  TCanvas *cvs_electron_kin = new TCanvas("cvs_electron_kin","cvs_electron_kin",pixelx,pixely);

  cvs_electron_kin->Divide(2,3);
  cvs_electron_kin->cd(1);
  h_xB->Draw();
  cvs_electron_kin->cd(2);
  h_QSq->Draw();
  cvs_electron_kin->cd(3);
  h_WSq->Draw();
  cvs_electron_kin->cd(4);
  h_xB_QSq->Draw("colz");
  cvs_electron_kin->cd(5);
  h_xB_WSq->Draw("colz");
  cvs_electron_kin->cd(6);
  h_QSq_WSq->Draw("colz");

  //  cvs_electron_kin->SaveAs("electron_kinematics.pdf");

  TCanvas *cvs_proton = new TCanvas("cvs_proton","cvs_proton",pixelx,pixely);
  cvs_proton->Divide(2,2);
  cvs_proton->cd(1);
  h_theta_L->Draw("colz");
  cvs_proton->cd(2);
  h_theta_Lq->Draw("colz");
  cvs_proton->cd(3);

  cvs_proton->cd(4);



  TCanvas *cvs_lead_proton_1 = new TCanvas("cvs_lead_proton_1","cvs_lead_proton_1",pixelx,pixely);
  cvs_lead_proton_1->Divide(2,3);

  cvs_lead_proton_1->cd(1);
  h_theta_L_FTOF->Draw();
  cvs_lead_proton_1->cd(2);
  h_theta_Lq_FTOF->Draw();
  cvs_lead_proton_1->cd(3);
  h_phi_e_L->Draw();
  cvs_lead_proton_1->cd(4);
  h_mmiss_FTOF->Draw();
  cvs_lead_proton_1->cd(5);
  h_mmiss_phi_e_L->Draw("colz");
  cvs_lead_proton_1->cd(6);
  h_xB_mmiss->Draw("colz");

  cvs_lead_proton_1->cd(7);
  h_pmiss_mmiss->Draw("colz");
  cvs_lead_proton_1->cd(8);
  h_xB_theta_1q->Draw("colz");
  cvs_lead_proton_1->cd(9);
  h_Loq_theta_1q->Draw("colz");
  cvs_lead_proton_1->cd(10);
  h_pmiss_theta_miss->Draw("colz");


  TCanvas *cvs_lead_proton_2 = new TCanvas("cvs_lead_proton_2","cvs_lead_proton_2",pixelx,pixely);
  cvs_lead_proton_2->Divide(2,3);

  cvs_lead_proton_2->cd(1);
  h_xB_theta_1q->Draw("colz");
  cvs_lead_proton_2->cd(2);
  h_Loq_theta_1q->Draw("colz");
  cvs_lead_proton_2->cd(3);
  h_pmiss_theta_miss->Draw("colz");


  TCanvas *cvs_leadsrc_proton = new TCanvas("cvs_leadsrc_proton","cvs_leadsrc_proton",pixelx,pixely);
  cvs_leadsrc_proton->Divide(2,3);

  cvs_leadsrc_proton->cd(1);
  h_pmiss_tight->Draw();
  cvs_leadsrc_proton->cd(2);
  h_mmiss_tight->Draw();
  cvs_leadsrc_proton->cd(3);
  h_pmiss_theta_miss_SRC_tight->Draw("colz");


  TCanvas *cvs_rec_proton = new TCanvas("cvs_rec_proton","cvs_rec_proton",pixelx,pixely);
  cvs_rec_proton->Divide(2,3);
  cvs_rec_proton->cd(1);
  h_num_R_tight->Draw("colz");
  cvs_rec_proton->cd(2);
  h_p_2_tight->Draw("colz");
  cvs_rec_proton->cd(3);

  cvs_proton->cd(4);


  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_p_2_high_tight->Draw();
  myCanvas->cd(2);
  h_p_rel_tight->Draw();
  myCanvas->cd(3);
  h_p_cm_tight->Draw();
  myCanvas->cd(4);
  h_p_t_cm_tight->Draw();
  myCanvas->cd(5);
  h_p_y_cm_tight->Draw();
  myCanvas->cd(6);
  h_p_x_cm_tight->Draw();

  //first pdf page needs ("
  //last pdf page needs )" becareful if you are haveing pdf display issues

  char fileName[100];
  sprintf(fileName,"%s[",argv[2]);

  TLatex text;
  text.SetTextSize(0.05);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",argv[2]);

  text.DrawLatex(0.2,0.9,"(e,e') Candidates:");
  text.DrawLatex(0.2,0.8,"No Cuts");
  myText->Print(fileName,"pdf");
  myText->Clear();

  cvs_electron->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e') Cuts:");
  text.DrawLatex(0.2,0.8,"V_{cal} and W_{cal} > 14 [cm]");
  text.DrawLatex(0.2,0.7,"0.18 < SF < 0.28");
  text.DrawLatex(0.2,0.6,"1 [GeV] < p_{e} < E_{beam}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  cvs_electron_kin->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e'p) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Proton Detected");
  myText->Print(fileName,"pdf");
  myText->Clear();


  cvs_proton->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  text.DrawLatex(0.2,0.7,"FTOF Proton");
  text.DrawLatex(0.2,0.6,"#theta_{p,q}<25^{o}");
  text.DrawLatex(0.2,0.5,"-3 < #chi^{2} PID<3 ");
  text.DrawLatex(0.2,0.4,"#theta_{p} <50^{o}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  cvs_lead_proton_1->Print(fileName,"pdf");
  cvs_lead_proton_2->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD}) Cuts");
  text.DrawLatex(0.2,0.7,"1.5 < Q^{2} [GeV]");
  text.DrawLatex(0.2,0.6,"0.25 [GeV] < p_{miss}");
  text.DrawLatex(0.2,0.5,"0.84 [GeV] < m_{mmiss} < 1.04 [GeV]");
  text.DrawLatex(0.2,0.4,"0.62 < |p|/|q| < 0.96");
  text.DrawLatex(0.2,0.3,"1.2 < x_{B}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  cvs_leadsrc_proton->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC}) Cuts");
  text.DrawLatex(0.2,0.7,"Second Proton Detected");
  myText->Print(fileName,"pdf");
  myText->Clear();

  cvs_rec_proton->Print(fileName,"pdf");

  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC},p_{Rec}) Cuts");
  text.DrawLatex(0.2,0.7,"0.25 [GeV] < p_{Rec}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Print(fileName,"pdf");

  myCanvas->Clear();
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_rel->Draw();
  myCanvas->cd(2);
  h_p_cm_theta_rel->Draw("colz");
  myCanvas->Print(fileName,"pdf");

  myCanvas->Clear();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mmiss_all->Draw();
  myCanvas->cd(2);
  h_mmiss_FTOF1->Draw();
  myCanvas->cd(3);
  h_mmiss_FTOF2->Draw();
  myCanvas->cd(4);
  h_mmiss_noFTOF->Draw();
  myCanvas->Print(fileName,"pdf");

  myCanvas->Clear();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mmiss_FTOF2_mag->Draw("colz");
  myCanvas->cd(2);
  h_mmiss_FTOF2_theta->Draw("colz");
  myCanvas->cd(3);
  h_mmiss_FTOF2_phi->Draw("colz");
  myCanvas->cd(4);
  h_mmiss_FTOF2_chi->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  /*
  myCanvas->Clear();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mmiss_FTOF2_PCAL->Draw("colz");
  myCanvas->cd(2);
  h_mmiss_FTOF2_ECIN->Draw("colz");
  myCanvas->cd(3);
  h_mmiss_FTOF2_ECOUT->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  */
  myCanvas->Clear();
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mmiss_noFTOF_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");

  myCanvas->Clear();
  sprintf(fileName,"%s]",argv[2]);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
  cout<<"Finished making file: "<< argv[1] <<"\n";



}


void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

