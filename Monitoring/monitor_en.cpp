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

  // CND::hits
  auto cnd_hits = config_c12->addBank("CND::hits");
  auto cnd_hit_id = config_c12->getBankOrder(cnd_hits,"id");
  auto cnd_hit_status = config_c12->getBankOrder(cnd_hits,"status");
  auto cnd_hit_trkID = config_c12->getBankOrder(cnd_hits,"trkID");
  auto cnd_hit_sector = config_c12->getBankOrder(cnd_hits,"sector");
  auto cnd_hit_layer = config_c12->getBankOrder(cnd_hits,"layer");
  auto cnd_hit_component = config_c12->getBankOrder(cnd_hits,"component");
  auto cnd_hit_energy = config_c12->getBankOrder(cnd_hits,"energy");
  auto cnd_hit_time = config_c12->getBankOrder(cnd_hits,"time");
  auto cnd_hit_energy_unc = config_c12->getBankOrder(cnd_hits,"energy_unc");
  auto cnd_hit_time_unc = config_c12->getBankOrder(cnd_hits,"time_unc");
  auto cnd_hit_x = config_c12->getBankOrder(cnd_hits,"x");
  auto cnd_hit_y = config_c12->getBankOrder(cnd_hits,"y");
  auto cnd_hit_z = config_c12->getBankOrder(cnd_hits,"z");
  auto cnd_hit_x_unc = config_c12->getBankOrder(cnd_hits,"x_unc");
  auto cnd_hit_y_unc = config_c12->getBankOrder(cnd_hits,"y_unc");
  auto cnd_hit_z_unc = config_c12->getBankOrder(cnd_hits,"z_unc");
  auto cnd_hit_tx = config_c12->getBankOrder(cnd_hits,"tx");
  auto cnd_hit_ty = config_c12->getBankOrder(cnd_hits,"ty");
  auto cnd_hit_tz = config_c12->getBankOrder(cnd_hits,"tz");
  auto cnd_hit_tlength = config_c12->getBankOrder(cnd_hits,"tlength");
  auto cnd_hit_pathlength = config_c12->getBankOrder(cnd_hits,"pathlength");
  auto cnd_hit_indexLadc = config_c12->getBankOrder(cnd_hits,"indexLadc");
  auto cnd_hit_indexRadc = config_c12->getBankOrder(cnd_hits,"indexRadc");
  auto cnd_hit_indexLtdc = config_c12->getBankOrder(cnd_hits,"indexLtdc");
  auto cnd_hit_indexRtdc = config_c12->getBankOrder(cnd_hits,"indexRtdc");

  // CND::clusters
  auto cnd_clusters = config_c12->addBank("CND::clusters");
  auto cnd_clusters_id = config_c12->getBankOrder(cnd_clusters,"id");
  auto cnd_clusters_sector = config_c12->getBankOrder(cnd_clusters,"sector");
  auto cnd_clusters_layer = config_c12->getBankOrder(cnd_clusters,"layer");
  auto cnd_clusters_component = config_c12->getBankOrder(cnd_clusters,"component");
  auto cnd_clusters_nhits = config_c12->getBankOrder(cnd_clusters,"nhits");
  auto cnd_clusters_energy = config_c12->getBankOrder(cnd_clusters,"energy");
  auto cnd_clusters_x = config_c12->getBankOrder(cnd_clusters,"x");
  auto cnd_clusters_y = config_c12->getBankOrder(cnd_clusters,"y");
  auto cnd_clusters_z = config_c12->getBankOrder(cnd_clusters,"z");
  auto cnd_clusters_time = config_c12->getBankOrder(cnd_clusters,"time");
  auto cnd_clusters_status = config_c12->getBankOrder(cnd_clusters,"status");




  // REC::ScintExtras
  

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


  /////////////////////////////////////
  //Electron Pid and Vertex
  /////////////////////////////////////
  TH1D * h_nphe = new TH1D("nphe","#Photo-electrons in HTCC;#Photo-electrons;Counts",40,0,40);
  hist_list_1.push_back(h_nphe);

  TH1D * h_vtz_e = new TH1D("vtz_e","Electron Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_e);


  
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
  //Neutral Particle Information
  /////////////////////////////////////
  TH2D * h_mom_beta_neutron_ECAL_zoom = new TH2D("mom_beta_neutron_ECAL","p vs. #beta ECAL (Neutrons);p;#beta",100,0,2.5,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_neutron_ECAL_zoom);
  TH2D * h_mom_beta_photon_ECAL_zoom = new TH2D("mom_beta_photon_ECAL","p vs. #beta ECAL (Photons);p;#beta",100,0,2.5,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_photon_ECAL_zoom);


  TH2D * h_mom_beta_neutron_CND_zoom = new TH2D("mom_beta_neutron_CND","p vs. #beta CND (Neutrons);p;#beta",100,0,2.5,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_neutron_CND_zoom);
  TH2D * h_mom_beta_photon_CND_zoom = new TH2D("mom_beta_photon_CND","p vs. #beta CND (Photons);p;#beta",100,0,2.5,100,0.0,1.2);
  hist_list_2.push_back(h_mom_beta_photon_CND_zoom);



  /////////////////////////////////////
  //ECAL Neutron Information
  /////////////////////////////////////
  TH2D * h_nsize_ECAL = new TH2D("nsize_ECAL","Number of Neutrons in ECAL;Number of PID 2112;Neutrons with Valid Momentum",10,0,10,10,0,10);
  hist_list_2.push_back(h_nsize_ECAL);

  TH1D * h_vtz_n_ECAL = new TH1D("vtz_n_ECAL","Neutron Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_n_ECAL);

  TH1D * h_phi_e_n_ECAL = new TH1D("phi_e_n_ECAL","|#phi_{e} - #phi_{n}|;|#phi_{e} - #phi_{n}|,Counts",180,0,180);
  hist_list_1.push_back(h_phi_e_n_ECAL);


  TH1D * h_theta_n_ECAL = new TH1D("theta_n_ECAL","#theta_{neutron};#theta_{neutron};Counts",100,0,60);
  hist_list_1.push_back(h_theta_n_ECAL);
  TH1D * h_theta_nq_ECAL = new TH1D("theta_nq_ECAL","#theta_{nq};#theta_{nq};Counts",100,0,80);
  hist_list_1.push_back(h_theta_nq_ECAL);
  TH2D * h_phi_theta_n_ECAL = new TH2D("phi_theta_n_ECAL","#phi_{n} vs. #theta_{n} ;#phi_{n};#theta_{n}",100,-180,180,100,3,40);
  hist_list_2.push_back(h_phi_theta_n_ECAL);
  TH2D * h_mom_beta_n_ECAL = new TH2D("mom_beta_n_ECAL","p_{n} vs. #beta_{n} ;p_{n};#beta_{n}",100,0.1,3,100,0.1,1);
  hist_list_2.push_back(h_mom_beta_n_ECAL);
  TH1D * h_timediff_n_ECAL = new TH1D("timediff_n_ECAL","ToF-ToF_{|n|} ;ToF-ToF_{|n|};Counts",100,-0.1,0.1);
  hist_list_1.push_back(h_timediff_n_ECAL);

  TH2D * h_mom_theta_n_ECAL[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_theta_n_%d",i+1);
    sprintf(temp_title,"p_{n} vs. #theta_{n} Sector=%d;Momentum;Theta",i+1);
    h_mom_theta_n_ECAL[i] = new TH2D(temp_name,temp_title,100,0.1,2,100,3,40);
    hist_list_2.push_back(h_mom_theta_n_ECAL[i]);
  }


  // extra bank info
  TH2D * h_dedx_ECAL = new TH2D("dedx_ECAL","Energy Deposition;#beta#gamma = p/m;dE/dx (MeV/cm)",100,0.1,1.5,100,1,20);
  hist_list_2.push_back(h_dedx_ECAL);
  TH2D * h_cluster_ECAL = new TH2D("cluster_ECAL","Cluster Size and Multiplicity;Cluster Size;Cluster Multiplicity",10,0,10,10,0,10);
  hist_list_2.push_back(h_cluster_ECAL);

  /////////////////////////////////////
  //CND Neutron Information
  /////////////////////////////////////
  TH2D * h_nsize_CND = new TH2D("nsize_CND","Number of Neutrons in ECAL;Number of Neutrons;Neutrons with Valid Momentum",10,0,10,10,0,10);
  hist_list_2.push_back(h_nsize_CND);

  TH1D * h_vtz_n_CND = new TH1D("vtz_n_CND","Neutron Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_n_CND);

  TH1D * h_phi_e_n_CND = new TH1D("phi_e_n_CND","|#phi_{e} - #phi_{n}|;|#phi_{e} - #phi_{n}|,Counts",180,0,180);
  hist_list_1.push_back(h_phi_e_n_CND);


  TH1D * h_theta_n_CND = new TH1D("theta_n_CND","#theta_{neutron};#theta_{neutron};Counts",100,20,140);
  hist_list_1.push_back(h_theta_n_CND);
  TH1D * h_theta_nq_CND = new TH1D("theta_nq_CND","#theta_{nq};#theta_{nq};Counts",100,0,80);
  hist_list_1.push_back(h_theta_nq_CND);
  TH2D * h_phi_theta_n_CND = new TH2D("phi_theta_n_CND","#phi_{n} vs. #theta_{n} ;#phi_{n};#theta_{n}",48,-180,180,100,35,145);
  hist_list_2.push_back(h_phi_theta_n_CND);
  TH2D * h_mom_beta_n_CND = new TH2D("mom_beta_n_CND","p_{n} vs. #beta_{n} ;p_{n};#beta_{n}",100,0.1,2,100,0.1,1);
  hist_list_2.push_back(h_mom_beta_n_CND);
  TH1D * h_timediff_n_CND = new TH1D("timediff_n_CND","ToF-ToF_{|n|} ;ToF-ToF_{|n|};Counts",100,-0.1,0.1);
  hist_list_1.push_back(h_timediff_n_CND);
  TH2D * h_mom_theta_n_CND = new TH2D("mom_theta_n_CND","p_{n} vs. #theta_{n} ;p_{n};#theta_{n}",100,0.2,1.5,100,35,145);
  hist_list_2.push_back(h_mom_theta_n_CND);
  TH2D * h_phi_mom_n_CND = new TH2D("phi_mom_n_CND","#phi_{n} vs p_{n};#phi_{n};p_{n}",48,-180,180,100,0.2,1.5);
  hist_list_2.push_back(h_phi_mom_n_CND);


  // extra bank info REC::ScintExtras
  TH2D * h_dedx_CND = new TH2D("dedx_CND","Energy Deposition;#beta#gamma = p/m;dE/dx (MeV/cm)",100,0.1,1.5,100,1,20);
  hist_list_2.push_back(h_dedx_CND);
  TH2D * h_cluster_CND = new TH2D("cluster_CND","Cluster Size and Multiplicity;Cluster Size;Cluster Multiplicity",10,0,10,10,0,10);
  hist_list_2.push_back(h_cluster_CND);



  // extra bank info CND::hits
  TH1D * h_sec = new TH1D("hit_sector","Sector of CND hit",48,0,48);
  hist_list_1.push_back(h_sec);
  TH1D * h_lay = new TH1D("hit_layer","Layer of CND hit",48,0,48);
  hist_list_1.push_back(h_lay);
  TH2D * h_sec_lay = new TH2D("hit_sector_layer","Hit Sector vs Layer;Sector;Layer",24,0,24,5,0,4);
  hist_list_2.push_back(h_sec_lay);
  //TH1D * h_Eunc = new TH1D("Eunc","Energy deposition uncertainty of hit (MeV)",100,0,50);
  //hist_list_1.push_back(h_Eunc);
  //TH1D * h_tunc = new TH1D("tunc","Hit time uncertainty (ns)",100,0,50);
  //hist_list_1.push_back(h_tunc);
  TH1D * h_x = new TH1D("hit_x","Hit position x;x (cm);Counts",50,-50,50);
  hist_list_1.push_back(h_x);
  TH1D * h_x_unc = new TH1D("hit_x_unc","Hit position x uncertainty;x uncertainty (cm);Counts",50,0,5);
  hist_list_1.push_back(h_x_unc);
  TH1D * h_y = new TH1D("hit_y","Hit position y;y (cm);Counts",50,-50,50);
  hist_list_1.push_back(h_y);
  TH1D * h_y_unc = new TH1D("hit_y_unc","Hit position y uncertainty;y uncertainty (cm);Counts",50,0,5);
  hist_list_1.push_back(h_y_unc);
  TH1D * h_z = new TH1D("hit_z","Hit position z;z (cm);Counts",50,-50,50);
  hist_list_1.push_back(h_z);
  TH1D * h_z_unc = new TH1D("hit_z_unc","Hit position z uncertainty;z uncertainty (cm);Counts",50,0,5);
  hist_list_1.push_back(h_z_unc);


  TH2D * h_hits_xy = new TH2D("hits_xy","Hit positions in CND;x (cm);y (cm)",100,-50,50,100,-50,50);
  hist_list_2.push_back(h_hits_xy);
  TH1D * h_hitE = new TH1D("hitE","CND Hit Energy Deposition;Hit Energy (MeV);Counts",100,0,100);
  hist_list_1.push_back(h_hitE);
  TH2D * h_xx = new TH2D("xx","X position from CND vs CVT;X coordinate of CND hit;X coordinate of hit from CVT info",100,-50,50,100,-50,50);
  hist_list_2.push_back(h_xx);
  TH2D * h_path_t = new TH2D("patht","Track Path: Entrance to Exit vs. Vertex;Path length from entrance to exit through hit;Path length from vertex",100,0,30,100,0,120);
  hist_list_2.push_back(h_path_t);
  //TH2D * h_adc_LR = new TH2D("adc_LR","ADC Index: Left vs Right;L adc index;R adc index",100,0,50,100,0,50);
  //hist_list_2.push_back(h_adc_LR);
  //TH2D * h_tdc_LR = new TH2D("tdc_LR","ADC Index: Left vs Right;L tdc index;R tdc index",100,0,50,100,0,50);
  //hist_list_2.push_back(h_tdc_LR);

  
  // extra bank info CND::clusters
  TH1D * h_cl_nhits = new TH1D("cl_nhits","Number of Hits in CND Cluster;Number of Hits;Counts",15,0,15);
  hist_list_1.push_back(h_cl_nhits);


  /////////////////////////////////////
  //CND Neutron Compared to Pmiss
  /////////////////////////////////////
  TH2D * h_pn_angles = new TH2D("pn_angles","Angular Separation between Proton and Neutron (before cut);#Delta#phi = #phi_{p} - #phi_{n};#Delta#theta = #theta_{p} - #theta_{n}",360,-180,180,360,-180,180);
  hist_list_2.push_back(h_pn_angles);
  TH2D * h_pn_pmiss = new TH2D("pn_pmiss","Neutron Momentum vs Missing Momentum;Missing Momentum (GeV/c);Neutron Momentum (GeV/c)",100,0,2,100,0,2);
  hist_list_2.push_back(h_pn_pmiss);
  TH1D * h_cos0 = new TH1D("cos0","Angle between p_{n} and p_{miss};cos #theta_{pn,pmiss};Counts",100,-1.1,1.1);
  hist_list_1.push_back(h_cos0);
  TH2D * h_theta_phi_p = new TH2D("theta_phi_p","Proton Theta vs Phi;#phi_{p};#theta_{n}",360,-180,180,180,0,180);
  hist_list_2.push_back(h_theta_phi_p);



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


  /////////////////////////////////////
  //Electron Pid
  /////////////////////////////////////
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
      h_xB_QSq->Fill(xB,QSq,weight);
      h_xB_WSq->Fill(xB,WSq,weight);
      h_QSq_WSq->Fill(QSq,WSq,weight);
      h_mom_theta[esector-1]->Fill(p_e.Mag(),theta_e,weight);

  /////////////////////////////////////
  //Neutral Particle Information
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

	if(PCAL || ECin || ECout){

          int PID = allParticles[j]->getPid();
          if(PID==2112){
            h_mom_beta_neutron_ECAL_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
          }  
          if(PID==22){
            h_mom_beta_photon_ECAL_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);
          }
	}
	if(CND1 || CND2 || CND3){
 	  
	  int PID = allParticles[j]->getPid();
	  if(PID==2112){
	    h_mom_beta_neutron_CND_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }
	  if(PID==22){
	    h_mom_beta_photon_CND_zoom->Fill(allParticles[j]->getP(),allParticles[j]->par()->getBeta(),weight);  
	  }
	}

      }
  /////////////////////////////////////
  //ECAL and CND Neutron Information
  /////////////////////////////////////
      int nn_ECAL = 0; int nn_CND = 0; int nn_ECAL_good = 0; int nn_CND_good = 0;

      for(int j = 0; j < neutrons.size(); j++){
	TVector3 p_n;
	p_n.SetMagThetaPhi(neutrons[j]->getP(),neutrons[j]->getTheta(),neutrons[j]->getPhi());
	double E_n = sqrt(mN*mN + p_n.Mag2());
	double theta_n = p_n.Theta() * 180 / M_PI;
	double phi_n = p_n.Phi() * 180 / M_PI;
	double theta_nq = p_n.Angle(p_q) * 180 / M_PI;
	double beta_n = neutrons[j]->par()->getBeta();
	double phi_diff = get_phi_diff(p_e,p_n);

	double path_n = neutrons[j]->getPath();
	double beta_frommom_n = p_n.Mag()/E_n;
	double time_frommom_n = path_n / (c*beta_frommom_n);
	double time_frombeta_n = path_n / (c*beta_n);

	bool PCAL = (neutrons[j]->cal(clas12::PCAL)->getDetector() == 7);
	bool ECin = (neutrons[j]->cal(clas12::ECIN)->getDetector() == 7);
        bool ECout = (neutrons[j]->cal(clas12::ECOUT)->getDetector() == 7);
        bool CND1 = (neutrons[j]->sci(clas12::CND1)->getDetector() == 3);
        bool CND2 = (neutrons[j]->sci(clas12::CND2)->getDetector() == 3);
        bool CND3 = (neutrons[j]->sci(clas12::CND3)->getDetector() == 3);

        if(PCAL || ECin || ECout)  { nn_ECAL += 1; }
        if(CND)                    { nn_CND += 1; }

        // "good neutrons" start here
        if (beta_frommom_n==0) {continue;}
        if (theta_n==0) {continue;}
        if (phi_n==0) {continue;}
        if (p_n.Mag()==0) {continue;}

	if(PCAL || ECin || ECout){
          nn_ECAL_good += 1;
	  h_phi_e_n_ECAL->Fill(phi_diff,weight);

	  h_theta_n_ECAL->Fill(theta_n,weight);
	  h_theta_nq_ECAL->Fill(theta_nq,weight);	
	  h_phi_theta_n_ECAL->Fill(phi_n,theta_n,weight);
	  h_mom_beta_n_ECAL->Fill(p_n.Mag(),beta_n,weight);
	  h_timediff_n_ECAL->Fill(time_frommom_n-time_frombeta_n,weight);

	  h_mom_theta_n_ECAL[neutrons[j]->getSector()-1]->Fill(p_n.Mag(),theta_n,weight);
	}

	if(CND){
          if (theta_n<40 || theta_n>140) {continue;}

          nn_CND_good += 1;
	  h_phi_e_n_CND->Fill(phi_diff,weight);

	  h_theta_n_CND->Fill(theta_n,weight);
	  h_theta_nq_CND->Fill(theta_nq,weight);	
	  h_phi_theta_n_CND->Fill(phi_n,theta_n,weight);
	  h_mom_beta_n_CND->Fill(p_n.Mag(),beta_n,weight);
	  h_timediff_n_CND->Fill(time_frommom_n-time_frombeta_n,weight);
	  
	  h_mom_theta_n_CND->Fill(p_n.Mag(),theta_n,weight);
          h_phi_mom_n_CND->Fill(phi_n,p_n.Mag(),weight);


          
          // fill histos with CND::hits bank info

          // extra bank info
          double dedx = neutrons[j]->sci(clas12::CND1)->getDedx() + neutrons[j]->sci(clas12::CND2)->getDedx() + neutrons[j]->sci(clas12::CND3)->getDedx();
          double csize = neutrons[j]->sci(clas12::CND1)->getSize() + neutrons[j]->sci(clas12::CND2)->getSize() + neutrons[j]->sci(clas12::CND3)->getSize();
          double layermult = neutrons[j]->sci(clas12::CND1)->getLayermulti() + neutrons[j]->sci(clas12::CND2)->getLayermulti() + neutrons[j]->sci(clas12::CND3)->getLayermulti();


          // fill histos with extra bank info
          h_dedx_CND->Fill(p_n.Mag()/mN,dedx,weight);
          h_cluster_CND->Fill(csize,layermult,weight);
	}


      }
      h_nsize_ECAL->Fill(nn_ECAL,nn_ECAL_good,weight);
      h_nsize_CND->Fill(nn_CND,nn_CND_good,weight);


      // CND::hits banks info
      for (auto iRow=0; iRow < c12->getBank(cnd_hits)->getRows(); iRow++){
        // define variables
        int sector = c12->getBank(cnd_hits)->getInt(cnd_hit_sector,iRow);
        int layer = c12->getBank(cnd_hits)->getInt(cnd_hit_layer,iRow);
        double E = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy,iRow);
        double time = c12->getBank(cnd_hits)->getFloat(cnd_hit_time,iRow);
        double E_unc = c12->getBank(cnd_hits)->getFloat(cnd_hit_energy_unc,iRow); // zero
        double t_unc = c12->getBank(cnd_hits)->getFloat(cnd_hit_time_unc,iRow); // zero
        double x = c12->getBank(cnd_hits)->getFloat(cnd_hit_x,iRow);
        double y = c12->getBank(cnd_hits)->getFloat(cnd_hit_y,iRow);
        double z = c12->getBank(cnd_hits)->getFloat(cnd_hit_z,iRow);
        double x_unc = c12->getBank(cnd_hits)->getFloat(cnd_hit_x_unc,iRow);
        double y_unc = c12->getBank(cnd_hits)->getFloat(cnd_hit_y_unc,iRow);
        double z_unc = c12->getBank(cnd_hits)->getFloat(cnd_hit_z_unc,iRow);
        double x_cvt = c12->getBank(cnd_hits)->getFloat(cnd_hit_tx,iRow);
        double y_cvt = c12->getBank(cnd_hits)->getFloat(cnd_hit_ty,iRow);
        double z_cvt = c12->getBank(cnd_hits)->getFloat(cnd_hit_tz,iRow);
        double tlength = c12->getBank(cnd_hits)->getFloat(cnd_hit_tlength,iRow);
        double pathlength = c12->getBank(cnd_hits)->getFloat(cnd_hit_pathlength,iRow);
        int indexLadc = c12->getBank(cnd_hits)->getFloat(cnd_hit_indexLadc,iRow); // zero
        int indexRadc = c12->getBank(cnd_hits)->getFloat(cnd_hit_indexRadc,iRow); // zero
        int indexLtdc = c12->getBank(cnd_hits)->getFloat(cnd_hit_indexLtdc,iRow); // zero
        int indexRtdc = c12->getBank(cnd_hits)->getFloat(cnd_hit_indexRtdc,iRow); // zero


        // fill histograms
        h_sec_lay->Fill(sector,layer,weight);

        h_hits_xy->Fill(x,y,weight);
        h_hitE->Fill(E,weight);
        h_xx->Fill(x,x_cvt,weight);
        h_x->Fill(x,weight);
        h_y->Fill(y,weight);
        h_z->Fill(z,weight);
        h_x_unc->Fill(x_unc,weight);
        h_y_unc->Fill(y_unc,weight);
        h_z_unc->Fill(z_unc,weight);
        //h_Eunc->Fill(E_unc,weight);
        //h_tunc->Fill(t_unc,weight);
        h_path_t->Fill(tlength,pathlength,weight);
        //h_adc_LR->Fill(indexLadc,indexRadc,weight);
        //h_tdc_LR->Fill(indexLtdc,indexRtdc,weight);
      }

      // CND::clusters banks info
      for (auto iRow=0; iRow < c12->getBank(cnd_clusters)->getRows(); iRow++){
        int nhits = c12->getBank(cnd_clusters)->getInt(cnd_clusters_nhits,iRow);
        h_cl_nhits->Fill(nhits,weight);
      }




  /////////////////////////////////////
  //Compare CND Neutrons to Pmiss
  /////////////////////////////////////

  // only events with 1 proton and 1 neutron
  if (protons.size()!=1) {continue;}
  if (neutrons.size()!=1) {continue;}

  // skip events with particles other than proton, neutrons, electron, photon
  bool trash = 0;
  for(int j = 0; j < allParticles.size(); j ++){
    int PID = allParticles[j]->getPid();
    if (PID!=2212 && PID!=2112 && PID!=22 && PID!=11){
      trash = 1;
    }
  }
  if (trash) {continue;}

  TVector3 p_p;
  p_p.SetMagThetaPhi(protons[0]->getP(),protons[0]->getTheta(),protons[0]->getPhi());
  TVector3 p_miss = p_p - p_q;
  double theta_pm = p_miss.Theta()*180 / M_PI;


  for(int j = 0; j < neutrons.size(); j++){
	TVector3 p_n;
	p_n.SetMagThetaPhi(neutrons[j]->getP(),neutrons[j]->getTheta(),neutrons[j]->getPhi());
	double E_n = sqrt(mN*mN + p_n.Mag2());
	double theta_n = p_n.Theta() * 180 / M_PI;
	double phi_n = p_n.Phi() * 180 / M_PI;
	double beta_frommom_n = p_n.Mag()/E_n;
        double cos0 = p_miss.Dot(p_n) / (p_miss.Mag()*p_n.Mag());

        double dphi_pn = protons[0]->getPhi()*180./M_PI - phi_n;
        double dtheta_pn = protons[0]->getTheta()*180./M_PI - theta_n;

        bool CND1 = (neutrons[j]->sci(clas12::CND1)->getDetector() == 3);
        bool CND2 = (neutrons[j]->sci(clas12::CND2)->getDetector() == 3);
        bool CND3 = (neutrons[j]->sci(clas12::CND3)->getDetector() == 3);

        // cuts to get good neutrons
        if (beta_frommom_n==0) {continue;}
        if (theta_n==0) {continue;}
        if (phi_n==0) {continue;}
        if (p_n.Mag()<0.2) {continue;}
        if (theta_n<40 || theta_n>140) {continue;}

        // require pmiss to be as expected
        if (p_miss.Mag()<0.2) {continue;}
        if (theta_pm<40 || theta_pm>140) {continue;}

        // "good neutrons" start here
        if(CND){
          h_pn_angles->Fill( dphi_pn , dtheta_pn , weight);
          if (dtheta_pn>-15) {continue;}
          h_pn_pmiss->Fill(p_miss.Mag(),p_n.Mag(),weight);
          h_cos0->Fill(cos0,weight);
          h_theta_phi_p->Fill(p_p.Theta()*180./M_PI,p_p.Phi()*180./M_PI,weight);
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
  
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_phi_theta->Draw("colz");
  myCanvas->cd(2);
  h_sector->Draw();
  myCanvas->cd(3);
  h_nphe->Draw();
  myCanvas->cd(4);
  h_vtz_e->Draw();
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

  /////////////////////////////////////
  //Neutral Particle Detected
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts with neutral particle");
  myText->Print(fileName,"pdf");  
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mom_beta_neutron_ECAL_zoom->Draw("colz");
  myCanvas->cd(2);
  h_mom_beta_neutron_CND_zoom->Draw("colz");
  myCanvas->cd(3);
  h_mom_beta_photon_ECAL_zoom->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_photon_CND_zoom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //ECAL Neutron Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'n) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in ECAL");
  text.DrawLatex(0.2,0.6,"Neutron #beta, #theta, #phi, p all nonzero");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_nsize_ECAL->Draw("colz");
  myCanvas->cd(2);
  h_phi_e_n_ECAL->Draw();
  myCanvas->cd(3);
  h_timediff_n_ECAL->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_n_ECAL->Draw();
  myCanvas->cd(2);
  h_theta_nq_ECAL->Draw();
  myCanvas->cd(3);
  h_phi_theta_n_ECAL->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_n_ECAL[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //CND Neutron Information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'n) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in CND");
  text.DrawLatex(0.2,0.6,"Neutron #beta, #theta, #phi, p all nonzero");
  text.DrawLatex(0.2,0.5,"Neuron #theta: min=40 deg, max=140 deg");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_nsize_CND->Draw("colz");
  myCanvas->cd(2);
  h_phi_e_n_CND->Draw();
  myCanvas->cd(3);
  h_timediff_n_CND->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_n_CND->Draw();
  myCanvas->cd(2);
  h_theta_nq_CND->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_phi_theta_n_CND->Draw("colz");
  myCanvas->cd(2);
  h_mom_theta_n_CND->Draw("colz");
  myCanvas->cd(3);
  h_phi_mom_n_CND->Draw("colz");
  myCanvas->cd(4);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();






  /////////////////////////////////////
  //CND::hits
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'n) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in CND");
  text.DrawLatex(0.2,0.6,"Data from CND::hits");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_sec_lay->Draw("colz");
  myCanvas->cd(2);
  h_hits_xy->Draw("colz");
  myCanvas->cd(3);
  h_hitE->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_xx->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1),
  h_x->Draw();
  myCanvas->cd(2);
  h_x_unc->Draw();
  myCanvas->cd(3);
  h_y->Draw();
  myCanvas->cd(4);
  h_y_unc->Draw();
  myCanvas->cd(5);
  h_z->Draw();
  myCanvas->cd(6);
  h_z_unc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_path_t->Draw("colz");
  myCanvas->cd(2);
  myCanvas->cd(3);
  //h_adc_LR->Draw("colz");
  myCanvas->cd(4);
  //h_tdc_LR->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  /////////////////////////////////////
  //REC::ScintExtras
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'n) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in CND");
  text.DrawLatex(0.2,0.6,"Data from REC::ScintExtras");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  myCanvas->SetLogy();
  h_dedx_CND->Draw("colz");
  myCanvas->cd(2);
  h_cluster_CND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  /////////////////////////////////////
  //CND Neutron Compared to Pmiss
  /////////////////////////////////////
    myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'n) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e') Cuts");
  text.DrawLatex(0.2,0.7,"Neutron Detected in CND");
  text.DrawLatex(0.2,0.6,"Neutron #beta, #theta, #phi, p all nonzero");
  text.DrawLatex(0.2,0.5,"Neutron #theta: min=40 deg, max=140 deg");
  text.DrawLatex(0.2,0.4,"p_{n}, p_{miss} > 0.2 GeV/c");
  text.DrawLatex(0.2,0.3,"#theta_{pmiss}: min=40 deg, max=140 deg");
  text.DrawLatex(0.2,0.2,"Neutron Momentum Compared to p_{miss}");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pn_pmiss->Draw("colz");
  myCanvas->cd(2);
  h_cos0->Draw();
  myCanvas->cd(3);
  h_pn_angles->Draw("colz");
  myCanvas->cd(4);
  h_theta_phi_p->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


/*  /////////////////////////////////////
  //CND::clusters
  /////////////////////////////////////
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_cl_nhits->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/


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

