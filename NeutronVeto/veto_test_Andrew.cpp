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
#include <TDatabasePDG.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "clas12reader.h"
#include "HipoChain.h"
#include "neutron-veto/veto_functions.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;
using namespace TMVA;
auto db=TDatabasePDG::Instance();

const double c = 29.9792458;
const double mN = 0.939;
const double mP = 0.938;
const double mD = 1.8756;

void printProgress(double percentage);

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/input.hipo> \n";
}


int main(int argc, char ** argv)
{

  if(argc < 6)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  
  bool isMC = false;
  if(atoi(argv[1]) == 1){isMC=true;}

  clas12ana clasAna;
  clasAna.printParams();

  double Ebeam = atof(argv[2]);
  
  TFile * outFile = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];

  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});




  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  chain.db()->turnOffQADB();
  
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());

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
  TH2D * h_eangles = new TH2D("eangles","Electron Angular Distribution",360,-180,180,180,0,180);
  hist_list_2.push_back(h_eangles);
  TH1D * h_mmiss = new TH1D("mmiss","m_{miss};m_{miss}",100,0.0,1.8);
  hist_list_1.push_back(h_mmiss);
  TH2D * h_xB_mmiss = new TH2D("xB_mmiss","x_{B} vs. m_{miss};x_{B};m_{miss}",100,1.2,2.0,100,0.0,1.8);
  hist_list_1.push_back(h_xB_mmiss);



  TH2D * h_mom_theta[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_theta_%d",i+1);
    sprintf(temp_title,"p_{e} vs. #theta_{e} Sector=%d;Momentum;Theta",i+1);
    h_mom_theta[i] = new TH2D(temp_name,temp_title,100,0,7,100,5,40);
    hist_list_2.push_back(h_mom_theta[i]);
  }


  /////////////////////////////////////
  //CND Neutron Information
  /////////////////////////////////////
  TH2D * h_nsize_CND = new TH2D("nsize_CND","Number of Neutrons in CND;Number of Neutrons;Neutrons with Valid Momentum",10,0,10,10,0,10);
  hist_list_2.push_back(h_nsize_CND);



  TH1D * h_theta_n_CND = new TH1D("theta_n_CND","#theta_{neutron};#theta_{neutron};Counts",100,20,140);
  hist_list_1.push_back(h_theta_n_CND);
  TH2D * h_phi_theta_n_CND = new TH2D("phi_theta_n_CND","#phi_{n} vs. #theta_{n} ;#phi_{n};#theta_{n}",48,-180,180,100,35,145);
  hist_list_2.push_back(h_phi_theta_n_CND);
  TH2D * h_mom_beta_n_CND = new TH2D("mom_beta_n_CND","p_{n} vs. #beta_{n} ;p_{n};#beta_{n}",100,0.1,2,100,0.1,1);
  hist_list_2.push_back(h_mom_beta_n_CND);

  TH2D * h_mom_theta_n_CND = new TH2D("mom_theta_n_CND","p_{n} vs. #theta_{n} ;p_{n};#theta_{n}",100,0.2,1.5,100,35,145);
  hist_list_2.push_back(h_mom_theta_n_CND);
  TH2D * h_phi_mom_n_CND = new TH2D("phi_mom_n_CND","#phi_{n} vs p_{n};#phi_{n};p_{n}",48,-180,180,100,0.2,1.5);
  hist_list_2.push_back(h_phi_mom_n_CND);


  TH1D * h_mva_n_CND = new TH1D("mva_n_CND","MVA Value of Neutron;MVA;Counts",100,0,1);
  hist_list_1.push_back(h_mva_n_CND);



  /////////////////////////////////////
  //CND Neutron Signal
  /////////////////////////////////////
  TH2D * h_res_mom_theta_sig = new TH2D("res_mom_theta_sig","Momentum Difference vs Angular Difference (Signal);(p_{miss}-p_{n})/p_{miss};#theta_{miss,n}",100,-3,1,100,0,180);
  hist_list_2.push_back(h_res_mom_theta_sig);
  TH2D * h_pn_angles_sig = new TH2D("pn_angles_sig","Angular Separation between Proton and Neutron (Signal);#Delta#phi = #phi_{p} - #phi_{n};#Delta#theta = #theta_{p} - #theta_{n}",360,-180,180,360,-180,180);
  hist_list_2.push_back(h_pn_angles_sig);
  TH2D * h_n_thetaphi_sig = new TH2D("n_thetaphi_sig","Angular Distribution of Neutrons (Signal);#phi (deg);#theta (deg)",48,-180,180,100,40,140);
  hist_list_2.push_back(h_n_thetaphi_sig);
  TH2D * h_pn_pmiss_sig = new TH2D("pn_pmiss_sig","Neutron Measured Momentum vs Predicted Momentum (Signal);Neutron Predicted Momentum (GeV/c);Neutron Momentum (GeV/c)",100,0,1.2,100,0,1.2);
  hist_list_2.push_back(h_pn_pmiss_sig);
  TH1D * h_cos0_sig = new TH1D("cos0_sig","Angle between p_{n} and p_{miss} (Signal);cos #theta_{pn,pmiss};Counts",20,-1.1,1.1);
  hist_list_1.push_back(h_cos0_sig);
  TH1D * h_pmiss_sig = new TH1D("pmiss_sig","p_{miss} (Signal);p_{miss};Counts",20,0.3,1.0);
  hist_list_1.push_back(h_pmiss_sig);
  TH1D * h_tof_sig = new TH1D("tof_sig","Time-of-Flight (Signal);TOF;Counts",60,-10,50);
  hist_list_1.push_back(h_tof_sig);
  TH1D * h_mom_sig = new TH1D("mom_sig","Measured Neutron Momentum (Signal);Momentum (GeV/c);Counts",100,0,1.5);
  hist_list_1.push_back(h_mom_sig);
  TH2D * h_ptheta_sig = new TH2D("ptheta_sig","Momentum vs #theta (Signal);Neutron Momentum (GeV/c);#theta_{n} (deg)",100,0.2,1.2,10,40,140);
  hist_list_1.push_back(h_ptheta_sig);
  TH1D * h_Edep_sig = new TH1D("edep_sig","Energy Deposition (Signal);Energy Deposition (MeVee);Counts",100,0,100);
  hist_list_1.push_back(h_Edep_sig);
  TH2D * h_Edep_p_sig = new TH2D("edep_p_sig","Energy Deposition vs Momentum;Momentum (GeV/c);Energy Deposition (MeVee)",100,0,1.5,100,0,100);
  hist_list_2.push_back(h_Edep_p_sig);
  TH2D * h_Edep_theta_sig = new TH2D("edep_theta_sig","Energy Deposition vs #theta_{n};#theta_{n};Energy Deposition (MeVee)",180,0,180,100,0,100);
  hist_list_2.push_back(h_Edep_theta_sig);
  TH2D * h_Edep_phi_sig = new TH2D("edep_phi_sig","Energy Deposition vs #phi;#phi;Energy Deposition (MeVee)",48,-180,180,100,0,100);
  hist_list_2.push_back(h_Edep_phi_sig);
  TH2D * h_Edep_thetapn_sig = new TH2D("edep_thetapn_sig","Energy Deposition vs #theta_{pn};theta_{pn};Energy Deposition (MeVee)",180,0,180,100,0,100);
  hist_list_2.push_back(h_Edep_thetapn_sig);


  /////////////////////////////////////
  //CND Neutron Background
  /////////////////////////////////////
  TH2D * h_res_mom_theta_back = new TH2D("res_mom_theta_back","Momentum Difference vs Angular Difference (Background);(p_{miss}-p_{n})/p_{miss};#theta_{miss,n}",100,-3,1,100,0,180);
  hist_list_2.push_back(h_res_mom_theta_back);
  TH2D * h_pn_angles_back = new TH2D("pn_angles_back","Angular Separation between Proton and Neutron (Background);#Delta#phi = #phi_{p} - #phi_{n};#Delta#theta = #theta_{p} - #theta_{n}",360,-180,180,360,-180,180);
  hist_list_2.push_back(h_pn_angles_back);
  TH2D * h_n_thetaphi_back = new TH2D("n_thetaphi_back","Angular Distribution of Neutrons (Background);#phi (deg);#theta (deg)",48,-180,180,100,40,140);
  hist_list_2.push_back(h_n_thetaphi_back);
  TH2D * h_pn_pmiss_back = new TH2D("pn_pmiss_back","Neutron Measured Momentum vs Predicted Momentum (Background);Neutron Predicted Momentum (GeV/c);Neutron Momentum (GeV/c)",100,0,1.2,100,0,1.2);
  hist_list_2.push_back(h_pn_pmiss_back);
  TH1D * h_cos0_back = new TH1D("cos0_back","Angle between p_{n} and p_{miss} (Background);cos #theta_{pn,pmiss};Counts",20,-1.1,1.1);
  hist_list_1.push_back(h_cos0_back);
  TH1D * h_pmiss_back = new TH1D("pmiss_back","p_{miss} (Background);p_{miss};Counts",20,0.3,1.0);
  hist_list_1.push_back(h_pmiss_back);
  TH1D * h_tof_back = new TH1D("tof_back","Time-of-Flight (Background);TOF;Counts",60,-10,50);
  hist_list_1.push_back(h_tof_back);
  TH1D * h_mom_back = new TH1D("mom_back","Measured Neutron Momentum (Background);Momentum (GeV/c);Counts",100,0,1.5);
  hist_list_1.push_back(h_mom_back);
  TH2D * h_ptheta_back = new TH2D("ptheta_back","Momentum vs #theta (Background);Neutron Momentum (GeV/c);#theta_{n} (deg)",100,0.2,1.2,10,40,140);
  hist_list_1.push_back(h_ptheta_back);
  TH1D * h_Edep_back = new TH1D("edep_back","Energy Deposition (Background);Energy Deposition (MeVee);Counts",100,0,100);
  hist_list_1.push_back(h_Edep_back);
  TH2D * h_Edep_p_back = new TH2D("edep_p_back","Energy Deposition vs Momentum;Momentum (GeV/c);Energy Deposition (MeVee)",100,0,1.5,100,0,100);
  hist_list_2.push_back(h_Edep_p_back);
  TH2D * h_Edep_theta_back = new TH2D("edep_theta_back","Energy Deposition vs #theta_{n};#theta_{n};Energy Deposition (MeVee)",180,0,180,100,0,100);
  hist_list_2.push_back(h_Edep_theta_back);
  TH2D * h_Edep_phi_back = new TH2D("edep_phi_back","Energy Deposition vs #phi;#phi;Energy Deposition (MeVee)",48,-180,180,100,0,100);
  hist_list_2.push_back(h_Edep_phi_back);
  TH2D * h_Edep_thetapn_back = new TH2D("edep_thetapn_back","Energy Deposition vs #theta_{pn};theta_{pn};Energy Deposition (MeVee)",180,0,180,100,0,100);
  hist_list_2.push_back(h_Edep_thetapn_back);


  /////////////////////////////////////
  //ML Features
  /////////////////////////////////////
  TH1D * h_energy_s = new TH1D("f_energy_s","Neutron Energy",100,0,100);
    hist_list_1.push_back(h_energy_s);
  TH1D * h_layermult_s = new TH1D("f_layermult_s","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_s);
  TH1D * h_size_s = new TH1D("f_size_s","Cluster Size",5,0,5);
    hist_list_1.push_back(h_size_s);
  TH1D * h_cnd_hits_s = new TH1D("f_cnd_hits_s","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_s);
  TH1D * h_cnd_energy_s = new TH1D("f_cnd_energy_s","Nearby CND Energy",100,0,100);
    hist_list_1.push_back(h_cnd_energy_s);
  TH1D * h_ctof_energy_s = new TH1D("f_ctof_energy_s","Nearby CTOF Energy",100,0,100);
    hist_list_1.push_back(h_ctof_energy_s);
  TH1D * h_ctof_hits_s = new TH1D("f_ctof_hits_s","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_s);
  TH1D * h_anglediff_s = new TH1D("f_anglediff_s","CVT Angle Diff",200,0,200);
    hist_list_1.push_back(h_anglediff_s);

  TH1D * h_energy_b = new TH1D("f_energy_b","Neutron Energy",100,0,100);
    hist_list_1.push_back(h_energy_b);
  TH1D * h_layermult_b = new TH1D("f_layermult_b","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_b);
  TH1D * h_size_b = new TH1D("f_size_b","Cluster Size",5,0,5);
    hist_list_1.push_back(h_size_b);
  TH1D * h_cnd_hits_b = new TH1D("f_cnd_hits_b","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_b);
  TH1D * h_cnd_energy_b = new TH1D("f_cnd_energy_b","Nearby CND Energy",100,0,100);
    hist_list_1.push_back(h_cnd_energy_b);
  TH1D * h_ctof_energy_b = new TH1D("f_ctof_energy_b","Nearby CTOF Energy",100,0,100);
    hist_list_1.push_back(h_ctof_energy_b);
  TH1D * h_ctof_hits_b = new TH1D("f_ctof_hits_b","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_b);
  TH1D * h_anglediff_b = new TH1D("f_anglediff_b","CVT Angle Diff",200,0,200);
    hist_list_1.push_back(h_anglediff_b);




  /////////////////////////////////////
  //Signal and Background Histograms
  /////////////////////////////////////
  



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




  // USE THIS SECTION TO ADD THE TMVA READER AND FEATURES
  // TMVA stuff
  //TMVA::Tools::Instance();
  TMVA::Reader * reader = new TMVA::Reader("!Color:!Silent");
  // define features for ML model
  Float_t energy, cnd_energy, ctof_energy, angle_diff, momentum;
  Float_t layermult, size, cnd_hits, ctof_hits;
  reader->AddVariable("energy", &energy);
  reader->AddVariable("layermult", &layermult);
  reader->AddVariable("size", &size);
  reader->AddVariable("cnd_hits", &cnd_hits);
  reader->AddVariable("cnd_energy", &cnd_energy);
  reader->AddVariable("ctof_energy", &ctof_energy);
  reader->AddVariable("ctof_hits", &ctof_hits);
  reader->AddVariable("angle_diff", &angle_diff);
  // spectator variable(s)
  reader->AddSpectator("momentum", &momentum);
  reader->BookMVA("MLP", "/w/hallb-scshelf2102/clas12/users/awild/RGM/rgm/NeutronVeto/dataset_6gev_pCD/weights/TrainNeutronVeto_TMVA_MLP.weights.xml");




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

      clasAna.Run(c12);
      // get particles by type
      auto allParticles = c12->getDetParticles();
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto neutrons=c12->getByID(2112);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);
      if(electrons.size()!=1){ continue;}

      double starttime = c12->event()->getStartTime();


  /////////////////////////////////////
  //Electron Kinematics  
  /////////////////////////////////////
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double EoP_e =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / p_e.Mag();
      int nphe = electrons[0]->che(HTCC)->getNphe();
      double vtz_e = electrons[0]->par()->getVz(); 
      h_phi_theta->Fill(p_e.Phi()*180/M_PI,p_e.Theta()*180/M_PI,weight);
      if(EoP_e<=0){ continue; }
      h_nphe->Fill(nphe,weight);
      h_vtz_e->Fill(vtz_e,weight);


      TVector3	p_q = p_b - p_e;
      double theta_q =  p_q.Theta() * 180 / M_PI;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;
      //if(WSq>1.25){continue;}

      if(xB<1.2){continue;}
      h_xB->Fill(xB,weight);
      h_QSq->Fill(QSq,weight);
      h_WSq->Fill(WSq,weight);
      h_xB_QSq->Fill(xB,QSq,weight);
      h_xB_WSq->Fill(xB,WSq,weight);
      h_QSq_WSq->Fill(QSq,WSq,weight);
      h_eangles->Fill(p_e.Phi()*180./M_PI,p_e.Theta()*180./M_PI,weight);

  /////////////////////////////////////
  //Proton Kinematics  
  /////////////////////////////////////
      SetLorentzVector(el,electrons[0]);
      TLorentzVector q = beam - el;
      for(auto p = protons.begin(); p != protons.end();++p){
	SetLorentzVector(proton_ptr,(*p));
	TLorentzVector miss = q + deut_ptr - proton_ptr;
	h_mmiss->Fill(miss.M(),weight);
	h_xB_mmiss->Fill(xB,miss.M(),weight);
      }

      clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
      auto lead    = clasAna.getLeadSRC();
      if(lead.size() != 1){continue;}
      SetLorentzVector(lead_ptr,lead[0]);
      TVector3 p_p = lead_ptr.Vect();
      /*
      int p_index = -1;
      for (int i=0; i<protons.size(); i++)
      {
        p_p.SetMagThetaPhi(protons[i]->getP(),protons[i]->getTheta(),protons[i]->getPhi());
        double dbeta = protons[i]->par()->getBeta() - p_p.Mag()/sqrt(p_p.Mag2()+mP*mP);
        double vzp = protons[i]->par()->getVz();
        double chipid = protons[i]->par()->getChi2Pid();
        if ((vzp-vtz_e)<-4. || (vzp-vtz_e)>4.) {continue;}
        if (chipid<-3. || chipid>3.) {continue;}
        if (dbeta<-0.05 || dbeta>0.05) {continue;}
        // require proton to be in correct angle and momentum range for the requested detector
        double p_theta = p_p.Theta()*180./M_PI;
        // forward proton
        //if (p_theta>40 || (p_p.Mag()<0.5 || p_p.Mag()>3.0)) {continue;}
        // central proton
        if ((p_theta<40 || p_theta>140) || (p_p.Mag()<0.2 || p_p.Mag()>1.2)) {continue;}
        p_index = i;
      }
      
      if (p_index<0) {continue;}
      p_p.SetMagThetaPhi(protons[p_index]->getP(),protons[p_index]->getTheta(),protons[p_index]->getPhi());
      */


  /////////////////////////////////////
  //Missing Momentum 
  /////////////////////////////////////
      TVector3 p_pred = p_q - p_p;
      /*if (p_pred.Mag()<0.2 || p_pred.Mag()>1.2) {continue;}
      if (p_pred.Theta()*180./M_PI<40 || p_pred.Theta()*180./M_PI>135) {continue;}
      double Ep = sqrt(mN*mN + p_p.Mag2());
      double Emiss = Ebeam + mD - p_e.Mag() - Ep;
      double mmiss = sqrt((Emiss*Emiss) - p_pred.Mag2());
      if (mmiss<0.7 || mmiss>1.05) {continue;}*/


  /////////////////////////////////////
  //CND Neutron Information
  /////////////////////////////////////
      int nn_CND = 0; int nn_CND_good = 0;
      if (neutrons.size()<1) {continue;}

      for(int j = 0; j < neutrons.size(); j++){
	TVector3 p_n;
	p_n.SetMagThetaPhi(neutrons[j]->getP(),neutrons[j]->getTheta(),neutrons[j]->getPhi());
	double E_n = sqrt(mN*mN + p_n.Mag2());
	double theta_n = p_n.Theta() * 180 / M_PI;
	double phi_n = p_n.Phi() * 180 / M_PI;
	double theta_nq = p_n.Angle(p_q) * 180 / M_PI;
	double beta_n = neutrons[j]->par()->getBeta();
	//double phi_diff = get_phi_diff(p_e,p_n);
	double cos0 = p_pred.Dot(p_n) / (p_pred.Mag()*p_n.Mag());

	double path_n = neutrons[j]->getPath();
	double beta_frommom_n = p_n.Mag()/E_n;
	double time_frommom_n = path_n / (c*beta_frommom_n);
	double time_frombeta_n = path_n / (c*beta_n);

        bool CND1 = (neutrons[j]->sci(clas12::CND1)->getDetector() == 3);
        bool CND2 = (neutrons[j]->sci(clas12::CND2)->getDetector() == 3);
        bool CND3 = (neutrons[j]->sci(clas12::CND3)->getDetector() == 3);
        bool CND = (CND1 || CND2 || CND3);
        bool CTOF = (neutrons[j]->sci(clas12::CTOF)->getDetector() == 4);

        double Edep = 0; double tof = 0;

        if (CND1){
          Edep = neutrons[j]->sci(clas12::CND1)->getEnergy();
          tof = neutrons[j]->sci(clas12::CND1)->getTime() - starttime;
        }
        else if (CND2){
          Edep = neutrons[j]->sci(clas12::CND2)->getEnergy();
          tof = neutrons[j]->sci(clas12::CND2)->getTime() - starttime;
        }
        else if (CND3){
          Edep = neutrons[j]->sci(clas12::CND3)->getEnergy();
          tof = neutrons[j]->sci(clas12::CND3)->getTime() - starttime;
        }
        else if (CTOF){
          Edep = neutrons[j]->sci(clas12::CTOF)->getEnergy();
          tof = neutrons[j]->sci(clas12::CTOF)->getTime() - starttime;
        }

        if (Edep<3) {continue;}

        //if (!!CTOF) {continue;}
        //if (!CND && !CTOF) {continue;}
        if (!(neutrons[j]->getRegion()==CD)) {continue;}

        nn_CND += 1;

        // cut out obviously bad neutrons
        if (beta_frommom_n==0) {continue;}
        if (theta_n==0) {continue;}
        if (phi_n==0) {continue;}
        if (p_n.Mag()<0.25) {continue;}

        // good neutron candidates start here
        if (theta_n<40 || theta_n>140) {continue;}
        //if (p_n.Mag()<0.25 || p_n.Mag()>1.2) {continue;}

        nn_CND_good += 1;

	h_theta_n_CND->Fill(theta_n,weight);
        h_phi_theta_n_CND->Fill(phi_n,theta_n,weight);
	h_mom_beta_n_CND->Fill(p_n.Mag(),beta_n,weight);
	h_mom_theta_n_CND->Fill(p_n.Mag(),theta_n,weight);
        h_phi_mom_n_CND->Fill(phi_n,p_n.Mag(),weight);



        // GET TMVA ML MODEL FEATURES FOR THIS NEUTRON HERE
        Struct ninfo = getFeatures(neutrons, allParticles, j);
        cnd_hits = ninfo.cnd_hits;
        cnd_energy = ninfo.cnd_energy;
        ctof_hits = ninfo.ctof_hits;
        ctof_energy = ninfo.ctof_energy;
        layermult = ninfo.layermult;
        energy = ninfo.energy;
        size = ninfo.size;
        angle_diff = ninfo.angle_diff;
        // spectator variable
        momentum = p_n.Mag();
        
        // GET TMVA VALUE - PLACE CUT TO ASSIGN NEUTRON AS SIGNAL/BACKGROUND
        double mvaValue = reader->EvaluateMVA("MLP");

        // BDT: signal cut>0.1, background cut<0.1
        // MLP: signal cut>0.5, background<0.5
	h_mva_n_CND->Fill(mvaValue,weight);
        //if ((mvaValue>0.9) && (ninfo.ctof_hits>0)) // signal
	if ((mvaValue>0.5)) // signal
        {
	  h_res_mom_theta_sig->Fill((p_pred.Mag()-p_n.Mag())/p_pred.Mag(),p_pred.Angle(p_n)*180/M_PI,weight);
          h_pn_angles_sig->Fill(p_p.Phi()*180./M_PI-p_n.Phi()*180./M_PI, p_p.Theta()*180./M_PI-p_n.Theta()*180./M_PI, weight);
          h_pn_pmiss_sig->Fill(p_pred.Mag(),p_n.Mag(),weight);
          h_cos0_sig->Fill(cos0,weight);
          h_pmiss_sig->Fill(p_pred.Mag(),weight);
	  h_n_thetaphi_sig->Fill(p_n.Phi()*180./M_PI, p_n.Theta()*180./M_PI, weight);
          h_tof_sig->Fill(tof,weight);
          h_mom_sig->Fill(p_n.Mag(),weight);
          h_ptheta_sig->Fill(p_n.Mag(),p_n.Theta()*180./M_PI,weight);
          h_Edep_sig->Fill(Edep,weight);
          h_Edep_p_sig->Fill(p_n.Mag(),Edep,weight);
          h_Edep_theta_sig->Fill(p_n.Theta()*180./M_PI,Edep,weight);
          h_Edep_phi_sig->Fill(p_n.Phi()*180./M_PI,Edep,weight);
          h_Edep_thetapn_sig->Fill(p_n.Angle(p_p)*180./M_PI,Edep,weight);
          // ML features
          h_energy_s->Fill(energy,weight);
          h_layermult_s->Fill(layermult,weight);
          h_size_s->Fill(size,weight);
          h_cnd_hits_s->Fill(cnd_hits,weight);
          h_cnd_energy_s->Fill(cnd_energy,weight);
          h_ctof_energy_s->Fill(ctof_energy,weight);
          h_ctof_hits_s->Fill(ctof_hits,weight);
          h_anglediff_s->Fill(angle_diff,weight);
        }
        else // background
        {
	  h_res_mom_theta_back->Fill((p_pred.Mag()-p_n.Mag())/p_pred.Mag(),p_pred.Angle(p_n)*180/M_PI,weight);
          h_pn_angles_back->Fill(p_p.Phi()*180./M_PI-p_n.Phi()*180./M_PI, p_p.Theta()*180./M_PI-p_n.Theta()*180./M_PI, weight);
          h_pn_pmiss_back->Fill(p_pred.Mag(),p_n.Mag(),weight);
          h_cos0_back->Fill(cos0,weight);
          h_pmiss_back->Fill(p_pred.Mag(),weight);
          h_n_thetaphi_back->Fill(p_n.Phi()*180./M_PI, p_n.Theta()*180./M_PI, weight);
          h_tof_back->Fill(tof,weight);
          h_mom_back->Fill(p_n.Mag(),weight);
          h_ptheta_back->Fill(p_n.Mag(),p_n.Theta()*180./M_PI,weight);
          h_Edep_back->Fill(Edep,weight);
          h_Edep_p_back->Fill(p_n.Mag(),Edep,weight);
          h_Edep_theta_back->Fill(p_n.Theta()*180./M_PI,Edep,weight);
          h_Edep_phi_back->Fill(p_n.Phi()*180./M_PI,Edep,weight);
          h_Edep_thetapn_back->Fill(p_n.Angle(p_p)*180./M_PI,Edep,weight);
          // ML features
          h_energy_b->Fill(energy,weight);
          h_layermult_b->Fill(layermult,weight);
          h_size_b->Fill(size,weight);
          h_cnd_hits_b->Fill(cnd_hits,weight);
          h_cnd_energy_b->Fill(cnd_energy,weight);
          h_ctof_energy_b->Fill(ctof_energy,weight);
          h_ctof_hits_b->Fill(ctof_hits,weight);
          h_anglediff_b->Fill(angle_diff,weight);
        }


        /*if (mvaValue<0.5 && Edep>18 && Edep<23) // background
        {
          std::cout << Edep << '\t' << p_n.Theta()*180./M_PI << '\t' << p_n.Mag() << '\t' << p_n.Phi()*180./M_PI << '\n';
        }*/
        /*if (mvaValue<0.5)
        {
          std::cout << p_e.Theta()*180./M_PI << '\t' << p_n.Angle(p_e)*180/M_PI << '\n';
        }*/


      }

      h_nsize_CND->Fill(nn_CND,nn_CND_good,weight);




  
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
  //Electron Kinematics  
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts:");
  double line = 0.8;

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

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_eangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mmiss->Draw();
  myCanvas->cd(2);
  h_xB_mmiss->Draw("colz");
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
  text.DrawLatex(0.2,0.5,"Neutron #theta: min=40 deg, max=140 deg");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_nsize_CND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_n_CND->Draw();
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

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_mva_n_CND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  /////////////////////////////////////
  //CND Neutron Signal & Background
  /////////////////////////////////////
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_n_thetaphi_sig->Draw("colz");
  myCanvas->cd(2);
  h_n_thetaphi_back->Draw("colz");
  myCanvas->cd(3);
  h_pn_angles_sig->Draw("colz");
  myCanvas->cd(4);
  h_pn_angles_back->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_res_mom_theta_sig->Draw("colz");
  myCanvas->cd(2);
  h_res_mom_theta_back->Draw("colz");
  myCanvas->cd(3);
  h_pmiss_sig->Draw("colz");
  myCanvas->cd(4);
  h_pmiss_back->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pn_pmiss_sig->Draw("colz");
  myCanvas->cd(2);
  h_pn_pmiss_back->Draw("colz");
  myCanvas->cd(3);
  h_cos0_sig->Draw();
  h_cos0_sig->SetStats(0);
  myCanvas->cd(4);
  h_cos0_back->Draw();
  h_cos0_back->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // ML features start here
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_energy_s->Draw();
  myCanvas->cd(2);
  h_energy_b->Draw();
  myCanvas->cd(3);
  h_layermult_s->Draw();
  myCanvas->cd(4);
  h_layermult_b->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_size_s->Draw();
  myCanvas->cd(2);
  h_size_b->Draw();
  myCanvas->cd(3);
  h_anglediff_s->Draw();
  myCanvas->cd(4);
  h_anglediff_b->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_cnd_hits_s->Draw();
  myCanvas->cd(2);
  h_cnd_hits_b->Draw();
  myCanvas->cd(3);
  h_cnd_energy_s->Draw();
  myCanvas->cd(4);
  h_cnd_energy_b->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_ctof_hits_s->Draw();
  myCanvas->cd(2);
  h_ctof_hits_b->Draw();
  myCanvas->cd(3);
  h_ctof_energy_s->Draw();
  myCanvas->cd(4);
  h_ctof_energy_b->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  // ML features end here

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_tof_sig->Draw();
  myCanvas->cd(2);
  h_tof_back->Draw();
  myCanvas->cd(3);
  h_mom_sig->Draw();
  myCanvas->cd(4);
  h_mom_back->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_ptheta_sig->Draw("colz");
  myCanvas->cd(2);
  h_ptheta_back->Draw("colz");
  myCanvas->cd(3);
  h_Edep_sig->Draw();
  myCanvas->cd(4);
  h_Edep_back->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_Edep_p_sig->Draw("colz");
  myCanvas->cd(2);
  h_Edep_p_back->Draw("colz");
  myCanvas->cd(3);
  h_Edep_theta_sig->Draw("colz");
  myCanvas->cd(4);
  h_Edep_theta_back->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_Edep_phi_sig->Draw("colz");
  myCanvas->cd(2);
  h_Edep_phi_back->Draw("colz");
  myCanvas->cd(3);
  h_Edep_thetapn_sig->Draw("colz");
  myCanvas->cd(4);
  h_Edep_thetapn_back->Draw("colz");
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

