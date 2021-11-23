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
bool pointsToBand(double theta,double phi,double z_m);

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
  //Lead Proton Checks
  /////////////////////////////////////
  TH1D * h_xB_Lead = new TH1D("xB_Lead","x_{B} Lead;x_{B};Counts",100,0.0,2.0);
  hist_list_1.push_back(h_xB_Lead);
  TH1D * h_Q2_Lead = new TH1D("Q2_Lead","Q^{2} Lead;Q^{2};Counts",100,0.0,5.0);
  hist_list_1.push_back(h_Q2_Lead);
  TH2D * h_xB_Q2_Lead = new TH2D("xB_Q2_Lead","x_{B} vs. Q^{2} Lead;x_{B};Q^{2};Counts",100,0.0,2.0,100,0.0,5.0);
  hist_list_2.push_back(h_xB_Q2_Lead);

  TH1D * h_theta_p_Lead = new TH1D("theta_p_Lead","#theta_{p,Lead};#theta_{p,Lead};Counts",180,0,180);
  hist_list_1.push_back(h_theta_p_Lead);
  TH1D * h_mom_p_Lead = new TH1D("mom_p_Lead","p_{p,Lead};p_{p,Lead};Counts",100,0,4);
  hist_list_1.push_back(h_mom_p_Lead);
  TH1D * h_phi_p_Lead = new TH1D("phi_p_Lead","#phi_{p,Lead};#phi_{p,Lead};Counts",180,-180,180);
  hist_list_1.push_back(h_phi_p_Lead);
  TH1D * h_theta_pq_Lead = new TH1D("theta_pq_Lead","#theta_{pq} Lead;#theta_{pq};Counts",180,0,90);
  hist_list_1.push_back(h_theta_pq_Lead);
  TH2D * h_mom_theta_p_Lead = new TH2D("mom_theta_p_Lead","#p_{p,Lead} vs. #theta_{p,Lead} ;#p_{p,Lead};#theta_{p,Lead}",100,0,4,100,0,135);
  hist_list_2.push_back(h_mom_theta_p_Lead);
  TH1D * h_phi_e_p_Lead = new TH1D("phi_e_p_Lead","|#phi_{e} - #phi_{p,Lead}|;|#phi_{e} - #phi_{p,Lead}|,Counts",100,120,180);
  hist_list_1.push_back(h_phi_e_p_Lead);
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
  //Points to BAND information
  /////////////////////////////////////
  TH1D * h_pmiss_BAND = new TH1D("pmiss_BAND","p_{miss} BAND;p_{miss};Counts",100,0,1.5);
  hist_list_1.push_back(h_pmiss_BAND);
  TH1D * h_thetamiss_BAND = new TH1D("thetamiss_BAND","#theta_{miss} BAND;#theta_{miss};Counts",100,145,180);
  hist_list_1.push_back(h_thetamiss_BAND);
  TH1D * h_mmiss_BAND = new TH1D("mmiss_BAND","m_{miss} BAND;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_BAND);  

  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  TH1D * h_xB_SRC = new TH1D("xB_SRC","x_{B} SRC;x_{B};Counts",100,1.0,2.0);
  hist_list_1.push_back(h_xB_SRC);
  TH1D * h_Q2_SRC = new TH1D("Q2_SRC","Q^{2} SRC;Q^{2};Counts",100,0.0,7.0);
  hist_list_1.push_back(h_Q2_SRC);
  TH1D * h_W_SRC = new TH1D("W_SRC","W SRC;W;Counts",100,0.0,1.0);
  hist_list_1.push_back(h_W_SRC);
  TH1D * h_pmiss_SRC = new TH1D("pmiss_SRC","p_{miss} SRC;p_{miss};Counts",50,0,1.5);
  hist_list_1.push_back(h_pmiss_SRC);
  TH1D * h_mmiss_SRC = new TH1D("mmiss_SRC","m_{miss} SRC;m_{miss};Counts",100,0.4,1.4);
  hist_list_1.push_back(h_mmiss_SRC);
  TH2D * h_pmiss_theta_miss_SRC = new TH2D("pmiss_theta_miss_SRC","p_{miss} vs. #theta_{miss};p_{miss};#theta_{miss};Counts",100,0,1.5,180,0,180);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_pmiss_xB_SRC = new TH2D("pmiss_xB_SRC","p_{miss} vs. x_{B};p_{miss};x_{B};Counts",100,0,1.5,100,1.0,2.0);
  hist_list_2.push_back(h_pmiss_xB_SRC);
  TH2D * h_pmiss_theta_L_SRC = new TH2D("pmiss_theta_L_SRC","p_{miss} vs. #theta_{L};p_{miss};#theta_{L};Counts",100,0,1.5,100,5,90);
  hist_list_2.push_back(h_pmiss_theta_miss_SRC);
  TH2D * h_xB_Loq_SRC = new TH2D("xB_Loq","x_{B} vs |p|/|q|;x_{B};|p|/|q|",100,1,2,100,0.5,1.2);
  hist_list_2.push_back(h_xB_Loq_SRC);
  TH2D * h_mmiss_xB_SRC = new TH2D("mmiss_xB_SRC","m_{miss} vs. x_{B};m_{miss};x_{B};Counts",100,0.4,1.4,100,1.0,2.0);
  hist_list_2.push_back(h_mmiss_xB_SRC);
  TH2D * h_mmiss_mom_p_SRC = new TH2D("mmiss_mom_p_SRC","m_{miss} vs. p_{p,Lead};m_{miss};p_{p,Lead};Counts",100,0.4,1.4,100,0,3);
  hist_list_2.push_back(h_mmiss_mom_p_SRC);
  TH2D * h_mmiss_theta_p_SRC = new TH2D("mmiss_theta_p_SRC","m_{miss} vs. #theta_{p,Lead};m_{miss};#theta_{p,Lead};Counts",100,0.4,1.4,100,0,90);
  hist_list_2.push_back(h_mmiss_theta_p_SRC);

  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
  TH1D * h_count_AllRec = new TH1D("count_AllRec","Number of Recoils;Multiplicity",5,-0.5,4.5);
  hist_list_1.push_back(h_count_AllRec);
  TH1D * h_p_rec_AllRec = new TH1D("p_rec_AllRec","p All Recoils;p_{rec}",50,0,3);
  hist_list_1.push_back(h_p_rec_AllRec);
  TH1D * h_theta_rec_AllRec = new TH1D("theta_rec_AllRec","Theta All Recoils;#theta_{rec}",100,0,135);
  hist_list_1.push_back(h_theta_rec_AllRec);
  TH1D * h_phi_rec_AllRec = new TH1D("phi_rec_AllRec","Phi All Recoils;#phi_{rec}",90,-180,180);
  hist_list_1.push_back(h_phi_rec_AllRec);
  TH2D * h_phi_theta_rec_AllRec = new TH2D("phi_theta_rec_AllRec","#phi_{rec} vs. #theta_{rec} ;#phi_{rec};#theta_{rec}",90,-180,180,50,0,135);
  hist_list_2.push_back(h_phi_theta_rec_AllRec);
  TH2D * h_mom_theta_rec_AllRec = new TH2D("mom_theta_rec_AllRec","p_{rec} vs. #theta_{rec} ;p_{rec};#theta_{rec}",50,0,3,50,0,135);
  hist_list_2.push_back(h_mom_theta_rec_AllRec);

  TH1D * h_vtz_rec_AllRec = new TH1D("vtz_rec_AllRec","Proton Z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_rec_AllRec);
  TH1D * h_vtz_erec_delta_AllRec = new TH1D("vtz_erec_delta_AllRec","#Delta Vertex (e^{-} and p_{rec});#Delta Vertex;Counts",50,-5,5);
  hist_list_1.push_back(h_vtz_erec_delta_AllRec);
  TH2D * h_vtz_e_vtz_rec_AllRec = new TH2D("vtz_e_vtz_rec_AllRec","Electron Z Vertex vs. Proton Recoil Z Vertex;vertex e;vertex p",100,-10,10,100,-10,10);
  hist_list_2.push_back(h_vtz_e_vtz_rec_AllRec);
  TH1D * h_chiSq_rec_AllRec = new TH1D("chiSq_rec_AllRec","#chi^{2}_{rec} All Recoils;#chi^{2}_{rec}",50,-10,10);
  hist_list_1.push_back(h_chiSq_rec_AllRec);
  TH1D * h_timediff_rec_AllRec = new TH1D("timediff_rec_AllRec","ToF-ToF_{|p|} Recoils;ToF-ToF_{|p|};Counts",50,-2,2);
  hist_list_1.push_back(h_timediff_rec_AllRec);
  TH2D * h_mom_beta_rec_AllRec = new TH2D("mom_beta_rec_AllRec","p_{rec} vs. #beta_{rec} ;p_{rec};#beta_{rec}",50,0,4,100,0.0,1);
  hist_list_2.push_back(h_mom_beta_rec_AllRec);

  TH2D * h_chiSq_rec_p_rec_AllRec = new TH2D("chiSq_rec_p_rec_AllRec","#chi^{2}_{rec} vs. p_{rec} All Recoils;#chi^{2}_{rec};p_{rec}",100,-20,20,100,0,4);
  hist_list_2.push_back(h_chiSq_rec_p_rec_AllRec);
  TH2D * h_timediff_rec_p_rec_AllRec = new TH2D("timediff_rec_p_rec_AllRec","ToF-ToF_{|p|} vs. p_{rec} All Recoils;ToF-ToF_{|p|};p_{rec}",100,-20,20,100,0,4);
  hist_list_2.push_back(h_timediff_rec_p_rec_AllRec);
  TH2D * h_chiSq_rec_vtzdiff_rec_AllRec = new TH2D("chiSq_rec_vtzdiff_rec_AllRec","#chi^{2}_{rec} vs. #Delta Vertex All Recoils;#chi^{2}_{rec};#Delta Vertex",100,-20,20,100,-15,15);
  hist_list_2.push_back(h_chiSq_rec_vtzdiff_rec_AllRec);
  TH2D * h_chiSq_rec_timediff_rec_AllRec = new TH2D("chiSq_rec_timediff_rec_AllRec","#chi^{2}_{rec} vs. ToF-ToF_{|p|} All Recoils;#chi^{2}_{rec};ToF-ToF_{|p|}",100,-20,20,100,-10,10);
  hist_list_2.push_back(h_chiSq_rec_vtzdiff_rec_AllRec);
  TH2D * h_p_rec_vtzdiff_rec_AllRec = new TH2D("p_rec_vtzdiff_rec_AllRec","p_{rec} vs. #Delta Vertex All Recoils;p_{rec};#Delta Vertex",100,0,4,100,-15,15);
  hist_list_2.push_back(h_p_rec_vtzdiff_rec_AllRec);


  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////
  TH1D * h_p_2_Rec = new TH1D("p_2_Rec","p_{rec};p_{rec};Counts",50,0,1.5);
  hist_list_1.push_back(h_p_2_Rec);
  TH1D * h_pmiss_Rec = new TH1D("pmiss_Rec","p_{miss} Rec;p_{miss};Counts",50,0,1.5);
  hist_list_1.push_back(h_pmiss_Rec);
  TH1D * h_p_rel_Rec = new TH1D("p_rel_Rec","p_{rel};p_{rel};Counts",50,0,1.5);
  hist_list_1.push_back(h_p_rel_Rec);
  TH1D * h_p_cm_Rec = new TH1D("p_cm_Rec","p_{C.M.};p_{C.M.};Counts",50,0,2.5);
  hist_list_1.push_back(h_p_cm_Rec);
  TH1D * h_p_t_cm_Rec = new TH1D("p_t_cm_Rec","p_{t,C.M.};p_{t,C.M.};Counts",50,-1,1);
  hist_list_1.push_back(h_p_t_cm_Rec);
  TH1D * h_p_y_cm_Rec = new TH1D("p_y_cm_Rec","p_{y,C.M.};p_{y,C.M.};Counts",50,-1,1);
  hist_list_1.push_back(h_p_y_cm_Rec);
  TH1D * h_p_x_cm_Rec = new TH1D("p_x_cm_Rec","p_{x,C.M.};p_{x,C.M.};Counts",50,-1,1);
  hist_list_1.push_back(h_p_x_cm_Rec);
  TH1D * h_theta_rel_Rec = new TH1D("theta_rel_Rec","#theta_{rel};#theta_{rel};Counts",90,0,180);
  hist_list_1.push_back(h_theta_rel_Rec);
  TH2D * h_p_cm_theta_rel_Rec = new TH2D("p_cm_theta_rel_Rec","p_{C.M.} vs. #theta_{rel};p_{C.M.};#theta_{rel}",50,0,2.5,90,0,180);
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

      if(!myCut.electroncut(c12)){continue;}      

      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double vtz_e = electrons[0]->par()->getVz();
      TVector3	p_q = p_b - p_e;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;

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
      h_mom_p_Lead->Fill(p_L.Mag(),weight);
      h_phi_p_Lead->Fill(phi_L,weight);
      h_theta_pq_Lead->Fill(theta_Lq,weight);
      h_mom_theta_p_Lead->Fill(p_L.Mag(),theta_L,weight);
      h_phi_e_p_Lead->Fill(phi_diff,weight);
      h_xB_Lead->Fill(xB,weight);
      h_Q2_Lead->Fill(QSq,weight);
      h_xB_Q2_Lead->Fill(xB,QSq,weight);
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

      
      if((p_miss.Theta()>(M_PI/2)) && pointsToBand(p_miss.Theta(),p_miss.Phi(),vtz_p)){
	if(p_miss.Mag()>0.2){
	  h_pmiss_BAND->Fill(p_miss.Mag(),weight);
	  h_thetamiss_BAND->Fill(theta_miss,weight);
	  h_mmiss_BAND->Fill(mmiss,weight);
	}
      }
        
  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
      if(!myCut.leadSRCnucleoncut(c12,index_L)){continue;}

      h_xB_SRC->Fill(xB,weight);
      h_Q2_SRC->Fill(QSq,weight);
      h_W_SRC->Fill(sqrt(WSq),weight);
      h_pmiss_SRC->Fill(p_miss.Mag(),weight);
      h_mmiss_SRC->Fill(mmiss,weight);

      h_pmiss_theta_miss_SRC->Fill(p_miss.Mag(),theta_miss,weight);
      h_pmiss_xB_SRC->Fill(p_miss.Mag(),xB,weight);
      h_pmiss_theta_L_SRC->Fill(p_miss.Mag(),theta_L,weight);
      h_xB_Loq_SRC->Fill(xB,Loq,weight);
      h_mmiss_xB_SRC->Fill(mmiss,xB,weight);
      h_mmiss_mom_p_SRC->Fill(mmiss,p_L.Mag(),weight);
      h_mmiss_theta_p_SRC->Fill(mmiss,theta_L,weight);


  /////////////////////////////////////
  //Recoil Nucleons
  /////////////////////////////////////
      for(int j = 0; j < protons.size(); j++){
	if(j==index_L){continue;}
	double mom = protons[j]->getP();
	double theta = protons[j]->getTheta()*180/M_PI;
	double phi = protons[j]->getPhi()*180/M_PI;
	double vtz_p = protons[j]->par()->getVz();
	double beta_p = protons[j]->par()->getBeta();
	double path_p = protons[j]->getPath();
	double beta_frommom_p = mom/sqrt(mom*mom + mN*mN);
	double time_frommom_p = path_p / (c*beta_frommom_p);

	double chi = protons[j]->par()->getChi2Pid();

	double time_frombeta_p = path_p / (c*beta_p);
	double time_diff = time_frombeta_p-time_frommom_p;

	h_p_rec_AllRec->Fill(mom,weight);
	h_theta_rec_AllRec->Fill(theta,weight);
	h_phi_rec_AllRec->Fill(phi,weight);
	h_phi_theta_rec_AllRec->Fill(phi,theta,weight);
	h_mom_theta_rec_AllRec->Fill(mom,theta,weight);

	h_vtz_rec_AllRec->Fill(vtz_p,weight);
	h_vtz_erec_delta_AllRec->Fill(vtz_e-vtz_p,weight);
	h_vtz_e_vtz_rec_AllRec->Fill(vtz_e,vtz_p,weight);
	h_chiSq_rec_AllRec->Fill(chi,weight);
	h_timediff_rec_AllRec->Fill(time_diff,weight);
	h_mom_beta_rec_AllRec->Fill(protons[j]->getP(),protons[j]->par()->getBeta(),weight);


	h_chiSq_rec_p_rec_AllRec->Fill(chi,mom,weight);
	h_timediff_rec_p_rec_AllRec->Fill(time_diff,mom,weight);
	h_chiSq_rec_vtzdiff_rec_AllRec->Fill(chi,vtz_e-vtz_p,weight);
	h_chiSq_rec_timediff_rec_AllRec->Fill(chi,time_diff,weight);
	h_p_rec_vtzdiff_rec_AllRec->Fill(mom,vtz_e-vtz_p,weight);
	
      }
      h_count_AllRec->Fill(protons.size()-1,weight);

  /////////////////////////////////////
  //Recoil SRC Proton Checks
  /////////////////////////////////////
      int index_R = myCut.recoilSRCnucleoncut(c12,index_L);
      if(index_R < 0){ continue; }
      //h_vtz_e_vtz_rec_AllRec->Fill(vtz_e,protons[index_R]->par()->getVz(),weight);      
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
      h_pmiss_Rec->Fill(p_miss.Mag(),weight);
      h_p_rel_Rec->Fill(p_rel.Mag(),weight);
      h_p_cm_Rec->Fill(p_cm.Mag(),weight);
      h_p_t_cm_Rec->Fill(p_cm.Dot(vt),weight);
      h_p_y_cm_Rec->Fill(p_cm.Dot(vy),weight);
      h_p_x_cm_Rec->Fill(p_cm.Dot(vx),weight);
      h_theta_rel_Rec->Fill(theta_rel,weight);
      h_p_cm_theta_rel_Rec->Fill(p_cm.Mag(),theta_rel,weight);
      
  }
  cout<<counter<<endl;


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
  //Lead Proton Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p) Cuts");
  double line = 0.7;
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
  h_xB_Lead->Draw();
  myCanvas->cd(2);
  h_Q2_Lead->Draw();
  myCanvas->cd(3);
  h_xB_Q2_Lead->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_theta_p_Lead->Draw();
  myCanvas->cd(2);
  h_mom_p_Lead->Draw();
  myCanvas->cd(3);
  h_phi_p_Lead->Draw();
  myCanvas->cd(4);
  h_theta_pq_Lead->Draw();
  myCanvas->cd(5);
  h_mom_theta_p_Lead->Draw("colz");
  myCanvas->cd(6);
  h_phi_e_p_Lead->Draw();
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
  myCanvas->cd(5);
  h_vtz_e_vtz_p_Lead->Draw("colz");
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
  //Points to BAND information
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead}->BAND) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{Lead}) Cuts");
  text.DrawLatex(0.2,0.7,"Lead Points to BAND");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_BAND->Draw();
  myCanvas->cd(2);
  h_thetamiss_BAND->Draw();
  myCanvas->cd(3);
  h_mmiss_BAND->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Lead SRC Proton Checks
  /////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD}) Cuts");
  line = 0.7;
  if(myCut.getDoCut(lsrc_cuts)){
    myCut.print_cut_onPDF(text,lsrc_Q2,line);
    myCut.print_cut_onPDF(text,lsrc_xB,line);
    myCut.print_cut_onPDF(text,lsrc_pmiss,line);
    myCut.print_cut_onPDF(text,lsrc_mmiss,line);
    myCut.print_cut_onPDF(text,lsrc_loq,line);
  }
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_xB_SRC->Draw();
  myCanvas->cd(2);
  h_Q2_SRC->Draw();
  myCanvas->cd(3);
  h_mmiss_SRC->Draw();
  myCanvas->cd(4);
  h_pmiss_SRC->Draw();
  myCanvas->cd(5);
  h_W_SRC->Draw();
  myCanvas->cd(6);
  h_pmiss_theta_miss_SRC->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_pmiss_xB_SRC->Draw("colz");
  myCanvas->cd(2);
  h_pmiss_theta_L_SRC->Draw("colz");
  myCanvas->cd(3);
  h_xB_Loq_SRC->Draw("colz");
  myCanvas->cd(4);
  h_mmiss_xB_SRC->Draw("colz");
  myCanvas->cd(5);
  h_mmiss_mom_p_SRC->Draw("colz");
  myCanvas->cd(6);
  h_mmiss_theta_p_SRC->Draw("colz");
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
  h_count_AllRec->Draw();
  myCanvas->cd(2);
  h_p_rec_AllRec->Draw();
  myCanvas->cd(3);
  h_theta_rec_AllRec->Draw();
  myCanvas->cd(4);
  h_phi_rec_AllRec->Draw();
  myCanvas->cd(5);
  h_phi_theta_rec_AllRec->Draw("colz");
  myCanvas->cd(6);
  h_mom_theta_rec_AllRec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_vtz_rec_AllRec->Draw();
  myCanvas->cd(2);
  h_vtz_erec_delta_AllRec->Draw();
  myCanvas->cd(3);
  h_vtz_e_vtz_rec_AllRec->Draw("colz");
  myCanvas->cd(4);
  h_chiSq_rec_AllRec->Draw("colz");
  myCanvas->cd(5);
  h_timediff_rec_AllRec->Draw();
  myCanvas->cd(6);
  h_mom_beta_rec_AllRec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  myCanvas->cd(1);
  h_chiSq_rec_p_rec_AllRec->Draw("colz");
  myCanvas->cd(2);
  h_timediff_rec_p_rec_AllRec->Draw("colz");
  myCanvas->cd(3);
  h_chiSq_rec_vtzdiff_rec_AllRec->Draw("colz");
  myCanvas->cd(4);
  h_chiSq_rec_timediff_rec_AllRec->Draw("colz");
  myCanvas->cd(5);
  h_p_rec_vtzdiff_rec_AllRec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  /////////////////////////////////////
  //Recoil SRC Nucleons
  /////////////////////////////////////

  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e'p_{Lead,SRC}p_{Rec,SRC}) Cuts:");
  text.DrawLatex(0.2,0.8,"(e,e'p_{LEAD,SRC},p_{Rec}) Cuts");
  line = 0.7;
  if(myCut.getDoCut(rsrc_cuts)){
    myCut.print_cut_onPDF(text,rsrc_pid,line);
    myCut.print_cut_onPDF(text,rsrc_mom,line);
    myCut.print_cut_onPDF(text,rsrc_chipid,line);
    myCut.print_cut_onPDF(text,rsrc_vtzdiff,line);
  }
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
  myCanvas->cd(3);
  h_pmiss_Rec->Draw();
  myCanvas->cd(4);
  TH1D * num = (TH1D*)h_pmiss_Rec->Clone("epp/ep");
  TH1D * den = (TH1D*)h_pmiss_SRC->Clone("clone");
  num->Divide(den);
  hist_list_1.push_back(num);
  num->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }
  outFile->Close();
}


void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

bool pointsToBand(double theta,double phi,double z_m){

  double inset = 0;//inset distance from edges of BAND in [cm]
  double z = z_m;

  // Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java  
  double thickness  = 7.3;                                // thickness of each bar (cm)
  double layerGap[] = {7.94, 7.62, 7.94, 7.62, 7.3};      // gap between center of neighbouring layers (cm), 1-2, 2-3, 3-4, 4-5, 5-6                                                 

  // Numbers taken from clas-band-calib/bin/src/org/clas/fcmon/band/BANDConstants.java  
  double bandlen[]  = {163.7,201.9,51.2,51.2,201.9};

  // Distance from ideal target to upstream end of BAND                                 
  // (from BAND survey report, 02/18/2019)                                              
  double zUpst = (-302.69-302.69-302.57-302.64)/4.; // [cm]                             

  // Distance from ideal target to downstream end of layer 5                            
  double zDown = (zUpst + 5*thickness) - z_m;

  double rho   = zDown/cos(theta);
  double xDown = rho*sin(theta)*cos(phi);
  double yDown = rho*sin(theta)*sin(phi);

  double globalX = (-240.5-240.5+241.0+243.7)/4./10.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)                                         
  double globalY = (-211.0+228.1-210.6+228.1)/4./10.; // [cm]

  // Sector boundaries                                                                  
  double topSec1  = globalY + 13*thickness - inset;
  double topSec2  = globalY + 10*thickness;
  double topSec34 = globalY +  3*thickness;
  double topSec5  = globalY -  3*thickness;
  double downSec5 = globalY -  5*thickness + inset;

  if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

  if(
     (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. - inset) ||
     ( (yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. - inset) &&   !(yDown >= topSec2 - inset && yDown < topSec2  && fabs(xDown) > (bandlen[0]/2.  - inset) ) && !(yDown >= topSec34 && yDown < topSec34 + inset && abs(xDown) < bandlen[1]/2 - bandlen[2] + inset) ) ||
     (yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2.- inset && fabs(xDown) > (bandlen[1]/2.-bandlen[2] + inset)) ||
     ( (yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. - inset) && !(yDown >= topSec5 -inset && yDown < topSec5 + inset && abs(xDown) < bandlen[1]/2 - bandlen[2] + inset) )

      ){
    return 1;
    }
  return 0;
}

