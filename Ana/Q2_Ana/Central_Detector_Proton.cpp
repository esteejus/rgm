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
#include <TGraph.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <TLine.h>
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLatex.h"
#include "HipoChain.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}


double Cut_Params[] = {0.0152222,
		       0.816844,
		       -0.0950375,
		       0.255628,
		       0.0760525,
		       0.240862,
		       -0.000276433,
		       0.229085};

double sq(double x){return x*x;}

double cut_func(double x, double a, double b, double c, double d){
  return a * (1 + (b/(x-d)) + (c/sq(x-d))); 
}

bool pass_cut(double mom, double DT, double w){
  double mu = cut_func(mom,Cut_Params[0],Cut_Params[1],Cut_Params[2],Cut_Params[3]);
  double sigma = cut_func(mom,Cut_Params[4],Cut_Params[5],Cut_Params[6],Cut_Params[7]);
  double upper = mu + w * sigma;
  double lower = mu - w * sigma;
  if((DT<upper)&&(DT>lower)){
    return true;
  }
  return false;
}


bool CD_fiducial(double phi, double theta, double momT){
  bool pass_fiducial = true;
  double fiducial_phi_width = 3;
  double fiducial_phi_shift = 0;
  double fiducial_momT_start = 0.15;
  double fiducial_phi_central = (-asin(fiducial_momT_start/momT) - (M_PI/2)) * 180/M_PI;
  if( (fabs(phi-fiducial_phi_central-fiducial_phi_shift)<fiducial_phi_width) ||
      (fabs(phi-fiducial_phi_central-fiducial_phi_shift-120)<fiducial_phi_width) || 
      (fabs(phi-fiducial_phi_central-fiducial_phi_shift-240)<fiducial_phi_width) || 
      (theta<40) ||
      (theta>125)){
    pass_fiducial = false;
  }
  return pass_fiducial;
}

void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
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

  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  //now get reference to (unique)ptr for accessing data in loop
  //this will point to the correct place when file changes
  //  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;

  double beam_E = 5.98;

  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;

  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  
  //Fid
  TH1D * h_theta_CD_bc = new TH1D("theta_CD_bc","#theta vs. Counts; #theta [degrees]; Counts",100,0,180);
  hist_list.push_back(h_theta_CD_bc);
  TH2D * h_phi_momT_CD_bc = new TH2D("phi_momT_CD_bc","#phi vs. p_{T}; #phi; Transverse Momentum [GeV]",100,-180,180,100,0,2);
  hist_list.push_back(h_phi_momT_CD_bc);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_bc = new TH2D("mom_ToFToF_d_ToFMom_CD_bc","Momentum vs. ToF - ToF_{expected} [ns];Momentum [GeV];ToF - ToF_{expected} [ns]",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_bc);
  TH1D * h_ToFToF_d_ToFMom_CD_bc_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_bc_bin_%d",i+1);
    sprintf(temp_title,"ToF - ToF_{expected}[ns] Bin=%d;ToF - ToF_{expected} [ns];Counts",i+1);
    h_ToFToF_d_ToFMom_CD_bc_bin[i] = new TH1D(temp_name,temp_title,100,-1,1);
    hist_list.push_back(h_ToFToF_d_ToFMom_CD_bc_bin[i]);
  }

  TH1D * h_edge_first_CD_bc = new TH1D("edge_first_CD_bc","First Edge; Distance to Edge [cm]; Counts",100,-5,25);
  hist_list.push_back(h_edge_first_CD_bc);
  TH1D * h_edge_last_CD_bc = new TH1D("edge_last_CD_bc","Last Edge; Distance to Edge [cm]; Counts",100,-5,25);
  hist_list.push_back(h_edge_last_CD_bc);

  TH1D * h_theta_CD_ac = new TH1D("theta_CD_ac","#theta vs. Counts; #theta [degrees]; Counts",100,0,180);
  hist_list.push_back(h_theta_CD_ac);
  TH2D * h_phi_momT_CD_ac = new TH2D("phi_momT_CD_ac","#phi vs. p_{T}; #phi; Transverse Momentum [GeV]",100,-180,180,100,0,2);
  hist_list.push_back(h_phi_momT_CD_ac);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_ac = new TH2D("mom_ToFToF_d_ToFMom_CD_ac","Momentum vs. ToF - ToF_{expected} [ns];Momentum [GeV];ToF - ToF_{expected} [ns]",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_ac);
  TH1D * h_ToFToF_d_ToFMom_CD_ac_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_ac_bin_%d",i+1);
    sprintf(temp_title,"ToF - ToF_{expected}[ns] Bin=%d;ToF - ToF_{expected} [ns];Counts",i+1);
    h_ToFToF_d_ToFMom_CD_ac_bin[i] = new TH1D(temp_name,temp_title,100,-1,1);
    hist_list.push_back(h_ToFToF_d_ToFMom_CD_ac_bin[i]);
  }

  TH2D * h_mom_ToFToF_d_ToFMom_CD_bad = new TH2D("mom_ToFToF_d_ToFMom_CD_bad","Momentum vs. ToF - ToF_{expected} [ns];Momentum [GeV];ToF - ToF_{expected} [ns]",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_bad);
  TH1D * h_ToFToF_d_ToFMom_CD_bad_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_bad_bin_%d",i+1);
    sprintf(temp_title,"ToF - ToF_{expected}[ns] Bin=%d;ToF - ToF_{expected} [ns];Counts",i+1);
    h_ToFToF_d_ToFMom_CD_bad_bin[i] = new TH1D(temp_name,temp_title,100,-1,1);
    hist_list.push_back(h_ToFToF_d_ToFMom_CD_bad_bin[i]);
  }

  //Vertex
  TH1D * h_vtz_CD_bc = new TH1D("vtz_CD_bc","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  hist_list.push_back(h_vtz_CD_bc);
  TH1D * h_diffvtz_CD_bc = new TH1D("diffvtz_CD_bc","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  hist_list.push_back(h_diffvtz_CD_bc);
  TH2D * h_vtz_e_p_CD_bc = new TH2D("vtz_e_p_CD_bc","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);
  hist_list.push_back(h_vtz_e_p_CD_bc);

  TH1D * h_vtz_CD_ac = new TH1D("vtz_CD_ac","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  hist_list.push_back(h_vtz_CD_ac);
  TH1D * h_diffvtz_CD_ac = new TH1D("diffvtz_CD_ac","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  hist_list.push_back(h_diffvtz_CD_ac);
  TH2D * h_vtz_e_p_CD_ac = new TH2D("vtz_e_p_CD_ac","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);
  hist_list.push_back(h_vtz_e_p_CD_ac);
  
  //PID
  TH1D * h_StartTime = new TH1D("StartTime","StartTime;StartTime;Counts",100,80,105);
  hist_list.push_back(h_StartTime);
  TH1D * h_HitTime = new TH1D("HitTime","HitTime;HitTime;Counts",100,80,105);
  hist_list.push_back(h_HitTime);
  TH1D * h_ToF = new TH1D("ToF","ToF;ToF;Counts",100,-1,6);
  hist_list.push_back(h_ToF);
  TH1D * h_Path = new TH1D("Path","Path;Path;Counts",100,0,50);
  hist_list.push_back(h_Path);
  TH2D * h_ToF_Path = new TH2D("ToF_Path","Path vs. ToF;ToF;Path",100,-1,6,100,0,50);
  hist_list.push_back(h_ToF_Path);
  TH2D * h_mom_DT = new TH2D("mom_DT","#Delta ToF vs. Momentum;p [GeV];#Delta ToF [ns]",100,0,3,200,-2,2);
  hist_list.push_back(h_mom_DT);
  TH2D * h_mom_beta = new TH2D("mom_beta","#beta vs. Momentum;Momentum [GeV];#beta",100,0,3,100,0,1.1);
  hist_list.push_back(h_mom_beta);
  TH1D * h_mom125_beta = new TH1D("mom125_beta","#beta;#beta;Counts",100,0.6,1.0);
  hist_list.push_back(h_mom125_beta);
  TH1D * h_DT_mom_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"DT_mom_bin_%d",i+1);
    h_DT_mom_bin[i] = new TH1D(temp_name,"#Delta Time;#Delta ToF [ns];Counts",200,-2,2);
    hist_list.push_back(h_DT_mom_bin[i]);
  }
  TH2D * h_mom_DT_theta_phi_bin[6][3];
  for(int i = 0; i < 6; i++){
    for(int j = 0; j < 3; j++){
      sprintf(temp_name,"mom_DT_theta_phi_bin_%d_%d",i+1,j+1);
      h_mom_DT_theta_phi_bin[i][j] = new TH2D(temp_name,"#Delta ToF vs. Momentum;p [GeV];#Delta ToF [ns]",100,0,3,200,-2,2);
      hist_list.push_back(h_mom_DT_theta_phi_bin[i][j]);
    }
  }

  TH2D * h_mom_beta_2212 = new TH2D("mom_beta_2212","#beta vs. Momentum;Momentum [GeV];#beta",100,0,3,100,0,1.1);
  hist_list.push_back(h_mom_beta_2212);
  TH1D * h_mom125_beta_2212 = new TH1D("mom125_beta_2212","#beta;#beta;Counts",100,0.6,1.0);
  hist_list.push_back(h_mom125_beta_2212);

  TH2D * h_mom_DT_wPID = new TH2D("mom_DT_wPID","#Delta ToF vs. Momentum;p [GeV];#Delta ToF [ns]",100,0,3,200,-2,2);
  hist_list.push_back(h_mom_DT_wPID);
  TH2D * h_mom_beta_wPID = new TH2D("mom_beta_wPID","#beta vs. Momentum;Momentum [GeV];#beta",100,0,3,100,0,1.1);
  hist_list.push_back(h_mom_beta_wPID);
  TH1D * h_mom125_beta_wPID = new TH1D("mom125_beta_wPID","#beta;#beta;Counts",100,0.6,1.0);
  hist_list.push_back(h_mom125_beta_wPID);
  TH2D * h_mom_Chi2PID_wPID = new TH2D("mom_Chi2PID_wPID","#chi^{2}_{PID} vs. Momentum [GeV];Momentum [GeV];#chi^{2}_{PID}",100,0,3,100,-10,10);
  hist_list.push_back(h_mom_Chi2PID_wPID);
  TH1D * h_Chi2PID_wPID = new TH1D("Chi2PID_wPID","#chi^{2}_{PID};#chi^{2}_{PID};Counts",100,-10,10);
  hist_list.push_back(h_Chi2PID_wPID);
  

  while(chain.Next())
  //while(chain.Next() && counter<100)
    {
      double weight = 1;
      if(isMC==1){
	weight = c12->mcevent()->getWeight(); //used if MC events have a weight
      }

      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

 clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pionplus = clasAna.getByPid(211);
      auto kaonplus = clasAna.getByPid(321);
      auto deuteronplus = clasAna.getByPid(45);
      auto particles = c12->getDetParticles(); //particles is now 
      if(electrons.size() == 1)
	{
	  SetLorentzVector(el,electrons[0]);
	  //	  SetLorentzVector(ptr,protons[0]);

	  TLorentzVector q = beam - el; //photon  4-vector            
          double Q2        = -q.M2(); // Q^2
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) ); //x-borken
	  double vtz_e = electrons[0]->par()->getVz();
	  
	  ///////////////////////////////
	  //Before cuts
	  ///////////////////////////////
	  for(auto p = particles.begin(); p != particles.end();++p){
	    if((*p)->par()->getCharge()<1){continue;}
	    int hpid = (*p)->getPid();

	    //Momenta
	    SetLorentzVector(lead_ptr,(*p));
	    double mom = lead_ptr.P();
	    double momT = lead_ptr.Perp();
	    double theta = lead_ptr.Theta() * 180 / M_PI;
	    double phi = lead_ptr.Phi() * 180 / M_PI;

	    double beta = (*p)->par()->getBeta();
	    double gamma = 1/sqrt(1-(beta*beta));
	    double path = (*p)->getPath();
	    double hadron_mass = (hpid==45)? mD : db->GetParticle(hpid)->Mass();
	    TLorentzVector hadron_ptr(0,0,0,hadron_mass);
	    SetLorentzVector(hadron_ptr,(*p));	    
	    double DT_proton = (path / (c*beta)) - (path / (c*lead_ptr.Beta()));
	    double DT_hadron = (path / (c*beta)) - (path / (c*hadron_ptr.Beta()));
	    double Dbeta_proton = beta - lead_ptr.Beta();
	    double Chi2PID = (*p)->par()->getChi2Pid();
	    Chi2PID *= DT_proton/DT_hadron;

	    double vtz_p = (*p)->par()->getVz();	    
	    int mom_bin = (mom<0.5)?0:(mom<1.0)?1:(mom<1.5)?2:3;
	    double pbg = lead_ptr.Rho()/(beta*gamma);

	    if(beta<0.2){continue;}
	    if(mom<0.3){continue;}
	    if((*p)->getRegion() != CD){continue;}

	    double edge_first = (*p)->traj(CVT,7)->getEdge();
	    TVector3 hit_first((*p)->traj(CVT,7)->getX(),(*p)->traj(CVT,7)->getY(),(*p)->traj(CVT,7)->getZ());
	    double hp_first = hit_first.Phi()*180/M_PI;
	    int hit_reg_first = hp_first<-90?1:hp_first<30?2:hp_first<150?3:1;

	    double edge_last = (*p)->traj(CVT,12)->getEdge();
	    TVector3 hit_last((*p)->traj(CVT,12)->getX(),(*p)->traj(CVT,12)->getY(),(*p)->traj(CVT,12)->getZ());
	    double hp_last = hit_last.Phi()*180/M_PI;
	    int hit_reg_last = hp_last<-90?1:hp_last<30?2:hp_last<150?3:1;

	    bool pass_fid = false;
	    if((edge_first>0.5) && (edge_last>0.5) && (hit_reg_first == hit_reg_last)){
	      pass_fid = true;
	    }

	    //vertex
	    h_vtz_CD_bc->Fill(vtz_p,weight);
	    h_diffvtz_CD_bc->Fill(vtz_e-vtz_p,weight);
	    h_vtz_e_p_CD_bc->Fill(vtz_e,vtz_p,weight);
	    if(fabs(vtz_e-vtz_p-0.62)>(2*0.86)){continue;}
	    if((vtz_p<-5.5) || (vtz_p>-0.5)){continue;}
	    h_vtz_CD_ac->Fill(vtz_p,weight);
	    h_diffvtz_CD_ac->Fill(vtz_e-vtz_p,weight);
	    h_vtz_e_p_CD_ac->Fill(vtz_e,vtz_p,weight);


	    //fid
	    h_mom_ToFToF_d_ToFMom_CD_bc->Fill(mom,DT_proton,weight);
	    h_theta_CD_bc->Fill(theta,weight);
	    h_phi_momT_CD_bc->Fill(phi,momT,weight);
	    h_ToFToF_d_ToFMom_CD_bc_bin[mom_bin]->Fill(DT_proton,weight);

	    h_edge_first_CD_bc->Fill(edge_first,weight);
	    h_edge_last_CD_bc->Fill(edge_last,weight);

	    if(!pass_fid){
	      h_mom_ToFToF_d_ToFMom_CD_bad->Fill(mom,DT_proton,weight);  		
	      h_ToFToF_d_ToFMom_CD_bad_bin[mom_bin]->Fill(DT_proton,weight);
	      continue;
	    }
	    
	    h_mom_ToFToF_d_ToFMom_CD_ac->Fill(mom,DT_proton,weight);
	    h_ToFToF_d_ToFMom_CD_ac_bin[mom_bin]->Fill(DT_proton,weight);
	    h_theta_CD_ac->Fill(theta,weight);		  
	    h_phi_momT_CD_ac->Fill(phi,momT,weight);


	    //pid
	    h_StartTime->Fill(c12->event()->getStartTime(),weight);
	    h_HitTime->Fill((*p)->getTime(),weight);
	    h_ToF->Fill((*p)->getTime()-c12->event()->getStartTime(),weight);
	    h_Path->Fill((*p)->getPath(),weight);
	    h_ToF_Path->Fill((*p)->getTime(),(*p)->getPath(),weight);
	    h_mom_DT->Fill(mom,DT_proton,weight);
	    h_mom_beta->Fill(mom,beta,weight);
	    if((mom>1.27) && (mom<1.3)){
	      h_mom125_beta->Fill(beta,weight);}
	    h_DT_mom_bin[mom_bin]->Fill(DT_proton,weight);
	    int theta_bin = theta<60?0:theta<80?1:2;
	    int phi_bin = phi<-120?0:phi<-60?1:phi<0?2:phi<60?3:phi<120?4:5;
	    h_mom_DT_theta_phi_bin[phi_bin][theta_bin]->Fill(mom,DT_proton,weight);

	    if(hpid==2212){
	      h_mom_beta_2212->Fill(mom,beta,weight);
	      if((mom>1.27) && (mom<1.3)){
		h_mom125_beta_2212->Fill(beta,weight);}
	    }

	    if(pass_cut(mom,DT_proton,2) && (DT_proton>-0.75)){
	      h_mom_DT_wPID->Fill(mom,DT_proton,weight);
	      h_mom_beta_wPID->Fill(mom,beta,weight);
	      if((mom>1.27) && (mom<1.3)){
		h_mom125_beta_wPID->Fill(beta,weight);}
	      h_mom_Chi2PID_wPID->Fill(mom,Chi2PID,weight);
	      h_Chi2PID_wPID->Fill(Chi2PID,weight);
	    }
	    

	    
	  }	  
	}
    }

  //clasAna.WriteDebugPlots();
  char temp[100];
  TGraph * g_phi[3];
  for(int i = 0; i<3; i++){
    //Make the graphs     
    sprintf(temp,"g_phi_%d",i);
    g_phi[i] = new TGraph();
    g_phi[i]->SetName(temp);
    g_phi[i]->SetLineColor(2);
    for(int j = 0; j < 100; j++){
      double shift = (-60 + i * 120) * (M_PI/180);
      double phi = (M_PI/2) * ((double)j/100) - shift;
      double momT = 0.15 * (1/sin(phi + (M_PI/2) + shift));
      g_phi[i]->SetPoint(g_phi[i]->GetN(),phi*180/M_PI,momT);
    }
    g_phi[i]->Write();
  }

  /*
  TF1 * f_vertex= new TF1("f_vertex","gaus(0)",-1,2);
  f_vertex->SetParameter(0,h_diffvtz_CD_bc->GetMaximum());
  f_vertex->SetParameter(1,0);
  f_vertex->SetParameter(2,1);
  TFitResultPtr p_vertex = h_diffvtz_CD_bc->Fit(f_vertex,"qesrn","",-1,2);

  TGraph * g_vertex = new TGraph();
  g_vertex->SetName("g_vertex");
  g_vertex->SetLineColor(2);
  for(int j = 0; j < 100; j++){
    double x = 0.62 + 0.86 *(-3 + 6*(double)j);
    double y = (39082/(0.86*sqrt(2*M_PI))) * exp(-0.5 * (x-0.62)*(x-0.62)/(0.86*0.86));
    g_vertex->SetPoint(g_vertex->GetN(),x,y);
  }
  */
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i = 0; i < hist_list.size(); i++){
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->SetTitleOffset(1.0);
    hist_list[i]->GetYaxis()->CenterTitle();

    hist_list[i]->GetXaxis()->SetTitleSize(0.08);
    hist_list[i]->GetYaxis()->SetTitleSize(0.05);
    hist_list[i]->GetYaxis()->SetTitleOffset(1.0);
    hist_list[i]->Write();
  }

  //Plot on pdf
  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

  // from ROOT plain style
  myStyle->SetPalette(1,0);
  myStyle->SetOptStat(0);
  // myStyle->SetOptTitle(0);
  myStyle->SetOptDate(0);
  myStyle->SetLabelSize(0.35, "xyz");
  myStyle->SetTitleSize(1.07, "xyz");
  //myStyle->SetTitleOffset(1.0, "xyz");
  //myStyle->SetTitleOffset(4.0, "y");
  myStyle->SetTitleFont(132, "xyz");
  myStyle->SetLabelFont(82, "xyz"); // size of axis values
  myStyle->SetNdivisions(4, "xyz");
  // Some canvas borders and stuff
  myStyle->SetCanvasColor(0);
  myStyle->SetCanvasBorderMode(0);
  myStyle->SetCanvasBorderSize(0);
  myStyle->SetPadBottomMargin(0.2);
  myStyle->SetPadTopMargin(0.07);
  myStyle->SetPadLeftMargin(0.1);
  myStyle->SetPadRightMargin(0.1);
  myStyle->SetPadGridX(0);
  myStyle->SetPadGridY(0);
  myStyle->SetPadTickX(1);
  myStyle->SetPadTickY(1);
  myStyle->SetFrameBorderMode(0);
  myStyle->SetPaperSize(20, 24);
  myStyle->SetPadBorderMode(0);
  // Title styles and histogram styles
  //myStyle->SetTitleStyle(0000);
  myStyle->SetHistFillStyle(3001);
  myStyle->SetHistLineColor(kBlack);
  myStyle->SetHistLineWidth(2); //Style option to make the plots look a certain way
  myStyle->SetHistFillColor(kYellow);
  myStyle->SetTitleSize(0.0002, "t");
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleTextColor(0);
  myStyle->SetTitleAlign(18);
  myStyle->SetPalette(kBird);
  //myStyle->SetLabelSize(3.5);
  myStyle->cd();

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  int pixelx = 2200;
  int pixely = 1600;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);
  /////////////////////////////////////

  /////////////////////////////////////
  //Fiducial Cuts
  /////////////////////////////////////
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_edge_first_CD_bc->Draw();
  TLine *line1 = new TLine(0.5,0,0.5,h_edge_first_CD_bc->GetMaximum());
  line1->SetLineColor(2);
  line1->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_edge_last_CD_bc->Draw();
  TLine *line2 = new TLine(0.5,0,0.5,h_edge_last_CD_bc->GetMaximum());
  line2->SetLineColor(2);
  line2->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phi_momT_CD_bc->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phi_momT_CD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_theta_CD_bc->Draw();
  h_theta_CD_ac->SetLineColor(2);
  h_theta_CD_ac->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  for(int i = 0; i < 4; i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    h_ToFToF_d_ToFMom_CD_ac_bin[i]->Draw("colz");
    h_ToFToF_d_ToFMom_CD_bad_bin[i]->Scale(h_ToFToF_d_ToFMom_CD_ac_bin[i]->Integral()/h_ToFToF_d_ToFMom_CD_bad_bin[i]->Integral());
    h_ToFToF_d_ToFMom_CD_bad_bin[i]->SetLineColor(2);
    h_ToFToF_d_ToFMom_CD_bad_bin[i]->Draw("same");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  /////////////////////////////////////
  //Vertex Cuts
  /////////////////////////////////////
  /*  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_vtz_CD_bc->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_diffvtz_CD_bc->Draw("colz");
  f_vertex->Draw("same");
  text.DrawLatex(-9.5,h_diffvtz_CD_bc->GetMaximum()*0.9,Form("#mu = %g #pm %g",p_vertex->Parameter(1),p_vertex->ParError(1)));
  text.DrawLatex(-9.5,h_diffvtz_CD_bc->GetMaximum()*0.8,Form("#sigma = %g #pm %g",p_vertex->Parameter(2),p_vertex->ParError(2)));
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_vtz_e_p_CD_bc->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_vtz_e_p_CD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */  
  /////////////////////////////////////
  //PID Cuts
  /////////////////////////////////////
  /*myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_StartTime->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_HitTime->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_ToF->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Path->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_ToF_Path->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_DT->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_beta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  for(int i = 0; i < 4; i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    h_DT_mom_bin[i]->Draw();
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_beta_2212->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_DT_wPID->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_beta_wPID->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom125_beta->SetLineColor(1);
  h_mom125_beta->Draw("SAME");
  h_mom125_beta_2212->SetLineColor(2);
  h_mom125_beta_2212->Draw("SAME");
  h_mom125_beta_wPID->SetLineColor(3);
  h_mom125_beta_wPID->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_Chi2PID_wPID->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Chi2PID_wPID->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}


/*
	    if(hpid==211){
	      h_mom_beta_CD_211_ac->Fill(lead_ptr.Rho(),beta);
	    }
	    if(hpid==321){
	      h_mom_beta_CD_321_ac->Fill(lead_ptr.Rho(),beta);
	    }
	    if(hpid==2212){
	      h_mom_beta_CD_2212_ac->Fill(lead_ptr.Rho(),beta);
	    }
	    if(hpid==45){
	      h_mom_beta_CD_45_ac->Fill(lead_ptr.Rho(),beta);
	      h_mom_Chi2PID_CD_45_DThad_ac->Fill(lead_ptr.Rho(),DT_hadron);
	      h_mom_Chi2PID_CD_45_DTpro_ac->Fill(lead_ptr.Rho(),DT_proton);

	      h_mom_Res_CD_45_ac->Fill(lead_ptr.Rho(),Chi2PID/DT_proton);
	      //cout<<Chi2PID/DT_proton<<endl;
	    }

  TH2D * h_mom_beta_CD_211_ac = new TH2D("mom_beta_CD_211_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_211_ac);
  TH2D * h_mom_beta_CD_321_ac = new TH2D("mom_beta_CD_321_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_321_ac);
  TH2D * h_mom_beta_CD_2212_ac = new TH2D("mom_beta_CD_2212_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_2212_ac);
  TH2D * h_mom_beta_CD_45_ac = new TH2D("mom_beta_CD_45_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_45_ac);
  TH2D * h_mom_Chi2PID_CD_45_DThad_ac = new TH2D("mom_Chi2PID_CD_45_DThad_ac","#Delta Time Deut CD",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_Chi2PID_CD_45_DThad_ac);
  TH2D * h_mom_Chi2PID_CD_45_DTpro_ac = new TH2D("mom_Chi2PID_CD_45_DTpro_ac","#Delta Time Deut CD",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_Chi2PID_CD_45_DTpro_ac);
  TH2D * h_mom_Res_CD_45_ac = new TH2D("mom_Res_CD_45_ac","#Delta Time Deut CD",100,0,3,100,-20,20);
  hist_list.push_back(h_mom_Res_CD_45_ac);
}


  //TH2D * h_mom_Dbeta_CD_bc = new TH2D("mom_Dbeta_CD_bc","#Delta #beta vs. Momentum;p;#beta",100,0,3,100,-1,1);
  //hist_list.push_back(h_mom_Dbeta_CD_bc);
  //TH2D * h_mom_DT_CD_bc = new TH2D("mom_DT_CD_bc","#Delta Time vs. Momentum;p;#Delta time",100,0,3,100,-1,1);
  //hist_list.push_back(h_mom_DT_CD_bc);
*/
