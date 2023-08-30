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

using namespace std;
using namespace clas12;

const double c = 29.9792458;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}

bool CD_fiducial(double phi, double theta, double momT){
  bool pass_fiducial = true;
  double fiducial_phi_width = 10;
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
  std::cerr << "Usage: ./code outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }



  TString outFile = argv[1];
  char * pdfFile = argv[2];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;

  //Read in target parameter files                                                                                                                                                           
  clasAna.readInputParam("../ana.par");
  clasAna.readEcalSFPar("../paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("../paramsPI_40Ca_x2.dat");
  clasAna.printParams();
    

  clas12root::HipoChain chain;
  for(int k = 3; k < argc; k++){
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
  TH1D * h_Q2_bc = new TH1D("Q2_bc","Q^{2} ",1000,0, 5);
  hist_list.push_back(h_Q2_bc);
  TH1D * h_xB_bc = new TH1D("xB_bc","x_{B} ",1000,0, 2);
  hist_list.push_back(h_xB_bc);
  TH2D * h_phi_theta_bc = new TH2D("phi_theta_bc","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);
  hist_list.push_back(h_phi_theta_bc);

  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  
  //Fid
  TH1D * h_theta_CD_bc = new TH1D("theta_CD_bc","#theta vs. Counts; #theta; Counts",100,0,180);
  hist_list.push_back(h_theta_CD_bc);
  TH2D * h_phi_momT_CD_bc = new TH2D("phi_momT_CD_bc","#phi vs. p_{T}; #phi; p_{T}",100,-180,180,100,0,2);
  hist_list.push_back(h_phi_momT_CD_bc);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_bc = new TH2D("mom_ToFToF_d_ToFMom_CD_bc","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_bc);
  TH1D * h_ToFToF_d_ToFMom_CD_bc_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_bc_bin_%d",i+1);
    sprintf(temp_title,"ToF_{ToF} - ToF_{p} Bin=%d;ToF_{ToF} - ToF_{p};Counts",i+1);
    h_ToFToF_d_ToFMom_CD_bc_bin[i] = new TH1D(temp_name,temp_title,100,-1,1);
    hist_list.push_back(h_ToFToF_d_ToFMom_CD_bc_bin[i]);
  }

  TH1D * h_theta_CD_ac = new TH1D("theta_CD_ac","#theta vs. Counts; #theta; Counts",100,0,180);
  hist_list.push_back(h_theta_CD_ac);
  TH2D * h_phi_momT_CD_ac = new TH2D("phi_momT_CD_ac","#phi vs. p_{T}; #phi; p_{T}",100,-180,180,100,0,2);
  hist_list.push_back(h_phi_momT_CD_ac);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_ac = new TH2D("mom_ToFToF_d_ToFMom_CD_ac","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_ac);
  TH1D * h_ToFToF_d_ToFMom_CD_ac_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_ac_bin_%d",i+1);
    sprintf(temp_title,"ToF_{ToF} - ToF_{p} Bin=%d;ToF_{ToF} - ToF_{p};Counts",i+1);
    h_ToFToF_d_ToFMom_CD_ac_bin[i] = new TH1D(temp_name,temp_title,100,-1,1);
    hist_list.push_back(h_ToFToF_d_ToFMom_CD_ac_bin[i]);
  }

  TH2D * h_mom_ToFToF_d_ToFMom_CD_bad = new TH2D("mom_ToFToF_d_ToFMom_CD_bad","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_ToFToF_d_ToFMom_CD_bad);
  TH1D * h_ToFToF_d_ToFMom_CD_bad_bin[4];
  for(int i = 0; i < 4; i++){
    sprintf(temp_name,"ToFToF_d_ToFMom_CD_bad_bin_%d",i+1);
    sprintf(temp_title,"ToF_{ToF} - ToF_{p} Bin=%d;ToF_{ToF} - ToF_{p};Counts",i+1);
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
  TH2D * h_mom_beta_CD_bc = new TH2D("mom_beta_CD_bc","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_bc);
  TH1D * h_Chi2PID_CD_bc = new TH1D("Chi2PID_CD_bc","#chi^{2} PID CD",100,-10,10);
  hist_list.push_back(h_Chi2PID_CD_bc);
  TH2D * h_mom_DT_CD_bc = new TH2D("mom_DT_CD_bc","#Delta Time CD",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_DT_CD_bc);
  TH2D * h_mom_Chi2PID_CD_bc = new TH2D("mom_Chi2PID_CD_bc","#chi^2 PID CD",100,0,3,100,-10,10);
  hist_list.push_back(h_mom_Chi2PID_CD_bc);

  TH2D * h_mom_beta_CD_ac = new TH2D("mom_beta_CD_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  hist_list.push_back(h_mom_beta_CD_ac);
  TH1D * h_Chi2PID_CD_ac = new TH1D("Chi2PID_CD_ac","#chi^{2} PID CD",100,-10,10);
  hist_list.push_back(h_Chi2PID_CD_ac);
  TH2D * h_mom_DT_CD_ac = new TH2D("mom_DT_CD_ac","#Delta Time CD",100,0,3,100,-1,1);
  hist_list.push_back(h_mom_DT_CD_ac);
  TH2D * h_mom_Chi2PID_CD_ac = new TH2D("mom_Chi2PID_CD_ac","#chi^2 PID CD",100,0,3,100,-10,10);
  hist_list.push_back(h_mom_Chi2PID_CD_ac);
  

  clasAna.setEcalSFCuts();
  clasAna.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();

  clasAna.setVertexCuts();
  //clasAna.setVertexCorrCuts();
  clasAna.setDCEdgeCuts();
  
  //clasAna.setVzcuts(-6,1);
  //clasAna.setVertexCorrCuts(-3,1);

  //while(chain.Next())
  while(chain.Next() && counter<10000)
    {

      double weight = c12->mcevent()->getWeight(); //used if MC events have a weight 

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
	  h_Q2_bc->Fill(Q2);
	  h_xB_bc->Fill(xB);
	  h_phi_theta_bc->Fill(el.Phi()*180/M_PI,el.Theta()*180/M_PI);
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
	    double path = (*p)->getPath();
	    double hadron_mass = (hpid==45)? mD : db->GetParticle(hpid)->Mass();
	    TLorentzVector hadron_ptr(0,0,0,hadron_mass);
	    SetLorentzVector(hadron_ptr,(*p));	    
	    double DT_proton = (path / (c*beta)) - (path / (c*lead_ptr.Beta()));
	    double DT_hadron = (path / (c*beta)) - (path / (c*hadron_ptr.Beta()));
	    double Chi2PID = (*p)->par()->getChi2Pid();
	    Chi2PID *= DT_proton/DT_hadron;

	    double vtz_p = (*p)->par()->getVz();	    
	    int mom_bin = (mom<0.5)?0:(mom<1.0)?1:(mom<1.5)?2:3;
	    
	    if(beta<0.2){continue;}

	    if((*p)->getRegion() != CD){continue;}

	    //fid
	    h_mom_ToFToF_d_ToFMom_CD_bc->Fill(mom,DT_proton);  
	    h_theta_CD_bc->Fill(theta);
	    h_phi_momT_CD_bc->Fill(phi,momT);
	    h_ToFToF_d_ToFMom_CD_bc_bin[mom_bin]->Fill(DT_proton);
	    if(!CD_fiducial(phi,theta,momT)){
	      h_mom_ToFToF_d_ToFMom_CD_bad->Fill(mom,DT_proton);  		
	      h_ToFToF_d_ToFMom_CD_bad_bin[mom_bin]->Fill(DT_proton);
	      continue;
	    }
	    
	    h_mom_ToFToF_d_ToFMom_CD_ac->Fill(mom,DT_proton);
	    h_ToFToF_d_ToFMom_CD_ac_bin[mom_bin]->Fill(DT_proton);
	    h_theta_CD_ac->Fill(theta);		  
	    h_phi_momT_CD_ac->Fill(phi,momT);

	    //vertex
	    h_vtz_CD_bc->Fill(vtz_p);
	    h_diffvtz_CD_bc->Fill(vtz_e-vtz_p);
	    h_vtz_e_p_CD_bc->Fill(vtz_e,vtz_p);
	    if(fabs(vtz_e-vtz_p)>2){continue;}
	    h_vtz_CD_ac->Fill(vtz_p);
	    h_diffvtz_CD_ac->Fill(vtz_e-vtz_p);
	    h_vtz_e_p_CD_ac->Fill(vtz_e,vtz_p);
	    
	    //pid
	    h_mom_beta_CD_bc->Fill(lead_ptr.Rho(),beta);
	    h_Chi2PID_CD_bc->Fill(Chi2PID);
	    h_mom_DT_CD_bc->Fill(lead_ptr.Rho(),DT_proton);
	    h_mom_Chi2PID_CD_bc->Fill(lead_ptr.Rho(),Chi2PID);
	    if(fabs(DT_proton)>0.4){continue;}
	    if(fabs(Chi2PID)>3.5){continue;}


	    h_mom_beta_CD_ac->Fill(lead_ptr.Rho(),beta);
	    h_Chi2PID_CD_ac->Fill(Chi2PID);
	    h_mom_DT_CD_ac->Fill(lead_ptr.Rho(),DT_proton);
	    h_mom_Chi2PID_CD_ac->Fill(lead_ptr.Rho(),Chi2PID);

	    
	  }
	  
	}
    }

  //clasAna.WriteDebugPlots();

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  h_Q2_bc->Write();
  h_xB_bc->Write();


  TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

  // from ROOT plain style
  myStyle->SetPalette(1,0);
  myStyle->SetOptStat(0);
  // myStyle->SetOptTitle(0);
  myStyle->SetOptDate(0);
  myStyle->SetLabelSize(0.35, "xyz");
  myStyle->SetTitleSize(0.07, "xyz");
  myStyle->SetTitleOffset(1.0, "xyz");
  myStyle->SetTitleOffset(0.0, "y");
  myStyle->SetTitleFont(132, "xyz");
  myStyle->SetLabelFont(82, "xyz"); // size of axis values
  myStyle->SetNdivisions(4, "xyz");
  // Some canvas borders and stuff
  myStyle->SetCanvasDefW(2000);
  myStyle->SetCanvasDefH(1500);
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
  //myStyle->SetTitleFont(132, "t");
  myStyle->SetTitleBorderSize(0);
  myStyle->SetTitleFillColor(0);
  myStyle->SetTitleTextColor(0);
  //myStyle->SetTitleH(0.12);
  //myStyle->SetTitleX(0.16);
  myStyle->SetTitleAlign(18);
  myStyle->SetPalette(kBird);
  //myStyle->SetLabelSize(3.5);
  myStyle->cd();

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  int pixelx = 2000;
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
  
  for(int i = 0; i < hist_list.size(); i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    hist_list[i]->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }
  
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
*/
