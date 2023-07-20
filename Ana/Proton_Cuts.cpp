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
  clasAna.readInputParam("ana.par");
  clasAna.readEcalSFPar("paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("paramsPI_40Ca_x2.dat");
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

  TH1D * h_Q2_bc = new TH1D("Q2_bc","Q^{2} ",1000,0, 5);
  TH1D * h_xB_bc = new TH1D("xB_bc","x_{B} ",1000,0, 2);
  TH2D * h_phi_theta_bc = new TH2D("phi_theta_bc","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);

  ///////////////////////////////////////////////////////
  //Forward Detector
  ///////////////////////////////////////////////////////  
  //Fid
  TH1D *h_DCedge_weight_FD_bc[3][6];
  TH1D *h_DCedge_FD_bc[3][6];
  for(int j=0; j<3; j++){
    for(int i=0; i<6; i++){
      sprintf(temp_name,"DCedge_weight_FD_bc_%d_%d",j+1,i+1);
      sprintf(temp_title,"Distance from DC Edge Sector=%d Region=%d;Distance [cm];Average #chi^{2}/DoF",j+1,i+1);
      h_DCedge_weight_FD_bc[j][i] = new TH1D(temp_name,temp_title,50,0,50);
      h_DCedge_weight_FD_bc[j][i]->Sumw2();
    
      sprintf(temp_name,"DCedge_FD_bc_%d_%d",j+1,i+1);
      h_DCedge_FD_bc[j][i] = new TH1D(temp_name,temp_title,50,0,50);
      h_DCedge_FD_bc[j][i]->Sumw2();
    }
  }

  //Vertex
  TH1D * h_vtz_FD_bc = new TH1D("vtz_FD_bc","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  TH1D * h_diffvtz_FD_bc = new TH1D("diffvtz_FD_bc","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  TH2D * h_vtz_e_p_FD_bc = new TH2D("vtz_e_p_FD_bc","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);

  TH1D * h_vtz_FD_ac = new TH1D("vtz_FD_ac","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  TH1D * h_diffvtz_FD_ac = new TH1D("diffvtz_FD_ac","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  TH2D * h_vtz_e_p_FD_ac = new TH2D("vtz_e_p_FD_ac","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);
  
  //PID
  TH2D * h_mom_beta_FD_bc = new TH2D("mom_beta_FD_bc","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  TH1D * h_Chi2PID_FD_bc = new TH1D("Chi2PID_FD_bc","#chi^{2} PID FD",100,-10,10);
  TH2D * h_mom_Chi2PID_FD_bc = new TH2D("mom_Chi2PID_FD_bc","#Delta Time FD",100,0,3,100,-10,10);

  TH2D * h_mom_beta_FD_ac = new TH2D("mom_beta_FD_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  TH1D * h_Chi2PID_FD_ac = new TH1D("Chi2PID_FD_ac","#chi^{2} PID FD",100,-10,10);
  TH2D * h_mom_Chi2PID_FD_ac = new TH2D("mom_Chi2PID_FD_ac","#Delta Time FD",100,0,3,100,-10,10);

  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  
  //Fid
  TH1D * h_theta_CD_bc = new TH1D("theta_CD_bc","#theta vs. Counts; #theta; Counts",100,0,180);
  TH2D * h_phi_momT_CD_bc = new TH2D("phi_momT_CD_bc","#phi vs. p_{T}; #phi; p_{T}",100,-180,180,100,0,2);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_bc = new TH2D("mom_ToFToF_d_ToFMom_CD_bc","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);

  TH1D * h_theta_CD_ac = new TH1D("theta_CD_ac","#theta vs. Counts; #theta; Counts",100,0,180);
  TH2D * h_phi_momT_CD_ac = new TH2D("phi_momT_CD_ac","#phi vs. p_{T}; #phi; p_{T}",100,-180,180,100,0,2);
  TH2D * h_mom_ToFToF_d_ToFMom_CD_ac = new TH2D("mom_ToFToF_d_ToFMom_CD_ac","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);

  TH2D * h_mom_ToFToF_d_ToFMom_CD_bad = new TH2D("mom_ToFToF_d_ToFMom_CD_bad","Momentum vs. ToF_{ToF} - ToF_{p};Momentum;ToF_{ToF} - ToF_{p}",100,0,3,100,-1,1);

  //Vertex
  TH1D * h_vtz_CD_bc = new TH1D("vtz_CD_bc","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  TH1D * h_diffvtz_CD_bc = new TH1D("diffvtz_CD_bc","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  TH2D * h_vtz_e_p_CD_bc = new TH2D("vtz_e_p_CD_bc","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);

  TH1D * h_vtz_CD_ac = new TH1D("vtz_CD_ac","Proton Z Vertex;Vertex_{p} [cm];Counts",100,-10,10);
  TH1D * h_diffvtz_CD_ac = new TH1D("diffvtz_CD_ac","Electron Minus Proton Z Vertex;Vertex_{e} - Vertex_{p} [cm];Counts",100,-10,10);
  TH2D * h_vtz_e_p_CD_ac = new TH2D("vtz_e_p_CD_ac","Electron Vertex vs. Proton Vertex;Vertex_{e} [cm];Vertex_{p} [cm]",100,-10,10,100,-10,10);
  
  //PID
  TH2D * h_mom_beta_CD_bc = new TH2D("mom_beta_CD_bc","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  TH1D * h_Chi2PID_CD_bc = new TH1D("Chi2PID_CD_bc","#chi^{2} PID CD",100,-10,10);
  TH2D * h_mom_Chi2PID_CD_bc = new TH2D("mom_Chi2PID_CD_bc","#Delta Time CD",100,0,3,100,-10,10);

  TH2D * h_mom_beta_CD_ac = new TH2D("mom_beta_CD_ac","Momentum vs. #beta ;p;#beta",100,0,2,100,0,1.1);
  TH1D * h_Chi2PID_CD_ac = new TH1D("Chi2PID_CD_ac","#chi^{2} PID CD",100,-10,10);
  TH2D * h_mom_Chi2PID_CD_ac = new TH2D("mom_Chi2PID_CD_ac","#Delta Time CD",100,0,3,100,-10,10);
  

  clasAna.setEcalSFCuts();
  clasAna.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();

  clasAna.setVertexCuts();
  //clasAna.setVertexCorrCuts();
  clasAna.setDCEdgeCuts();
  
  //clasAna.setVzcuts(-6,1);
  //clasAna.setVertexCorrCuts(-3,1);

  while(chain.Next())
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

	    if(beta<0.2){continue;}

	    if((*p)->getRegion() == FD){
	      int psector = (*p)->getSector();
	      double DCedge[3];
	      DCedge[0]  = (*p)->traj(DC,6 )->getFloat("edge",(*p)->traj(DC,6 )->getIndex());
	      DCedge[1]  = (*p)->traj(DC,18)->getFloat("edge",(*p)->traj(DC,18)->getIndex());
	      DCedge[2]  = (*p)->traj(DC,36)->getFloat("edge",(*p)->traj(DC,36)->getIndex());
	      double Chi2DoF = (*p)->trk(DC)->getChi2()/(*p)->trk(DC)->getNDF();
	      //fid
	      bool pass_fiducial = true;
	      for(int k=0; k<3; k++){
		h_DCedge_FD_bc[k][psector-1]->Fill(DCedge[k]);
		h_DCedge_weight_FD_bc[k][psector-1]->Fill(DCedge[k],Chi2DoF);
		if(DCedge[k]<10){pass_fiducial=false;}
	      }
	      if(!pass_fiducial){continue;}
	      //vertex
	      h_vtz_FD_bc->Fill(vtz_p);
	      h_diffvtz_FD_bc->Fill(vtz_e-vtz_p);
	      h_vtz_e_p_FD_bc->Fill(vtz_e,vtz_p);
	      if(fabs(vtz_e-vtz_p)>2){continue;}
	      h_vtz_FD_ac->Fill(vtz_p);
	      h_diffvtz_FD_ac->Fill(vtz_e-vtz_p);
	      h_vtz_e_p_FD_ac->Fill(vtz_e,vtz_p);

	      //pid
	      h_mom_beta_FD_bc->Fill(lead_ptr.Rho(),beta);
	      h_Chi2PID_FD_bc->Fill(Chi2PID);
	      h_mom_Chi2PID_FD_bc->Fill(lead_ptr.Rho(),Chi2PID);
	      if(fabs(Chi2PID)>4){continue;}
	      h_mom_beta_FD_ac->Fill(lead_ptr.Rho(),beta);
	      h_Chi2PID_FD_ac->Fill(Chi2PID);
	      h_mom_Chi2PID_FD_ac->Fill(lead_ptr.Rho(),Chi2PID);


	    }
	    else if((*p)->getRegion() == CD){

	      //fid
	      h_mom_ToFToF_d_ToFMom_CD_bc->Fill(mom,DT_proton);  
	      h_theta_CD_bc->Fill(theta);
	      h_phi_momT_CD_bc->Fill(phi,momT);
	      if(!CD_fiducial(phi,theta,momT)){
		h_mom_ToFToF_d_ToFMom_CD_bad->Fill(mom,DT_proton);  		
		continue;
	      }
	      	      
	      h_mom_ToFToF_d_ToFMom_CD_ac->Fill(mom,DT_proton);
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
	      h_mom_Chi2PID_CD_bc->Fill(lead_ptr.Rho(),Chi2PID);
	      if(fabs(Chi2PID)>4){continue;}
	      h_mom_beta_CD_ac->Fill(lead_ptr.Rho(),beta);
	      h_Chi2PID_CD_ac->Fill(Chi2PID);
	      h_mom_Chi2PID_CD_ac->Fill(lead_ptr.Rho(),Chi2PID);

	    }
	    else{
	      cout<<"Not Either"<<endl;
	    }	    
	  }
	  
	  ///////////////////////////////
	  //After cuts
	  ///////////////////////////////

	  /*
	  for(auto p = protons.begin(); p != protons.end();++p){
	    if((*p)->par()->getCharge()<1){continue;}

	    double Chi2PID = (*p)->par()->getChi2Pid();
	    SetLorentzVector(lead_ptr,(*p));
	    double path = (*p)->getPath();
	    double beta = (*p)->par()->getBeta();
	    double beta_frommom = lead_ptr.Beta();
	    double time_frommom = path / (c*beta_frommom);
	    double td = ((*p)->getTime()-c12->event()->getStartTime()) - time_frommom;

	    if((*p)->getRegion() == FD){
	      h_mom_beta_FD_ac->Fill(lead_ptr.Rho(),beta);
	      h_Chi2PID_FD_ac->Fill(Chi2PID);
	      h_TimeDiff_FD_ac->Fill(td);
	    }
	    else if((*p)->getRegion() == CD){
	      h_mom_beta_CD_ac->Fill(lead_ptr.Rho(),beta);
	      h_Chi2PID_CD_ac->Fill(Chi2PID);
	      h_TimeDiff_CD_ac->Fill(td);
	    }
	    else{
	      cout<<"Not Either"<<endl;
	    }	    
	  }
	  */
	}
    }

  //clasAna.WriteDebugPlots();

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  h_Q2_bc->Write();
  h_xB_bc->Write();


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

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Q2_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xB_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phi_theta_bc->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  ///////////////////////////////////////////////////////
  //Forward Detector
  ///////////////////////////////////////////////////////  
  //fid
  for(int j = 0; j < 3; j++){
    myCanvas->Divide(2,3);
    for(int i = 0; i < 6; i++){
      myCanvas->cd(i+1);
      h_DCedge_weight_FD_bc[j][i]->Divide(h_DCedge_FD_bc[j][i]);
      h_DCedge_weight_FD_bc[j][i]->Draw();  
      h_DCedge_weight_FD_bc[j][i]->SetMaximum(100);
      h_DCedge_weight_FD_bc[j][i]->SetMinimum(0);
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  //vertex
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_diffvtz_FD_bc->Draw("colz");
  myCanvas->cd(2);
  h_vtz_e_p_FD_bc->Draw("colz");
  myCanvas->cd(3);
  h_vtz_FD_bc->Draw("colz");
  myCanvas->cd(4);
  h_vtz_FD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  //pid
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_Chi2PID_FD_bc->Draw("colz");
  myCanvas->cd(2);
  h_mom_Chi2PID_FD_bc->Draw("colz");
  myCanvas->cd(3);
  h_mom_beta_FD_bc->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_FD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  
  //fid
  myCanvas->Divide(2,1);
  myCanvas->cd(1);
  h_mom_ToFToF_d_ToFMom_CD_ac->Draw("colz");
  myCanvas->cd(2);
  h_mom_ToFToF_d_ToFMom_CD_bad->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_CD_bc->Draw("colz");
  myCanvas->cd(2);
  h_theta_CD_ac->Draw("colz");
  myCanvas->cd(3);
  h_phi_momT_CD_bc->Draw("colz");
  myCanvas->cd(4);
  h_phi_momT_CD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  //vertex
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_diffvtz_CD_bc->Draw("colz");
  myCanvas->cd(2);
  h_vtz_e_p_CD_bc->Draw("colz");
  myCanvas->cd(3);
  h_vtz_CD_bc->Draw("colz");
  myCanvas->cd(4);
  h_vtz_CD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  //pid
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_Chi2PID_CD_bc->Draw("colz");
  myCanvas->cd(2);
  h_mom_Chi2PID_CD_bc->Draw("colz");
  myCanvas->cd(3);
  h_mom_beta_CD_bc->Draw("colz");
  myCanvas->cd(4);
  h_mom_beta_CD_ac->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

