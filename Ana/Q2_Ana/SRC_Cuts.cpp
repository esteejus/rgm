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
const double mN = 0.938272;

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
  //clasAna.readEcalSFPar("../paramsSF_40Ca_x2.dat");
  //clasAna.readEcalPPar("../paramsPI_40Ca_x2.dat");
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

  TH1D * h_xb_fd = new TH1D("xb_fd","x-Bjorken x_{B};x_{B};Counts",100,0,3);
  TH1D * h_pmiss_fd = new TH1D("pmiss_fd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0,1.2);
  TH1D * h_q2_fd = new TH1D("q2_fd","Q^{2};Q^{2} (GeV^{2});Counts",100,0,4);
  TH2D * h_thetapq_pq_fd = new TH2D("thetapq_pq_fd","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_fd = new TH1D("mmiss_fd","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  TH2D * h_p_theta_fd = new TH2D("p_theta_fd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p (GeV/c)",180,0,180,100,0,2.5);


  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  

  TH1D * h_xb_cd = new TH1D("xb_cd","x-Bjorken x_{B};x_{B};Counts",100,0,3);
  TH1D * h_pmiss_cd = new TH1D("pmiss_cd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0,1.2);
  TH1D * h_q2_cd = new TH1D("q2_cd","Q^{2};Q^{2} (GeV^{2});Counts",100,0,4);
  TH2D * h_thetapq_pq_cd = new TH2D("thetapq_pq_cd","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_cd = new TH1D("mmiss_cd","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  TH2D * h_p_theta_cd = new TH2D("p_theta_cd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p (GeV/c)",180,0,180,100,0,2.5);



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

	    if(beta<0.2){continue;} // proton cut

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
		if(DCedge[k]<10){pass_fiducial=false;}
	      }
	      if(!pass_fiducial){continue;} // proton cut

              // calculate SRC kinematics
              TVector3 pmiss = lead_ptr.Vect() - q.Vect();
              double thetapq = lead_ptr.Vect().Angle(q.Vect())*180./M_PI;
              double pq = (lead_ptr.Vect().Mag()) / (q.Vect().Mag());
              double mmiss = (q + TLorentzVector(TVector3(0.,0.,0.),2*mN) - lead_ptr).Mag();

              // SRC histograms here
              h_xb_fd->Fill(xB);
              if (xB<1.2) {continue;}
              h_pmiss_fd->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.3 || pmiss.Mag()>1.0) {continue;}
              h_q2_fd->Fill(Q2);
              if (Q2<1.5) {continue;}
              h_thetapq_pq_fd->Fill(pq,thetapq);
              if (thetapq>25) {continue;}
              if (pq<0.62 || pq>0.96) {continue;}
              h_mmiss_fd->Fill(mmiss);
              if (mmiss>1.1) {continue;}
              h_p_theta_fd->Fill(theta,mom);


	    }
	    else if((*p)->getRegion() == CD){

	      if(!CD_fiducial(phi,theta,momT)){
		continue;
	      }

              // calculate SRC kinematics
              TVector3 pmiss = lead_ptr.Vect() - q.Vect();
              double thetapq = lead_ptr.Vect().Angle(q.Vect())*180./M_PI;
              double pq = (lead_ptr.Vect().Mag()) / (q.Vect().Mag());
              double mmiss = (q + TLorentzVector(TVector3(0.,0.,0.),2*mN) - lead_ptr).Mag();

              // SRC histograms here
              h_xb_cd->Fill(xB);
              if (xB<1.2) {continue;}
              h_pmiss_cd->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.3 || pmiss.Mag()>1.0) {continue;}
              h_q2_cd->Fill(Q2);
              if (Q2<1.5) {continue;}
              h_thetapq_pq_cd->Fill(pq,thetapq);
              if (thetapq>25) {continue;}
              if (pq<0.62 || pq>0.96) {continue;}
              h_mmiss_cd->Fill(mmiss);
              if (mmiss>1.1) {continue;}
              h_p_theta_cd->Fill(theta,mom);
	    }
	    else{
	      cout<<"Not Either"<<endl;
	    }	    
	  }
	  

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

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xb_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_theta_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 


  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xb_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_cd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_theta_cd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 


  
  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

