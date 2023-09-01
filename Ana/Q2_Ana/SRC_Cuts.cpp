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
#include <TLine.h>
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



  // ELECTRON AND PROTON CUTS
  clas12ana clasAna;

  //Read in target parameter files
  // before run 15542 - D
  clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/Ana/cutFiles/paramsSF_LD2_x2.dat");
  clasAna.readEcalPPar("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/Ana/cutFiles/paramsPI_LD2_x2.dat");
  // after run 15542 - 40Ca
  //clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/Ana/cutFiles/paramsSF_40Ca_x2.dat");
  //clasAna.readEcalPPar("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/Ana/cutFiles/paramsPI_40Ca_x2.dat");

  //clasAna.printParams();

  clasAna.setVzcuts(-6,1);
  clasAna.setVertexCorrCuts(-3,1);


    

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
  TH2D * h_qmag_qtheta = new TH2D("qmag_qtheta","Q Magnitude vs Theta;q Magnitude (GeV);#theta_{q} (deg)",100,0,4,100,0,100);


  ///////////////////////////////////////////////////////
  //Forward Detector
  ///////////////////////////////////////////////////////  

  TH2D * h_thetae_q2_fd = new TH2D("thetae_q2_fd","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_mmiss_xb_fd = new TH2D("mmiss_xb_fd","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_fd = new TH2D("mmiss_q2_fd","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_fd = new TH1D("mmiss_nocuts_fd","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2});Counts",100,0,2);
  TH1D * h_xb_fd = new TH1D("xb_fd","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_pmiss_fd = new TH1D("pmiss_fd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0,1.2);
  TH1D * h_q2_fd = new TH1D("q2_fd","Q^{2};Q^{2} (GeV^{2});Counts",100,0,4);
  TH2D * h_mmiss_thetapq_fd = new TH2D("mmiss_thetapq_fd","Missing Mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_fd = new TH2D("mmiss_pq_fd","Missing Mass vs p/q;p/q;M_{miss} (GeV/c^{2})",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_fd = new TH2D("thetapq_pq_fd","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_fd = new TH1D("mmiss_fd","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  TH2D * h_p_theta_fd = new TH2D("p_theta_fd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p (GeV/c)",180,0,180,100,0,2.5);






  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  

  TH2D * h_thetae_q2_cd = new TH2D("thetae_q2_cd","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_mmiss_xb_cd = new TH2D("mmiss_xb_cd","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_cd = new TH2D("mmiss_q2_cd","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_cd = new TH1D("mmiss_nocuts_cd","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2});Counts",100,0,2);
  TH1D * h_xb_cd = new TH1D("xb_cd","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_pmiss_cd = new TH1D("pmiss_cd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0,1.2);
  TH1D * h_q2_cd = new TH1D("q2_cd","Q^{2};Q^{2} (GeV^{2});Counts",100,0,4);
  TH2D * h_mmiss_thetapq_cd = new TH2D("mmiss_thetapq_cd","Missing mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_cd = new TH2D("mmiss_pq_cd","Mmiss vs pq;p/q;Mmiss",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_cd = new TH2D("thetapq_pq_cd","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_cd = new TH1D("mmiss_cd","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  TH2D * h_p_theta_cd = new TH2D("p_theta_cd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p (GeV/c)",180,0,180,100,0,2.5);

  ///////////////////////////////////////////////////////
  //Both Forward Detector & Central Detector
  ///////////////////////////////////////////////////////  

  TH2D * h_thetae_q2_all = new TH2D("thetae_q2_all","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_thetap_q2_all = new TH2D("thetap_q2_all","Proton Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{p} (deg)",100,0,4,60,0,60);
  TH1D * h_xb_all = new TH1D("xb_all","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH2D * h_mmiss_xb_all = new TH2D("mmiss_xb_all","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_all = new TH2D("mmiss_q2_all","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_all = new TH1D("mmiss_nocuts_all","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2});Counts",100,0,2);
  TH1D * h_pmiss_all = new TH1D("pmiss_all","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0,1.2);
  TH1D * h_q2_all = new TH1D("q2_all","Q^{2};Q^{2} (GeV^{2});Counts",100,0,4);
  TH2D * h_mmiss_thetapq_all = new TH2D("mmiss_thetapq_all","Missing Mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_all = new TH2D("mmiss_pq_all","Missing Mass vs p/q;p/q;M_{miss} (GeV/c^{2})",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_all = new TH2D("thetapq_pq_all","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_all = new TH1D("mmiss_all","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);

  TH2D * h_p_theta_all = new TH2D("p_theta_all","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p (GeV/c)",180,0,180,100,0,2.5);
  TH2D * h_emiss_omega_all = new TH2D("emiss_omega_all","Missing Energy vs #omega;#omega (GeV);E_{miss} (GeV)",100,0,2.5,100,0,0.5);
  
  TH1D * h_pmiss_src_all = new TH1D("pmiss_src_all","Missing Momentum of Lead SRC Protons;p_{miss} (GeV/c);Counts",50,0.35,1);
  TH2D * h_xb_q2_src_all = new TH2D("xb_q2_src_all","x_{B} vs Q^{2} of Lead SRC Protons;Q^{2} (GeV^{2});x_{B}",50,1.5,4,50,1.2,2.5);
  /*TH2D * h_q2_omega_src_all = new TH2D("emiss_omega","Missing Energy vs #omega;#omega (GeV);E_{miss}"
  TH2D * h_mmiss_emiss_src_all*/
  


  
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
          h_qmag_qtheta->Fill(-1*q.Mag(),q.Vect().Theta()*180./M_PI);
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
	    double vtz_p = (*p)->par()->getVz();

            // calculate SRC kinematics
            TVector3 pmiss = lead_ptr.Vect() - q.Vect();
            double thetapq = lead_ptr.Vect().Angle(q.Vect())*180./M_PI;
            double pq = (lead_ptr.Vect().Mag()) / (q.Vect().Mag());
            double mmiss = (q + TLorentzVector(TVector3(0.,0.,0.),2*mN) - lead_ptr).Mag();

	    if(beta<0.2){continue;} // proton cut



            // FORWARD DETECTOR PROTONS
	    if((*p)->getRegion() == FD){

              h_thetae_q2_fd->Fill(Q2,el.Theta()*180./M_PI);
              h_thetae_q2_all->Fill(Q2,el.Theta()*180./M_PI);
              h_thetap_q2_all->Fill(Q2,theta);

              // SRC histograms here - FD
              h_mmiss_nocuts_fd->Fill(mmiss);
              h_mmiss_xb_fd->Fill(xB,mmiss);
              h_mmiss_nocuts_all->Fill(mmiss);
              h_mmiss_xb_all->Fill(xB,mmiss);
              h_xb_fd->Fill(xB);
            h_xb_all->Fill(xB);
              if (xB<1.2) {continue;}


              h_mmiss_q2_fd->Fill(Q2,mmiss);
            h_mmiss_q2_all->Fill(Q2,mmiss);

              h_q2_fd->Fill(Q2);
            h_q2_all->Fill(Q2);
              if (Q2<1.5) {continue;}


              h_mmiss_thetapq_fd->Fill(thetapq,mmiss);
              h_mmiss_pq_fd->Fill(pq,mmiss);
              h_thetapq_pq_fd->Fill(pq,thetapq);
            h_thetapq_pq_all->Fill(pq,thetapq);
            h_mmiss_thetapq_all->Fill(thetapq,mmiss);
            h_mmiss_pq_all->Fill(pq,mmiss);
              if (thetapq>25) {continue;}
              if (pq<0.62 || pq>0.96) {continue;}


              h_pmiss_fd->Fill(pmiss.Mag());
            h_pmiss_all->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}


              h_mmiss_fd->Fill(mmiss);
            h_mmiss_all->Fill(mmiss);
              if (mmiss>1.1) {continue;}
              h_p_theta_fd->Fill(theta,mom);
            h_p_theta_all->Fill(theta,mom);
            h_pmiss_src_all->Fill(pmiss.Mag());
            h_xb_q2_src_all->Fill(Q2,xB);


              // recoil selection here


	    }
            // CENTRAL DETECTOR PROTONS
	    else if((*p)->getRegion() == CD){

	      if(!CD_fiducial(phi,theta,momT)){
		continue;
	      }

              h_thetae_q2_cd->Fill(Q2,el.Theta()*180./M_PI);
              h_thetae_q2_all->Fill(Q2,el.Theta()*180./M_PI);

              // SRC histograms here - CD
              h_mmiss_nocuts_cd->Fill(mmiss);
              h_mmiss_xb_cd->Fill(xB,mmiss);  // define
              h_mmiss_nocuts_all->Fill(mmiss);
              h_mmiss_xb_all->Fill(xB,mmiss);
              h_xb_cd->Fill(xB);
              h_xb_all->Fill(xB);
              if (xB<1.2) {continue;}


              h_mmiss_q2_cd->Fill(Q2,mmiss);  // define
            h_mmiss_q2_all->Fill(Q2,mmiss);

              h_q2_cd->Fill(Q2);
            h_q2_all->Fill(Q2);
              if (Q2<1.5) {continue;}

              h_mmiss_thetapq_cd->Fill(thetapq,mmiss);
              h_mmiss_pq_cd->Fill(pq,mmiss);
              h_thetapq_pq_cd->Fill(pq,thetapq);
            h_thetapq_pq_all->Fill(pq,thetapq);
            h_mmiss_thetapq_all->Fill(thetapq,mmiss);
            h_mmiss_pq_all->Fill(pq,mmiss);
              if (thetapq>25) {continue;}
              if (pq<0.62 || pq>0.96) {continue;}


              h_pmiss_cd->Fill(pmiss.Mag());
            h_pmiss_all->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}


              h_mmiss_cd->Fill(mmiss);
            h_mmiss_all->Fill(mmiss);

              if (mmiss>1.1) {continue;}
              h_p_theta_cd->Fill(theta,mom);
            h_p_theta_all->Fill(theta,mom);
            h_pmiss_src_all->Fill(pmiss.Mag());
            h_xb_q2_src_all->Fill(Q2,xB);

	    }

	    else{
	      cout<<"Not Either"<<endl;
	    }	   

            // ALL PROTONS - FD and CD
            if ((*p)->getRegion()==CD && !CD_fiducial(phi,theta,momT)) {continue;} // CD fiducial cut



            if (xB<1.2) {continue;}



            if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}



            if (Q2<1.5) {continue;}

            if (thetapq>25) {continue;}
            if (pq<0.62 || pq>0.96) {continue;}
            if (mmiss>1.1) {continue;}
 
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
  myCanvas->cd(1)->SetLogy();
  h_xB_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  //h_phi_theta_bc->Draw("colz");
  h_thetap_q2_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_qmag_qtheta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  



  ///////////////////////////////////////////////////////
  //Both Forward Detector & Central Detector
  ///////////////////////////////////////////////////////  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_nocuts_all->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetae_q2_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_xb_all->Draw("colz");
  TLine * line_mmiss_xb_all = new TLine(1.2,0,1.2,2);
  line_mmiss_xb_all->SetLineColor(kRed);
  line_mmiss_xb_all->SetLineWidth(3);
  line_mmiss_xb_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogy();
  h_xb_all->Draw();
  TLine * line_xb_all = new TLine(1.2,0,1.2,h_xb_all->GetMaximum());
  line_xb_all->SetLineColor(kRed);
  line_xb_all->SetLineWidth(3);
  line_xb_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_q2_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_all->Draw();
  TLine * line_q2_all = new TLine(1.5,0,1.5,h_q2_all->GetMaximum());
  line_q2_all->SetLineColor(kRed);
  line_q2_all->SetLineWidth(3);
  line_q2_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_thetapq_all->Draw("colz");
  TLine * line_mmiss_thetapq_all = new TLine(25,0,25,2);
  line_mmiss_thetapq_all->SetLineColor(kRed);
  line_mmiss_thetapq_all->SetLineWidth(3);
  line_mmiss_thetapq_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pq_all->Draw("colz");
  TLine * line_mmiss_pq_all_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_all_1->SetLineColor(kRed);
  line_mmiss_pq_all_1->SetLineWidth(3);
  line_mmiss_pq_all_1->Draw("same");
  TLine * line_mmiss_pq_all_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_all_2->SetLineColor(kRed);
  line_mmiss_pq_all_2->SetLineWidth(3);
  line_mmiss_pq_all_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_all->Draw("colz");
  TLine * line_pq_all_1 = new TLine(0.62,0,0.62,60);
  line_pq_all_1->SetLineColor(kRed);
  line_pq_all_1->SetLineWidth(3);
  line_pq_all_1->Draw("same");
  TLine * line_pq_all_2 = new TLine(0.96,0,0.96,60);
  line_pq_all_2->SetLineColor(kRed);
  line_pq_all_2->SetLineWidth(3);
  line_pq_all_2->Draw("same");
  TLine * line_pq_all_3 = new TLine(0,25,1.2,25);
  line_pq_all_3->SetLineColor(kRed);
  line_pq_all_3->SetLineWidth(3);
  line_pq_all_3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_all->Draw();
  TLine * line_pmiss_all_1 = new TLine(0.35,0,0.35,h_pmiss_all->GetMaximum());
  TLine * line_pmiss_all_2 = new TLine(1,0,1,h_pmiss_all->GetMaximum());
  line_pmiss_all_1->SetLineColor(kRed);
  line_pmiss_all_1->SetLineWidth(3);
  line_pmiss_all_1->Draw("same");
  line_pmiss_all_2->SetLineColor(kRed);
  line_pmiss_all_2->SetLineWidth(3);
  line_pmiss_all_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_all->Draw();
  TLine * line_mmiss_all = new TLine(1.1,0,1.1,h_mmiss_all->GetMaximum());
  line_mmiss_all->SetLineColor(kRed);
  line_mmiss_all->SetLineWidth(3);
  line_mmiss_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_theta_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_src_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xb_q2_src_all->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
 

  ///////////////////////////////////////////////////////
  //Forward Detector
  ///////////////////////////////////////////////////////  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_nocuts_fd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetae_q2_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_xb_fd->Draw("colz");
  TLine * line_mmiss_xb_fd = new TLine(1.2,0,1.2,2);
  line_mmiss_xb_fd->SetLineColor(kRed);
  line_mmiss_xb_fd->SetLineWidth(3);
  line_mmiss_xb_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogy();
  h_xb_fd->Draw();
  TLine * line_xb_fd = new TLine(1.2,0,1.2,h_xb_fd->GetMaximum());
  line_xb_fd->SetLineColor(kRed);
  line_xb_fd->SetLineWidth(3);
  line_xb_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_q2_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_fd->Draw();
  TLine * line_q2_fd = new TLine(1.5,0,1.5,h_q2_fd->GetMaximum());
  line_q2_fd->SetLineColor(kRed);
  line_q2_fd->SetLineWidth(3);
  line_q2_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_thetapq_fd->Draw("colz");
  TLine * line_mmiss_thetapq_fd = new TLine(25,0,25,2);
  line_mmiss_thetapq_fd->SetLineColor(kRed);
  line_mmiss_thetapq_fd->SetLineWidth(3);
  line_mmiss_thetapq_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pq_fd->Draw("colz");
  TLine * line_mmiss_pq_fd_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_fd_1->SetLineColor(kRed);
  line_mmiss_pq_fd_1->SetLineWidth(3);
  line_mmiss_pq_fd_1->Draw("same");
  TLine * line_mmiss_pq_fd_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_fd_2->SetLineColor(kRed);
  line_mmiss_pq_fd_2->SetLineWidth(3);
  line_mmiss_pq_fd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_fd->Draw("colz");
  TLine * line_pq_fd_1 = new TLine(0.62,0,0.62,60);
  line_pq_fd_1->SetLineColor(kRed);
  line_pq_fd_1->SetLineWidth(3);
  line_pq_fd_1->Draw("same");
  TLine * line_pq_fd_2 = new TLine(0.96,0,0.96,60);
  line_pq_fd_2->SetLineColor(kRed);
  line_pq_fd_2->SetLineWidth(3);
  line_pq_fd_2->Draw("same");
  TLine * line_pq_fd_3 = new TLine(0,25,1.2,25);
  line_pq_fd_3->SetLineColor(kRed);
  line_pq_fd_3->SetLineWidth(3);
  line_pq_fd_3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_fd->Draw();
  TLine * line_pmiss_fd_1 = new TLine(0.35,0,0.35,h_pmiss_fd->GetMaximum());
  TLine * line_pmiss_fd_2 = new TLine(1,0,1,h_pmiss_fd->GetMaximum());
  line_pmiss_fd_1->SetLineColor(kRed);
  line_pmiss_fd_1->SetLineWidth(3);
  line_pmiss_fd_1->Draw("same");
  line_pmiss_fd_2->SetLineColor(kRed);
  line_pmiss_fd_2->SetLineWidth(3);
  line_pmiss_fd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_fd->Draw();
  TLine * line_mmiss_fd = new TLine(1.1,0,1.1,h_mmiss_fd->GetMaximum());
  line_mmiss_fd->SetLineColor(kRed);
  line_mmiss_fd->SetLineWidth(3);
  line_mmiss_fd->Draw("same");
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
  h_mmiss_nocuts_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetae_q2_cd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_xb_cd->Draw("colz");
  TLine * line_mmiss_xb_cd = new TLine(1.2,0,1.2,2);
  line_mmiss_xb_cd->SetLineColor(kRed);
  line_mmiss_xb_cd->SetLineWidth(3);
  line_mmiss_xb_cd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogy();
  h_xb_cd->Draw();
  TLine * line_xb_cd = new TLine(1.2,0,1.2,h_xb_cd->GetMaximum());
  line_xb_cd->SetLineColor(kRed);
  line_xb_cd->SetLineWidth(3);
  line_xb_cd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_cd->Draw();
  TLine * line_q2_cd = new TLine(1.5,0,1.5,h_q2_cd->GetMaximum());
  line_q2_cd->SetLineColor(kRed);
  line_q2_cd->SetLineWidth(3);
  line_q2_cd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_thetapq_cd->Draw("colz");
  TLine * line_mmiss_thetapq_cd = new TLine(25,0,25,2);
  line_mmiss_thetapq_cd->SetLineColor(kRed);
  line_mmiss_thetapq_cd->SetLineWidth(3);
  line_mmiss_thetapq_cd->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pq_cd->Draw("colz");
  TLine * line_mmiss_pq_cd_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_cd_1->SetLineColor(kRed);
  line_mmiss_pq_cd_1->SetLineWidth(3);
  line_mmiss_pq_cd_1->Draw("same");
  TLine * line_mmiss_pq_cd_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_cd_2->SetLineColor(kRed);
  line_mmiss_pq_cd_2->SetLineWidth(3);
  line_mmiss_pq_cd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_cd->Draw("colz");
  TLine * line_pq_cd_1 = new TLine(0.62,0,0.62,60);
  line_pq_cd_1->SetLineColor(kRed);
  line_pq_cd_1->SetLineWidth(3);
  line_pq_cd_1->Draw("same");
  TLine * line_pq_cd_2 = new TLine(0.96,0,0.96,60);
  line_pq_cd_2->SetLineColor(kRed);
  line_pq_cd_2->SetLineWidth(3);
  line_pq_cd_2->Draw("same");
  TLine * line_pq_cd_3 = new TLine(0,25,1.2,25);
  line_pq_cd_3->SetLineColor(kRed);
  line_pq_cd_3->SetLineWidth(3);
  line_pq_cd_3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_cd->Draw();
  TLine * line_pmiss_cd_1 = new TLine(0.35,0,0.35,h_pmiss_cd->GetMaximum());
  TLine * line_pmiss_cd_2 = new TLine(1,0,1,h_pmiss_cd->GetMaximum());
  line_pmiss_cd_1->SetLineColor(kRed);
  line_pmiss_cd_1->SetLineWidth(3);
  line_pmiss_cd_1->Draw("same");
  line_pmiss_cd_2->SetLineColor(kRed);
  line_pmiss_cd_2->SetLineWidth(3);
  line_pmiss_cd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_cd->Draw();
  TLine * line_mmiss_cd = new TLine(1.1,0,1.1,h_mmiss_cd->GetMaximum());
  line_mmiss_cd->SetLineColor(kRed);
  line_mmiss_cd->SetLineWidth(3);
  line_mmiss_cd->Draw("same");
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

