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
  std::cerr << "Usage: ./code A outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }


  int A = atoi(argv[1]);


  // set up output root file (and tree)
  TString outFile = argv[2];
  TFile *f = new TFile(outFile,"RECREATE");
  TTree *tree = new TTree("T","myTree");
  double e_mom, e_theta, e_phi, mom, theta, phi, xB, Q2; 
  double thetapq, pq, mmiss, pmiss_mom, pmiss_theta, pmiss_phi, theta_pmq;
  tree->Branch("e_mom",&e_mom,"e_mom/D");
  tree->Branch("e_theta",&e_theta,"e_theta/D");
  tree->Branch("e_phi",&e_phi,"e_phi/D");
  tree->Branch("xB",&xB,"xB/D");
  tree->Branch("Q2",&Q2,"Q2/D");
  tree->Branch("mom",&mom,"mom/D");
  tree->Branch("theta",&theta,"theta/D");
  tree->Branch("phi",&phi,"phi/D");
  tree->Branch("thetapq",&thetapq,"thetapq/D");
  tree->Branch("pq",&pq,"pq/D");
  tree->Branch("mmiss",&mmiss,"mmiss/D");
  tree->Branch("pmiss_mom",&pmiss_mom,"pmiss_mom/D");
  tree->Branch("pmiss_theta",&pmiss_theta,"pmiss_theta/D");
  tree->Branch("pmiss_phi",&pmiss_phi,"pmiss_phi/D");
  tree->Branch("theta_pmq",&theta_pmq,"theta_pmq/D");


  // set up output pdf file
  char * pdfFile = argv[3];

  cout<<"Output file "<< outFile <<endl;
  cout<<"Output PDF file "<< pdfFile <<endl;



  // ELECTRON AND PROTON CUTS
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
  double mU = 0.9315;
  double me = 0.000510998;
  double mA = 0;
  switch (A) {
    case 4: // He-4
      mA = 4.002603*mU - 2*me;
    case 12: // C-12
      mA = 12*mU - 6*me;
    case 40: // Ca-40
      mA = 39.962591*mU - 20*me;
    case 48: // Ca-48
      mA = 47.952523*mU - 20*me;
  }

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
  TH1D * h_mmiss_nocuts_fd = new TH1D("mmiss_nocuts_fd","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2})",100,0,2);
  TH1D * h_xb_fd = new TH1D("xb_fd","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_pmiss_fd = new TH1D("pmiss_fd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_fd = new TH1D("q2_fd","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_fd = new TH2D("mmiss_thetapq_fd","Missing Mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_fd = new TH2D("mmiss_pq_fd","Missing Mass vs p_{L}/q;p_{L}/q;M_{miss} (GeV/c^{2})",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_fd = new TH2D("thetapq_pq_fd","#theta_{pq} vs p_{L}/q;p_{L}/q;#theta_{pq} (degrees)",50,0,1.2,50,0,60);
  TH1D * h_mmiss_fd = new TH1D("mmiss_fd","Missing Mass;M_{miss} (GeV/c^{2});Counts",50,0,2);
  TH2D * h_p_theta_fd = new TH2D("p_theta_fd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p_{L} (GeV/c)",45,0,180,50,0,2.5);






  ///////////////////////////////////////////////////////
  //Central Detector
  ///////////////////////////////////////////////////////  

  TH2D * h_thetae_q2_cd = new TH2D("thetae_q2_cd","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_mmiss_xb_cd = new TH2D("mmiss_xb_cd","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_cd = new TH2D("mmiss_q2_cd","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_cd = new TH1D("mmiss_nocuts_cd","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2})",100,0,2);
  TH1D * h_xb_cd = new TH1D("xb_cd","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH1D * h_pmiss_cd = new TH1D("pmiss_cd","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_cd = new TH1D("q2_cd","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_cd = new TH2D("mmiss_thetapq_cd","Missing mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_cd = new TH2D("mmiss_pq_cd","Mmiss vs p_{L}/q;p_{L}/q;Mmiss",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_cd = new TH2D("thetapq_pq_cd","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",50,0,1.2,50,0,60);
  TH1D * h_mmiss_cd = new TH1D("mmiss_cd","Missing Mass;M_{miss} (GeV/c^{2});Counts",50,0,2);
  TH2D * h_p_theta_cd = new TH2D("p_theta_cd","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p_{L} (GeV/c)",45,0,180,50,0,2.5);




  ///////////////////////////////////////////////////////
  //Both Forward Detector & Central Detector
  ///////////////////////////////////////////////////////  

  TH2D * h_thetae_q2_all = new TH2D("thetae_q2_all","Electron Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{e} (deg)",100,0,4,60,0,60);
  TH2D * h_thetap_q2_all = new TH2D("thetap_q2_all","Proton Angle vs Q^{2};Q^{2} (GeV^{2});#theta_{p} (deg)",100,0,4,60,0,60);
  TH1D * h_xb_all = new TH1D("xb_all","x-Bjorken x_{B};x_{B};Counts",100,0.5,2.5);
  TH2D * h_mmiss_xb_all = new TH2D("mmiss_xb_all","Missing Mass vs x_{B};x_{B};M_{miss} (GeV/c^{2})",50,0,3,50,0,2);
  TH2D * h_mmiss_q2_all = new TH2D("mmiss_q2_all","Missing Mass vs Q^{2};Q^{2} (GeV^{2});M_{miss} (GeV/c^{2})",50,0,4,50,0,2);
  TH1D * h_mmiss_nocuts_all = new TH1D("mmiss_nocuts_all","Missing Mass (before SRC cuts);Missing Mass (GeV/c^{2})",100,0,2);
  TH1D * h_pmiss_all = new TH1D("pmiss_all","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",25,0,1.2);
  TH1D * h_q2_all = new TH1D("q2_all","Q^{2};Q^{2} (GeV^{2});Counts",25,0,4);
  TH2D * h_mmiss_thetapq_all = new TH2D("mmiss_thetapq_all","Missing Mass vs #theta_{pq};#theta_{pq} (deg);M_{miss} (GeV/c^{2})",60,0,60,60,0,2);
  TH2D * h_mmiss_pq_all = new TH2D("mmiss_pq_all","Missing Mass vs p_{L}/q;p_{L}/q;M_{miss} (GeV/c^{2})",60,0,1.2,60,0,2);
  TH2D * h_thetapq_pq_all = new TH2D("thetapq_pq_all","#theta_{pq} vs p_{L}/q;p+{L}/q;#theta_{pq} (degrees)",50,0,1.2,50,0,60);
  TH1D * h_mmiss_all = new TH1D("mmiss_all","Missing Mass;M_{miss} (GeV/c^{2});Counts",50,0,2);

  TH2D * h_p_theta_all = new TH2D("p_theta_all","Phase Space Distribution of Leading Protons;#theta_{p} (Degrees);Momentum p_{L} (GeV/c)",45,0,180,50,0,2.5);
  TH2D * h_emiss_omega_all = new TH2D("emiss_omega_all","Missing Energy vs #omega;#omega (GeV);E_{miss} (GeV)",100,0,2.5,100,0,0.5);
  
  TH1D * h_pmiss_src_all = new TH1D("pmiss_src_all","Missing Momentum of Lead SRC Protons;p_{miss} (GeV/c);Counts",50,0.35,1);
  TH2D * h_xb_q2_src_all = new TH2D("xb_q2_src_all","x_{B} vs Q^{2} of Lead SRC Protons;Q^{2} (GeV^{2});x_{B}",50,1.5,4,50,1.2,2.5);
  /*TH2D * h_q2_omega_src_all = new TH2D("emiss_omega","Missing Energy vs #omega;#omega (GeV);E_{miss}"
  TH2D * h_mmiss_emiss_src_all*/

  ///////////////////////////////////////////////////////
  //Final Distributions after Cuts
  ///////////////////////////////////////////////////////  
  TH1D * h_xb_final = new TH1D("xb_final","x-Bjorken x_{B};x_{B};Counts",100,1.2,2.5);
  TH1D * h_pmiss_final = new TH1D("pmiss_final","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0.35,1.2);
  TH1D * h_q2_final = new TH1D("q2_final","Q^{2};Q^{2} (GeV^{2});Counts",100,1.4,4);
  TH2D * h_thetapq_pq_final = new TH2D("thetapq_pq_final","#theta_{pq} vs p/q;p/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_final = new TH1D("mmiss_final","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);

  ///////////////////////////////////////////////////////
  //(e,e'p) Kinematics
  ///////////////////////////////////////////////////////  
  TH1D * h_emiss = new TH1D("emiss","Missing Energy;E_{miss} (GeV);Counts",50,-0.2,0.5);
  TH2D * h_emiss_omega = new TH2D("emiss_omega","Missing Energy vs #omega;#omega (GeV);E_{miss} (GeV)",50,0,2.5,50,-0.2,0.6);
  TH2D * h_q2_omega = new TH2D("q2_omega","Q^{2} vs #omega;#omega (GeV);Q^{2} (GeV/c^{2})",50,0,2.5,50,1,6);
  TH2D * h_mmiss_emiss = new TH2D("mmiss_emiss","Missing Mass vs Missing Energy;E_{miss} (GeV);M_{miss}",50,-0.3,0.6,50,0.4,1.3);
  


  ///////////////////////////////////////////////////////
  //Recoil Selection
  /////////////////////////////////////////////////////// 
  TH2D * h_lead_rec = new TH2D("lead_rec","Lead Theta vs Recoil Theta;#theta_{rec} (deg);#theta_{L} (deg)",90,0,180,90,0,180); 
  TH1D * h_prec = new TH1D("prec","Recoil Proton Momentum;p_{rec} (GeV/c);Counts",100,0,1.5);
  TH1D * h_cos0 = new TH1D("cos0","Angle between Lead and Recoil Protons",100,-1,1);
  

  ///////////////////////////////////////////////////////
  //(e,e'pp) Kinematics
  /////////////////////////////////////////////////////// 
  TH1D * h_emiss_pp = new TH1D("emiss_pp","Missing Energy;E_{miss} (GeV)",100,-0.1,0.5);
  TH2D * h_emiss_omega_pp = new TH2D("emiss_omega_pp","Missing Energy vs #omega;#omega (GeV);E_{miss} (GeV)",50,0,2.5,50,-0.1,0.5);
  TH1D * h_xb_pp = new TH1D("xb_pp","x-Bjorken x_{B};x_{B};Counts",100,1.2,2.5);
  TH1D * h_pmiss_pp = new TH1D("pmiss_pp","Missing Momentum p_{miss};p_{miss} (GeV/c);Counts",120,0.3,1.2);
  TH1D * h_q2_pp = new TH1D("q2_pp","Q^{2};Q^{2} (GeV^{2});Counts",100,1,4);
  TH2D * h_thetapq_pq_pp = new TH2D("thetapq_pq_pp","#theta_{pq} vs p_{L}/q;p_{L}/q;#theta_{pq} (degrees)",100,0,1.2,100,0,60);
  TH1D * h_mmiss_pp = new TH1D("mmiss_pp","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,1.2);
  TH2D * h_plead_prec = new TH2D("plead_prec","Lead Proton Momentum vs Recoil Momentum;p_{rec} (GeV/c);p_{L} (GeV/c)",50,0,1,50,0,3);
  TH2D * h_tlead_trec = new TH2D("tlead_trec","Lead Proton #theta vs Recoil #theta;#theta_{rec} (deg);#theta_{L} (deg)",90,0,180,90,0,180);
  TH2D * h_pmiss_prec = new TH2D("pmiss_prec","Missing Momentum vs Recoil Momentum;p_{rec} (GeV/c);p_{miss} (GeV/c)",50,0.3,1,50,0.3,1.2);
  TH1D * h_thetapmq_pp = new TH1D("thetapmq_pp","Angle between Missing Momentum and q;#theta_{pmiss,q};Counts",50,90,180);


  ///////////////////////////////////////////////////////
  //Additional Histos
  /////////////////////////////////////////////////////// 
  TH1D * h_plead_all = new TH1D("plead_all","Momentum of Lead Proton Candidate;p_{L} (GeV/c)",100,0,3);
  TH1D * h_plead_fd = new TH1D("plead_fd","Momentum of Lead Proton Candidate;p_{L} (GeV/c)",100,0,3);
  TH1D * h_plead_cd = new TH1D("plead_cd","Momentum of Lead Proton Candidate;p_{L} (GeV/c)",100,0,3);
  TH2D * h_thetapmq_pq_all = new TH2D("thetapmq_pq_all","Angle between Missing Momentum and q;p/q;#theta_{pmiss,q}",100,0,1.2,90,0,180);
  TH2D * h_thetapmq_pq_fd = new TH2D("thetapmq_pq_fd","Angle between Missing Momentum and q;p/q;#theta_{pmiss,q}",100,0,1.2,90,0,180);
  TH2D * h_thetapmq_pq_cd = new TH2D("thetapmq_pq_cd","Angle between Missing Momentum and q;p/q;#theta_{pmiss,q}",100,0,1.2,90,0,180);
  TH1D * h_doublelead = new TH1D("doublelead","Number of Protons Passing SRC Cuts;Proton Number",5,1,6);


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
	  e_mom = el.P();
          e_theta = el.Theta()*180./M_PI;
          e_phi = el.Phi()*180./M_PI;

	  TLorentzVector q = beam - el; //photon  4-vector            
          Q2        = -q.M2(); // Q^2
          xB       = Q2/(2 * mass_p * (beam.E() - el.E()) ); //x-bjorken
          double omega = beam.E() - el.E();
	  h_Q2_bc->Fill(Q2);
	  h_xB_bc->Fill(xB);
	  h_phi_theta_bc->Fill(el.Phi()*180/M_PI,el.Theta()*180/M_PI);
          h_qmag_qtheta->Fill(-1*q.Mag(),q.Vect().Theta()*180./M_PI);
	  double vtz_e = electrons[0]->par()->getVz();

          // define final lead kinematics - use these after first proton loop (see warning below)
          TVector3 lead_pmiss(0,0,0);
          double lead_mom = 0; double lead_theta = 0;
          double lead_thetapq = 0; double lead_pq = 0;
          double lead_mmiss = 0;
          double Tp = 0; double Tb = 0; double lead_emiss = 0;
	  
	  ///////////////////////////////
	  //Before cuts
	  ///////////////////////////////
	  int pindex = -1;
          int num_lead = 0;

	  for(auto p = protons.begin(); p != protons.end();++p){
	    if((*p)->par()->getCharge()<1){continue;}

            // Warning: these proton kinematics apply to the current proton in the loop. If the lead proton is not the last proton in the proton array, these values will not be the correct values for the lead.

	    //Momenta
	    SetLorentzVector(lead_ptr,(*p));
	    mom = lead_ptr.P();
	    double momT = lead_ptr.Perp();
	    theta = lead_ptr.Theta() * 180 / M_PI;
	    phi = lead_ptr.Phi() * 180 / M_PI;

	    double beta = (*p)->par()->getBeta();
	    double path = (*p)->getPath();
	    double vtz_p = (*p)->par()->getVz();

            // calculate SRC kinematics
            TVector3 pmiss = lead_ptr.Vect() - q.Vect();
            thetapq = lead_ptr.Vect().Angle(q.Vect())*180./M_PI;
            pq = (lead_ptr.Vect().Mag()) / (q.Vect().Mag());
            mmiss = (q + TLorentzVector(TVector3(0.,0.,0.),2*mN) - lead_ptr).Mag();
            pmiss_mom = pmiss.Mag();
            pmiss_theta = pmiss.Theta()*180./M_PI;
            pmiss_phi = pmiss.Phi()*180./M_PI;
            theta_pmq = pmiss.Angle(q.Vect())*180./M_PI;

            // end warning

	    if(beta<0.2){continue;} // proton cut

            tree->Fill();


            // FORWARD DETECTOR PROTONS
	    if((*p)->getRegion() == FD){

              h_thetae_q2_fd->Fill(Q2,el.Theta()*180./M_PI);
              h_thetae_q2_all->Fill(Q2,el.Theta()*180./M_PI);
              h_thetap_q2_all->Fill(Q2,theta);

              // plead cut
              h_plead_all->Fill(mom);
              h_plead_fd->Fill(mom);
              if (mom<1.0) {continue;}

              // SRC histograms here - FD
              h_mmiss_nocuts_fd->Fill(mmiss);
              h_mmiss_xb_fd->Fill(xB,mmiss);
              h_mmiss_nocuts_all->Fill(mmiss);
              h_mmiss_xb_all->Fill(xB,mmiss);
              h_xb_fd->Fill(xB);
            h_xb_all->Fill(xB);
              if (xB<1.2) {continue;}


              h_pmiss_fd->Fill(pmiss.Mag());
            h_pmiss_all->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}

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
              h_thetapmq_pq_fd->Fill(pq,theta_pmq);
            h_thetapmq_pq_all->Fill(pq,theta_pmq);

              //if (thetapq>25) {continue;}
              //if (pq<0.62) {continue;}
              if (pq>0.96) {continue;}


              h_mmiss_fd->Fill(mmiss);
            h_mmiss_all->Fill(mmiss);
              if (mmiss>1.1) {continue;}
              h_p_theta_fd->Fill(theta,mom);
            h_p_theta_all->Fill(theta,mom);
            h_pmiss_src_all->Fill(pmiss.Mag());
            h_xb_q2_src_all->Fill(Q2,xB);


              // define kinematics for lead proton (valid past this proton loop)
              pindex = (*p)->trk(DC)->getPindex();
              lead_mom = mom;
              lead_theta = theta;
              lead_pmiss = pmiss;
              lead_thetapq = thetapq;
              lead_pq = pq;
              lead_mmiss = mmiss;
              Tp = lead_ptr.E() - mN;
              Tb = omega + mA - lead_ptr.E() - pow( pow( (omega + mA - lead_ptr.E()), 2.) - lead_pmiss*lead_pmiss, 0.5 );
              lead_emiss = omega - Tp - Tb;


	    }
            // CENTRAL DETECTOR PROTONS
	    else if((*p)->getRegion() == CD){

	      if(!CD_fiducial(phi,theta,momT)){
		continue;
	      }

              h_thetae_q2_cd->Fill(Q2,el.Theta()*180./M_PI);
              h_thetae_q2_all->Fill(Q2,el.Theta()*180./M_PI);

              // plead cut
              h_plead_all->Fill(mom);
              h_plead_cd->Fill(mom);
              if (mom<1.0) {continue;}

              // SRC histograms here - CD
              h_mmiss_nocuts_cd->Fill(mmiss);
              h_mmiss_xb_cd->Fill(xB,mmiss);  // define
              h_mmiss_nocuts_all->Fill(mmiss);
              h_mmiss_xb_all->Fill(xB,mmiss);
              h_xb_cd->Fill(xB);
              h_xb_all->Fill(xB);
              if (xB<1.2) {continue;}


              h_pmiss_cd->Fill(pmiss.Mag());
            h_pmiss_all->Fill(pmiss.Mag());
              if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}


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
              h_thetapmq_pq_cd->Fill(pq,theta_pmq);
            h_thetapmq_pq_all->Fill(pq,theta_pmq);

              //if (thetapq>25) {continue;}
              //if (pq<0.62) {continue;}
              if (pq>0.96) {continue;}

              h_mmiss_cd->Fill(mmiss);
            h_mmiss_all->Fill(mmiss);

              if (mmiss>1.1) {continue;}
              h_p_theta_cd->Fill(theta,mom);
            h_p_theta_all->Fill(theta,mom);
            h_pmiss_src_all->Fill(pmiss.Mag());
            h_xb_q2_src_all->Fill(Q2,xB);


              // define kinematics for lead proton (valid past this proton loop)
              pindex = (*p)->trk(CVT)->getPindex();
              lead_mom = mom;
              lead_theta = theta;
              lead_pmiss = pmiss;
              lead_thetapq = thetapq;
              lead_pq = pq;
              lead_mmiss = mmiss;
              Tp = lead_ptr.E() - mN;
              Tb = omega + mA - lead_ptr.E() - pow( pow( (omega + mA - lead_ptr.E()), 2.) - lead_pmiss*lead_pmiss, 0.5 );
              lead_emiss = omega - Tp - Tb;

	    }

	    else{
	      cout<<"Not Either"<<endl;
	    }  // end if statement deciding whether lead proton is in FD or CD




            // ALL PROTONS - FD and CD
            if ((*p)->getRegion()==CD && !CD_fiducial(phi,theta,momT)) {continue;} // CD fiducial cut if applicable

            // SRC cuts
            if (mom<1.0) {continue;}
            if (xB<1.2) {continue;}
            if (pmiss.Mag()<0.35 || pmiss.Mag()>1.0) {continue;}
            if (Q2<1.5) {continue;}
            //if (thetapq>25) {continue;}
            //if (pq<0.62) {continue;}
            if (pq>0.96) {continue;}
            if (mmiss>1.1) {continue;}


            // count protons passing SRC cuts
            num_lead = num_lead + 1;


            // final e'p distributions after SRC cuts
            h_xb_final->Fill(xB);
            h_pmiss_final->Fill(pmiss.Mag());
            h_q2_final->Fill(Q2);
            h_thetapq_pq_final->Fill(pq,thetapq);
            h_mmiss_final->Fill(mmiss);
            h_emiss->Fill(lead_emiss);
            h_emiss_omega->Fill(omega,lead_emiss);
            h_q2_omega->Fill(omega,Q2);
            h_mmiss_emiss->Fill(lead_emiss,mmiss);
 
	  } // end first loop over protons

          if (num_lead>0) {h_doublelead->Fill(num_lead);}

          if (pindex==-1) {continue;} // if no lead proton, skip to next event (don't look for recoil)


          // e'pp starts here - look for CD recoils only


	  for(auto p = protons.begin(); p != protons.end();++p){

	    if((*p)->par()->getCharge()<1){continue;} // look only at charged particles in the CVT

            if ((*p)->trk(CVT)->getPindex()==pindex) {continue;} // make sure the recoil is not the same as the lead


              //Momenta
              SetLorentzVector(recoil_ptr,(*p));
              double mom = recoil_ptr.P();
              double momT = recoil_ptr.Perp();
              double theta = recoil_ptr.Theta() * 180 / M_PI;
              double phi = recoil_ptr.Phi() * 180 / M_PI;

              double beta = (*p)->par()->getBeta();
              double path = (*p)->getPath();
              double vtz_p = (*p)->par()->getVz();


              // calculate SRC kinematics
              TVector3 pmiss = recoil_ptr.Vect() - q.Vect();
              double thetapq = recoil_ptr.Vect().Angle(q.Vect())*180./M_PI;
              double pq = (recoil_ptr.Vect().Mag()) / (q.Vect().Mag());
              mmiss = (q + TLorentzVector(TVector3(0.,0.,0.),2*mN) - recoil_ptr).Mag();
              double recoil_emiss = lead_emiss - recoil_ptr.E();

	      //if(mom==0){continue;} // proton cut
	      h_lead_rec->Fill(theta,lead_theta);

            if((*p)->getRegion() == CD){ // look in CD only


              // cut on momentum
              h_prec->Fill(mom);
              if (mom<0.35) {continue;}
              
              h_cos0->Fill(cos(lead_pmiss.Angle(recoil_ptr.Vect())));

              // histos for e'pp kinematics
              h_emiss_pp->Fill(lead_emiss);
              h_emiss_omega_pp->Fill(omega,lead_emiss);
              h_xb_pp->Fill(xB);
              h_pmiss_pp->Fill(lead_pmiss.Mag());
              h_q2_pp->Fill(Q2);
              h_thetapq_pq_pp->Fill(lead_pq,lead_thetapq);
              h_mmiss_pp->Fill(lead_mmiss);
              h_plead_prec->Fill(mom,lead_mom);
              h_tlead_trec->Fill(theta,lead_theta);
              h_pmiss_prec->Fill(mom,lead_pmiss.Mag());
              h_thetapmq_pp->Fill(lead_pmiss.Angle(q.Vect())*180./M_PI);
              

            } // end requirement for recoil to be in CD
          

          } // end second loop over protons


	} // end requirement for one electron in event
    } // end event loop

  //clasAna.WriteDebugPlots();



  f->cd();
  tree->Write();

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
  /*TLine * line_mmiss_pq_all_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_all_1->SetLineColor(kRed);
  line_mmiss_pq_all_1->SetLineWidth(3);
  line_mmiss_pq_all_1->Draw("same");*/
  TLine * line_mmiss_pq_all_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_all_2->SetLineColor(kRed);
  line_mmiss_pq_all_2->SetLineWidth(3);
  line_mmiss_pq_all_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_all->Draw("colz");
  /*TLine * line_pq_all_1 = new TLine(0.62,0,0.62,60);
  line_pq_all_1->SetLineColor(kRed);
  line_pq_all_1->SetLineWidth(3);
  line_pq_all_1->Draw("same");*/
  TLine * line_pq_all_2 = new TLine(0.96,0,0.96,60);
  line_pq_all_2->SetLineColor(kRed);
  line_pq_all_2->SetLineWidth(3);
  line_pq_all_2->Draw("same");
  /*TLine * line_pq_all_3 = new TLine(0,25,1.2,25);
  line_pq_all_3->SetLineColor(kRed);
  line_pq_all_3->SetLineWidth(3);
  line_pq_all_3->Draw("same");*/
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
  /*TLine * line_mmiss_pq_fd_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_fd_1->SetLineColor(kRed);
  line_mmiss_pq_fd_1->SetLineWidth(3);
  line_mmiss_pq_fd_1->Draw("same");*/
  TLine * line_mmiss_pq_fd_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_fd_2->SetLineColor(kRed);
  line_mmiss_pq_fd_2->SetLineWidth(3);
  line_mmiss_pq_fd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_fd->Draw("colz");
  /*TLine * line_pq_fd_1 = new TLine(0.62,0,0.62,60);
  line_pq_fd_1->SetLineColor(kRed);
  line_pq_fd_1->SetLineWidth(3);
  line_pq_fd_1->Draw("same");*/
  TLine * line_pq_fd_2 = new TLine(0.96,0,0.96,60);
  line_pq_fd_2->SetLineColor(kRed);
  line_pq_fd_2->SetLineWidth(3);
  line_pq_fd_2->Draw("same");
  /*TLine * line_pq_fd_3 = new TLine(0,25,1.2,25);
  line_pq_fd_3->SetLineColor(kRed);
  line_pq_fd_3->SetLineWidth(3);
  line_pq_fd_3->Draw("same");*/
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
  /*TLine * line_mmiss_pq_cd_1 = new TLine(0.62,0,0.62,2);
  line_mmiss_pq_cd_1->SetLineColor(kRed);
  line_mmiss_pq_cd_1->SetLineWidth(3);
  line_mmiss_pq_cd_1->Draw("same");*/
  TLine * line_mmiss_pq_cd_2 = new TLine(0.96,0,0.96,2);
  line_mmiss_pq_cd_2->SetLineColor(kRed);
  line_mmiss_pq_cd_2->SetLineWidth(3);
  line_mmiss_pq_cd_2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_cd->Draw("colz");
  /*TLine * line_pq_cd_1 = new TLine(0.62,0,0.62,60);
  line_pq_cd_1->SetLineColor(kRed);
  line_pq_cd_1->SetLineWidth(3);
  line_pq_cd_1->Draw("same");*/
  TLine * line_pq_cd_2 = new TLine(0.96,0,0.96,60);
  line_pq_cd_2->SetLineColor(kRed);
  line_pq_cd_2->SetLineWidth(3);
  line_pq_cd_2->Draw("same");
  /*TLine * line_pq_cd_3 = new TLine(0,25,1.2,25);
  line_pq_cd_3->SetLineColor(kRed);
  line_pq_cd_3->SetLineWidth(3);
  line_pq_cd_3->Draw("same");*/
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



  ///////////////////////////////////////////////////////
  //Final Distributions after SRC Cuts
  ///////////////////////////////////////////////////////  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xb_final->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_final->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_final->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_final->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_final->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  ///////////////////////////////////////////////////////
  //(e,e'p) Kinematics
  ///////////////////////////////////////////////////////  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_emiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_emiss_omega->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_omega->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_emiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  ///////////////////////////////////////////////////////
  //Lead and Recoil (e,e'pp)
  ///////////////////////////////////////////////////////  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_lead_rec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_prec->Draw();
  TLine * line_prec = new TLine(0.35,0,0.35,h_prec->GetMaximum());
  line_prec->SetLineColor(kRed);
  line_prec->SetLineWidth(3);
  line_prec->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  ///////////////////////////////////////////////////////
  //(e,e'pp) Kinematics
  ///////////////////////////////////////////////////////  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_emiss_pp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_emiss_omega_pp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xb_pp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_pp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_q2_pp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapq_pq_pp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_plead_prec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_tlead_trec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_prec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapmq_pp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  ///////////////////////////////////////////////////////
  //Additional Histos
  ///////////////////////////////////////////////////////  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_plead_all->Draw();
  TLine * line_plead_all = new TLine(1,0,1,h_plead_all->GetMaximum());
  line_plead_all->SetLineColor(kRed);
  line_plead_all->SetLineWidth(3);
  line_plead_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_plead_fd->Draw();
  TLine * line_plead_fd = new TLine(1,0,1,h_plead_fd->GetMaximum());
  line_plead_fd->SetLineColor(kRed);
  line_plead_fd->SetLineWidth(3);
  line_plead_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_plead_cd->Draw();
  TLine * line_plead_cd = new TLine(1,0,1,h_plead_cd->GetMaximum());
  line_plead_cd->SetLineColor(kRed);
  line_plead_cd->SetLineWidth(3);
  line_plead_cd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapmq_pq_all->Draw("colz");
  TLine * line_thetapmq_pq_all = new TLine(0.96,0,0.96,180);
  line_thetapmq_pq_all->SetLineColor(kRed);
  line_thetapmq_pq_all->SetLineWidth(3);
  line_thetapmq_pq_all->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapmq_pq_fd->Draw("colz");
  TLine * line_thetapmq_pq_fd = new TLine(0.96,0,0.96,180);
  line_thetapmq_pq_fd->SetLineColor(kRed);
  line_thetapmq_pq_fd->SetLineWidth(3);
  line_thetapmq_pq_fd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetapmq_pq_cd->Draw("colz");
  TLine * line_thetapmq_pq_cd = new TLine(0.96,0,0.96,180);
  line_thetapmq_pq_cd->SetLineColor(kRed);
  line_thetapmq_pq_cd->SetLineWidth(3);
  line_thetapmq_pq_cd->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_doublelead->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

