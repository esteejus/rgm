#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TCutG.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include <mysql.h>
#include "RCDB/Connection.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}

void LowEnergyReader(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!
  TString inputFile;
  TString outputFile;

  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".hipo"))){
      inputFile=opt(5,opt.Sizeof());
    }
  }
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////
  //recon_qe_he4_6gev_9_torus-1.0.hipo
  outputFile = inputFile(inputFile.Index("recon_qe_")+9,inputFile.Index(".hipo"));

  cout<<"Analysing hipo file "<<inputFile<<endl;
  cout<<"Outfile "<<outputFile<<endl;
  //  clas12databases::SetCCDBLocalConnection( string(gSystem->Getenv("PWD")) + "/db/ccdb.sqlite");
  //  clas12databases::SetQADBConnection( string(gSystem->Getenv("PWD")) + "/db/qaDB.json");
  //  clas12databases::SetRCDBRootConnection( string(gSystem->Getenv("PWD")) + "rcdb_LowEnergy.root");

  //  clas12root::HipoChain chain;
  //  chain.Add(inputFile.Data());

  auto pdg_db=TDatabasePDG::Instance();
  double mass_e = pdg_db->GetParticle(11)->Mass();
  double mass_p = pdg_db->GetParticle(2212)->Mass();
  double mass_n = pdg_db->GetParticle(2112)->Mass();

  TLorentzVector beam(0,0,6,6);
  //  TLorentzVector beam(0,0,2,2);
  TLorentzVector target(0,0,0,11.1779);
  TLorentzVector el(0,0,0,mass_e);
  TLorentzVector pr(0,0,0,mass_p);
  TLorentzVector lead_pr(0,0,0,mass_p);
  TLorentzVector recoil_pr(0,0,0,mass_p);

  TLorentzVector ntr(0,0,0,mass_n);

  TLorentzVector el_t(0,0,0,0);
  TLorentzVector lead_t(0,0,0,0);
  TLorentzVector recoil_t(0,0,0,0);

  int counter=0;
  Double_t x_b = 0;
  Double_t x_prime = 0;
  Double_t theta = 0;
  Double_t p_q = 0;
  Double_t miss_m = 0;
  Double_t miss_p = 0;
  Double_t miss_pz = 0;
  Double_t q2 = 0;
  Double_t vz_diff = 0;
  Double_t vz = 0;

  Bool_t ecal = false;
  Bool_t p_cd = false;
  Bool_t p_fd = false;

  /*
  TFile *tree_file = new TFile("/work/clas12/users/esteejus/output/LowEnergy/tree_"+outputFile+".root","RECREATE");

  TTree *tree = new TTree("tree","Low Energy Data");
  tree->Branch("xb",&x_b,"xb/D");
  tree->Branch("x_prime",&x_prime,"x_prime/D");
  tree->Branch("theta",&theta,"theta/D");
  tree->Branch("p_q",&p_q,"p_q/D");
  tree->Branch("miss_m",&miss_m,"miss_m/D");
  tree->Branch("miss_p",&miss_p,"miss_p/D");
  tree->Branch("miss_pz",&miss_pz,"miss_pz/D");
  tree->Branch("q2",&q2,"q2/D");
  tree->Branch("vz_diff",&vz_diff,"vz_diff/D");
  tree->Branch("vz",&vz,"vz/D");

  tree->Branch("counter",&counter,"counter/I");
  tree->Branch("ecal",&ecal,"ecal/B");
  tree->Branch("p_cd",&p_cd,"p_cd/B");
  tree->Branch("p_fd",&p_fd,"p_fd/B");
  */

  auto *p_minus_a = new TH1D("p_minus_a","PID -proton",180,0,180);
  auto *p_a = new TH1D("p_a","PID proton",180,0,180);

  auto *pid_p = new TH1D("pid_p","PID of proton",3,-1,2);

  auto *pid_e_ec = new TH2D("pid_p_ec","EC/p vs e electron",100,0,5,100,0,.6);
  auto *pid_p_fd = new TH2D("pid_p_fd","FD TOF vs p PID proton",500,0,4,500,0,1.3);
  auto *pid_p_cd = new TH2D("pid_p_cd","CD TOF vs p PID proton",200,0,4,200,0,1.3);

  auto *chi2_p = new TH1D("chi2_p","Chi2 proton",100,-20,20);
  auto *chi2_e = new TH1D("chi2_e","Chi2 electron",100,-20,20);

  auto *hbeam_e = new TH1F("beam_e","Beam Energy",100,3,5);
  auto *hmiss = new TH1F("missM","missM",200,-1,4);
  auto *missm_pmiss = new TH2F("missm_pmiss","Missing mass vs. pmiss",600,-3,3,100,0,4);
  auto *missm_xb = new TH2F("missm_xb","Missing mass vs. x_{B}",600,-3,3,100,0,4);


  auto *q_xb_dist_truth=new TH2F("q_xb_dist_truth","Q xb distribution truth",1000,-1,4,100,-1,4);
  auto *q_xb_dist_pp =new TH2F("q_xb_dist_pp","Q xb distribution",1000,-1,4,100,-1,4); 
  auto *q_xb_dist_p =new TH2F("q_xb_dist_p","Q xb distribution",1000,-1,4,100,-1,4); 

  auto *xb_dist=new TH1F("xb_dist","xb distribution",100,-1,4);
  auto *q_dist=new TH1F("q_dist","#Q^2 distribution",1000,-1,4);
  auto *theta_pq_ratio = new TH2D("theta_pq_ratio","Theata pq vs p/q",1000,0,2,1000,0,50);
  auto *theta_pmiss = new TH2D("theta_pmiss","Pmiss vs Theta",180,0,180,100,0,4);
  auto *theta_phi = new TH2D("theta_phi","Phi vs theta",180,0,180,360,-180,180);

  auto *vz_corr = new TH1D("vz_corr","e_vz - p_vz",1000,-30,30);

  auto *vz_el_p = new TH2D("vz_el_p","vz_el_p",1000,-30,30,1000,-30,30);

  auto *xbor=new TH1F("xBor","xBor",200,-1,10);   




  //A(e,e'p) vs pmiss
  //A(e,e'pp) vs pmiss
  //A(e,e'pn) vs pmiss

  int bin_num = 100;

  auto *p_pmiss = new TH1F("p_pmiss","A(e,e'p) vs. pmiss",bin_num,0,4);
  auto *pp_pmiss = new TH1F("pp_pmiss","A(e,e'pp) vs. pmiss",bin_num,0,4);
  auto *pn_pmiss = new TH1F("pn_pmiss","A(e,e'pn) vs. pmiss",bin_num,0,4);

  auto *p_pmiss_w = new TH1F("p_pmiss_w","A(e,e'p) vs. pmiss",bin_num,0,4);
  auto *pp_pmiss_w = new TH1F("pp_pmiss_w","A(e,e'pp) vs. pmiss",bin_num,0,4);
  auto *pn_pmiss_w = new TH1F("pn_pmiss_w","A(e,e'pn) vs. pmiss",bin_num,0,4);

  auto *pmiss_resol_cd_w = new TH2F("pmiss_resol_cd_w","Proton resolution weight",bin_num,0,4,200,-.2,.2);
  auto *pmiss_resol_cd = new TH2F("pmiss_resol_cd","Proton resolution",bin_num,0,4,200,-.2,.2);
  auto *p_resol_cd_w = new TH2F("p_resol_cd_w","Proton resolution weight",bin_num,0,4,200,-.2,.2);
  auto *p_resol_cd = new TH2F("p_resol_cd","Proton resolution",bin_num,0,4,200,-.2,.2);

  auto *pmiss_resol_fd_w = new TH2F("pmiss_resol_fd_w","Proton resolution weight",bin_num,0,4,200,-.2,.2);
  auto *pmiss_resol_fd = new TH2F("pmiss_resol_fd","Proton resolution",bin_num,0,4,200,-.2,.2);
  auto *p_resol_fd_w = new TH2F("p_resol_fd_w","Proton resolution weight",bin_num,0,4,200,-.2,.2);
  auto *p_resol_fd = new TH2F("p_resol_fd","Proton resolution",bin_num,0,4,200,-.2,.2);

  //Lead proton anglular distribution
  //Recoil proton angular distribution
  //Multiplicity of proton, neutron, electron

  auto *p_pmiss_truth = new TH1F("p_pmiss_truth","A(e,e'p) vs. pmiss_truth",bin_num,0,4);
  auto *pp_pmiss_truth = new TH1F("pp_pmiss_truth","A(e,e'pp) vs. pmiss_truth",bin_num,0,4);
  auto *pn_pmiss_truth = new TH1F("pn_pmiss_truth","A(e,e'pn) vs. pmiss_truth",bin_num,0,4);

  auto *pq_truth_lead = new TH2F("pq_truth_lead","#theta vs |p|/|q| Truth lead",100,0,1,100,0,100);
  auto *pq_truth_recoil = new TH2F("pq_truth_recoil","#theta vs |p|/|q| Truth recoil",100,0,1,100,0,100);

  auto *p_pmiss_truth_w = new TH1F("p_pmiss_truth_w","A(e,e'p) vs. pmiss_truth",bin_num,0,4);
  auto *pp_pmiss_truth_w = new TH1F("pp_pmiss_truth_w","A(e,e'pp) vs. pmiss_truth",bin_num,0,4);
  auto *pn_pmiss_truth_w = new TH1F("pn_pmiss_truth_w","A(e,e'pn) vs. pmiss_truth",bin_num,0,4);

  auto *pq_truth_lead_w = new TH2F("pq_truth_lead_w","#theta vs |p|/|q| Truth lead",100,0,1,100,0,100);
  auto *pq_truth_recoil_w = new TH2F("pq_truth_recoil_w","#theta vs |p|/|q| Truth recoil",100,0,1,100,0,100);

  gBenchmark->Start("timer");
 
  TFile *cutf = TFile::Open("cuts_v1.root");
  TCutG *elec_cut = (TCutG *)cutf->Get("elec");
  TCutG *ftof_p = (TCutG *)cutf->Get("ftof_p");
  TCutG *ctof_p = (TCutG *)cutf->Get("ctof_p");
  
  clas12reader c12(inputFile.Data());

    while(c12.next()==true){

      auto mceve=c12.mcevent(); 

      x_b = 0;
      x_prime = 0;
      theta = 0;
      p_q = 0;
      miss_m = 0;
      miss_p = 0;
      miss_pz = 0;
      q2 = 0;
      vz_diff = 0;
      vz = 0;
      
      ecal = false;
      p_cd = false;
      p_fd = false;

      //can get an estimate of the beam current to this event
      //c12.getCurrApproxCharge();//if called c12.scalerReader();
	
      //c12.event()->getStartTime();

	
      //Loop over all particles to see how to access detector info.
    
      for(auto& p : c12.getDetParticles()){
	//  get predefined selected information
	p->getTime();
	p->getDetEnergy();
	p->getDeltaEnergy();

	//check trigger bits
	//	 if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
	//	 else cout<<"NOT"<<endl;

	// get any detector information (if exists for this particle)
	// there should be a get function for any entry in the bank
	switch(p->getRegion()) {
	case FD :
	  p->cal(PCAL)->getEnergy();
	  p->cal(PCAL)->getLu();
	  p->cal(PCAL)->getLv();
	  p->cal(PCAL)->getLw();
	  p->cal(ECIN)->getEnergy();
	  p->cal(ECOUT)->getEnergy();
	  p->sci(FTOF1A)->getEnergy();
	  p->sci(FTOF1B)->getEnergy();
	  p->sci(FTOF2)->getEnergy();
	  p->trk(DC)->getSector();
	  p->trk(DC)->getChi2();
	  p->che(HTCC)->getNphe();
	  p->che(LTCC)->getNphe();
	  //trajectories
	  p->traj(LTCC)->getX();
	  // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
	  break;
	case FT :
	  p->ft(FTCAL)->getEnergy();
	  p->ft(FTHODO)->getEnergy();
	  break;
	case CD:
	  p->sci(CTOF)->getEnergy();
	  p->sci(CND)->getEnergy();
	  break;
	}
	//   covariance matrix (comment in to see!)
	// p->covmat()->print();
	p->cmat();
      }

      // MC::Lund
      // For the jth particle in the lund file
      //first particle is electron, lead, recoil

      double mc_weight = mceve->getWeight();
      cout<<"Weight "<<mc_weight<<endl;
      mc_weight = mc_weight/10000; //Need to change in future
      int mc_p = 0;
      int mc_n = 0;

      el_t.SetXYZM( c12.mcparts()->getPx(0),c12.mcparts()->getPy(0),c12.mcparts()->getPz(0),c12.mcparts()->getMass(0));
      lead_t.SetXYZM( c12.mcparts()->getPx(1),c12.mcparts()->getPy(1),c12.mcparts()->getPz(1),c12.mcparts()->getMass(1));
      recoil_t.SetXYZM( c12.mcparts()->getPx(2),c12.mcparts()->getPy(2),c12.mcparts()->getPz(2),c12.mcparts()->getMass(2));


      bool lead_is_p = (c12.mcparts()->getPid(1)==2212);

      for(int iMC = 1; iMC < 3; iMC++)
	{
	  if( c12.mcparts()->getPid(iMC) == 2212)
	    mc_p++;
	  else if( c12.mcparts()->getPid(iMC) == 2112)
	    mc_n++;
	}

      TLorentzVector miss_ex_truth = beam + target - el_t - lead_t; //missing 4-vector
      TLorentzVector q_truth = beam - el_t; //missing 4-vector
      TVector3 miss_p_truth = miss_ex_truth.Vect();
      double x_prime_truth = -q_truth.M2()/(2 * miss_p_truth.Dot(q_truth.Vect()) * (beam.E() - el_t.E()) ); //x-borken
      double  x_b_truth = -q_truth.M2()/(2 * mass_p * (beam.E() - el_t.E()) ); //x-borken


      double theta_truth_recoil = recoil_t.Vect().Angle(q_truth.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q
      double p_q_truth_recoil = recoil_t.Vect().Mag()/q_truth.Vect().Mag();

      double theta_truth_lead = lead_t.Vect().Angle(q_truth.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q
      double p_q_truth_lead = lead_t.Vect().Mag()/q_truth.Vect().Mag();

      xb_dist -> Fill(x_b_truth,mc_weight);

      //      if(x_b_truth > 1.2 && p_q_truth_lead > .6 && p_q_truth_lead < .9 && theta_truth_lead < 25)
      if(x_b_truth > .5 && p_q_truth_lead > .6 && p_q_truth_lead < .9 && theta_truth_lead < 25)
	{


	  if(mc_p == 2)
	    {
	      pp_pmiss_truth_w -> Fill(miss_ex_truth.P(),mc_weight);
	      pq_truth_recoil_w->Fill(p_q_truth_recoil,theta_truth_recoil,mc_weight);
	      pq_truth_lead_w->Fill(p_q_truth_lead,theta_truth_lead,mc_weight);
	      
	      pp_pmiss_truth -> Fill(miss_ex_truth.P());
	      pq_truth_recoil->Fill(p_q_truth_recoil,theta_truth_recoil);
	      pq_truth_lead->Fill(p_q_truth_lead,theta_truth_lead);

	      q_xb_dist_truth->Fill(x_b_truth,-q_truth.M2(),mc_weight);

	      
	    }
	  else if(mc_p == 1 && mc_n == 1)
	    {
	      pn_pmiss_truth_w -> Fill(miss_ex_truth.P(),mc_weight);
	      pn_pmiss_truth -> Fill(miss_ex_truth.P());
	      cout<<"pmiss "<<miss_ex_truth.P()<<endl;
	    }

	  if(lead_is_p)
	    {
	      p_pmiss_truth_w -> Fill(miss_ex_truth.P(),mc_weight);
	      p_pmiss_truth -> Fill(miss_ex_truth.P());
	    }

	}

      // get particles by type
      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      auto protons_m=c12.getByID(-2212);
      auto neutrons=c12.getByID(2112);

      pid_p->Fill(-1,protons_m.size());
      pid_p->Fill(1,protons.size());

      /*
      if(protons_m.size() == 1)
	{
	  SetLorentzVector(pr1,protons_m[0]);
	  p_minus_a->Fill(pr1.Theta()*TMath::RadToDeg());

	}
      if(protons.size() == 1)
	{
	  SetLorentzVector(pr1,protons[0]);
	  p_a->Fill(pr1.Theta()*TMath::RadToDeg());

	}
      */

      if(electrons.size()==1)
	{

	  double energy =  electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy();

	  bool ecal = (electrons[0]->getRegion()==FD && (electrons[0]->par()->getP() < beam.E() && electrons[0]->par()->getP() > 1) && ( energy/electrons[0]->par()->getP() < .275 && energy/electrons[0]->par()->getP() > .175) );


	  if(!ecal)
	    continue;


	  if (protons.size()==2)
	    {

	      bool p1_fd = (protons[0]->getRegion()==FD && abs(protons[0]->par()->getChi2Pid()) < 3);
	      bool p2_fd = (protons[1]->getRegion()==FD && abs(protons[1]->par()->getChi2Pid()) < 3);

	      bool p1_cd = (protons[0]->getRegion()==CD && abs(protons[0]->par()->getChi2Pid()) < 3);
	      bool p2_cd = (protons[1]->getRegion()==CD && abs(protons[1]->par()->getChi2Pid()) < 3);

	      if( !((p1_fd || p1_cd) && (p2_fd || p2_cd)) )
		continue;


	  vz_diff = electrons[0]->par()->getVz() - protons[0]->par()->getVz();
	  vz = electrons[0]->par()->getVz(); 

	  double vz_diff_p = electrons[0]->par()->getVz() + protons[0]->par()->getVz() + 5.5;
	  vz_corr->Fill(vz_diff);
	  

	  if( !(abs(vz_diff) < 2 && (vz < 6 && vz > -6) ) )
	    continue;

	  SetLorentzVector(el,electrons[0]);

	  double lead_miss = 0;
	  double lead_pr_p = 0;
	  int lead_prs = 0;
	  bool is_lead_fd = false;

	  for(int iPr = 0; iPr < 2; iPr++)
	    {
	      SetLorentzVector(pr,protons[iPr]);
	      TLorentzVector miss=beam+target-el-pr; //missing 4-vector
	      TLorentzVector q = beam - el;          //photon  4-vector
	      miss_p = miss.P();
	      miss_pz = miss.Pz();
	      miss_m = miss.M();
	      q2 = -q.M2();
	      x_prime = -q.M2()/(2 * miss.Dot(q) * (beam.E() - el.E()) ); //x-borken
	      x_b = -q.M2()/(2 * mass_p * (beam.E() - el.E()) ); //x-borken
	      
	      theta = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q
	      p_q = pr.Vect().Mag()/q.Vect().Mag();


	      if(theta < 25 && p_q < 1 && p_q > .6 && p_q < .9 && x_b > .5)
	      //	      if(theta < 25 && p_q < 1 && p_q > .6 && p_q < .9 && x_b > 1.2)
		{
		  cout<<"here ?"<<endl;
		  lead_miss = miss.P();
		  lead_prs++;
		  lead_pr_p = pr.P();

		  is_lead_fd = (protons[iPr]->getRegion()==FD && abs(protons[iPr]->par()->getChi2Pid()) < 3);
		}
	    }
	  
	  //	  theta_pmiss->Fill(miss.Theta()*TMath::RadToDeg(),miss.P());
	  if(lead_prs==1)
	    {
	      pp_pmiss_w->Fill(lead_miss,mc_weight);
	      pp_pmiss->Fill(lead_miss);

	      q_xb_dist_p->Fill(x_b,q2,mc_weight);

	      cout<<"lead miss "<<lead_miss<<endl;
	      if(lead_is_p)
		{

		  if(is_lead_fd)
		    {
		      pmiss_resol_fd_w->Fill( miss_ex_truth.P(),(lead_miss - miss_ex_truth.P())/miss_ex_truth.P(),mc_weight);
		      pmiss_resol_fd->Fill( miss_ex_truth.P(), (lead_miss - miss_ex_truth.P())/miss_ex_truth.P());

		      p_resol_fd_w->Fill( lead_t.P(),(lead_pr_p - lead_t.P())/lead_t.P(),mc_weight);
		      p_resol_fd->Fill( lead_t.P(), (lead_pr_p - lead_t.P())/lead_t.P() );
		    }
		  else
		    {
		      pmiss_resol_cd_w->Fill( miss_ex_truth.P(),(lead_miss - miss_ex_truth.P())/miss_ex_truth.P() ,mc_weight);
		      pmiss_resol_cd->Fill( miss_ex_truth.P(), (lead_miss - miss_ex_truth.P())/miss_ex_truth.P() );

		      p_resol_cd_w->Fill( lead_t.P(),(lead_pr_p - lead_t.P())/lead_t.P() ,mc_weight);
		      p_resol_cd->Fill( lead_t.P(), (lead_pr_p - lead_t.P())/lead_t.P() );
		    }

		}
	    }

	}//(e,e'pp) events
      


	  if (protons.size()==1)
	    {
	      


	      bool p1_fd = (protons[0]->getRegion()==FD && abs(protons[0]->par()->getChi2Pid()) < 3);
	      bool p1_cd = (protons[0]->getRegion()==CD && abs(protons[0]->par()->getChi2Pid()) < 3);

	      if( !(p1_fd || p1_cd) )
		continue;
	      
	      vz_diff = electrons[0]->par()->getVz() - protons[0]->par()->getVz();
	      vz = electrons[0]->par()->getVz(); 
	      
	      double vz_diff_p = electrons[0]->par()->getVz() + protons[0]->par()->getVz() + 5.5;
	      vz_corr->Fill(vz_diff);
	      
	      
	      if( !(abs(vz_diff) < 2 && (vz < 6 && vz > -6) ) )
		continue;
	      
	      SetLorentzVector(el,electrons[0]);
	      SetLorentzVector(pr,protons[0]);
	      TLorentzVector miss=beam+target-el-pr; //missing 4-vector
	      TLorentzVector q = beam - el;          //photon  4-vector
	  
	      miss_p = miss.P();
	      miss_pz = miss.Pz();
	      miss_m = miss.M();
	      q2 = -q.M2();
	      x_prime = -q.M2()/(2 * miss.Dot(q) * (beam.E() - el.E()) ); //x-borken
	      x_b = -q.M2()/(2 * mass_p * (beam.E() - el.E()) ); //x-borken


	      theta = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q
	      p_q = pr.Vect().Mag()/q.Vect().Mag();
	  
	      
	      //	      if(theta < 25 && p_q < 1 && p_q > .6 && p_q < .9 && x_b > 1.2)
	      if(theta < 25 && p_q < 1 && p_q > .6 && p_q < .9 && x_b > .5)
		{
		  p_pmiss->Fill(miss.P());
		  p_pmiss_w->Fill(miss.P(),mc_weight);

		  q_xb_dist_pp->Fill(x_b,q2,mc_weight);
		}	      
	      
	    }//(e,e'p) events
	  


	}

      counter++;
    }
    
    gBenchmark->Stop("timer");
    gBenchmark->Print("timer");
    
    
    TCanvas* can = new TCanvas();
    xb_dist->Draw();
    //    pid_p->Draw();
    //    pp_pmiss_truth->Draw();
  
    TCanvas* can1 = new TCanvas();
    //    pn_pmiss_truth->Draw();
    pn_pmiss->SetLineColor(2);
    //    pn_pmiss->Draw();
    //    pid_p->Draw();
    p_minus_a->SetLineColor(2);

    //    p_a->Draw();
    //    p_minus_a->Draw("same");

    TFile *out = new TFile("/volatile/clas12/users/esteejus/histoutput/hist_"+outputFile+".root","RECREATE");
    pn_pmiss_truth->Write();
    pp_pmiss_truth->Write();
    p_pmiss_truth->Write();

    pmiss_resol_fd_w->Write();
    pmiss_resol_fd->Write();
    p_resol_fd_w->Write();
    p_resol_fd->Write();

    pmiss_resol_cd_w->Write();
    pmiss_resol_cd->Write();
    p_resol_cd_w->Write();
    p_resol_cd->Write();

    p_pmiss_truth_w->Write();
    pn_pmiss_truth_w->Write();
    pp_pmiss_truth_w->Write();

    pq_truth_recoil_w->Write();
    pq_truth_lead_w->Write();

    q_xb_dist_truth->Write();
    q_xb_dist_p->Write();
    q_xb_dist_pp->Write();

    pn_pmiss->Write();
    pp_pmiss->Write();
    p_pmiss->Write();
    theta_pmiss->Write();

    pn_pmiss_w->Write();
    pp_pmiss_w->Write();
    p_pmiss_w->Write();

    pq_truth_recoil_w->Write();
    pq_truth_lead_w->Write();
    pid_p->Write();
    xb_dist->Write();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

    }
