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

using namespace clas12;


bool pointsToBand(double theta,double phi,double z_m){
  double z = z_m*100; // from m to cm

  // Numbers taken from band/src/main/java/org/jlab/rec/band/constants/Parameters.java
  double thickness  = 7.2;                                // thickness of each bar (cm)
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

  double globalX = (-240.5-240.5+241.0+243.7)/4.; // [cm] --> Not using this yet (need to make sure we have the right coordinate system)
  double globalY = (-211.0+228.1-210.6+228.1)/4.; // [cm]

  // Sector boundaries
  double topSec1  = globalY + 13*thickness;
  double topSec2  = globalY + 10*thickness;
  double topSec34 = globalY +  3*thickness;
  double topSec5  = globalY -  3*thickness;
  double downSec5 = globalY -  5*thickness;

  if( yDown >= topSec1 || yDown <= downSec5 ) return 0;

  if(             (yDown < topSec1  && yDown >= topSec2  && fabs(xDown) < bandlen[0]/2. )||
		  (yDown < topSec2  && yDown >= topSec34 && fabs(xDown) < bandlen[1]/2. )||
		  (yDown < topSec34 && yDown >= topSec5  && fabs(xDown) < bandlen[1]/2. && fabs(xDown) > bandlen[1]/2.-bandlen[2])||
		  (yDown < topSec5  && yDown >= downSec5 && fabs(xDown) < bandlen[4]/2. )
		  )
    return 1;

  return 0;
}


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	     rp->par()->getPz(),p4.M());

}

void LowEnergyReader_MC(){
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
  outputFile = inputFile(inputFile.Index("recon_qe"),inputFile.Index(".hipo"));

  cout<<"Analysing hipo file "<<inputFile<<endl;

  //  clas12databases::SetCCDBLocalConnection( string(gSystem->Getenv("PWD")) + "/db/ccdb.sqlite");
  //  clas12databases::SetQADBConnection( string(gSystem->Getenv("PWD")) + "/db/qaDB.json");
  //  clas12databases::SetRCDBRootConnection( string(gSystem->Getenv("PWD")) + "rcdb_LowEnergy.root");

  //  clas12root::HipoChain chain;
  //  chain.Add(inputFile.Data());

  auto pdg_db=TDatabasePDG::Instance();
  double mass_e = pdg_db->GetParticle(11)->Mass();
  double mass_p = pdg_db->GetParticle(2212)->Mass();
  double mass_n = pdg_db->GetParticle(2112)->Mass();

  TLorentzVector beam(0,0,4.1,4.1);
  TLorentzVector target(0,0,0,1.8756129);
  TLorentzVector el(0,0,0,mass_e);
  TLorentzVector pr(0,0,0,mass_p);


  Double_t x_b = 0;
  Double_t x_prime = 0;
  Double_t theta = 0;
  Double_t p_q = 0;
  Double_t miss_m = 0;
  Double_t miss_p = 0;
  Double_t miss_m_MC = 0;
  Double_t miss_p_MC = 0;
  Double_t miss_pz = 0;
  Double_t q2 = 0;
  Double_t vz_diff = 0;
  Double_t vz = 0;

  Bool_t ecal = false;
  Bool_t p_cd = false;
  Bool_t p_fd = false;
  Bool_t point_band = false;


  TFile *tree_file = new TFile("/work/clas12/users/esteejus/output/LowEnergyMC/MCtree_"+outputFile+".root","RECREATE");

  TTree *tree = new TTree("tree","Low Energy Data");
  tree->Branch("xb",&x_b,"xb/D");
  tree->Branch("x_prime",&x_prime,"x_prime/D");
  tree->Branch("theta",&theta,"theta/D");
  tree->Branch("p_q",&p_q,"p_q/D");
  tree->Branch("miss_m",&miss_m,"miss_m/D");
  tree->Branch("miss_p",&miss_p,"miss_p/D");
  tree->Branch("miss_m_MC",&miss_m_MC,"miss_m_MC/D");
  tree->Branch("miss_p_MC",&miss_p_MC,"miss_p_MC/D");
  tree->Branch("miss_pz",&miss_pz,"miss_pz/D");
  tree->Branch("q2",&q2,"q2/D");
  tree->Branch("vz_diff",&vz_diff,"vz_diff/D");
  tree->Branch("vz",&vz,"vz/D");

  tree->Branch("ecal",&ecal,"ecal/B");
  tree->Branch("p_cd",&p_cd,"p_cd/B");
  tree->Branch("p_fd",&p_fd,"p_fd/B");
  tree->Branch("point_band",&point_band,"point_band/B");

  auto *pid_e_ec = new TH2D("pid_p_ec","EC/p vs e electron",100,0,5,100,0,.6);
  auto *pid_p_fd = new TH2D("pid_p_fd","FD TOF vs p PID proton",500,0,4,500,0,1.3);
  auto *pid_p_cd = new TH2D("pid_p_cd","CD TOF vs p PID proton",200,0,4,200,0,1.3);

  auto *chi2_p = new TH1D("chi2_p","Chi2 proton",100,-20,20);
  auto *chi2_e = new TH1D("chi2_e","Chi2 electron",100,-20,20);

  auto *hbeam_e = new TH1F("beam_e","Beam Energy",100,3,5);
  auto *hmiss = new TH1F("missM","missM",200,-1,4);
  auto *missm_pmiss = new TH2F("missm_pmiss","Missing mass vs. pmiss",600,-3,3,100,0,4);
  auto *missm_xb = new TH2F("missm_xb","Missing mass vs. x_{B}",600,-3,3,100,0,4);

  auto *q_dist=new TH1F("q_dist","#Q^2 distribution",1000,-1,4);
  auto *theta_pq_ratio = new TH2D("theta_pq_ratio","Theata pq vs p/q",1000,0,2,1000,0,50);
  auto *vz_corr = new TH1D("vz_corr","e_vz - p_vz",1000,-30,30);
  auto *theta_pmiss = new TH2D("theta_pmiss","Pmiss vs Theta",180,0,180,100,0,4);
  auto *theta_phi = new TH2D("theta_phi","Phi vs theta",180,0,180,360,-180,180);

  auto *vz_el_p = new TH2D("vz_el_p","vz_el_p",1000,-30,30,1000,-30,30);

  auto *xbor=new TH1F("xBor","xBor",200,-1,10);   

  gBenchmark->Start("timer");
  int counter=0;
 
  TFile *cutf = TFile::Open("cuts_v1.root");
  TCutG *elec_cut = (TCutG *)cutf->Get("elec");
  TCutG *ftof_p = (TCutG *)cutf->Get("ftof_p");
  TCutG *ctof_p = (TCutG *)cutf->Get("ctof_p");
  
  clas12reader c12(inputFile.Data());
  //  clas12databases db;
  //  c12.connectDataBases(&db);
  //  for(Int_t iFile = 0; iFile < chain.GetNFiles(); iFile++) {
      //create the event reader
  //      clas12reader c12(chain.GetFileName(iFile).Data());
      //rcdb info
  /*
  c12.queryRcdb();
  auto rcdbData= c12.getRcdbVals();
      //The following run conditions can be returned directly by c12
      cout<<"Event count: "<<rcdbData.event_count<<endl;
      cout<<"Beam energy: "<<rcdbData.beam_energy<<endl;
      cout<<"Beam current: "<<rcdbData.beam_current<<endl;
  */
      //      beam.SetE(rcdbData.beam_energy/1000);
      //      beam.SetPz(rcdbData.beam_energy/1000);
      //      auto beam_e = rcdbData.beam_energy;
      //      cout<<"Beam "<<beam_e<<endl;
      //      hbeam_e->Fill(rcdbData.beam_energy/1000);
      //    clas12databases db;
      //    c12.connectDataBases(&db);
    
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}
      
    //Add some event Pid based selections
    //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
    //c12.addExactPid(11,1);    //exactly 1 electron
    //c12.addExactPid(211,1);    //exactly 1 pi+
    //c12.addExactPid(-211,1);    //exactly 1 pi-
    //c12.addExactPid(2212,1);    //exactly 1 proton
    //c12.addExactPid(22,2);    //exactly 2 gamma
    //////c12.addZeroOfRestPid();  //nothing else
    //////c12.useFTBased(); //and use the Pids from RECFT

    //can also access the integrated current at this point
    //c12.scalerReader();//must call this first
    //c12.getRunBeamCharge();
    //    c12.getRcdbVals();

    //c12.setEntries(1E3); //only process 1E3 events per file
    while(c12.next()==true){

      x_b = 0;
      x_prime = 0;
      theta = 0;
      p_q = 0;
      miss_m = 0;
      miss_p = 0;
      miss_m_MC = 0;
      miss_p_MC = 0;
      miss_pz = 0;
      q2 = 0;
      vz_diff = 0;
      vz = 0;

      ecal = false;
      p_cd = false;
      p_fd = false;
      point_band = false;


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

      //MC::Lund
      // For the jth particle in the lund file
      // 0:e , 1:p ,2:n in this Lund file
      /*
      c12.mcparts()->setEntry(2);
      c12.mcparts()->getPid();
      c12.mcparts()->getPx();
      c12.mcparts()->getPy();
      c12.mcparts()->getPz();
      c12.mcparts()->getVx();
      c12.mcparts()->getVy();
      c12.mcparts()->getVz();
      c12.mcparts()->getMass();
      TVector3 mc_neutron(c12.mcparts()->getPx(), c12.mcparts()->getPy(),c12.mcparts()->getPz());
      */

      TLorentzVector el_MC;
      el_MC.SetXYZM(c12.mcparts()->getPx(0),c12.mcparts()->getPy(0),c12.mcparts()->getPz(0), c12.mcparts()->getMass(0) );

      TLorentzVector pr_MC;
      pr_MC.SetXYZM(c12.mcparts()->getPx(1),c12.mcparts()->getPy(1),c12.mcparts()->getPz(1), c12.mcparts()->getMass(1) );

      TLorentzVector miss_MC=beam+target-el_MC-pr_MC; //missing 4-vector

      if(counter < 10)
	{

	  cout<<"el "<<c12.mcparts()->getPx(0)<<" "<< c12.mcparts()->getPy(0)<<" " << c12.mcparts()->getPz(0)<<" "<< c12.mcparts()->getMass(0)<<endl;

	  cout<<"pr "<<c12.mcparts()->getPx(1)<<" "<< c12.mcparts()->getPy(1)<<" " << c12.mcparts()->getPz(1)<<" "<< c12.mcparts()->getMass(1)<<endl;

	  cout<<"Lorentz "<< miss_MC.Px()<<" "<< miss_MC.Py()<<" "<< miss_MC.Pz()<<" "<< miss_MC.E()<<endl;

	cout<<"mis mass "<<miss_MC.M()<<endl;

	}
      // get particles by type
      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      //       auto gammas=c12.getByID(22);
      //       auto pips=c12.getByID(211);
      //       auto pims=c12.getByID(-211);
       
      //      if(protons.size()==1)
	//	cout<< protons[0]->getRegion()<<endl;

      if(electrons.size()==1 && protons.size()==1 ){
       

	miss_m_MC = miss_MC.M();
	miss_p_MC = miss_MC.P();

	//       	if(protons[0]->getRegion()==FD && !ftof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta()) )
	//	  continue;

	//	if( protons[0]->getRegion()==CD && ctof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta() )
	   //	    continue;
	    
	//	if( electronss[0]->getRegion()==FD && elec_cut ->IsInside(electrons[0]->par()->getP(),electrons[0]->cal(PCAL)->getEnergy()/electrons[0]-par()->getP() )
	   //	    continue;

	double energy =  electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy();

	if(protons[0]->getRegion()==FD)
	  pid_p_fd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );

	if(protons[0]->getRegion()==CD)
	  pid_p_cd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );

	if(electrons[0]->getRegion()==FD)
	  pid_e_ec -> Fill(electrons[0]->par()->getP(), electrons[0]->cal(PCAL)->getEnergy()/electrons[0]->par()->getP() );
	  
	//	chi2_p->Fill(protons[0]->par()->getChi2Pid() );
	//	chi2_e->Fill(electrons[0]->par()->getChi2Pid() );
	
	// set the particle momentum
	SetLorentzVector(el,electrons[0]);
	SetLorentzVector(pr,protons[0]);
	
	TLorentzVector miss=beam+target-el-pr; //missing 4-vector
       	TLorentzVector q = beam - el;          //photon  4-vector

	miss_p = miss.P();
	miss_m = miss.M();
	miss_pz = miss.Pz();
	q2 = -q.M2();
	x_prime = -q.M2()/(2 * miss.Dot(q) * (beam.E() - el.E()) ); //x-borken              
	x_b = -q.M2()/(2 * mass_p * (beam.E() - el.E()) ); //x-borken                        
	theta = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg();  //angle between vectors p_miss and q                                                                                  
	p_q = pr.Vect().Mag()/q.Vect().Mag();

	point_band = pointsToBand(miss.Vect().Theta(), miss.Vect().Phi(), electrons[0]->par()->getVz());
	p_fd = (protons[0]->getRegion()==FD && ftof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta()) );
	p_cd = (protons[0]->getRegion()==CD && ctof_p ->IsInside(protons[0]->par()->getP(),protons[0]->par()->getBeta() ) );
	//	ecal = (electrons[0]->getRegion()==FD && elec_cut ->IsInside(electrons[0]->par()->getP(),electrons[0]->cal(PCAL)->getEnergy()/electrons[0]->par()->getP()) );

	ecal = (electrons[0]->getRegion()==FD && (electrons[0]->par()->getP() < 4 && electrons[0]->par()->getP() > 1) && ( energy/electrons[0]->par()->getP() < .275 && energy/electrons[0]->par()->getP() > .175) );

	vz_diff = electrons[0]->par()->getVz() - protons[0]->par()->getVz();
	vz = electrons[0]->par()->getVz();
	double vz_diff_p = electrons[0]->par()->getVz() + protons[0]->par()->getVz() + 5.5;
	vz_corr->Fill(vz_diff);

	tree->Fill();

	//	if( (pcut_fd || pcut_cd) && ecut && x_b > 1.2 && abs(vz_diff) < 3 && (vz_diff_p > -10 && vz_diff_p < 10) )
	if( p_fd && ecal && point_band && miss.Pz() < 0   )
	  {
	    hmiss->Fill(miss.M() - mass_n);
	    missm_pmiss->Fill(miss.M() - mass_n,miss.P());
	    missm_xb->Fill(miss.M() - mass_n,x_b);
	    q_dist->Fill(-q.M2());
	    vz_el_p->Fill(protons[0]->par()->getVz(),electrons[0]->par()->getVz());
	    theta_phi->Fill(miss.Theta()*TMath::RadToDeg(), miss.Phi()*TMath::RadToDeg());
	    theta_pmiss->Fill(miss.Vect().Theta()*TMath::RadToDeg(), miss.P());
	    theta_pq_ratio->Fill(p_q,theta);
	  }

	//Double_t eTime=electrons[0]->sci(FTOF1A)->getTime();
      }
    
       
      counter++;
    }

    //}
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");

  theta_pq_ratio->GetYaxis()->SetTitle("#theta_{pq}");
  theta_pq_ratio->GetXaxis()->SetTitle("|p|/|q|");
  missm_pmiss->GetYaxis()->SetTitle("Missing mass (GeV)");
  missm_pmiss->GetXaxis()->SetTitle("P_{miss} (GeV)");
  pid_e_ec->GetYaxis()->SetTitle("E_{loss}/p");
  pid_e_ec->GetXaxis()->SetTitle("p (GeV)");				
  pid_p_fd->GetYaxis()->SetTitle("FD TOF");
  pid_p_fd->GetXaxis()->SetTitle("p_{proton} (GeV)");				
  pid_p_cd->GetYaxis()->SetTitle("CD TOF");
  pid_p_cd->GetXaxis()->SetTitle("p_{proton} (GeV)");				

  theta_pq_ratio->GetYaxis()->CenterTitle();
  theta_pq_ratio->GetXaxis()->CenterTitle();
  missm_pmiss->GetYaxis()->CenterTitle();
  missm_pmiss->GetXaxis()->CenterTitle();
  pid_e_ec->GetYaxis()->CenterTitle();
  pid_e_ec->GetXaxis()->CenterTitle();
  pid_p_fd->GetYaxis()->CenterTitle();
  pid_p_fd->GetXaxis()->CenterTitle();
  pid_p_cd->GetYaxis()->CenterTitle();
  pid_p_cd->GetXaxis()->CenterTitle();

  TCanvas* can = new TCanvas();
  can->Divide(2,2);
  can->cd(1);
  missm_pmiss->Draw("colz");
  can->cd(2);
  missm_xb->Draw("colz");
  can->cd(3);
  //  theta_pq_ratio->Draw("colz");
  q_dist->Draw();

  TCanvas* can1 = new TCanvas();
  theta_phi->Draw("colz");
  //  theta_pmiss->Draw("colz");
  //  pid_e_ec->Draw("colz");

  tree_file->cd();
  tree->Write();

  TFile *out = new TFile("/work/clas12/users/esteejus/output/LowEnergyMC/hist_"+outputFile+".root","RECREATE");
  hmiss->Write();
  pid_e_ec->Write();
  pid_p_fd->Write();
  pid_p_cd->Write();
  theta_pq_ratio->Write();
  missm_pmiss->Write();
  missm_xb->Write();
  q_dist->Write();
  vz_el_p->Write();
  vz_corr->Write();
  hbeam_e->Write();
  theta_pmiss->Write();
  theta_phi->Write();

  /*
  can1->Divide(2,2);
  can1->cd(1);
  pid_p_fd->Draw("colz");
  can1->cd(2);
  pid_p_cd->Draw("colz");
  can1->cd(3);
  pid_e_ec->Draw("colz");
  */
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
