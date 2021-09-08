#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"
using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void CLAS12SkimmerTree(TString inFile = "", TString outputFile = "", double beamE = 0){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  Double_t mD = 1.8756;


  Double_t q2 = 0;
  Double_t w2 = 0;
  Double_t x_b = 0;
  Double_t x_prime = 0;

  Double_t theta_pq = 0;
  Double_t p_q = 0;

  Double_t miss_m = 0;
  Double_t miss_p_l = 0; //lead
  Double_t miss_p_r = 0; //recoil 

  Double_t vz_e_p = 0;
  Double_t vz_e = 0;

  Int_t p_lead_det = -1; // -1 - No Det 0- CD   1 - FD
  Int_t p_recoil_det = -1; // -1- No Det 0- CD   1 - FD

  Int_t n_lead = 0;
  Double_t p_lead_cd_chi2;
  Double_t p_lead_fd_chi2;

  Double_t p_recoil_cd_chi2;
  Double_t p_recoil_fd_chi2;

  Bool_t ecal = false;



  /////////////////////////////////////
  TFile *tree_file = new TFile(outputFile+".root","RECREATE");

  TTree *tree = new TTree("tree","Low Energy Data");
  tree->Branch("q2",&q2,"q2/D");
  tree->Branch("xb",&x_b,"xb/D");
  tree->Branch("w2",&w2,"w2/D");
  tree->Branch("x_prime",&x_prime,"x_prime/D");

  tree->Branch("theta_pq",&theta_pq,"theta_pq/D");
  tree->Branch("p_q",&p_q,"p_q/D");

  tree->Branch("miss_m",&miss_m,"miss_m/D");
  tree->Branch("miss_p_l",&miss_p_l,"miss_p_l/D");
  tree->Branch("miss_p_r",&miss_p_r,"miss_p_r/D");

  tree->Branch("vz_e_p",&vz_e_p,"vz_e_p/D");
  tree->Branch("vz_e",&vz_e,"vz_e/D");

  tree->Branch("ecal",&ecal,"ecal/B");

  tree->Branch("n_lead",&n_lead,"n_lead/I");
  tree->Branch("p_lead_det",&p_lead_det,"p_lead_det/B");
  tree->Branch("p_recoil_det",&p_recoil_det,"p_recoil_det/B");

  tree->Branch("p_lead_cd_chi2",&p_lead_cd_chi2,"p_lead_cd_chi2/D");
  tree->Branch("p_lead_fd_chi2",&p_lead_fd_chi2,"p_lead_fd_chi2/D");

  tree->Branch("p_recoil_cd_chi2",&p_recoil_cd_chi2,"p_recoil_cd_chi2/D");
  tree->Branch("p_recoil_fd_chi2",&p_recoil_fd_chi2,"p_recoil_fd_chi2/D");


  //initialising clas12writer with path to output file
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");

  clas12root::HipoChain chain;
  chain.Add(inFile);
  chain.SetReaderTags({0});
  auto config_c12=chain.GetC12Reader();

  auto& c12=chain.C12ref();
  auto& rcdbData= config_c12->rcdb()->current();//struct with all relevent rcdb values        
  auto beam_energy = rcdbData.beam_energy/1000;
  cout << "Beam energy is " << beam_energy << endl;

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();

  //some particles
  TLorentzVector beam(0,0,beam_energy, beam_energy);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   
  gBenchmark->Start("timer");
 
  int counter = 0;
  while(chain.Next())
    {

      beam_energy = rcdbData.beam_energy/1000;
      beam.SetE(beam_energy);
      beam.SetPz(beam_energy);

      if( beamE > 1e-4)
	{
	  if( counter == 0)
	    cout<<"Beam energy set manually: "<<beamE<<endl;
	  beam.SetE(beamE);
	  beam.SetPz(beamE);
	}

    q2 = 0;
    w2 = 0.;
    x_b = 0;
    x_prime = 0;

    theta_pq = 0;
    p_q = 0;

    miss_m = 0;
    miss_p_l = 0;
    miss_p_r = 0;

    vz_e_p = 0;
    vz_e = 0;

    ecal = false;
    p_lead_det = -1;
    p_recoil_det = -1;

    n_lead = 0;

    p_lead_cd_chi2 = 0.;
    p_lead_fd_chi2 = 0.;
    p_recoil_cd_chi2 = 0.;
    p_recoil_fd_chi2 = 0.;


    // get particles by type
    auto electrons=c12->getByID(11);
    auto protons=c12->getByID(2212);
       

    if(electrons.size()==1 && protons.size()==2 )
      {

	double energy =  electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy();
	
       ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 1 && electrons[0]->getP() < 10 );


       TLorentzVector q = beam - el;                  //photon  4-vector                     
       q2       = -q.M2();
       x_b      = q2/(2 * mass_p * (beam.E() - el.E()) ); //x-borken 
       
       vz_e   = electrons[0]->par()->getVz();
       vz_e_p = electrons[0]->par()->getVz() - protons[0]->par()->getVz();
	   


       int idx_lead = -1;
       for(int iPr = 0; iPr < 2; iPr++)
	 {
	   // set the particle momentum
	   SetLorentzVector(el,electrons[0]);
	   SetLorentzVector(pr,protons[iPr]);
	   
	   TLorentzVector miss = beam + target - el - pr; //missing 4-vector

	   //angle between vectors p_miss and q                                         
	   theta_pq = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); 
	   p_q      = pr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|

	   if( theta_pq < 25 && p_q < 0.96 && p_q > 0.62)
	     {
	       idx_lead = iPr;
	       n_lead++;
	     }
	   
	 }
	   
       if(idx_lead != -1)
	 {
	   if(protons[idx_lead]->getRegion()==CD)
	     {
	       p_lead_cd_chi2 = protons[idx_lead]->par()->getChi2Pid();
	       p_lead_det = 0;
	     }
	   if(protons[idx_lead]->getRegion()==FD)
	     {
	       p_lead_fd_chi2 = protons[idx_lead]->par()->getChi2Pid();
	       p_lead_det = 1;
	     }

	   //Determine recoil particle
	   if(idx_lead == 0)
	     {
	       if(protons[1]->getRegion()==CD)
		 {
		   p_recoil_cd_chi2 = protons[1]->par()->getChi2Pid();
		   p_recoil_det = 0;
		 }	       

	       if(protons[1]->getRegion()==FD)
		 {
		   p_recoil_fd_chi2 = protons[1]->par()->getChi2Pid();
		   p_recoil_det = 1;
		 }
	     }
	   else
	     {
	       if(protons[0]->getRegion()==CD)
		 {
		   p_recoil_cd_chi2 = protons[0]->par()->getChi2Pid();
		   p_recoil_det = 0;
		 }	       
	       if(protons[0]->getRegion()==FD)
		 {
		   p_recoil_fd_chi2 = protons[0]->par()->getChi2Pid();
		   p_recoil_det = 1;
		 }
	     }
	   
	   // set the particle momentum
	   SetLorentzVector(el,electrons[0]);
	   SetLorentzVector(pr,protons[idx_lead]);
	   
	   TLorentzVector miss = beam + target - el - pr; //missing 4-vector
	   miss_p_l = miss.P();
	   miss_m   = miss.M2();
	   
	   theta_pq = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); //angle between vectors p_miss and q                                                                              
	   p_q      = pr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|
	   x_prime  = q2/(2 * miss.Dot(q) * (beam.E() - el.E()) ); //x-borken prime             
	 }

       tree->Fill();
      }

    counter++;

    }
  
  tree_file->cd();
  tree->Write();
  tree_file->Close();
  
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " read events = "<<counter<<endl;
  }
