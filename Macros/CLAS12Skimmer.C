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

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void CLAS12Skimmer(TString inFile = "", TString outputFile = "",double beamE = 0){

  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

   
  cout<<"Analysing hipo file "<<inFile<<endl;

  TChain fake("hipo");
  fake.Add(inFile.Data());
  auto files=fake.GetListOfFiles();

  //initialising clas12writer with path to output file
  clas12writer c12writer(outputFile.Data());
  
  //can as writer not to write certain banks
  //c12writer.skipBank("REC::Cherenkov");
  //c12writer.skipBank("REC::Scintillator");

  cout<<"Analysing hipo file "<<inFile<<endl;

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;

  //some particles
  TLorentzVector beam(0,0,0,0);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
   
  gBenchmark->Start("timer");
 
  int counter = 0;
  int writeCounter = 0;
   

  //create the event reader

  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");


  for(Int_t i=0;i<files->GetEntries();i++)
    {
      
      clas12reader c12(files->At(i)->GetTitle(),{0});
      
      clas12databases rcdb;
      c12.connectDataBases(&rcdb);
      auto& rcdbData= c12.rcdb()->current();//using standalone clas12reader object
      
      
      //assign a reader to the writer
      c12writer.assignReader(c12);
      
      
      while(c12.next()){
	
	if(counter == 0)
	  cout<<"Beam energy: "<<rcdbData.beam_energy<<endl;                                      
	beam.SetE(rcdbData.beam_energy/1000);                                           
	beam.SetPz(rcdbData.beam_energy/1000);     
	
	if(beamE > 1e-4)
	  {
	    if(counter == 0)
	      cout<<"Beam energy manually set: "<<beamE<<endl;
	    
	    beam.SetE(beamE);
	    beam.SetPz(beamE);
	  }
	
	
	// get particles by type
	auto electrons=c12.getByID(11);
	auto protons=c12.getByID(2212);
	auto pips=c12.getByID(211);
	auto pims=c12.getByID(-211);
    
    
	if(electrons.size()==1 && protons.size() >= 1)
	  {
	    double energy_sf =  (electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy()) / electrons[0]->getP();
	    
	    bool electron_cut = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy_sf > 0.18 && energy_sf < 0.28) && electrons[0]->getP() > 1 && electrons[0]->getP() < 10 );
	    
	    // set the particle momentum
	    SetLorentzVector(el,electrons[0]);
	    SetLorentzVector(pr,protons[0]);
	    
	    TLorentzVector miss = beam + target - el - pr; //missing 4-vector
	    TLorentzVector q = beam - el;                  //photon  4-vector                     
	    
	    double q2       = -q.M2();
	    double x_prime  = q2/(2 * miss.Dot(q) * (beam.E() - el.E()) ); //x-borken prime             
	    double x_b      = q2/(2 * mass_p * (beam.E() - el.E()) ); //x-borken 
	    double theta_pq = pr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); //angle between vectors p_miss and q                                                                              
	    double p_q      = pr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|
	    
	    
	    if( electron_cut && x_b > 1.2 && q2 > 1.5)
	      {
		//write out an event
		c12writer.writeEvent(); 
		writeCounter++;
	      }
	  }
	
	counter++;

      }
    }
  //close writer
  c12writer.closeWriter();
  
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " read events = "<<counter<<" wrote events = "<<writeCounter<<" s\n";
  
  gROOT->ProcessLine(".q");

}
