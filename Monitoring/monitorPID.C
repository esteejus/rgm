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

void monitorPID(){
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
  outputFile = inputFile(inputFile.Index("inc_")+4,inputFile.Index(".hipo"));

  cout<<"Analysing hipo file "<<inputFile<<endl;


  //=======Quality Track Cuts============
  auto *pid_e_ec = new TH2D("pid_p_ec","EC/p vs e electron",100,0,8,100,0,.4);
  auto *pid_p_fd = new TH2D("pid_p_fd","FD TOF vs p PID proton",500,0,4,500,0,1.3);
  auto *pid_p_cd = new TH2D("pid_p_cd","CD TOF vs p PID proton",200,0,4,200,0,1.3);

  auto *chi2_p_cd = new TH1D("chi2_p_cd","Chi2 CD proton",100,-10,10);
  auto *chi2_p_fd = new TH1D("chi2_p_fd","Chi2 FD proton",100,-10,10);


  auto *vz_corr = new TH1D("vz_corr","e_vz - p_vz",1000,-10,10);
  auto *vz_el_p = new TH2D("vz_el_p","vz_el_p",1000,-30,30,1000,-30,30);
  auto *vz_el = new TH1D("vz_el","vz_el",1000,-10,4);
  //====================================


  //=====(e,e')N_l N_r===========


  //===========



  gBenchmark->Start("timer");
  int counter = 0;
  
  clas12reader c12(inputFile.Data());

    while(c12.next()==true){
    
      for(auto& p : c12.getDetParticles()){
	// get particles by type
      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      auto neutrons=c12.getByID(2112);

      if(electrons.size()==1 && protons.size()==1)
	{

	  double sampling_frac =  electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy();

	  if(protons[0]->getRegion()==FD)
	    {
	      pid_p_fd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );
	      chi2_p_fd->Fill(protons[0]->par()->getChi2Pid() );
	    }	  
	  if(protons[0]->getRegion()==CD)
	    {
	    pid_p_cd -> Fill(protons[0]->par()->getP(),protons[0]->par()->getBeta() );
	    chi2_p_cd->Fill(protons[0]->par()->getChi2Pid() );
	    }	  

	  if(electrons[0]->getRegion()==FD)
	    pid_e_ec -> Fill(electrons[0]->par()->getP(), sampling_frac/electrons[0]->par()->getP());
	
	  

	  double vz_diff = electrons[0]->par()->getVz() - protons[0]->par()->getVz();
	  double vz = electrons[0]->par()->getVz(); 
	  double vz_diff_p = electrons[0]->par()->getVz() + protons[0]->par()->getVz() + 5.5;

	  vz_el->Fill(vz);
	  vz_corr->Fill(vz_diff);

	}

      }
      
       
      counter++;
    }
    
    
    gBenchmark->Stop("timer");
    gBenchmark->Print("timer");
    
    
    TCanvas* cvs = new TCanvas("cvs","cvs",1200,1200);
    cvs->Divide(3,3);

    cvs->cd(1);
    pid_p_fd->Draw("colz");

    cvs->cd(2);
    pid_p_cd->Draw("colz");

    cvs->cd(3);
    pid_e_ec->Draw("colz");

    cvs->cd(4);
    chi2_p_fd->Draw();

    cvs->cd(5);
    chi2_p_cd->Draw();

    cvs->cd(6);
    vz_el->Draw();

    cvs->cd(7);
    vz_corr->Draw();


    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";

}
