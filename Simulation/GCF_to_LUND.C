#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include <fstream>
#include <iostream>
#include "targets.h"
using namespace std;


void GCF_to_LUND(TString inputFile = "", TString outputFile = "", string target = "liquid", int A = 1, int Z = 1)
{
  //Read in target parameter files
  cout << "Converting file " << inputFile << endl;
  TFile* inFile = new TFile(inputFile);
  cout << "Making LUND file " << outputFile <<endl;
  
  TTree* T = (TTree*)inFile->Get("genT");
  
  const int nParticles = 3;
  const int nMom = 3;

  //Not useful usually
  double targP = 0.; // polarization
  double beamP = 0.; // polarization
  int interactN = 1;
  int beamType = 11;
  
  double beamE = -99;	// GeV
  
  int lead_pid, rec_pid;
  double pLead[nMom];
  double pRec[nMom];
  double pe[nMom];
  double weight;
  
  T->SetBranchAddress("pLead", &pLead);
  T->SetBranchAddress("pRec", &pRec);
  T->SetBranchAddress("pe", &pe);
  T->SetBranchAddress("lead_type", &lead_pid);
  T->SetBranchAddress("rec_type", &rec_pid);
  T->SetBranchAddress("weight", &weight);
  
  int nEvents = T->GetEntries();
  cout<<"Number of events "<<nEvents<<endl;
  
  ofstream outfile;
  outfile.open(outputFile); 
  if(!outfile.is_open())
    cout<<"Output file cannot be created"<<endl;

  TString formatstring, outstring;
  
  for (int i = 0; i < nEvents; i++)
    {
      T->GetEntry(i);
      
      // LUND header for the event:
      formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
      outstring = Form(formatstring, nParticles, A, Z, targP, beamP, beamType, beamE, interactN, i, weight);
      outfile << outstring; 

      //Add particles in event below
      double mass_lead = 0;
      double mass_rec = 0;
      
      if(lead_pid == 2212)
	  mass_lead = mass_p;
      else if(lead_pid == 2112)
	mass_lead = mass_n;

      if(rec_pid == 2212)
	mass_rec = mass_p;
      else if(rec_pid == 2112)
	mass_rec = mass_n;
      
      auto vtx = randomVertex(target); //get vertex of event 
      //electron
      outfile << addParticle(1,11,TVector3(pe[0],pe[1],pe[2]),mass_e,vtx);
      //lead SRC
      outfile << addParticle(2,lead_pid,TVector3(pLead[0],pLead[1],pLead[2]),mass_lead,vtx);
      //rec SRC
      outfile << addParticle(3,rec_pid,TVector3(pRec[0],pRec[1],pRec[2]),mass_rec,vtx);
    }
  
  outfile.close();
  //  gSystem->Exec(".q");

}
