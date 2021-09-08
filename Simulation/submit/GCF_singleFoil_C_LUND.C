#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include <fstream>
#include <iostream>

using namespace std;


//void LUND_tagged(int fileNum,TString prefix);

void LUND_tagged(int iFile, TString prefix);


void GCF_singleFoil_C_LUND(int nFiles = 300, TString prefix = "") {

  for(int i = 1; i < nFiles+1; i++) 
    {
      LUND_tagged(i,prefix);
      cout << i <<endl;
    }
}


void LUND_tagged(int iFile, TString prefix) {
//void LUND_genQESuite(int iFile, TString prefix) {

	TString genPath = "../../rootfiles";
	TString inFileName = Form("%s/%s_%d.root",genPath.Data(),prefix.Data(),iFile);
	cout << "Converting file " << inFileName << endl;
	TFile* inFile = new TFile(inFileName);
	TString lundPath = "../../lundfiles";
	TString outfilename = Form("%s/lund_qe_%s_%d.dat",lundPath.Data(),prefix.Data(),iFile);
	cout << "Making LUND file " << Form("%s/lund_qe_%s_%d.dat",lundPath.Data(),prefix.Data(),iFile) <<endl;

	TTree* T = (TTree*)inFile->Get("genT");

	double targetMin = 2.5;	   //vz min [cm]
	double targetMax = 2.500001; //vz max [cm]
	double rasterX = 0.04;	   // cm
	double rasterY = 0.04; 	   // cm

	const int nParticles = 3;
	const int nMom = 3;
	double mass_e = 0.511e-3;
	double mass_p = 0.938272;	
	double mass_n = 0.93957;

	int targA = 12;
	int targZ = 6;
	double targP = 0.; // polarization
	double beamP = 0.; // polarization
	int interactN = 1;
	int beamType = 11;

	double beamE = -99;	// GeV
	TRandom3 *myRand = new TRandom3(0);

	int lead_type, rec_type;
	double pLead[nMom];
	double pRec[nMom];
	double pe[nMom];
	double weight;

	T->SetBranchAddress("pLead", &pLead);
	T->SetBranchAddress("pRec", &pRec);
	T->SetBranchAddress("pe", &pe);
	T->SetBranchAddress("lead_type", &lead_type);
	T->SetBranchAddress("rec_type", &rec_type);
	T->SetBranchAddress("weight", &weight);


	int nEvents = T->GetEntries();
	cout<<"Number of events "<<nEvents<<endl;

	ofstream outfile;
	outfile.open(outfilename); 

	TString formatstring, outstring;

	for (int i = 0; i<nEvents; i++)
	  {
	    T->GetEntry(i);
	    
	    //	    weight = 1./weight;
	    // LUND header for the event:
	    formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
	    outstring = Form(formatstring, nParticles, targA, targZ, targP, beamP, beamType, beamE, interactN, i, weight);
	    
	    outfile << outstring; 
	    
	    double vx = myRand->Uniform(-rasterX/2., rasterX/2.);	
	    double vy = myRand->Uniform(-rasterY/2., rasterY/2.);
	    double vz = myRand->Uniform(targetMin, targetMax);	
	    
	    // LUND info for each particle in the event
		  
	    formatstring = "%i \t %.3f \t %i \t %i \t %i \t %i \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \n";

	    double mass_lead = 0;
	    double mass_rec = 0;
	    
	    if(lead_type == 2212)
	      mass_lead = mass_p;
	    else if(lead_type == 2112)
	      mass_lead = mass_n;
	    
	    if(rec_type == 2212)
	      mass_rec = mass_p;
	    else if(rec_type == 2112)
	      mass_rec = mass_n;
	    
	    double px = 0;
	    double py = 0;
	    double pz = 0;
	    double energy = 0;
	    
	    px = pe[0];
	    py = pe[1];
	    pz = pe[2];
	    energy = sqrt(px*px + py*py + pz*pz + mass_e*mass_e );			
	    outstring = Form(formatstring, 1, 0.0, 1, 11, 0, 0, px, py, pz, energy, mass_e, vx, vy, vz);  
	    outfile << outstring;


	    px = pLead[0];
	    py = pLead[1];
	    pz = pLead[2];
	    energy = sqrt(px*px + py*py + pz*pz + mass_lead*mass_lead );			
	    outstring = Form(formatstring, 2, 0.0, 1, lead_type, 0, 0, px, py, pz, energy, mass_lead, vx, vy, vz);  
	    outfile << outstring;
	    
	    px = pRec[0];
	    py = pRec[1];
	    pz = pRec[2];
	    energy = sqrt(px*px + py*py + pz*pz + mass_rec*mass_rec );			
	    outstring = Form(formatstring, 3, 0.0, 1, rec_type, 0, 0, px, py, pz, energy, mass_rec, vx, vy, vz);
	    outfile << outstring;
	    
	  }
	
	outfile.close();
	
}
