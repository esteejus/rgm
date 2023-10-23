#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include <fstream>
#include <iostream>

using namespace std;

TString addParticle(int part_idx,int pid, TVector3 momentum, double mass, TVector3 vtx)
{
  // LUND info for each particle in the event
  TString formatstring = "%i \t %.3f \t %i \t %i \t %i \t %i \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \n";
  double energy = sqrt(momentum.Mag2() + mass*mass );			
  TString outstring = Form(formatstring, part_idx, 0.0, 1, pid, 0, 0, momentum.Px(), momentum.Py(), momentum.Pz(), energy, mass, vtx.X(), vtx.Y(), vtx.Z());  

  return outstring;
}

void generate_protons(TString outputFile = "", int nEvents = 10000, double theta_min = 35, double theta_max = 135)
{
  TRandom3 ran(0);
  const int nParticles = 2;
  double mass_e = 0.511e-3;
  double mass_p = 0.938272;	
  double mass_n = 0.93957;
  
  //Not useful usually
  double targP = 0.; // polarization
  double beamP = 0.; // polarization
  int interactN = 1;
  int beamType = 11;
  
  double beamE = -99;	// GeV
  
  double weight = 1;
  ofstream outfile;
  outfile.open(outputFile);
  if(!outfile.is_open())
    cout<<"Output file cannot be created"<<endl;

  TString formatstring, outstring;

  for (int i = 0; i < nEvents; i++)
    {
      
      // LUND header for the event:
      formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
      outstring = Form(formatstring, nParticles, 1, 1, targP, beamP, beamType, beamE, interactN, i, weight);
      outfile << outstring; 
      
      //      double theta_min = 35;
      //      double theta_max = 140;

      //Add particles in event below
      TVector3 vtx(0,0,-3);

      //electron
      double mom = 6.;
      double phi = 0;
      double theta = 25 * TMath::RadToDeg();
      double px = mom * sin(theta)*cos(phi); 
      double py = mom * sin(theta)*sin(phi); 
      double pz = mom * cos(theta);

      outfile << addParticle(1,11,TVector3(px,py,pz),mass_e,vtx);

      mom = ran.Uniform(.2,1);
      phi = ran.Uniform(-TMath::Pi(),TMath::Pi());
      theta = acos(2*ran.Uniform() - 1); //isotropic

      while( !(theta*TMath::RadToDeg() > theta_min && theta*TMath::RadToDeg() < theta_max) )
	theta = acos(2*ran.Uniform() - 1);

      px = mom * sin(theta)*cos(phi); 
      py = mom * sin(theta)*sin(phi); 
      pz = mom * cos(theta);
      outfile << addParticle(2,2212,TVector3(px,py,pz),mass_p,vtx);
      
    }
  
  outfile.close();

}
