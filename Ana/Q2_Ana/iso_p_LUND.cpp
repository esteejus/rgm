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
#include <TRandom3.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

//Target parameters
const double mass_e = 0.511e-3;
const double mass_p = 0.938272;
const double mass_n = 0.93957;
const double mass_pi = 0.13957;
const double mu_p = 2.79;
const double GeVfm  =0.1973;
const double alpha = 0.0072973525664;
const double cmSqGeVSq = GeVfm*GeVfm*1.E-26;
const double nbGeVSq = cmSqGeVSq*1.E33;
TRandom3 * myRand = new TRandom3(0);

double beamspot_x = 0.;
double beamspot_y = 0.;

double beamspread_x = 0.04;
double beamspread_y = 0.04;

double init_Ebeam = 5.98636;// GeV
TVector3 init_vBeam(0,0,init_Ebeam);
double mN = 0.938;

//Targets: liquid, 4-foil, 1-foil, Ar, Ca
double global_z = -3;  // center of hallB in GEMC in cm

std::map <std::string, std::vector<double> > targets {
  {"4-foil", {-1.875 + global_z,-0.625 + global_z,0.625 + global_z, 1.875 + global_z}},
  {"1-foil", {global_z + 2.5}},
  {"Ar",     {global_z - 2.5}},
  {"Ca",     {global_z}},
  {"liquid", {global_z - 2.5, global_z + 2.5}}
};

TVector3 randomVertex(string target);
TString addParticle(int part_idx, int active, int pid, TVector3 momentum, double mass, TVector3 vtx);
void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp);
void Usage();
void generate_event(TVector3& ve, TVector3 &vp);

int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }

  int nEvents = atoi(argv[1]);
  ofstream outfile;
  outfile.open(argv[2]); 
  TFile * root_outfile = new TFile(argv[3],"RECREATE");

  cout<<"Number of events "<<nEvents<<endl;
  cout << "Making LUND file " << argv[2] <<endl;
  cout << "Making ROOT file " << argv[3] <<endl;

  const int nParticles = 2;
  const int nMom = 3;

  //Not useful usually
  double targP = 0.; // polarization
  double beamP = 0.; // polarization
  int interactN = 1;
  int beamType = 11;
      
  TString formatstring, outstring;

  if(!outfile.is_open()){
    cout<<"Output file cannot be created"<<endl;}

  double pe[3], pp[3];
  root_outfile->cd();
  TTree * outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pp",pp,"pp[3]/D");


  for (int i = 0; i < nEvents; i++)
    {
      TVector3 ve;
      TVector3 vp;
      generate_event(ve,vp);
      pe[0]=ve.X();
      pe[1]=ve.Y();
      pe[2]=ve.Z();
      pp[0]=vp.X();
      pp[1]=vp.Y();
      pp[2]=vp.Z();
      outtree->Fill();

      // LUND header for the event:
      formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
      outstring = Form(formatstring, nParticles, 1, 1, targP, beamP, beamType, init_Ebeam, interactN, i, 1);
      outfile << outstring; 
      
      auto vtx = randomVertex("liquid"); //get vertex of event 
      //electron
      outfile << addParticle(1,1,11,ve,mass_e,vtx);
      //p
      outfile << addParticle(2,1,2212,vp,mN,vtx);
    }
  
  outfile.close();
  outtree->SetName("genT");
  outtree->Write();
  root_outfile->Delete("genTbuffer;*");
  root_outfile->Close();

  return 0;
}

TVector3 randomVertex(string target)
{
  double x=-999,y=-999,z=-999;

  if(targets.find(target) != targets.end())
    {

      if(target == "liquid")
	{
	  x = myRand->Gaus(beamspot_x,beamspread_x);
	  y = myRand->Gaus(beamspot_y,beamspread_y);
	  z = myRand->Uniform(targets[target].at(0),targets[target].at(1));
	  return (TVector3(x,y,z));
	}
      else
	{
	  // foil targets                                                                     
	  x = myRand->Gaus(beamspot_x,beamspread_x);
	  y = myRand->Gaus(beamspot_y,beamspread_y);
	  z = targets[target].at(myRand->Integer(targets[target].size()) ); 
	  return (TVector3(x,y,z));
	}

    }

  cout<<"Error in target initialization. Target not found in map. Vertex will be set to -999,-999,-999 "<<endl;

  return (TVector3(-999,-999,-999));
}


TString addParticle(int part_idx, int active, int pid, TVector3 momentum, double mass, TVector3 vtx)
{
  // LUND info for each particle in the event                                                                                                                                                 
  TString formatstring = "%i \t %.3f \t %i \t %i \t %i \t %i \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \t %.5f \n";
  double energy = sqrt(momentum.Mag2() + mass*mass );
  TString outstring = Form(formatstring, part_idx, 0.0, active, pid, 0, 0, momentum.Px(), momentum.Py(), momentum.Pz(), energy, mass, vtx.X(), vtx.Y(), vtx.Z());

  return outstring;
}


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}

void Usage()
{
  std::cerr << "Usage: ./code nEvents outputfile \n\n\n";

}

void generate_event(TVector3& ve, TVector3 &vp){

  double costhetae = (myRand->Uniform()*0.3)+0.7;
  double thetae = acos(costhetae);
  double phie = (myRand->Uniform()*2*M_PI)-M_PI;
  double mome = (myRand->Uniform()*5)+1;
  ve.SetMagThetaPhi(mome,thetae,phie);

  double costhetap = (myRand->Uniform()*2)-1.0;
  double thetap = acos(costhetap);
  double phip = (myRand->Uniform()*2*M_PI)-M_PI;
  double momp = (myRand->Uniform()*3.3)+0.2;
  vp.SetMagThetaPhi(momp,thetap,phip);
}

