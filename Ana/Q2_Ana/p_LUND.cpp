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
TVector3 radiateElectron(TVector3 ve);
double radiationFactor(double Ebeam, double Ek, double QSq);
void generate_event(double &weight, TVector3& ve, TVector3 &vp);
double sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton);
double Gkelly(double QSq,double a1, double b1, double b2, double b3);
double deltaHard(double QSq);

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
  double W, weight;
  root_outfile->cd();
  TTree * outtree = new TTree("genTbuffer","Generator Tree");
  outtree->Branch("pe",pe,"pe[3]/D");
  outtree->Branch("pp",pp,"pp[3]/D");
  outtree->Branch("W",&W,"W/D");
  outtree->Branch("weight",&weight,"weight/D");


  for (int i = 0; i < nEvents; i++)
    {
      TVector3 ve;
      TVector3 vp;
      generate_event(weight,ve,vp);
      pe[0]=ve.X();
      pe[1]=ve.Y();
      pe[2]=ve.Z();
      pp[0]=vp.X();
      pp[1]=vp.Y();
      pp[2]=vp.Z();
      TVector3 q = init_vBeam - ve;
      double omega = init_vBeam.Mag() - ve.Mag();
      double Q2 = q.Mag2() - omega*omega;
      double WSq = (mN*mN) - Q2 + (2*omega*mN);
      W = sqrt(WSq);
      if(weight!=0){
	outtree->Fill();

	// LUND header for the event:
	formatstring = "%i \t %i \t %i \t %.3f \t %.3f \t %i \t %.1f \t %i \t %i \t %.3f \n";
	outstring = Form(formatstring, nParticles, 1, 1, targP, beamP, beamType, init_Ebeam, interactN, i, weight);
	outfile << outstring; 
	
	auto vtx = randomVertex("liquid"); //get vertex of event 
	//electron
	outfile << addParticle(1,1,11,ve,mass_e,vtx);
	//p
	outfile << addParticle(2,1,2212,vp,mN,vtx);
      }
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

TVector3 radiateElectron(TVector3 ve)
{
  double Ee = ve.Mag();
  double lambda_e = alpha/M_PI*(log(4*(Ee*Ee)/(mass_e*mass_e)) -1);
  double DeltaE = pow(myRand->Rndm(),1./lambda_e)*Ee;
  double Ee_rad = Ee - DeltaE;
  TVector3 ve_rad = ve;
  ve_rad.SetMag(Ee_rad);
  return ve_rad;
}

double radiationFactor(double Ebeam, double Ek, double QSq)
{
  double lambda_ei = alpha/M_PI*(log(4*(Ebeam*Ebeam)/(mass_e*mass_e)) -1);
  double lambda_ef = alpha/M_PI*(log(4*(Ek*Ek)/(mass_e*mass_e)) -1);

  return (1 - deltaHard(QSq)) * pow(Ebeam/sqrt(Ebeam*Ek),lambda_ei) * pow(Ek/sqrt(Ebeam*Ek),lambda_ef);
}

void generate_event(double &weight, TVector3& ve, TVector3 &vp){

  TVector3 vBeam(0,0,init_Ebeam);
  //First radiation
  vBeam = radiateElectron(vBeam);
  double Ebeam = vBeam.Mag();

  //Random Sample
  double Q2 = (myRand->Uniform()*5)+0.5;
  double Q2_max = 4*Ebeam*Ebeam/(1+2*Ebeam/mN);
  double phie = (myRand->Uniform()*2*M_PI)-M_PI;
  
  //Calculate Kinematics
  double omega = Q2 / (2 * mN);
  double Eprime = Ebeam - omega;
  double mome = Eprime;
  double thetae = 2 * asin(sqrt( Q2 / (4*Ebeam*Eprime) ));
  ve.SetMagThetaPhi(mome,thetae,phie);
  TVector3 init_ve = ve;
  vp = vBeam - ve;

  //Second Radiation
  ve = radiateElectron(ve);

  //Calculate weight
  //Values at interaction
  if(Q2>Q2_max){
    weight=0;
  }
  else{
    weight = radiationFactor(Ebeam, Eprime, Q2);
    weight*= sigma_onShell_by_Etheta(Ebeam,init_ve,true);
  }

  //Measured proton mass
  double Ep = Ebeam + mN - ve.Mag();
  double m = sqrt(Ep*Ep - vp.Mag2());
  
  if(weight!=0){
    if(!(m>0)){      
      cout<<"missed\n";
      weight=0;
    }
  }
}


double sigma_onShell_by_Etheta(double Ebeam, TVector3 k, bool isProton)
{
  double theta=k.Theta();
  double E3 = Ebeam * mN/ (mN + Ebeam*(1.-k.CosTheta()));
  double QSq = 2. * Ebeam * E3 * (1.-k.CosTheta());
  double tau = QSq/(4.*mN*mN);
  double GE = Gkelly(QSq,-0.24,10.98,12.82,21.97);//isProton ? GEp(QSq) : GEn(QSq);
  double GM = mu_p * Gkelly(QSq,0.12,10.97,18.86,6.55);//isProton ? GMp(QSq) : GMn(QSq);
  double epsilon = 1./(1.+2.*(1.+tau)*(tan(theta/2.)*tan(theta/2.)));

  double sigmaMott = nbGeVSq * (2.*alpha*E3 * cos(theta/2.)/QSq) * (2.*alpha*E3 * cos(theta/2.)/QSq) * (E3/Ebeam);

  double sigmaRosenbluth = sigmaMott * ((GE*GE) + tau/epsilon * (GM*GM))/(1. + tau);
  return sigmaRosenbluth * Ebeam / (E3 * (2.*tau + 1.));
}

double Gkelly(double QSq,double a1, double b1, double b2, double b3)
{
  double tau = QSq/(4.*mN*mN);
  double denom = 1. + b1*tau + b2*tau*tau + b3*tau*tau*tau;
  double numer = 1. + a1*tau;
  return numer/denom;
}

double deltaHard(double QSq)
{
  return 2.*alpha/M_PI * ( -13./12.*log(QSq/(mass_e*mass_e)) + 8./3.);
}


