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
#include <TGraph.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "many_plots.h"
#include "reweighter.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

double params_Theta_FD[6][4][4]={{{-0.191151,-0.669629,10,1.00023},
			    {0.113741,-0.0599251,19.0221,31.6623},
			    {350.715,460.172,17.1609,1},
			    {3.32022,2.7311,29.6054,25.1635}},
			   {{-0.0891744,-0.14685,7.50659,25.2638},
			    {0.0450352,-0.109521,6.0011,19.0888},
			    {350.357,140.077,6,15.8616},
			    {-3.25284,-1.24789,10.6854,24.956}},
			   {{-0.0849098,-0.0651043,11.7406,19.6928},
			    {0.0209555,-0.0273879,6.00308,20.1268},
			    {98.1241,20.8179,6,29.528},
			    {-5.06788,-4.9607,20.1145,1.00211}},
			   {{-0.0739777,0.0963433,6.00048,17.6548},
			    {0.014477,-0.156077,9.29576,12.1747},
			    {341.693,154.122,4.66104,22.2281},
			    {-6.38361,-1.95148,6,29.7998}},
			   {{-0.0528101,-0.100172,6.00982,12.4738},
			    {0.010489,-0.0339853,6.09838,16.0642},
			    {348.057,161.215,6,21.4269},
			    {-4.07125,-10,6,2.57348}},
			   {{-0.0623622,-0.179411,7.77701,9.66338},
			    {0.040617,-0.185013,6.00574,9.96181},
			    {199.345,74.073,6,28.6864},
			    {-3.83103,-1.41067,6,11.1274}}};
  
double params_Theta_CD[3][2]={{-0.19764,8.14947},
			      {-0.0207092,4.52151},
			      {-1.45427,0}};

double params_EnergyLoss_FD[3][3]={{-0.000695124,-0.000355869,0},
				   {0.00182181,7.77933e-05,0},
				   {0.000266541,0.424055,49.068}};

double params_EnergyLoss_CD[3][2]={{-0.00555986,-6.06201e-05},
				   {0.00695634,8.24535e-05},
				   {0.00155103,1.74283e-05}};

double params_Momentum_FD[6][4][4]={{{0.00999685,-0.0274693,97.8319,27.1649},
				{0.107308,0.0482255,6.00016,26.5809},
				{126.732,40.9353,4,15.7632},
				{-1.29807,-2.08763,5.6061,21.9216}},
			       {{0.0863814,0.0662239,78.8143,22.7513},
				{0.0819542,0.0353392,4.05495,26.0594},
				{123.313,39.0737,4,17.3371},
				{-1.23367,-2.10145,4.37322,21.3527}},
			       {{0.136301,0.0624967,4.06521,39.5042},
				{0.0573779,0.0258018,4.99407,24.1826},
				{119.322,74.2621,7.26117,10},
				{-4.36463,-2.59454,7.47078,36.4093}},
			       {{0.0273149,-0.009257,20.2085,24.9454},
				{0.0689756,0.0280117,4.01384,30.0473},
				{128.494,53.5087,4.35928,17.9939},
				{-1.59758,-1.5448,4.0001,29.9412}},
			       {{0.114001,0.0792582,37.7721,36.9068},
				{0.075125,0.0282635,4.00108,29.5526},
				{138.358,59.4015,4,10},
				{-0.996918,-1.03077,4.00001,27.4224}},
			       {{0.108582,0.0801866,79.5092,28.1964},
				{0.068627,0.0341774,16.8015,27.416},
				{130.348,48.6141,4,17.9232},
				{-0.926971,-0.686175,4,20.191}}};

double params_Momentum_CD[3][2]={{-0.107267,0.00193308},
				 {0.0267032,0.000545869},
				 {-2.31836,0.0177897}};

double Function_Erf(double x, double A, double B, double C, double D){
  return A - B*(1+erf(((-x+D)/C))); 
}

double Function_Trig(double x, double A, double B, double C, double D){
  return A + B*sin((x*2*M_PI/C)+D); 
}

double Function_TrigFixedPeriod(double x, double A, double B, double C){
  return A + B*sin((x*2*M_PI/180)+C); 
}

double Function_Algebraic1(double x, double A, double B){
  return A + B*x; 
}

double Function_Algebraic2(double x, double A, double B){
  return A + B/(x); 
}

double Function_Algebraic3(double x, double A, double B, double C){
  return A + B*x + (C/x); 
}

double Function_Algebraic4(double x, double A, double B, double C){
  return A - B/(x-C); 
}

void GetLorentzVector_ReconVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void SetLorentzVector_ThetaCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Theta_FD[sector-1][i][0],params_Theta_FD[sector-1][i][1],params_Theta_FD[sector-1][i][2],params_Theta_FD[sector-1][i][3]);
    }
    theta+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic2(theta,params_Theta_CD[i][0],params_Theta_CD[i][1]);
    }
    theta+=Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);        
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetTheta(theta*M_PI/180);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

void SetLorentzVector_EnergyLossCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    double params0 = Function_Algebraic1(theta,params_EnergyLoss_FD[0][0],params_EnergyLoss_FD[0][1]);
    double params1 = Function_Algebraic1(theta,params_EnergyLoss_FD[1][0],params_EnergyLoss_FD[1][1]);
    double params2 = Function_Algebraic4(theta,params_EnergyLoss_FD[2][0],params_EnergyLoss_FD[2][1],params_EnergyLoss_FD[2][2]);
    mom+=Function_Algebraic3(mom,params0,params1,params2);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic1(theta,params_EnergyLoss_CD[i][0],params_EnergyLoss_CD[i][1]);
    }
    mom+=Function_Algebraic3(mom,params[0],params[1],params[2]);    
      }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
} 

void SetLorentzVector_MomentumCorrection(TLorentzVector &p4, clas12::region_particle* p){
  TVector3 v3 = p4.Vect();
  double mom = v3.Mag();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Momentum_FD[sector-1][i][0],params_Momentum_FD[sector-1][i][1],params_Momentum_FD[sector-1][i][2],params_Momentum_FD[sector-1][i][3]);
    }    
    mom+=Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
    for(int i = 0; i < 3; i++){
      params[i] = Function_Algebraic1(theta,params_Momentum_CD[i][0],params_Momentum_CD[i][1]);
    }    
    mom+=Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);
  }
  else{
    cout<<"Problem\n\n\n\n\n\n\n\n\n";
    exit(-2);
  }
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
}

double GetThetaCorrected(clas12::region_particle* p){
  TVector3 v(p->par()->getPx(),p->par()->getPy(),p->par()->getPz());
  double theta = v.Theta()*180/M_PI;
  double phi = v.Phi()*180/M_PI;
  if(p->getRegion()==FD){
    int sector = p->getSector();
    double shift = 0;
    shift += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi>0)?180:(sector==4 && phi<0)?-180:(sector==5)?-120:(sector==6)?-60:0;
    phi-=shift;

    double params[4];
    for(int i = 0; i < 4; i++){
      params[i] = Function_Erf(theta,params_Theta_FD[sector-1][i][0],params_Theta_FD[sector-1][i][1],params_Theta_FD[sector-1][i][2],params_Theta_FD[sector-1][i][3]);
    }
    return Function_Trig(phi,params[0],params[1],params[2],params[3]);
  }
  else if(p->getRegion()==CD){
    double params[3];
      for(int i = 0; i < 3; i++){
	params[i] = Function_Algebraic2(theta,params_Theta_CD[i][0],params_Theta_CD[i][1]);
      }
    return Function_TrigFixedPeriod(phi,params[0],params[1],params[2]);    
  }
  cout<<"Problem\n\n\n\n\n\n\n\n\n";
  return -1000.00;
}
////////////////////////
///////////////////////
int tBinCD(double t)
{
  int b = (t - 35) / 5;
  if(b < 0){return -1;}
  if(b > 9){return -1;}
  return b;
}

int tBinFD(double t)
{
  int b = (t - 20) / 5;
  if(b < 0){return -1;}
  if(b > 4){return -1;}
  return b;
}

int tBinFDe(double t)
{
  int b = (t - 10) / 3;
  if(b < 0){return -1;}
  if(b > 8){return -1;}
  return b;
}

vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,45};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD = 1.8756;

//double beam_E = 5.98636;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double cotan(double x){ return cos(x)/sin(x);}

double E(double x, double N, double tau){
  return N * exp( x / tau) ; 
}

double ThetaE(double ThetaP){
  double ThetaP_rad = ThetaP*M_PI/180;
  return (180/M_PI) * 2*atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaP_rad) );
}


double ThetaP(double ThetaE){
  double ThetaE_rad = ThetaE*M_PI/180;
  return (180/M_PI) * atan( (1/((beam_E/mass_p)+1)) * cotan(ThetaE_rad/2) );
}

double SQ(double x){ return x*x;}

double ClosestPoint(double x, double ThetaEp, double ThetaPp){
  return sqrt(SQ(x-ThetaEp) + SQ(ThetaP(x)-ThetaPp));
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph){
  int ctr = 0;
  char temp[100];
  //Now project the histogram    
  for(int j = 0; j < h_myhist->GetXaxis()->GetNbins(); j++){
    //Define x and y(1D histogram)
    double x = h_myhist->GetXaxis()->GetBinCenter(j+1);
    ctr++;
    sprintf(temp,"Proj_num%d",ctr);
    TH1D * proj = h_myhist->ProjectionY(temp,j+1,j+1);
    //proj->Rebin(2);
    //Now preform a guassian fit
    if(proj->GetEntries()<50){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-1,1,3);
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,0.0);
    gFit->SetParLimits(1,-0.8,0.8);
    gFit->SetParameter(2,0.2);
    gFit->SetParLimits(2,0.001,1.0);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-1,1);
    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    proj->Write();
  }
}

void getFunction(TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint){  
  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-1,1);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,-1,1);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-1,1);
  f_myfunc->SetParameter(3,0);
  f_myfunc->SetParLimits(3,-1,1);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",-40,40);
      
}


void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  int isMC = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TVector3 vbeam(0,0,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector el_corrected(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr_corrected(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector lead_pion_ptr(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector stationary_proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;
  TH1D * h_Delta_Int_Eprime = new TH1D("Delta_Eprime","Delta_Eprime",100,-0.5,0.5);
  TH1D * h_Delta_Int_Eprime_Corrected = new TH1D("Delta_Eprime_Corrected","Delta_Eprime_Corrected",100,-0.5,0.5);
  TH1D * h_Delta_Int_pMomFD = new TH1D("Delta_pMomFD","Delta_pMomFD",100,-0.5,0.5);
  TH1D * h_Delta_Int_pMomFD_Corrected = new TH1D("Delta_pMomFD_Corrected","Delta_pMomFD_Corrected",100,-0.5,0.5);
  TH1D * h_Delta_Int_pMomCD = new TH1D("Delta_pMomCD","Delta_pMomCD",100,-0.5,0.5);
  TH1D * h_Delta_Int_pMomCD_Corrected = new TH1D("Delta_pMomCD_Corrected","Delta_pMomCD_Corrected",100,-0.5,0.5);


  
  TH1D * h_Delta_Eprime[9];
  TH1D * h_Delta_Eprime_Corrected[9];
  TH1D * h_Delta_Eprime_sectors[6][9];

  TH1D * h_Q2[9];
  TH1D * h_p_e[9];
  TH1D * h_theta_q[9];
  TH1D * h_phi_diff[9];
  TH1D * h_phi_diff_FD[9];
  TH1D * h_phi_diff_CD[9];
  TH1D * h_theta_p[9];
  int bins[9] = {100,100,100,100,100,100,50,50,50};
  for(int i=0; i<9; i++){
    int min = 10+(3*i);
    int max = 13+(3*i);

    sprintf(temp_name,"Delta_Eprime_%d",min);
    sprintf(temp_title,"#Delta E' (%d< #theta_{e} < %d);#Delta E';Counts",min,max);
    h_Delta_Eprime[i] = new TH1D(temp_name,temp_title,bins[i],-0.5,0.5);
    hist_list.push_back(h_Delta_Eprime[i]);

    sprintf(temp_name,"Delta_Eprime_Corrected_%d",min);
    sprintf(temp_title,"#Delta E' Corrected (%d< #theta_{e} < %d);#Delta E' Corrected;Counts",min,max);
    h_Delta_Eprime_Corrected[i] = new TH1D(temp_name,temp_title,bins[i],-0.5,0.5);
    hist_list.push_back(h_Delta_Eprime_Corrected[i]);

    for(int j=1; j<=6; j++){
      sprintf(temp_name,"Delta_Eprime_%d_sector_%d",min,j);
      sprintf(temp_title,"#Delta E' (%d< #theta_{e} < %d) Sector %d;#Delta E';Counts",min,max,j);
      h_Delta_Eprime_sectors[j-1][i] = new TH1D(temp_name,temp_title,bins[i],-0.5,0.5);
      hist_list.push_back(h_Delta_Eprime_sectors[j-1][i]);
    }
    
    sprintf(temp_name,"Q2_%d",min);
    sprintf(temp_title,"Q^{2} (%d< #theta_{e} < %d);Q^{2};Counts",min,max);
    h_Q2[i] = new TH1D(temp_name,temp_title,100,0,5);
    hist_list.push_back(h_Q2[i]);

    sprintf(temp_name,"p_e_%d",min);
    sprintf(temp_title,"p_{e} (%d< #theta_{e} < %d);p_{e};Counts",min,max);
    h_p_e[i] = new TH1D(temp_name,temp_title,100,0,6);
    hist_list.push_back(h_p_e[i]);

    sprintf(temp_name,"theta_q_%d",min);
    sprintf(temp_title,"#theta_{q} (%d< #theta_{e} < %d);#theta_{q};Counts",min,max);
    h_theta_q[i] = new TH1D(temp_name,temp_title,100,0,100);
    hist_list.push_back(h_theta_q[i]);


    sprintf(temp_name,"phi_diff_%d",min);
    sprintf(temp_title,"#Delta #phi (e,e'p) (%d< #theta_{e} < %d);#Delta #phi;Counts",min,max);
    h_phi_diff[i] = new TH1D(temp_name,temp_title,100,-10,10);
    hist_list.push_back(h_phi_diff[i]);

    sprintf(temp_name,"phi_diff_FD_%d",min);
    sprintf(temp_title,"#Delta #phi (e,e'p_{FD}) (%d< #theta_{e} < %d);#Delta #phi;Counts",min,max);
    h_phi_diff_FD[i] = new TH1D(temp_name,temp_title,100,-10,10);
    hist_list.push_back(h_phi_diff_FD[i]);

    sprintf(temp_name,"phi_diff_CD_%d",min);
    sprintf(temp_title,"#Delta #phi (e,e'p_{CD}) (%d< #theta_{e} < %d);#Delta #phi;Counts",min,max);
    h_phi_diff_CD[i] = new TH1D(temp_name,temp_title,100,-10,10);
    hist_list.push_back(h_phi_diff_CD[i]);

    sprintf(temp_name,"theta_p_%d",min);
    sprintf(temp_title,"#theta_{p} (%d< #theta_{e} < %d);#theta_{p};Counts",min,max);
    h_theta_p[i] = new TH1D(temp_name,temp_title,100,0,100);
    hist_list.push_back(h_theta_p[i]);
  }

  TH1D * h_EbeamRat = new TH1D("EbeamRat","E_{0}/E_{beam};E_{0}/E_{beam};Counts",100,0.9,1.1);
  hist_list.push_back(h_EbeamRat);
  TH1D * h_EbeamRatFD = new TH1D("EbeamRatFD","E_{0}/E_{beam} FD;E_{0}/E_{beam};Counts",100,0.9,1.1);
  hist_list.push_back(h_EbeamRatFD);
  TH1D * h_EbeamRatCD = new TH1D("EbeamRatCD","E_{0}/E_{beam} CD;E_{0}/E_{beam};Counts",100,0.9,1.1);
  hist_list.push_back(h_EbeamRatCD);
  TH1D * h_EbeamRatCorrFD = new TH1D("EbeamRatCorrFD","E_{0}/E_{beam} Corrected FD;E_{0}/E_{beam};Counts",100,0.9,1.1);
  hist_list.push_back(h_EbeamRatCorrFD);
  TH1D * h_EbeamRatCorrCD = new TH1D("EbeamRatCorrCD","E_{0}/E_{beam} Corrected CD;E_{0}/E_{beam};Counts",100,0.9,1.1);
  hist_list.push_back(h_EbeamRatCorrCD);


  
  TH2D * h_phie_EbeamRat_FD = new TH2D("phie_EbeamRat_FD","E_{0}/E_{beam} vs #phi_{e};#phi_{e};E_{0}/E_{beam};Counts",180,-180,180,100,0.9,1.1);
  hist_list.push_back(h_phie_EbeamRat_FD);
  TH2D * h_phie_EbeamRat_CD = new TH2D("phie_EbeamRat_CD","E_{0}/E_{beam} vs #phi_{e};#phi_{e};E_{0}/E_{beam};Counts",180,-180,180,100,0.9,1.1);
  hist_list.push_back(h_phie_EbeamRat_CD);
  TH2D * h_phip_EbeamRat_FD = new TH2D("phip_EbeamRat_FD","E_{0}/E_{beam} vs #phi_{p};#phi_{p};E_{0}/E_{beam};Counts",180,-180,180,100,0.9,1.1);
  hist_list.push_back(h_phip_EbeamRat_FD);
  TH2D * h_phip_EbeamRat_CD = new TH2D("phip_EbeamRat_CD","E_{0}/E_{beam} vs #phi_{p};#phi_{p};E_{0}/E_{beam};Counts",180,-180,180,100,0.9,1.1);
  hist_list.push_back(h_phip_EbeamRat_CD);

  
  TH1D * h_phie_FD = new TH1D("phie_FD","Counts vs #phi_{e};#phi_{e};Counts",180,-180,180);
  hist_list.push_back(h_phie_FD);	 
  TH1D * h_phie_CD = new TH1D("phie_CD","Counts vs #phi_{e};#phi_{e};Counts",180,-180,180);
  hist_list.push_back(h_phie_CD);	 
  TH1D * h_phip_FD = new TH1D("phip_FD","Counts vs #phi_{p};#phi_{p};Counts",180,-180,180);
  hist_list.push_back(h_phip_FD);	 
  TH1D * h_phip_CD = new TH1D("phip_CD","Counts vs #phi_{p};#phi_{p};Counts",180,-180,180);
  hist_list.push_back(h_phip_CD);


  TH2D * h_thetae_thetap_FD = new TH2D("thetae_thetap_FD","#theta_{p} vs. #theta_{e};#theta_{e};#theta_{p};Counts",100,15,40,100,20,45);
  hist_list.push_back(h_thetae_thetap_FD);
  TH2D * h_thetae_thetap_CD = new TH2D("thetae_thetap_CD","#theta_{p} vs. #theta_{e};#theta_{e};#theta_{p};Counts",100,5,25,100,35,65);
  hist_list.push_back(h_thetae_thetap_CD);


  TH2D * h_theta_corr_eFD_binSector[6];
  TH2D * h_theta_corr_eCD_binSector[6];
  TH2D * h_theta_corr_eCD_FixP_binSector[6];
  TH2D * h_theta_corr_pFD_binSector[6];
  TH2D * h_phi_corr_eFD_binSector[6];
  TH2D * h_phi_corr_eCD_binSector[6];
  TH2D * h_phi_corr_eCD_FixP_binSector[6];
  TH2D * h_phi_corr_pFD_binSector[6];

  TGraphErrors * g_theta_corr_eFD_binSector[6];
  TGraphErrors * g_theta_corr_eCD_binSector[6];
  TGraphErrors * g_theta_corr_eCD_FixP_binSector[6];
  TGraphErrors * g_theta_corr_pFD_binSector[6];
  TGraphErrors * g_phi_corr_eFD_binSector[6];
  TGraphErrors * g_phi_corr_eCD_binSector[6];
  TGraphErrors * g_phi_corr_eCD_FixP_binSector[6];
  TGraphErrors * g_phi_corr_pFD_binSector[6];
  for(int j=1; j<=6; j++){
    sprintf(temp_name,"theta_corr_eFD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #theta eFD Sector %d;#theta;Correction;Counts",j);
    h_theta_corr_eFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,5,45,100,-2.0,2.0);
    hist_list.push_back(h_theta_corr_eFD_binSector[j-1]);
    g_theta_corr_eFD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_theta_corr_eFD_sector_%d",j);
    g_theta_corr_eFD_binSector[j-1]->SetName(temp_name);
      
    sprintf(temp_name,"theta_corr_eCD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #theta eCD Sector %d;#theta;Correction;Counts",j);
    h_theta_corr_eCD_binSector[j-1] = new TH2D(temp_name,temp_title,100,5,45,100,-2.0,2.0);
    hist_list.push_back(h_theta_corr_eCD_binSector[j-1]);
    g_theta_corr_eCD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_theta_corr_eCD_sector_%d",j);
    g_theta_corr_eCD_binSector[j-1]->SetName(temp_name);

    sprintf(temp_name,"theta_corr_eCD_FixP_sector_%d",j);
    sprintf(temp_title,"Correction vs. #theta eCD_FixP Sector %d;#theta;Correction;Counts",j);
    h_theta_corr_eCD_FixP_binSector[j-1] = new TH2D(temp_name,temp_title,100,5,45,100,-2.0,2.0);
    hist_list.push_back(h_theta_corr_eCD_FixP_binSector[j-1]);
    g_theta_corr_eCD_FixP_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_theta_corr_eCD_FixP_sector_%d",j);
    g_theta_corr_eCD_FixP_binSector[j-1]->SetName(temp_name);

    sprintf(temp_name,"theta_corr_pFD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #theta pFD Sector %d;#theta;Correction;Counts",j);
    h_theta_corr_pFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,5,45,100,-2.0,2.0);
    hist_list.push_back(h_theta_corr_pFD_binSector[j-1]);
    g_theta_corr_pFD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_theta_corr_pFD_sector_%d",j);
    g_theta_corr_pFD_binSector[j-1]->SetName(temp_name);

    sprintf(temp_name,"phi_corr_eFD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #phi eFD Sector %d;#phi;Correction;Counts",j);
    h_phi_corr_eFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-35,35,100,-2.0,2.0);
    hist_list.push_back(h_phi_corr_eFD_binSector[j-1]);
    g_phi_corr_eFD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_eFD_sector_%d",j);
    g_phi_corr_eFD_binSector[j-1]->SetName(temp_name);
      
    sprintf(temp_name,"phi_corr_eCD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #phi eCD Sector %d;#phi;Correction;Counts",j);
    h_phi_corr_eCD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-35,35,100,-2.0,2.0);
    hist_list.push_back(h_phi_corr_eCD_binSector[j-1]);
    g_phi_corr_eCD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_eCD_sector_%d",j);
    g_phi_corr_eCD_binSector[j-1]->SetName(temp_name);

    sprintf(temp_name,"phi_corr_eCD_FixP_sector_%d",j);
    sprintf(temp_title,"Correction vs. #phi eCD_FixP Sector %d;#phi;Correction;Counts",j);
    h_phi_corr_eCD_FixP_binSector[j-1] = new TH2D(temp_name,temp_title,100,-35,35,100,-2.0,2.0);
    hist_list.push_back(h_phi_corr_eCD_FixP_binSector[j-1]);
    g_phi_corr_eCD_FixP_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_eCD_FixP_sector_%d",j);
    g_phi_corr_eCD_FixP_binSector[j-1]->SetName(temp_name);

    sprintf(temp_name,"phi_corr_pFD_sector_%d",j);
    sprintf(temp_title,"Correction vs. #phi pFD Sector %d;#phi;Correction;Counts",j);
    h_phi_corr_pFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-35,35,100,-2.0,2.0);
    hist_list.push_back(h_phi_corr_pFD_binSector[j-1]);
    g_phi_corr_pFD_binSector[j-1] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_pFD_sector_%d",j);
    g_phi_corr_pFD_binSector[j-1]->SetName(temp_name);

  }



  TH2D * h_phi_theta_eFD_binSector[6];
  TH2D * h_phi_theta_eCD_binSector[6];
  TH2D * h_phi_theta_pFD_binSector[6];
  TH2D * h_phi_theta_COMB_binSector[6];
  for(int j=1; j<=6; j++){
    sprintf(temp_name,"phi_theta_eFD_sector_%d",j);
    sprintf(temp_title,"#theta vs. #phi eFD Sector %d;#phi;#theta;Counts",j);
    h_phi_theta_eFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-40,40,100,5,45);
    hist_list.push_back(h_phi_theta_eFD_binSector[j-1]);
    
    sprintf(temp_name,"phi_theta_eCD_sector_%d",j);
    sprintf(temp_title,"#theta vs. #phi eCD Sector %d;#phi;#theta;Counts",j);
    h_phi_theta_eCD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-40,40,100,5,45);
    hist_list.push_back(h_phi_theta_eCD_binSector[j-1]);

    sprintf(temp_name,"phi_theta_pFD_sector_%d",j);
    sprintf(temp_title,"#theta vs. #phi pFD Sector %d;#phi;#theta;Counts",j);
    h_phi_theta_pFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,-40,40,100,5,45);
    hist_list.push_back(h_phi_theta_pFD_binSector[j-1]);

    sprintf(temp_name,"phi_theta_COMB_sector_%d",j);
    sprintf(temp_title,"#theta vs. #phi COMB Sector %d;#phi;#theta;Counts",j);
    h_phi_theta_COMB_binSector[j-1] = new TH2D(temp_name,temp_title,100,-40,40,100,5,45);
    hist_list.push_back(h_phi_theta_COMB_binSector[j-1]);
  }


  
  TH2D * h_theta_corr_binSector_binPhi[6][8];
  TGraphErrors * g_theta_corr_binSector_binPhi[6][8];
  for(int j=1; j<=6; j++){
    for(int i=0; i<8; i++){
      int min = bE_Phi[i];
      int max = bE_Phi[i+1];
      sprintf(temp_name,"theta_corr_sector_%d_phi_%d",j,i);
      sprintf(temp_title,"Correction vs. #theta Sector %d (%d< #phi < %d);#theta;Correction;Counts",j,min,max);
      h_theta_corr_binSector_binPhi[j-1][i] = new TH2D(temp_name,temp_title,40,5,45,100,-1.5,1.5);
      hist_list.push_back(h_theta_corr_binSector_binPhi[j-1][i]);

      g_theta_corr_binSector_binPhi[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_theta_corr_sector_%d_phi_%d",j,i);
      g_theta_corr_binSector_binPhi[j-1][i]->SetName(temp_name);

    }
  }

  TH2D * h_phi_corr_binSector_binTheta[6][12];
  TGraphErrors * g_phi_corr_binSector_binTheta[6][12];
  for(int j=1; j<=6; j++){
    for(int i=0; i<12; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction;Counts",j,min,max);
      h_phi_corr_binSector_binTheta[j-1][i] = new TH2D(temp_name,temp_title,40,-40,40,100,-1.5,1.5);
      hist_list.push_back(h_phi_corr_binSector_binTheta[j-1][i]);

      g_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_phi_corr_sector_%d_theta_%d",j,i);
      g_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

    }
  }

  TH2D * h_theta_corr_pCD = new TH2D("h_theta_corr_eFD","Correction vs. #theta pCD;#theta;Correction;Counts",100,30,90,100,-2.0,2.0);
  hist_list.push_back(h_theta_corr_pCD);
  TGraphErrors * g_theta_corr_pCD = new TGraphErrors(); 
  g_theta_corr_pCD->SetName("g_theta_corr_eFD");

  TH2D * h_phi_corr_pCD = new TH2D("h_phi_corr_eFD","Correction vs. #phi pCD;#phi;Correction;Counts",100,-180,180,100,-2.0,2.0);
  hist_list.push_back(h_phi_corr_pCD);
  TGraphErrors * g_phi_corr_pCD = new TGraphErrors(); 
  g_phi_corr_pCD->SetName("g_phi_corr_eFD");
  
  TH2D * h_phi_theta_pCD = new TH2D("h_phi_theta_eFD","#theta vs. #phi pCD;#phi;#theta;Counts",100,-180,180,100,30,90);
  hist_list.push_back(h_phi_theta_pCD);
  TGraphErrors * g_phi_theta_pCD = new TGraphErrors(); 
  g_phi_theta_pCD->SetName("g_phi_theta_eFD");

  TH2D * h_phi_corr_binThetaCD[6];
  TGraphErrors * g_phi_corr_binThetaCD[6];
  for(int i=0; i<6; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"phi_corr_theta_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< #theta < %d);#phi;Correction;Counts",min,max);
    h_phi_corr_binThetaCD[i] = new TH2D(temp_name,temp_title,180,-180,180,100,-1.5,1.5);
    hist_list.push_back(h_phi_corr_binThetaCD[i]);
    
    g_phi_corr_binThetaCD[i] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_theta_%d",i);
    g_phi_corr_binThetaCD[i]->SetName(temp_name);
  }
  
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  int counter = 0;
  //while(chain.Next())
  while(chain.Next() && (counter<100000000))
    {
      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight;
      }
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto particles = c12->getDetParticles(); //particles is now 
      if(electrons.size() == 1 && protons.size() >= 0)
	{

	  GetLorentzVector_ReconVector(el,electrons[0]);
	  TLorentzVector el_corrected = el;
	  SetLorentzVector_ThetaCorrection(el_corrected,electrons[0]);
	  SetLorentzVector_MomentumCorrection(el_corrected,electrons[0]);


	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vel_fromAngle;
	  vel_fromAngle.SetMagThetaPhi(Eprime,el.Theta(),el.Phi());
	  TVector3 vp_fromAngle = vbeam-vel_fromAngle;
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();


	  double Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN );
	  TVector3 vel_fromAngle_corrected;
	  vel_fromAngle_corrected.SetMagThetaPhi(Eprime_corrected,el_corrected.Theta(),el_corrected.Phi());
	  TVector3 vp_fromAngle_corrected = vbeam-vel_fromAngle_corrected;
	  double Delta_Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN ) - el_corrected.P();

		  
	  int sector = electrons[0]->getSector();
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
	  double omega_diff = omega - (Q2/(2*mN));

	  double vtz_e = electrons[0]->par()->getVz();
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double phi_e = el.Phi() * 180/M_PI;
	  double theta_e = el.Theta()*180/M_PI;
	  double theta_e_corrected = el_corrected.Theta()*180/M_PI;
	  double theta_q = q.Theta()*180/M_PI;
	  
	  if(tBinFDe(theta_e_corrected)!=-1){
	    h_Delta_Eprime_Corrected[tBinFDe(theta_e_corrected)]->Fill(Delta_Eprime_corrected,wep);
	  }
	  if(tBinFDe(theta_e)!=-1){
	    h_Delta_Eprime[tBinFDe(theta_e)]->Fill(Delta_Eprime,wep);
	    h_Delta_Eprime_sectors[sector-1][tBinFDe(theta_e)]->Fill(Delta_Eprime,wep);
	    if(Delta_Eprime<0.1){
	      h_Q2[tBinFDe(theta_e)]->Fill(Q2,wep);
	      h_p_e[tBinFDe(theta_e)]->Fill(el.P(),wep);
	      h_theta_q[tBinFDe(theta_e)]->Fill(theta_q,wep);
	    }
	  }
	  if(protons.size() <= 0){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  TLorentzVector proton_ptr_corrected = proton_ptr;
	  SetLorentzVector_ThetaCorrection(proton_ptr_corrected,protons[0]);
	  SetLorentzVector_EnergyLossCorrection(proton_ptr_corrected,protons[0]);
	  SetLorentzVector_MomentumCorrection(proton_ptr_corrected,protons[0]);


	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  
	  double mom_p = proton_ptr.P();
	  double mom_p_corrected = proton_ptr_corrected.P();
	  double Delta_Mom_FromAngle = vp_fromAngle.Mag()-mom_p;
	  double Delta_Mom_FromAngle_corrected = vp_fromAngle_corrected.Mag()-mom_p_corrected;

	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  
	  //double mom_t = beta * gamma * mN;
	  
	  double Delta_mom = (mom_p-mom_q)/mom_q;
	  //cout<<Delta_mom<<endl;
	  double Delta_theta = theta_p-theta_q;
	  double Delta_phi = (phi_p-phi_q);
	  if(Delta_phi<-180){Delta_phi+=360;}
	  else if(Delta_phi>180){Delta_phi-=360;}

	  if(fabs(Delta_Eprime)>0.15){continue;}
	  
	  if(tBinFDe(theta_e)!=-1){
	    h_phi_diff[tBinFDe(theta_e)]->Fill(Delta_phi,wep);
	    h_theta_p[tBinFDe(theta_e)]->Fill(theta_p,wep);
	    if(protons[0]->getRegion()==FD){
	      h_phi_diff_FD[tBinFDe(theta_e)]->Fill(Delta_phi,wep);
	    }
	    else if(protons[0]->getRegion()==CD){
	      h_phi_diff_CD[tBinFDe(theta_e)]->Fill(Delta_phi,wep);
	    }

	    if(fabs(Delta_phi)>3){continue;}
	    double theta_correction_e = (theta_e+GetThetaCorrected(electrons[0]))*M_PI/180;
	    double theta_correction_p = (theta_p+GetThetaCorrected(protons[0]))*M_PI/180;
	    double theta_e = theta_correction_e*180/M_PI;
	    double theta_p = theta_correction_p*180/M_PI;
	    //double E_Beam_Calc = mN*(cotan(1)*cotan(1) - 1);
	    double E_Beam_Calc = mN*(cotan((theta_e*M_PI/180)/2)*cotan(theta_p*M_PI/180) - 1);
	    TF1 * fDist = new TF1("fDist",[&](double *x, double *p){ return ClosestPoint(x[0],p[0],p[1]); },1,50,2);
	    fDist->SetParameter(0,theta_e);
	    fDist->SetParameter(1,theta_p);
	    double theta_e_corr = fDist->GetMinimumX(1,50);
	    double theta_p_corr = ThetaP(theta_e_corr);

	    double theta_e_corr_FixP = ThetaE(theta_p);

	    double shift_e = 0;
	    shift_e += (sector==1)?0:(sector==2)?60:(sector==3)?120:(sector==4 && phi_e>0)?180:(sector==4 && phi_e<0)?-180:(sector==5)?-120:(sector==6)?-60:0;

	    double E_Beam_Calc_Corrected = mN*(cotan((theta_correction_e)/2)*cotan(theta_correction_p) - 1);

	    h_EbeamRat->Fill(E_Beam_Calc/beam_E,wep);
	    h_Delta_Int_Eprime->Fill(Delta_Eprime,wep);
	    h_Delta_Int_Eprime_Corrected->Fill(Delta_Eprime_corrected,wep);
	    if(protons[0]->getRegion()==FD){
	      h_EbeamRatFD->Fill(E_Beam_Calc/beam_E,wep);
	      h_EbeamRatCorrFD->Fill(E_Beam_Calc_Corrected/beam_E,wep);
	      h_Delta_Int_pMomFD->Fill(Delta_Mom_FromAngle,wep);
	      h_Delta_Int_pMomFD_Corrected->Fill(Delta_Mom_FromAngle_corrected,wep);
	    }
	    else if(protons[0]->getRegion()==CD){
	      h_EbeamRatCD->Fill(E_Beam_Calc/beam_E,wep);
	      h_EbeamRatCorrCD->Fill(E_Beam_Calc_Corrected/beam_E,wep);
	      h_Delta_Int_pMomCD->Fill(Delta_Mom_FromAngle,wep);
	      h_Delta_Int_pMomCD_Corrected->Fill(Delta_Mom_FromAngle_corrected,wep);
	    }
	    
	    ////////
	    if(binX(bE_Phi,phi_e-shift_e)!=-1){
	      h_theta_corr_binSector_binPhi[sector-1][binX(bE_Phi,phi_e-shift_e)]->Fill(theta_e,theta_e_corr-theta_e,wep);
	    }
	    if(binX(bE_Theta,theta_e)!=-1){
	      h_phi_corr_binSector_binTheta[sector-1][binX(bE_Theta,theta_e)]->Fill(phi_e-shift_e,theta_e_corr-theta_e,wep);
	    }


	    if(protons[0]->getRegion()==FD){
	      int sector_p = protons[0]->getSector();
	      double shift_p = 0;
	      shift_p += (sector_p==1)?0:(sector_p==2)?60:(sector_p==3)?120:(sector_p==4 && phi_p>0)?180:(sector_p==4 && phi_p<0)?-180:(sector_p==5)?-120:(sector_p==6)?-60:0;

	      //////////////
	      if(binX(bE_Phi,phi_p-shift_p)!=-1){
		h_theta_corr_binSector_binPhi[sector_p-1][binX(bE_Phi,phi_p-shift_p)]->Fill(theta_p,theta_p_corr-theta_p,wep);
	      }
	      if(binX(bE_Theta,theta_p)!=-1){
		h_phi_corr_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(phi_p-shift_p,theta_p_corr-theta_p,wep);
	      }

	      
	      h_phie_EbeamRat_FD->Fill(phi_e,E_Beam_Calc/beam_E,wep);
	      h_phip_EbeamRat_FD->Fill(phi_p,E_Beam_Calc/beam_E,wep);
	      h_phie_FD->Fill(phi_e,wep);
	      h_phip_FD->Fill(phi_p,wep);

	      h_thetae_thetap_FD->Fill(theta_e,theta_p,wep);

	      //
	      h_theta_corr_eFD_binSector[sector-1]->Fill(theta_e,theta_e_corr-theta_e,wep);
	      h_theta_corr_pFD_binSector[sector_p-1]->Fill(theta_p,theta_p_corr-theta_p,wep);

	      h_phi_corr_eFD_binSector[sector-1]->Fill(phi_e-shift_e,theta_e_corr-theta_e,wep);
	      h_phi_corr_pFD_binSector[sector_p-1]->Fill(phi_p-shift_p,theta_p_corr-theta_p,wep);
	      
	      h_phi_theta_eFD_binSector[sector-1]->Fill(phi_e-shift_e,theta_e,wep);
	      h_phi_theta_pFD_binSector[sector_p-1]->Fill(phi_p-shift_p,theta_p,wep);

	      h_phi_theta_COMB_binSector[sector-1]->Fill(phi_e-shift_e,theta_e,wep);
	      h_phi_theta_COMB_binSector[sector_p-1]->Fill(phi_p-shift_p,theta_p,wep);

	    }
	    else if(protons[0]->getRegion()==CD){
	      h_phie_EbeamRat_CD->Fill(phi_e,E_Beam_Calc/beam_E,wep);
	      h_phip_EbeamRat_CD->Fill(phi_p,E_Beam_Calc/beam_E,wep);
	      h_phie_CD->Fill(phi_e,wep);
	      h_phip_CD->Fill(phi_p,wep);

	      
	      h_thetae_thetap_CD->Fill(theta_e,theta_p,wep);

	      //
	      h_theta_corr_eCD_binSector[sector-1]->Fill(theta_e,theta_e_corr-theta_e,wep);
	      h_theta_corr_eCD_FixP_binSector[sector-1]->Fill(theta_e,theta_e_corr_FixP-theta_e,wep);

	      h_phi_corr_eCD_binSector[sector-1]->Fill(phi_e-shift_e,theta_e_corr-theta_e,wep);
	      h_phi_corr_eCD_FixP_binSector[sector-1]->Fill(phi_e-shift_e,theta_e_corr_FixP-theta_e,wep);

	      h_phi_theta_eCD_binSector[sector-1]->Fill(phi_e-shift_e,theta_e,wep);

	      h_phi_theta_COMB_binSector[sector-1]->Fill(phi_e-shift_e,theta_e,wep);


	      //
	      h_theta_corr_pCD->Fill(theta_p,theta_p_corr-theta_p,wep);
	      h_phi_corr_pCD->Fill(phi_p,theta_p_corr-theta_p,wep);
	      h_phi_theta_pCD->Fill(phi_p,theta_p,wep);
	      if(binX(bE_ThetaCD,theta_p)!=-1){
		h_phi_corr_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(phi_p,theta_p_corr-theta_p,wep);
	      }

	    }
	    
	  }
	}
    }
  
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  }
  //g_omega_diff_mu->Write();
  for(int j=1; j<=6; j++){
    //g_omega_diff_mu_sectors[j-1]->Write();
  }
  int pixelx = 1980;
  int pixely = 1530;
  TCanvas * myCanvas = new TCanvas("myPage","myPage",pixelx,pixely);
  TCanvas * myText = new TCanvas("myText","myText",pixelx,pixely);
  TLatex text;
  text.SetTextSize(0.05);
  
  char fileName[100];
  sprintf(fileName,"%s[",pdfFile);
  myText->SaveAs(fileName);
  sprintf(fileName,"%s",pdfFile);

  //////////////////////////////////////////////////
  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_Delta_Eprime_Corrected[i-1]->SetLineColor(2);    
    h_Delta_Eprime_Corrected[i-1]->Draw();    
    h_Delta_Eprime[i-1]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);    
    for(int j=1; j<=6; j++){
      h_Delta_Eprime_sectors[j-1][i-1]->SetLineColor(j);    
      h_Delta_Eprime_sectors[j-1][i-1]->Draw("SAME");    
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  
  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_Q2[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_p_e[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_theta_q[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_phi_diff[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_phi_diff_FD[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_phi_diff_CD[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_theta_p[i-1]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_EbeamRat->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_EbeamRatCorrFD->SetLineColor(2);
  h_EbeamRatCorrFD->Draw();
  h_EbeamRatFD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_EbeamRatCorrCD->SetLineColor(2);
  h_EbeamRatCorrCD->Draw();
  h_EbeamRatCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_phie_EbeamRat_FD->Draw("colz");
  myCanvas->cd(2);
  h_phip_EbeamRat_FD->Draw("colz");
  myCanvas->cd(3);
  h_phie_EbeamRat_CD->Draw("colz");
  myCanvas->cd(4);
  h_phip_EbeamRat_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Delta_Int_Eprime_Corrected->SetLineColor(2);
  h_Delta_Int_Eprime_Corrected->Draw();
  h_Delta_Int_Eprime->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Delta_Int_pMomFD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomFD_Corrected->Draw();
  h_Delta_Int_pMomFD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Delta_Int_pMomCD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomCD_Corrected->Draw();
  h_Delta_Int_pMomCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_phie_FD->Draw("colz");
  myCanvas->cd(2);
  h_phip_FD->Draw("colz");
  myCanvas->cd(3);
  h_phie_CD->Draw("colz");
  myCanvas->cd(4);
  h_phip_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_thetae_thetap_FD->Draw("colz");
  myCanvas->cd(2);
  h_thetae_thetap_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_theta_corr_eFD_binSector[i]->Draw("colz");
    getGraph(h_theta_corr_eFD_binSector[i],g_theta_corr_eFD_binSector[i]);
    g_theta_corr_eFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_theta_corr_pFD_binSector[i]->Draw("colz");
    getGraph(h_theta_corr_pFD_binSector[i],g_theta_corr_pFD_binSector[i]);
    g_theta_corr_pFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_theta_corr_eCD_binSector[i]->Draw("colz");
    getGraph(h_theta_corr_eCD_binSector[i],g_theta_corr_eCD_binSector[i]);
    g_theta_corr_eCD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_theta_corr_eCD_FixP_binSector[i]->Draw("colz");
    getGraph(h_theta_corr_eCD_FixP_binSector[i],g_theta_corr_eCD_FixP_binSector[i]);
    g_theta_corr_eCD_FixP_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    g_theta_corr_eCD_binSector[i]->Draw();
    g_theta_corr_eCD_FixP_binSector[i]->SetLineColor(2);
    g_theta_corr_eCD_FixP_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //////////////////
  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_eFD_binSector[i]->Draw("colz");
    getGraph(h_phi_corr_eFD_binSector[i],g_phi_corr_eFD_binSector[i]);
    g_phi_corr_eFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_pFD_binSector[i]->Draw("colz");
    getGraph(h_phi_corr_pFD_binSector[i],g_phi_corr_pFD_binSector[i]);
    g_phi_corr_pFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_eCD_binSector[i]->Draw("colz");
    getGraph(h_phi_corr_eCD_binSector[i],g_phi_corr_eCD_binSector[i]);
    g_phi_corr_eCD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_eCD_FixP_binSector[i]->Draw("colz");
    getGraph(h_phi_corr_eCD_FixP_binSector[i],g_phi_corr_eCD_FixP_binSector[i]);
    g_phi_corr_eCD_FixP_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    g_phi_corr_eCD_binSector[i]->Draw();
    g_phi_corr_eCD_FixP_binSector[i]->SetLineColor(2);
    g_phi_corr_eCD_FixP_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  /////////////////////////////////

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_eFD_binSector[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_pFD_binSector[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_eCD_binSector[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_COMB_binSector[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  ///////////////////////////////////
  for(int j = 0; j < 6; j++){
    myCanvas->Divide(3,3);
    for(int i = 0; i < 8; i++){
      myCanvas->cd(i+1);
      h_theta_corr_binSector_binPhi[j][i]->Draw("colz");
      getGraph(h_theta_corr_binSector_binPhi[j][i],g_theta_corr_binSector_binPhi[j][i]);
      g_theta_corr_binSector_binPhi[j][i]->Draw("SAME");
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }


  ///////////////////////////////////
  for(int j = 0; j < 6; j++){
    myCanvas->Divide(3,4);
    for(int i = 0; i < 12; i++){
      myCanvas->cd(i+1);
      h_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      getGraph(h_phi_corr_binSector_binTheta[j][i],g_phi_corr_binSector_binTheta[j][i]);
      g_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_phi_corr_pCD->Draw("colz");
  getGraph(h_phi_corr_pCD,g_phi_corr_pCD);
  g_phi_corr_pCD->Draw("SAME");
  myCanvas->cd(2);
  h_theta_corr_pCD->Draw("colz");
  getGraph(h_theta_corr_pCD,g_theta_corr_pCD);
  g_theta_corr_pCD->Draw("SAME");
  myCanvas->cd(3);
  h_phi_theta_pCD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    getGraph(h_phi_corr_binThetaCD[i],g_phi_corr_binThetaCD[i]);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}

/*
	  if(Delta_Eprime<0.1){	  
	    if((protons[0]->getRegion()==CD) && (tBinCD(theta_p)!=-1)){
	      h_mom_res_CTOF_bin_thetaq[tBinCD(theta_p)]->Fill(Delta_mom,wep);
	      h_mom_CTOF_bin_thetaq[tBinCD(theta_p)]->Fill(mom_p,wep);
	      //h_momt_res_CTOF_bin_thetaq[tBinCD(theta_p)]->Fill((mom_t-mom_q)/mom_q,wep);
	      h_theta_res_CTOF_bin_thetaq[tBinCD(theta_p)]->Fill(Delta_theta,wep);
	      h_phi_res_CTOF_bin_thetaq[tBinCD(theta_p)]->Fill(Delta_phi,wep);
	    }
	  }
*/

  /*
  TH1D * h_mom_res_CTOF_bin_thetaq[10];
  TH1D * h_mom_CTOF_bin_thetaq[10];
  TH1D * h_momt_res_CTOF_bin_thetaq[10];
  TH1D * h_theta_res_CTOF_bin_thetaq[10];
  TH1D * h_phi_res_CTOF_bin_thetaq[10];
  for(int i=0; i<10; i++){
    int min = 35+(5*i);
    int max = 40+(5*i);

    sprintf(temp_name,"mom_res_CTOF_bin_thetaq_%d",min);
    sprintf(temp_title,"#Delta p/p (%d< #theta_{p} < %d);#Delta p/p;Counts",min,max);
    h_mom_res_CTOF_bin_thetaq[i] = new TH1D(temp_name,temp_title,100,-0.5,0.5);
    hist_list.push_back(h_mom_res_CTOF_bin_thetaq[i]);

    sprintf(temp_name,"mom_CTOF_bin_thetaq_%d",min);
    sprintf(temp_title,"p (%d< #theta_{p} < %d);p;Counts",min,max);
    h_mom_CTOF_bin_thetaq[i] = new TH1D(temp_name,temp_title,100,1,3);
    hist_list.push_back(h_mom_CTOF_bin_thetaq[i]);

    sprintf(temp_name,"momt_res_CTOF_bin_thetaq_%d",min);
    sprintf(temp_title,"#Delta p/p (%d< #theta_{p} < %d);#Delta p/p;Counts",min,max);
    h_momt_res_CTOF_bin_thetaq[i] = new TH1D(temp_name,temp_title,100,-0.5,0.5);
    hist_list.push_back(h_momt_res_CTOF_bin_thetaq[i]);

    sprintf(temp_name,"theta_res_CTOF_bin_thetaq_%d",min);
    sprintf(temp_title,"#Delta #theta (%d< #theta_{p} < %d);#Delta #theta;Counts",min,max);
    h_theta_res_CTOF_bin_thetaq[i] = new TH1D(temp_name,temp_title,100,-5,5);
    hist_list.push_back(h_theta_res_CTOF_bin_thetaq[i]);

    sprintf(temp_name,"phi_res_CTOF_bin_thetaq_%d",min);
    sprintf(temp_title,"cos(#theta) #Delta #phi (%d< #theta_{p} < %d);cos(#theta) #Delta #phi;Counts",min,max);
    h_phi_res_CTOF_bin_thetaq[i] = new TH1D(temp_name,temp_title,100,-5,5);
    hist_list.push_back(h_phi_res_CTOF_bin_thetaq[i]);
  }
  */

