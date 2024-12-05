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
#include "reweighter.h"
#include "Corrections.h"

using namespace std;
using namespace clas12;

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mass_n = db->GetParticle(2112)->Mass();
double mD = 1.8756;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

const double c = 29.9792458;

void SetMom(TLorentzVector &p4,double mom){
  TVector3 v3 = p4.Vect();
  double theta = v3.Theta()*180/M_PI;
  double phi = v3.Phi()*180/M_PI;
  double M = p4.M();        
  v3.SetMag(mom);
  p4.SetXYZM(v3.X(),v3.Y(),v3.Z(),M);  
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
vector<double> bE_MomFD = {0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0};
vector<double> bE_MomCD = {0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_ThetaCD = {40,45,50,55,60,70,90};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,45};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double cotan(double x){ return cos(x)/sin(x);}

double E(double x, double N, double tau){
  return N * exp( x / tau) ; 
}


double getExp(TLorentzVector balance_ptr, TLorentzVector par){
  double theta_bpar = balance_ptr.Vect().Angle(par.Vect());
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_n*mass_n) - balance_ptr.M2() - par.M2()) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*par.M2();
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}

double getExpProton(double eMom, double eTheta, double ePhi, double pTheta, double pPhi){
  TVector3 v3e;
  v3e.SetMagThetaPhi(eMom,eTheta,ePhi);
  TLorentzVector vLe;
  vLe.SetXYZM(v3e.X(),v3e.Y(),v3e.Z(),me);

  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector balance_ptr = beam + deut_ptr - vLe;

  TVector3 vp;
  vp.SetMagThetaPhi(1,pTheta,pPhi);
  
  double theta_bpar = balance_ptr.Vect().Angle(vp);
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_n*mass_n) - balance_ptr.M2() - mass_p*mass_p) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*mass_p*mass_p;
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
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

double ClosestPointMom(double x, double eMomP, double eThetaP, double ePhiP, double pMomP, double pThetaP, double pPhiP){
  return sqrt(SQ(x-eMomP) + SQ(getExpProton(x,eThetaP,ePhiP,pThetaP,pPhiP)-pMomP));
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
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",-0.75,0.75);
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
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());
  reweighter newWeight(beam_E,2,2);

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;

  TH1D * h_Delta_Eprime[9];
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
    h_phi_diff[i] = new TH1D(temp_name,temp_title,100,-30,30);
    hist_list.push_back(h_phi_diff[i]);

    sprintf(temp_name,"phi_diff_FD_%d",min);
    sprintf(temp_title,"#Delta #phi (e,e'p_{FD}) (%d< #theta_{e} < %d);#Delta #phi;Counts",min,max);
    h_phi_diff_FD[i] = new TH1D(temp_name,temp_title,100,-30,30);
    hist_list.push_back(h_phi_diff_FD[i]);

    sprintf(temp_name,"phi_diff_CD_%d",min);
    sprintf(temp_title,"#Delta #phi (e,e'p_{CD}) (%d< #theta_{e} < %d);#Delta #phi;Counts",min,max);
    h_phi_diff_CD[i] = new TH1D(temp_name,temp_title,100,-30,30);
    hist_list.push_back(h_phi_diff_CD[i]);

    sprintf(temp_name,"theta_p_%d",min);
    sprintf(temp_title,"#theta_{p} (%d< #theta_{e} < %d);#theta_{p};Counts",min,max);
    h_theta_p[i] = new TH1D(temp_name,temp_title,100,0,100);
    hist_list.push_back(h_theta_p[i]);
  }

  TH2D * h_pmiss_mmiss_FD = new TH2D("pmiss_mmiss_FD","pmiss mmiss;mmiss;pmiss;Counts",100,0.4,1.4,100,0.0,1.0);
  TH2D * h_pmiss_mmiss_CD = new TH2D("pmiss_mmiss_CD","pmiss mmiss;mmiss;pmiss;Counts",100,0.4,1.4,100,0.0,1.0);
  TH2D * h_kmiss_mmiss_FD = new TH2D("kmiss_mmiss_FD","kmiss mmiss;mmiss;kmiss;Counts",100,0.4,1.4,100,0.0,0.8);
  TH2D * h_kmiss_mmiss_CD = new TH2D("kmiss_mmiss_CD","kmiss mmiss;mmiss;kmiss;Counts",100,0.4,1.4,100,0.0,0.8);
  TH2D * h_pmissf_mmiss_FD = new TH2D("pmissf_mmiss_FD","pmissf mmiss;mmiss;pmissf;Counts",100,0.4,1.4,100,0.0,1.0);
  TH2D * h_pmissf_mmiss_CD = new TH2D("pmissf_mmiss_CD","pmissf mmiss;mmiss;pmissf;Counts",100,0.4,1.4,100,0.0,1.0);

  TH1D * h_mmissdiffFD = new TH1D("mmissdiffFD","m_{miss}-m_{n} FD;m_{miss}-m_{n};Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffFD);
  TH1D * h_mmissdiffCD = new TH1D("mmissdiffCD","m_{miss}-m_{n} CD;m_{miss}-m_{n};Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffCD);
  TH1D * h_mmissdiffcorrFD = new TH1D("mmissdiffcorrFD","m_{miss}-m_{n} FD;m_{miss}-m_{n};Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffcorrFD);
  TH1D * h_mmissdiffcorrCD = new TH1D("mmissdiffcorrCD","m_{miss}-m_{n} CD;m_{miss}-m_{n};Counts",100,-0.6,0.6);
  hist_list.push_back(h_mmissdiffcorrCD);

  TH1D * h_corr_eFD = new TH1D("corr_eFD","Correction eFD;Correction;Counts",100,-0.6,0.6);
  hist_list.push_back(h_corr_eFD);
  TH1D * h_corr_pFD = new TH1D("corr_pFD","Correction pFD;Correction;Counts",100,-0.6,0.6);
  hist_list.push_back(h_corr_pFD);
  TH1D * h_corr_eCD = new TH1D("corr_eCD","Correction eCD;Correction;Counts",100,-0.6,0.6);
  hist_list.push_back(h_corr_eCD);
  TH1D * h_corr_pCD = new TH1D("corr_pCD","Correction pCD;Correction;Counts",100,-0.6,0.6);
  hist_list.push_back(h_corr_pCD);

  TH2D * h_corr_phi_eFD = new TH2D("corr_phi_eFD","Correction eFD;#phi;Correction;Counts",100,-180,180,100,-0.6,0.6);
  hist_list.push_back(h_corr_phi_eFD);
  TH2D * h_corr_phi_pFD = new TH2D("corr_phi_pFD","Correction pFD;#phi;Correction;Counts",100,-180,180,100,-0.6,0.6);
  hist_list.push_back(h_corr_phi_pFD);
  TH2D * h_corr_phi_eCD = new TH2D("corr_phi_eCD","Correction eCD;#phi;Correction;Counts",100,-180,180,100,-0.6,0.6);
  hist_list.push_back(h_corr_phi_eCD);
  TH2D * h_corr_phi_pCD = new TH2D("corr_phi_pCD","Correction pCD;#phi;Correction;Counts",100,-180,180,100,-0.6,0.6);
  hist_list.push_back(h_corr_phi_pCD);

  TH2D * h_corr_theta_eFD = new TH2D("corr_theta_eFD","Correction eFD;#theta;Correction;Counts",100,8,40,100,-0.6,0.6);
  hist_list.push_back(h_corr_theta_eFD);
  TH2D * h_corr_theta_pFD = new TH2D("corr_theta_pFD","Correction pFD;#theta;Correction;Counts",100,8,40,100,-0.6,0.6);
  hist_list.push_back(h_corr_theta_pFD);
  TH2D * h_corr_theta_eCD = new TH2D("corr_theta_eCD","Correction eCD;#theta;Correction;Counts",100,8,40,100,-0.6,0.6);
  hist_list.push_back(h_corr_theta_eCD);
  TH2D * h_corr_theta_pCD = new TH2D("corr_theta_pCD","Correction pCD;#theta;Correction;Counts",100,30,100,100,-0.6,0.6);
  hist_list.push_back(h_corr_theta_pCD);

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
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
  TH2D * h_mom_theta_eFD_binSector[6];
  TH2D * h_mom_theta_eCD_binSector[6];
  TH2D * h_mom_theta_pFD_binSector[6];
  TH2D * h_mom_theta_COMB_binSector[6];
  for(int j=1; j<=6; j++){
    sprintf(temp_name,"mom_theta_eFD_sector_%d",j);
    sprintf(temp_title,"#theta vs. p eFD Sector %d;p;#theta;Counts",j);
    h_mom_theta_eFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,0,6,100,5,45);
    hist_list.push_back(h_mom_theta_eFD_binSector[j-1]);
    
    sprintf(temp_name,"mom_theta_eCD_sector_%d",j);
    sprintf(temp_title,"#theta vs. p eCD Sector %d;p;#theta;Counts",j);
    h_mom_theta_eCD_binSector[j-1] = new TH2D(temp_name,temp_title,100,0,6,100,5,45);
    hist_list.push_back(h_mom_theta_eCD_binSector[j-1]);

    sprintf(temp_name,"mom_theta_pFD_sector_%d",j);
    sprintf(temp_title,"#theta vs. p pFD Sector %d;p;#theta;Counts",j);
    h_mom_theta_pFD_binSector[j-1] = new TH2D(temp_name,temp_title,100,0,6,100,5,45);
    hist_list.push_back(h_mom_theta_pFD_binSector[j-1]);

    sprintf(temp_name,"mom_theta_COMB_sector_%d",j);
    sprintf(temp_title,"#theta vs. p COMB Sector %d;p;#theta;Counts",j);
    h_mom_theta_COMB_binSector[j-1] = new TH2D(temp_name,temp_title,100,0,6,100,5,45);
    hist_list.push_back(h_mom_theta_COMB_binSector[j-1]);
  }
  TH2D * h_phi_theta_pCD = new TH2D("h_phi_theta_eFD","#theta vs. #phi pCD;#phi;#theta;Counts",100,-180,180,100,30,90);
  hist_list.push_back(h_phi_theta_pCD);
  TH2D * h_mom_theta_pCD = new TH2D("h_mom_theta_eFD","#theta vs. p pCD;p;#theta;Counts",100,0,3,100,30,90);
  hist_list.push_back(h_mom_theta_pCD);
  TH2D * h_mom_theta_pFDCD = new TH2D("h_mom_theta_eFD","#theta vs. p pFDCD;p;#theta;Counts",100,0,4.5,100,20,90);
  hist_list.push_back(h_mom_theta_pFDCD);

  
  TH2D * h_phi_corr_binSector_binTheta[6][12];
  TGraphErrors * g_phi_corr_binSector_binTheta[6][12];
  for(int j=1; j<=6; j++){
    for(int i=0; i<12; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction;Counts",j,min,max);
      h_phi_corr_binSector_binTheta[j-1][i] = new TH2D(temp_name,temp_title,45,-45,45,100,-0.5,0.5);
      hist_list.push_back(h_phi_corr_binSector_binTheta[j-1][i]);

      g_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_phi_corr_sector_%d_theta_%d",j,i);
      g_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

    }
  }

  TH2D * h_phi_corr_binSector_binMom[6][11];
  TGraphErrors * g_phi_corr_binSector_binMom[6][11];
  for(int j=1; j<=6; j++){
    for(int i=0; i<11; i++){
      int min = bE_MomFD[i]*1000;
      int max = bE_MomFD[i+1]*1000;
      sprintf(temp_name,"phi_corr_sector_%d_mom_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< p < %d);#phi;Correction;Counts",j,min,max);
      h_phi_corr_binSector_binMom[j-1][i] = new TH2D(temp_name,temp_title,45,-45,45,100,-0.5,0.5);
      hist_list.push_back(h_phi_corr_binSector_binMom[j-1][i]);

      g_phi_corr_binSector_binMom[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_phi_corr_sector_%d_mom_%d",j,i);
      g_phi_corr_binSector_binMom[j-1][i]->SetName(temp_name);

    }
  }

  
  TH2D * h_phi_corr_binThetaCD[6];
  TGraphErrors * g_phi_corr_binThetaCD[6];
  for(int i=0; i<6; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"phi_corr_theta_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< #theta < %d);#phi;Correction;Counts",min,max);
    h_phi_corr_binThetaCD[i] = new TH2D(temp_name,temp_title,180,-180,180,100,-0.5,0.5);
    hist_list.push_back(h_phi_corr_binThetaCD[i]);
    
    g_phi_corr_binThetaCD[i] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_theta_%d",i);
    g_phi_corr_binThetaCD[i]->SetName(temp_name);
  }

  TH2D * h_phi_corr_binMomCD[6];
  TGraphErrors * g_phi_corr_binMomCD[6];
  for(int i=0; i<6; i++){
    int min = bE_MomCD[i]*1000;
    int max = bE_MomCD[i+1]*1000;
    sprintf(temp_name,"phi_corr_mom_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< p < %d);#phi;Correction;Counts",min,max);
    h_phi_corr_binMomCD[i] = new TH2D(temp_name,temp_title,180,-180,180,100,-0.5,0.5);
    hist_list.push_back(h_phi_corr_binMomCD[i]);
    
    g_phi_corr_binMomCD[i] = new TGraphErrors();
    sprintf(temp_name,"g_phi_corr_mom_%d",i);
    g_phi_corr_binMomCD[i]->SetName(temp_name);
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
	  SetLorentzVector_ThetaCorrection(el,electrons[0]);
	  TLorentzVector el_corrected = el;
	  SetLorentzVector_MomentumCorrection(el_corrected,electrons[0]);
	  TLorentzVector q_corrected = beam - el_corrected;

	  
	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
	  int sector_e = electrons[0]->getSector();
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
	  double omega_diff = omega - (Q2/(2*mN));

	  double vtz_e = electrons[0]->par()->getVz();
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double theta_e = el.Theta()*180/M_PI;
	  double phi_e = el.Phi() * 180/M_PI;
	  double shift_e = 0;
	  shift_e += (sector_e==1)?0:(sector_e==2)?60:(sector_e==3)?120:(sector_e==4 && phi_e>0)?180:(sector_e==4 && phi_e<0)?-180:(sector_e==5)?-120:(sector_e==6)?-60:0;
	  double theta_q = q.Theta()*180/M_PI;
	  	  
	  if(tBinFDe(theta_e)!=-1){
	    h_Delta_Eprime[tBinFDe(theta_e)]->Fill(Delta_Eprime,wep);
	    h_Delta_Eprime_sectors[sector_e-1][tBinFDe(theta_e)]->Fill(Delta_Eprime,wep);
	    if(Delta_Eprime<0.1){
	      h_Q2[tBinFDe(theta_e)]->Fill(Q2,wep);
	      h_p_e[tBinFDe(theta_e)]->Fill(el.P(),wep);
	      h_theta_q[tBinFDe(theta_e)]->Fill(theta_q,wep);
	    }
	  }
	  if(protons.size() <= 0){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  SetLorentzVector_ThetaCorrection(proton_ptr,protons[0]);
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);

	  TLorentzVector proton_ptr_fix = proton_ptr;
	  double pMom_fixed = getExpProton(el.P(),el.Theta(),el.Phi(),proton_ptr.Theta(),proton_ptr.Phi());
	  SetMom(proton_ptr_fix,pMom_fixed);
	  TLorentzVector miss_fixed = q + deut_ptr - proton_ptr_fix;

	  TLorentzVector proton_ptr_corrected = proton_ptr;
	  SetLorentzVector_MomentumCorrection(proton_ptr_corrected,protons[0]);
	  TLorentzVector miss_corrected = q_corrected + deut_ptr - proton_ptr_corrected;
	  
	  TLorentzVector miss = q + deut_ptr - proton_ptr;

	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  
	  double mom_p = proton_ptr.P();	  
	  double theta_p = proton_ptr.Theta()*180/M_PI;
	  double phi_p = proton_ptr.Phi()*180/M_PI;
	  double shift_p = 0;
	  
	  double Delta_mom = (mom_p-mom_q)/mom_q;
	  //cout<<Delta_mom<<endl;
	  double Delta_theta = theta_p-theta_q;
	  double Delta_phi = (phi_p-phi_q);

	  TLorentzVector miss_LC = proton_ptr - q;      
	  TVector3 u_ZQ = q.Vect().Unit();
	  double pmm_ZQ = miss_LC.E() - miss_LC.Vect().Dot(u_ZQ);
	  double pmp_ZQ = miss_LC.Vect().Perp(u_ZQ);
	  double kmiss_ZQ = sqrt(mN*mN*((pmp_ZQ*pmp_ZQ+mN*mN)/(pmm_ZQ*(2*mN-pmm_ZQ))) - mN*mN);


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
	  }


	  
	  if(fabs(Delta_phi)>7){continue;}

	  if(protons[0]->getRegion()==FD){
	    h_pmiss_mmiss_FD->Fill(miss.M(),miss.P(),wep);	  
	    h_kmiss_mmiss_FD->Fill(miss.M(),kmiss_ZQ,wep);
	    h_pmissf_mmiss_FD->Fill(miss.M(),miss_fixed.P(),wep);    
	  }
	  else if(protons[0]->getRegion()==CD){
	    h_pmiss_mmiss_CD->Fill(miss.M(),miss.P(),wep);	  
	    h_kmiss_mmiss_CD->Fill(miss.M(),kmiss_ZQ,wep);
	    h_pmissf_mmiss_CD->Fill(miss.M(),miss_fixed.P(),wep);    
	  }
	  //if(miss.M()>0.7){continue;}
	  TF1 * fDist = new TF1("fDist",[&](double *x, double *p){ return ClosestPointMom(x[0],p[0],p[1],p[2],p[3],p[4],p[5]); },0.1,6.1,6);
	  fDist->SetParameter(0,el.P());
	  fDist->SetParameter(1,el.Theta());
	  fDist->SetParameter(2,el.Phi());
	  fDist->SetParameter(3,proton_ptr.P());
	  fDist->SetParameter(4,proton_ptr.Theta());
	  fDist->SetParameter(5,proton_ptr.Phi());
	  double eMom = fDist->GetMinimumX(0.1,6.1);
	  double pMom = getExpProton(eMom,el.Theta(),el.Phi(),proton_ptr.Theta(),proton_ptr.Phi());
	  
	  h_mom_theta_pFDCD->Fill(proton_ptr.P(),theta_p,wep);

		    
	  /////////////
	  if(binX(bE_Theta,theta_e)!=-1){
	    h_phi_corr_binSector_binTheta[sector_e-1][binX(bE_Theta,theta_e)]->Fill(phi_e-shift_e,eMom-el.P(),wep);
	  }
	  if(binX(bE_MomFD,el.P())!=-1){
	    h_phi_corr_binSector_binMom[sector_e-1][binX(bE_MomFD,el.P())]->Fill(phi_e-shift_e,eMom-el.P(),wep);
	  }
	  /////////////

	  if(protons[0]->getRegion()==FD){
	    int sector_p = protons[0]->getSector();	  
	    shift_p += (sector_p==1)?0:(sector_p==2)?60:(sector_p==3)?120:(sector_p==4 && phi_p>0)?180:(sector_p==4 && phi_p<0)?-180:(sector_p==5)?-120:(sector_p==6)?-60:0;
	    

	    h_mmissdiffFD->Fill(miss.M()-mN,wep);
	    h_mmissdiffcorrFD->Fill(miss_corrected.M()-mN,wep);
	    h_corr_eFD->Fill(eMom-el.P(),wep);
	    h_corr_phi_eFD->Fill(el.Phi()*180/M_PI,eMom-el.P(),wep);
	    h_corr_theta_eFD->Fill(el.Theta()*180/M_PI,eMom-el.P(),wep);
	    h_corr_pFD->Fill(pMom-proton_ptr.P(),wep);
	    h_corr_phi_pFD->Fill(proton_ptr.Phi()*180/M_PI,pMom-proton_ptr.P(),wep);
	    h_corr_theta_pFD->Fill(proton_ptr.Theta()*180/M_PI,pMom-proton_ptr.P(),wep);

	    //
	    h_phi_theta_eFD_binSector[sector_e-1]->Fill(phi_e-shift_e,theta_e,wep);
	    h_phi_theta_pFD_binSector[sector_p-1]->Fill(phi_p-shift_p,theta_p,wep);
	    h_phi_theta_COMB_binSector[sector_e-1]->Fill(phi_e-shift_e,theta_e,wep);
	    h_phi_theta_COMB_binSector[sector_p-1]->Fill(phi_p-shift_p,theta_p,wep);

	    h_mom_theta_eFD_binSector[sector_e-1]->Fill(el.P(),theta_e,wep);
	    h_mom_theta_pFD_binSector[sector_p-1]->Fill(proton_ptr.P(),theta_p,wep);
	    h_mom_theta_COMB_binSector[sector_e-1]->Fill(el.P(),theta_e,wep);
	    h_mom_theta_COMB_binSector[sector_p-1]->Fill(proton_ptr.P(),theta_p,wep);

	    /////////////
	    if(binX(bE_Theta,theta_p)!=-1){
	      h_phi_corr_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(phi_p-shift_p,pMom-proton_ptr.P(),wep);
	    }
	    if(binX(bE_MomFD,mom_p)!=-1){
	      h_phi_corr_binSector_binMom[sector_p-1][binX(bE_MomFD,mom_p)]->Fill(phi_p-shift_p,pMom-proton_ptr.P(),wep);
	    }
	    /////////////
	  }
	  else if(protons[0]->getRegion()==CD){
	    h_mmissdiffCD->Fill(miss.M()-mN,wep);
	    h_mmissdiffcorrCD->Fill(miss_corrected.M()-mN,wep);
	    h_corr_eCD->Fill(eMom-el.P(),wep);
	    h_corr_phi_eCD->Fill(el.Phi()*180/M_PI,eMom-el.P(),wep);
	    h_corr_theta_eCD->Fill(el.Theta()*180/M_PI,eMom-el.P(),wep);
	    h_corr_pCD->Fill(pMom-proton_ptr.P(),wep);
	    h_corr_phi_pCD->Fill(proton_ptr.Phi()*180/M_PI,pMom-proton_ptr.P(),wep);
	    h_corr_theta_pCD->Fill(proton_ptr.Theta()*180/M_PI,pMom-proton_ptr.P(),wep);

	    //
	    h_phi_theta_eCD_binSector[sector_e-1]->Fill(phi_e-shift_e,theta_e,wep);
	    h_phi_theta_COMB_binSector[sector_e-1]->Fill(phi_e-shift_e,theta_e,wep);
	    h_phi_theta_pCD->Fill(phi_p-shift_p,theta_p,wep);

	    h_mom_theta_eCD_binSector[sector_e-1]->Fill(el.P(),theta_e,wep);
	    h_mom_theta_COMB_binSector[sector_e-1]->Fill(el.P(),theta_e,wep);
	    h_mom_theta_pCD->Fill(proton_ptr.P(),theta_p,wep);

	    /////////////
	    if(binX(bE_ThetaCD,theta_p)!=-1){
	      h_phi_corr_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(phi_p,pMom-proton_ptr.P(),wep);
	    }
	    if(binX(bE_MomCD,proton_ptr.P())!=-1){
	      h_phi_corr_binMomCD[binX(bE_MomCD,proton_ptr.P())]->Fill(phi_p,pMom-proton_ptr.P(),wep);
	    }
	    /////////////
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
    h_Delta_Eprime[i-1]->Draw();    
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
  h_pmiss_mmiss_FD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_mmiss_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_kmiss_mmiss_FD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_kmiss_mmiss_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmissf_mmiss_FD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmissf_mmiss_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmissdiffcorrFD->SetLineColor(2);
  h_mmissdiffcorrFD->Draw();
  h_mmissdiffFD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmissdiffcorrCD->SetLineColor(2);
  h_mmissdiffcorrCD->Draw();
  h_mmissdiffCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_corr_eFD->Draw();
  myCanvas->cd(2);
  h_corr_pFD->Draw();
  myCanvas->cd(3);
  h_corr_eCD->Draw();
  myCanvas->cd(4);
  h_corr_pCD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_corr_phi_eFD->Draw("colz");
  TGraphErrors * g_corr_phi_eFD = new TGraphErrors();
  getGraph(h_corr_phi_eFD,g_corr_phi_eFD);
  g_corr_phi_eFD->Draw("SAME");
  myCanvas->cd(2);
  h_corr_phi_pFD->Draw("colz");
  h_corr_phi_pFD->Draw("colz");
  TGraphErrors * g_corr_phi_pFD = new TGraphErrors();
  getGraph(h_corr_phi_pFD,g_corr_phi_pFD);
  g_corr_phi_pFD->Draw("SAME");
  myCanvas->cd(3);
  h_corr_phi_eCD->Draw("colz");
  h_corr_phi_eCD->Draw("colz");
  TGraphErrors * g_corr_phi_eCD = new TGraphErrors();
  getGraph(h_corr_phi_eCD,g_corr_phi_eCD);
  g_corr_phi_eCD->Draw("SAME");
  myCanvas->cd(4);
  h_corr_phi_pCD->Draw("colz");
  h_corr_phi_pCD->Draw("colz");
  TGraphErrors * g_corr_phi_pCD = new TGraphErrors();
  getGraph(h_corr_phi_pCD,g_corr_phi_pCD);
  g_corr_phi_pCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_corr_theta_eFD->Draw("colz");
  TGraphErrors * g_corr_theta_eFD = new TGraphErrors();
  getGraph(h_corr_theta_eFD,g_corr_theta_eFD);
  g_corr_theta_eFD->Draw("SAME");
  myCanvas->cd(2);
  h_corr_theta_pFD->Draw("colz");
  h_corr_theta_pFD->Draw("colz");
  TGraphErrors * g_corr_theta_pFD = new TGraphErrors();
  getGraph(h_corr_theta_pFD,g_corr_theta_pFD);
  g_corr_theta_pFD->Draw("SAME");
  myCanvas->cd(3);
  h_corr_theta_eCD->Draw("colz");
  h_corr_theta_eCD->Draw("colz");
  TGraphErrors * g_corr_theta_eCD = new TGraphErrors();
  getGraph(h_corr_theta_eCD,g_corr_theta_eCD);
  g_corr_theta_eCD->Draw("SAME");
  myCanvas->cd(4);
  h_corr_theta_pCD->Draw("colz");
  h_corr_theta_pCD->Draw("colz");
  TGraphErrors * g_corr_theta_pCD = new TGraphErrors();
  getGraph(h_corr_theta_pCD,g_corr_theta_pCD);
  g_corr_theta_pCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_eFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_eCD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_pFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_theta_COMB_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_eFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_eCD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_pFD_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,2);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_theta_COMB_binSector[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phi_theta_pCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_theta_pCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mom_theta_pFDCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

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

  for(int j = 0; j < 6; j++){
    myCanvas->Divide(3,4);
    for(int i = 0; i < 11; i++){
      myCanvas->cd(i+1);
      h_phi_corr_binSector_binMom[j][i]->Draw("colz");
      getGraph(h_phi_corr_binSector_binMom[j][i],g_phi_corr_binSector_binMom[j][i]);
      g_phi_corr_binSector_binMom[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    getGraph(h_phi_corr_binThetaCD[i],g_phi_corr_binThetaCD[i]);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binMomCD[i]->Draw("colz");
    getGraph(h_phi_corr_binMomCD[i],g_phi_corr_binMomCD[i]);
    g_phi_corr_binMomCD[i]->Draw("SAME");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}
