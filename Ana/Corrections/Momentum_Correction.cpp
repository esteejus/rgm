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

const double c = 29.9792458;
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

//vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_MomCD = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};//{0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
vector<double> bE_ThetaCD = {38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92};//{40,45,50,55,60,70,90};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};

vector<double> bE_ThetaE = {10,13,16,19,22,25,28,31,34,37};
vector<double> bE_ThetapFD = {19,22,25,28,31,34,37,40,43,46};

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

double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
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


void getMax(TH1D * hist, TF1 * func, double & mean_fit, double & sigma_fit){
  int N = hist->GetEntries();
  if(N<20000){hist->Rebin(2);}
  else if(N<2000){hist->Rebin(4);}
  double max = hist->GetMaximum();
  double mode = hist->GetBinCenter(hist->GetMaximumBin());
  double diff = fabs(mode - hist->GetMean());
  if(mode<0){mode=0;}
  double stddev = hist->GetStdDev();

  func->SetParameter(0,max);
  func->SetParameter(1,mode);
  func->SetParLimits(1,mode-diff,mode+diff);
  func->SetParameter(2,stddev);
  func->SetParLimits(2,0.001,stddev);
  TFitResultPtr point = hist->Fit(func,"SrBeqn","",mode-stddev,mode+0.9*stddev);
  mean_fit = point->Parameter(1);
  sigma_fit = point->Parameter(2);
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
  vector<double> x_rad = {11.5,14.5,17.5,20.5,23.5,26.5,29.5,32.5,35.5};  
  vector<double> y_rad = {12.24,14.7628,14.7636,15.2256,17.0205,22.3517,23.9472,29.5068};  
  TGraph * g_rad = new TGraph;
  g_rad->SetName("g_rad");
  g_rad->SetTitle("#Delta E' Shift Due to Radiation;#theta;#Delta E'");
  for(int i = 0; i < 8; i++){
    g_rad->SetPoint(g_rad->GetN(),x_rad[i],y_rad[i]);
  }
  TF1 * f_rad = new TF1("f_rad",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },10,40,3);
  f_rad->SetParameter(0,5.0);
  f_rad->SetParameter(1,1.0);
  f_rad->SetParameter(2,1.0);
  TFitResultPtr p_rad = g_rad->Fit(f_rad,"SrBeqn","",10,40);



  //new after redo
  vector<double> x_rad_new = {9,11,13,15,17,19,21.5,24.5,27.5,30.5,33.5};
  vector<double> y_rad_new = {6.81449,9.87645,11.6931,12.5859,12.0127,13.2819,14.7657,16.5664,15.8999,22.4731,92.194};
  vector<double> e_rad_new = {14.1635,15.2244,17.2681,19.2462,21.3319,21.9251,24.3228,31.3548,33.7588,34.8207,39.8671};
  TGraph * g_rad_new = new TGraph;
  g_rad_new->SetName("g_rad_new");
  g_rad_new->SetTitle("#Delta p_{e} Shift Due to Radiation;#theta;#mu_{#Delta E'} [MeV]");
  TGraph * g_rade_new = new TGraph;
  g_rade_new->SetName("g_rade_new");
  g_rade_new->SetTitle("#Delta E' Resolution;#theta;#sigma_{#Delta E'}");
  for(int i = 0; i < 10; i++){
    g_rad_new->SetPoint(g_rad_new->GetN(),x_rad_new[i],y_rad_new[i]);
    g_rade_new->SetPoint(g_rad_new->GetN(),x_rad_new[i],e_rad_new[i]);
  }
  TF1 * f_rad_new = new TF1("f_rad_new",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },8,40,3);
  f_rad_new->SetParameter(0,5.0);
  f_rad_new->SetParameter(1,1.0);
  f_rad_new->SetParameter(2,1.0);
  TFitResultPtr p_rad_new = g_rad_new->Fit(f_rad_new,"SrBeqn","",8,32);

  
  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;

  TH1D * h_Delta_Int_Eprime_beforeRad = new TH1D("Delta_Eprime_beforeRad","#Delta p_{e} (e,e'p);#Delta p_{e} [GeV];Counts",100,-0.2,0.2);
  TH1D * h_Delta_Int_Eprime = new TH1D("Delta_Eprime","#Delta p_{e} - #Delta p_{e}^{radiation} (e,e'p);#Delta p_{e} [GeV];Counts",100,-0.2,0.2);
  TH1D * h_Delta_Int_Eprime_Corrected = new TH1D("Delta_Eprime_Corrected","#Delta p_{e} (e,e'p);#Delta p_{e} [GeV];Counts",100,-0.2,0.2);
  TH1D * h_Delta_Int_pMomFD = new TH1D("Delta_pMomFD","#Delta p_{p} (e,e'p_{FD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  TH1D * h_Delta_Int_pMomFD_Corrected = new TH1D("Delta_pMomFD_Corrected","#Delta p_{p} (e,e'p_{FD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  TH1D * h_Delta_Int_pMomCD = new TH1D("Delta_pMomCD","#Delta p_{p} (e,e'p_{CD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);
  TH1D * h_Delta_Int_pMomCD_Corrected = new TH1D("Delta_pMomCD_Corrected","#Delta p_{p} (e,e'p_{CD});#Delta p_{p} [GeV]; Counts",100,-0.4,0.4);

  TH2D * h_phie_Dp_FD = new TH2D("phie_Dp_FD","#Delta p_{e} vs. #phi_{e} (e,e'p_{FD});#phi_{e}^{#circ};#Delta p_{e} [GeV];Counts",180,-180,180,100,-0.2,0.2);
  hist_list.push_back(h_phie_Dp_FD);
  TH2D * h_phie_Dp_CD = new TH2D("phie_Dp_CD","#Delta p_{e} vs. #phi_{e} (e,e'p_{CD});#phi_{e}^{#circ};#Delta p_{e} [GeV];Counts",180,-180,180,100,-0.2,0.2);
  hist_list.push_back(h_phie_Dp_CD);
  TH2D * h_phip_Dp_FD = new TH2D("phip_Dp_FD","#Delta p_{p} vs. #phi_{p} (e,e'p_{FD});#phi_{p}^{#circ};#Delta p_{p} [GeV];Counts",180,-180,180,100,-0.4,0.4);
  hist_list.push_back(h_phip_Dp_FD);
  TH2D * h_phip_Dp_CD = new TH2D("phip_Dp_CD","#Delta p_{p} vs. #phi_{p} (e,e'p_{CD});#phi_{p}^{#circ};#Delta p_{p} [GeV];Counts",180,-180,180,100,-0.4,0.4);
  hist_list.push_back(h_phip_Dp_CD);

  TH2D * h_thetae_Dp_FD = new TH2D("thetae_Dp_FD","#Delta p_{e} vs. #theta_{e} (e,e'p_{FD})^{#circ};#theta_{e};#Delta p_{e} [GeV];Counts",100,15,27,100,-0.2,0.2);
  hist_list.push_back(h_thetae_Dp_FD);
  TH2D * h_thetae_Dp_CD = new TH2D("thetae_Dp_CD","#Delta p_{e} vs. #theta_{e} (e,e'p_{CD})^{#circ};#theta_{e};#Delta p_{e} [GeV];Counts",100,7,23,100,-0.2,0.2);
  hist_list.push_back(h_thetae_Dp_CD);
  TH2D * h_thetap_Dp_FD = new TH2D("thetap_Dp_FD","#Delta p_{p} vs. #theta_{p} (e,e'p_{FD})^{#circ};#theta_{p};#Delta p_{p} [GeV];Counts",100,25,42,100,-0.4,0.4);
  hist_list.push_back(h_thetap_Dp_FD);
  TH2D * h_thetap_Dp_CD = new TH2D("thetap_Dp_CD","#Delta p_{p} vs. #theta_{p} (e,e'p_{CD})^{#circ};#theta_{p};#Delta p_{p} [GeV];Counts",100,35,65,100,-0.4,0.4);
  hist_list.push_back(h_thetap_Dp_CD);

  TH2D * h_phie_Dp_Corr_FD = new TH2D("phie_Dp_Corr_FD","#Delta p_{e} vs. #phi_{e} (e,e'p_{FD});#phi_{e}^{#circ};#Delta p_{e} [GeV];Counts",180,-180,180,100,-0.2,0.2);
  hist_list.push_back(h_phie_Dp_Corr_FD);
  TH2D * h_phie_Dp_Corr_CD = new TH2D("phie_Dp_Corr_CD","#Delta p_{e} vs. #phi_{e} (e,e'p_{CD});#phi_{e}^{#circ};#Delta p_{e} [GeV];Counts",180,-180,180,100,-0.2,0.2);
  hist_list.push_back(h_phie_Dp_Corr_CD);
  TH2D * h_phip_Dp_Corr_FD = new TH2D("phip_Dp_Corr_FD","#Delta p_{p} vs. #phi_{p} (e,e'p_{FD});#phi_{p}^{#circ};#Delta p_{p} [GeV];Counts",180,-180,180,100,-0.4,0.4);
  hist_list.push_back(h_phip_Dp_Corr_FD);
  TH2D * h_phip_Dp_Corr_CD = new TH2D("phip_Dp_Corr_CD","#Delta p_{p} vs. #phi_{p} (e,e'p_{CD});#phi_{p}^{#circ};#Delta p_{p} [GeV];Counts",180,-180,180,100,-0.4,0.4);
  hist_list.push_back(h_phip_Dp_Corr_CD);
  
  
  TH1D * h_Delta_Res_Eprime_Corrected = new TH1D("Res_Delta_Eprime_Corrected","#Delta E'/E';#Delta E'/E';Counts",100,-0.05,0.05);
  TH1D * h_Delta_Res_pMomFD_Corrected = new TH1D("Res_Delta_pMomFD_Corrected","#Delta p/p for FD Protons;#Delta p/p; Counts",100,-0.2,0.2);
  TH1D * h_Delta_Res_pMomCD_Corrected = new TH1D("Res_Delta_pMomCD_Corrected","#Delta p/p for CD Protons;#Delta p/p; Counts",100,-0.4,0.4);

  
  TH1D * h_E_Res[14];
  TH1D * h_E_Res_Corrected[14];
  for(int i=0; i<14; i++){
    int min = bE_Theta[i];
    int max = bE_Theta[i+1];
    sprintf(temp_name,"h_E_Res_%d",i);
    sprintf(temp_title,"Counts vs. #Delta p_{e} (%d< #theta < %d);#Delta p_{e} [GeV];Counts",min,max);
    h_E_Res[i] = new TH1D(temp_name,temp_title,100,-0.15,0.15);
    sprintf(temp_name,"h_E_Res_%d_Corrected",i);
    h_E_Res_Corrected[i] = new TH1D(temp_name,temp_title,100,-0.15,0.15);
  }
  
  TH1D * h_pFD_Res[9];
  TH1D * h_pFD_Res_Corrected[9];
  for(int i=0; i<9; i++){
    int min = bE_ThetapFD[i];
    int max = bE_ThetapFD[i+1];
    sprintf(temp_name,"h_pFD_Res_%d",i);
    sprintf(temp_title,"Counts vs. #Delta p FD Protons (%d< #theta < %d);#Delta E';Counts",min,max);
    h_pFD_Res[i] = new TH1D(temp_name,temp_title,100,-0.3,0.3);
    sprintf(temp_name,"h_pFD_Res_%d_Corrected",i);
    h_pFD_Res_Corrected[i] = new TH1D(temp_name,temp_title,100,-0.3,0.3);
  }
  
  
  TH2D * h_phi_corr_binSector_binTheta[6][14];
  TGraphErrors * g_phi_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
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
  
  TH2D * h_e_phi_corr_binSector_binTheta[6][14];
  TGraphErrors * g_e_phi_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"e_phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction;Counts",j,min,max);
      h_e_phi_corr_binSector_binTheta[j-1][i] = new TH2D(temp_name,temp_title,45,-45,45,100,-0.2,0.2);
      hist_list.push_back(h_e_phi_corr_binSector_binTheta[j-1][i]);

      g_e_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_e_phi_corr_sector_%d_theta_%d",j,i);
      g_e_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

    }
  }

  TH2D * h_p_phi_corr_binSector_binTheta[6][14];
  TGraphErrors * g_p_phi_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"p_phi_corr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Correction vs. #phi Sector %d (%d< #theta < %d);#phi;Correction;Counts",j,min,max);
      h_p_phi_corr_binSector_binTheta[j-1][i] = new TH2D(temp_name,temp_title,45,-55,35,100,-0.5,0.5);
      hist_list.push_back(h_p_phi_corr_binSector_binTheta[j-1][i]);

      g_p_phi_corr_binSector_binTheta[j-1][i] = new TGraphErrors();
      sprintf(temp_name,"g_p_phi_corr_sector_%d_theta_%d",j,i);
      g_p_phi_corr_binSector_binTheta[j-1][i]->SetName(temp_name);

    }
  }

  
  TH2D * h_phi_corr_binThetaCD[18];
  TGraphErrors * g_phi_corr_binThetaCD[18];
  for(int i=0; i<18; i++){
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

  TH2D * h_aftercorr_phi_corr_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"aftercorr_phi_corr_theta_%d",i);
    sprintf(temp_title,"Correction vs. #phi (%d< #theta < %d);#phi;Correction;Counts",min,max);
    h_aftercorr_phi_corr_binThetaCD[i] = new TH2D(temp_name,temp_title,180,-180,180,100,-0.5,0.5);
    hist_list.push_back(h_aftercorr_phi_corr_binThetaCD[i]);
  }

  
  TH2D * h_phi_corr_binMomCD[12];
  TGraphErrors * g_phi_corr_binMomCD[12];
  for(int i=0; i<12; i++){
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

	  
	  //SetLorentzVectorCorrectedTheta(el_corrected,electrons[0]);
	  //SetLorentzVectorCorrectedMomentum(el_corrected,electrons[0]);


	  double Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN );
	  TVector3 vel_fromAngle;
	  vel_fromAngle.SetMagThetaPhi(Eprime,el.Theta(),el.Phi());
	  TVector3 vp_fromAngle = vbeam-vel_fromAngle;
	  double Delta_Eprime = ( mN*beam_E )/( beam_E*(1-cos(el.Theta())) + mN ) - el.P();
	  ///
	  double Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN );
	  TVector3 vel_fromAngle_corrected;
	  vel_fromAngle_corrected.SetMagThetaPhi(Eprime_corrected,el_corrected.Theta(),el_corrected.Phi());
	  TVector3 vp_fromAngle_corrected = vbeam-vel_fromAngle_corrected;
	  double Delta_Eprime_corrected = ( mN*beam_E )/( beam_E*(1-cos(el_corrected.Theta())) + mN ) - el_corrected.P();
	  ///
	  
	  int sector_e = electrons[0]->getSector();
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
	  double Electron_Radiation = Quad(theta_e,p_rad_new->Parameter(0),p_rad_new->Parameter(1),p_rad_new->Parameter(2))/1000.0;

	  double theta_e_corrected = el_corrected.Theta()*180/M_PI;
	  double theta_q = q.Theta()*180/M_PI;
	  double shift_e = 0;
	  shift_e += (sector_e==1)?0:(sector_e==2)?60:(sector_e==3)?120:(sector_e==4 && phi_e>0)?180:(sector_e==4 && phi_e<0)?-180:(sector_e==5)?-120:(sector_e==6)?-60:0;

	  if(protons.size() <= 0){continue;}
	  GetLorentzVector_ReconVector(proton_ptr,protons[0]);
	  SetLorentzVector_ThetaCorrection(proton_ptr,protons[0]);
	  SetLorentzVector_EnergyLossCorrection(proton_ptr,protons[0]);
	  TLorentzVector proton_ptr_corrected = proton_ptr;
	  SetLorentzVector_MomentumCorrection(proton_ptr_corrected,protons[0]);

	  double mom_q = q.P();
	  double phi_q = q.Phi()*180/M_PI;
	  
	  double mom_p = proton_ptr.P();
	  double Delta_Mom_FromAngle = vp_fromAngle.Mag()-mom_p;
	  double mom_p_corrected = proton_ptr_corrected.P();
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
	  if(fabs(Delta_phi)>3){continue;}
	  
	  h_Delta_Int_Eprime_beforeRad->Fill(Delta_Eprime,wep);
	  h_Delta_Int_Eprime->Fill(Delta_Eprime-Electron_Radiation,wep);
	  h_Delta_Int_Eprime_Corrected->Fill(Delta_Eprime_corrected,wep);
	  h_Delta_Res_Eprime_Corrected->Fill(Delta_Eprime_corrected/Eprime_corrected,wep);

	  
	  if(protons[0]->getRegion()==FD){
	    h_phie_Dp_FD->Fill(phi_e,Delta_Eprime-Electron_Radiation,wep);
	    h_phip_Dp_FD->Fill(phi_p,Delta_Mom_FromAngle,wep);
	    h_phie_Dp_Corr_FD->Fill(phi_e,Delta_Eprime_corrected-Electron_Radiation,wep);
	    h_phip_Dp_Corr_FD->Fill(phi_p,Delta_Mom_FromAngle_corrected,wep);

	    h_thetae_Dp_FD->Fill(theta_e,Delta_Eprime-Electron_Radiation,wep);
	    h_thetap_Dp_FD->Fill(theta_p,Delta_Mom_FromAngle,wep);
	  }
	  else{
	    h_phie_Dp_CD->Fill(phi_e,Delta_Eprime-Electron_Radiation,wep);
	    h_phip_Dp_CD->Fill(phi_p,Delta_Mom_FromAngle,wep);
	    h_phie_Dp_Corr_CD->Fill(phi_e,Delta_Eprime_corrected-Electron_Radiation,wep);
	    h_phip_Dp_Corr_CD->Fill(phi_p,Delta_Mom_FromAngle_corrected,wep);
	    
	    h_thetae_Dp_CD->Fill(theta_e,Delta_Eprime-Electron_Radiation,wep);
	    h_thetap_Dp_CD->Fill(theta_p,Delta_Mom_FromAngle,wep);
	  }
	    
	  if(binX(bE_Theta,theta_e)!=-1){
	    h_E_Res[binX(bE_Theta,theta_e)]->Fill(Delta_Eprime,wep);
	    h_E_Res_Corrected[binX(bE_Theta,theta_e)]->Fill(Delta_Eprime_corrected,wep);
	  }

	  if(binX(bE_Theta,theta_e)!=-1){
	    h_phi_corr_binSector_binTheta[sector_e-1][binX(bE_Theta,theta_e)]->Fill(phi_e-shift_e,Delta_Eprime-Electron_Radiation,wep);
	    h_e_phi_corr_binSector_binTheta[sector_e-1][binX(bE_Theta,theta_e)]->Fill(phi_e-shift_e,Delta_Eprime-Electron_Radiation,wep);
	  }
	  if(protons[0]->getRegion()==FD){
	    int sector_p = protons[0]->getSector();
	    double shift_p = 0;
	    shift_p += (sector_p==1)?0:(sector_p==2)?60:(sector_p==3)?120:(sector_p==4 && phi_p>0)?180:(sector_p==4 && phi_p<0)?-180:(sector_p==5)?-120:(sector_p==6)?-60:0;

	    h_Delta_Int_pMomFD->Fill(Delta_Mom_FromAngle,wep);
	    h_Delta_Int_pMomFD_Corrected->Fill(Delta_Mom_FromAngle_corrected,wep);
	    h_Delta_Res_pMomFD_Corrected->Fill(Delta_Mom_FromAngle_corrected/mom_p_corrected,wep);
	    if(binX(bE_ThetapFD,theta_p)!=-1){
	      h_pFD_Res[binX(bE_ThetapFD,theta_p)]->Fill(Delta_Mom_FromAngle,wep);
	      h_pFD_Res_Corrected[binX(bE_ThetapFD,theta_p)]->Fill(Delta_Mom_FromAngle_corrected,wep);
	    }


	    if(binX(bE_Theta,theta_p)!=-1){
	      h_phi_corr_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(phi_p-shift_p,Delta_Mom_FromAngle,wep);
	      h_p_phi_corr_binSector_binTheta[sector_p-1][binX(bE_Theta,theta_p)]->Fill(phi_p-shift_p,Delta_Mom_FromAngle,wep);
	    }
	  }
	  else if(protons[0]->getRegion()==CD){
	    h_Delta_Int_pMomCD->Fill(Delta_Mom_FromAngle,wep);
	    h_Delta_Int_pMomCD_Corrected->Fill(Delta_Mom_FromAngle_corrected,wep);
	    h_Delta_Res_pMomCD_Corrected->Fill(Delta_Mom_FromAngle_corrected/mom_p_corrected,wep);
	    if(binX(bE_ThetaCD,theta_p)!=-1){
	      h_phi_corr_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(phi_p,Delta_Mom_FromAngle,wep);
	      h_aftercorr_phi_corr_binThetaCD[binX(bE_ThetaCD,theta_p)]->Fill(phi_p,Delta_Mom_FromAngle_corrected,wep);
	    }
	    if(binX(bE_MomCD,proton_ptr.P())!=-1){
	      h_phi_corr_binMomCD[binX(bE_MomCD,proton_ptr.P())]->Fill(phi_p,Delta_Mom_FromAngle,wep);
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

  ///////////////////////////////////
  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  g_rad->Draw();
  f_rad->Draw("SAME");
  text.DrawLatex(11,28,"#Delta E' = A + B #theta + C #theta^2");
  text.DrawLatex(11,26,Form("A=%g",p_rad->Parameter(0)));
  text.DrawLatex(11,24,Form("B=%g",p_rad->Parameter(1)));
  text.DrawLatex(11,22,Form("C=%g",p_rad->Parameter(2)));
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  
  
  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Delta_Int_Eprime_Corrected->SetLineColor(2);
  h_Delta_Int_Eprime_Corrected->Draw();
  h_Delta_Int_Eprime->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_Eprime_beforeRad->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_pMomFD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_pMomCD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  g_rad_new->Draw();
  f_rad_new->Draw("SAME");
  text.DrawLatex(11,22,"#Delta E' = A + B #theta + C #theta^2");
  text.DrawLatex(11,20,Form("A=%g",p_rad_new->Parameter(0)));
  text.DrawLatex(11,18,Form("B=%g",p_rad_new->Parameter(1)));
  text.DrawLatex(11,16,Form("C=%g",p_rad_new->Parameter(2)));
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_Eprime->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_phie_Dp_FD->Draw("colz");
  myCanvas->cd(2);
  myCanvas->GetPad(2)->SetLeftMargin(0.15);
  h_phip_Dp_FD->Draw("colz");
  myCanvas->cd(3);
  myCanvas->GetPad(3)->SetLeftMargin(0.15);
  h_phie_Dp_CD->Draw("colz");
  myCanvas->cd(4);
  myCanvas->GetPad(4)->SetLeftMargin(0.15);
  h_phip_Dp_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  /*
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_thetae_Dp_FD->Draw("colz");
  myCanvas->cd(2);
  myCanvas->GetPad(2)->SetLeftMargin(0.15);
  h_thetap_Dp_FD->Draw("colz");
  myCanvas->cd(3);
  myCanvas->GetPad(3)->SetLeftMargin(0.15);
  h_thetae_Dp_CD->Draw("colz");
  myCanvas->cd(4);
  myCanvas->GetPad(4)->SetLeftMargin(0.15);
  h_thetap_Dp_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  TGraph * g_mode = new TGraph;
  TGraph * g_sigma = new TGraph;
  TGraph * g_mode_corr = new TGraph;
  TGraph * g_sigma_corr = new TGraph;
  myCanvas->Divide(3,4);
  for(int i = 0; i < 11; i++){
    double mode, sigma, mode_corr, sigma_corr;
    
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    TF1 * f_thetabin_corr = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    getMax(h_E_Res[i],f_thetabin,mode,sigma);
    getMax(h_E_Res_Corrected[i],f_thetabin_corr,mode_corr,sigma_corr);    

    myCanvas->cd(i+1);
    h_E_Res_Corrected[i]->SetLineColor(2);
    h_E_Res_Corrected[i]->Draw();
    f_thetabin_corr->SetLineColor(2);
    f_thetabin_corr->Draw("SAME");
    h_E_Res[i]->SetLineColor(4);
    h_E_Res[i]->Draw("SAME");      
    f_thetabin->SetLineColor(4);
    f_thetabin->Draw("SAME");

    double x = (bE_Theta[i]+bE_Theta[i+1])/2;
    g_mode->SetPoint(g_mode->GetN(),x,1000*mode);
    g_mode_corr->SetPoint(g_mode_corr->GetN(),x,1000*mode_corr);
    g_sigma->SetPoint(g_sigma->GetN(),x,1000*sigma);
    g_sigma_corr->SetPoint(g_sigma_corr->GetN(),x,1000*sigma_corr);
    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  double x_ab1[2] = {10,35};
  double y_ab1[2] = {0,80};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetLineColor(0);
  r_ab1->SetTitle("#mu_{#Delta p_{e}} vs. #theta_{e} (e,e'p);#theta_{e}^{#circ};#mu_{#Delta p_{e}} [MeV]");
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  r_ab1->Draw();
  g_mode_corr->SetLineColor(2);
  g_mode_corr->Draw("SAME");
  g_mode->SetLineColor(4);
  g_mode->Draw("SAME");
  //g_rad->Draw("SAME");
  g_rad_new->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  /*
  double x_ab2[2] = {10,35};
  double y_ab2[2] = {0,60};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetLineColor(0);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab2->Draw();
  g_sigma_corr->SetLineColor(2);
  g_sigma_corr->Draw();
  g_sigma->SetLineColor(4);
  g_sigma->Draw("SAME");
  g_rade_new->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_Eprime_Corrected->SetLineColor(2);
  h_Delta_Int_Eprime_Corrected->Draw();
  h_Delta_Int_Eprime_beforeRad->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_pMomFD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomFD_Corrected->Draw();
  h_Delta_Int_pMomFD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_Delta_Int_pMomCD_Corrected->SetLineColor(2);
  h_Delta_Int_pMomCD_Corrected->Draw();
  h_Delta_Int_pMomCD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_phie_Dp_Corr_FD->Draw("colz");
  myCanvas->cd(2);
  myCanvas->GetPad(2)->SetLeftMargin(0.15);
  h_phip_Dp_Corr_FD->Draw("colz");
  myCanvas->cd(3);
  myCanvas->GetPad(3)->SetLeftMargin(0.15);
  h_phie_Dp_Corr_CD->Draw("colz");
  myCanvas->cd(4);
  myCanvas->GetPad(4)->SetLeftMargin(0.15);
  h_phip_Dp_Corr_CD->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    h_phi_corr_binThetaCD[i]->Draw("colz");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    h_phi_corr_binThetaCD[i]->Draw("colz");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    myCanvas->cd(i+1);
    h_aftercorr_phi_corr_binThetaCD[i]->Draw("colz");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    myCanvas->cd(i-8);
    h_aftercorr_phi_corr_binThetaCD[i]->Draw("colz");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();


  return 0;
}

  /*

  */  
  /*
  TGraph * g_pFD_mode = new TGraph;
  TGraph * g_pFD_sigma = new TGraph;
  TGraph * g_pFD_mode_corr = new TGraph;
  TGraph * g_pFD_sigma_corr = new TGraph;
  myCanvas->Divide(3,3);
  for(int i = 1; i < 8; i++){
    double mode, sigma, mode_corr, sigma_corr;

    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    TF1 * f_thetabin_corr = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    getMax(h_pFD_Res[i],f_thetabin,mode,sigma);
    getMax(h_pFD_Res_Corrected[i],f_thetabin_corr,mode_corr,sigma_corr);    

    myCanvas->cd(i+1);
    h_pFD_Res_Corrected[i]->SetLineColor(2);
    h_pFD_Res_Corrected[i]->Draw();
    f_thetabin_corr->SetLineColor(2);
    f_thetabin_corr->Draw("SAME");
    h_pFD_Res[i]->SetLineColor(4);
    h_pFD_Res[i]->Draw("SAME");      
    f_thetabin->SetLineColor(4);
    f_thetabin->Draw("SAME");
    
    g_pFD_mode->SetPoint(g_pFD_mode->GetN(),x_rad[i],1000*mode);
    g_pFD_mode_corr->SetPoint(g_pFD_mode_corr->GetN(),x_rad[i],1000*mode_corr);
    g_pFD_sigma->SetPoint(g_pFD_sigma->GetN(),x_rad[i],1000*sigma);
    g_pFD_sigma_corr->SetPoint(g_pFD_sigma_corr->GetN(),x_rad[i],1000*sigma_corr);
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  double x_ab3[2] = {22,45};
  double y_ab3[2] = {-60,60};  
  TGraph * r_ab3 = new TGraph(2,x_ab3,y_ab3);
  r_ab3->SetLineColor(0);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab3->Draw();
  g_pFD_mode_corr->SetLineColor(2);
  g_pFD_mode_corr->Draw("SAME");
  g_pFD_mode->SetLineColor(4);
  g_pFD_mode->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  double x_ab4[2] = {22,45};
  double y_ab4[2] = {0,60};  
  TGraph * r_ab4 = new TGraph(2,x_ab4,y_ab4);
  r_ab4->SetLineColor(0);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab4->Draw();
  g_pFD_sigma_corr->SetLineColor(2);
  g_pFD_sigma_corr->Draw();
  g_pFD_sigma->SetLineColor(4);
  g_pFD_sigma->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  

  */
  /*
  myCanvas->cd(1);
  h_Delta_Res_Eprime_Corrected->Draw();
  TF1 * f_Delta_Res_Eprime_Corrected = new TF1("Delta_Res_Eprime",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.1,0.1,3);
  f_Delta_Res_Eprime_Corrected->SetParameter(0,1000);
  f_Delta_Res_Eprime_Corrected->SetParameter(1,0);
  f_Delta_Res_Eprime_Corrected->SetParameter(2,0.005);
  TFitResultPtr p_Delta_Res_Eprime_Corrected = h_Delta_Res_Eprime_Corrected->Fit(f_Delta_Res_Eprime_Corrected,"qeSrn","",-0.02,0.02);
  text.DrawLatex(-0.04,h_Delta_Res_Eprime_Corrected->GetMaximum()*0.9,Form("#mu = %g %%",p_Delta_Res_Eprime_Corrected->Parameter(1)*100));
  text.DrawLatex(-0.04,h_Delta_Res_Eprime_Corrected->GetMaximum()*0.8,Form("#sigma = %g %%",p_Delta_Res_Eprime_Corrected->Parameter(2)*100));
  */
  /*
  text.DrawLatex(-0.18,h_Delta_Res_Eprime_Corrected->GetMaximum()*0.9,Form("#mu = %g #pm %g MeV",p_Delta_Res_Eprime_Corrected->Parameter(1)*1000,p_Delta_Res_Eprime_Corrected->ParError(1)*1000));
  text.DrawLatex(-0.18,h_Delta_Res_Eprime_Corrected->GetMaximum()*0.8,Form("#sigma = %g #pm %g MeV",p_Delta_Res_Eprime_Corrected->Parameter(2)*1000,p_Delta_Res_Eprime_Corrected->ParError(2)*1000));
  */

  /*
  f_Delta_Res_Eprime_Corrected->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->cd(1);
  h_Delta_Res_pMomFD_Corrected->Draw();
  TF1 * f_Delta_Res_pMomFD_Corrected = new TF1("Delta_Res_pMomFD",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
  f_Delta_Res_pMomFD_Corrected->SetParameter(0,1000);
  f_Delta_Res_pMomFD_Corrected->SetParameter(1,0);
  f_Delta_Res_pMomFD_Corrected->SetParameter(2,0.05);
  TFitResultPtr p_Delta_Res_pMomFD_Corrected = h_Delta_Res_pMomFD_Corrected->Fit(f_Delta_Res_pMomFD_Corrected,"qeSrn","",-0.1,0.05);
  text.DrawLatex(-0.18,h_Delta_Res_pMomFD_Corrected->GetMaximum()*0.9,Form("#mu = %g %%",p_Delta_Res_pMomFD_Corrected->Parameter(1)*100));
  text.DrawLatex(-0.18,h_Delta_Res_pMomFD_Corrected->GetMaximum()*0.8,Form("#sigma = %g %%",p_Delta_Res_pMomFD_Corrected->Parameter(2)*100));
  */
  /*
  text.DrawLatex(-0.38,h_Delta_Res_pMomFD_Corrected->GetMaximum()*0.9,Form("#mu = %g #pm %g MeV",p_Delta_Res_pMomFD_Corrected->Parameter(1)*1000,p_Delta_Res_pMomFD_Corrected->ParError(1)*1000));
  text.DrawLatex(-0.38,h_Delta_Res_pMomFD_Corrected->GetMaximum()*0.8,Form("#sigma = %g #pm %g MeV",p_Delta_Res_pMomFD_Corrected->Parameter(2)*1000,p_Delta_Res_pMomFD_Corrected->ParError(2)*1000));
  */
  /*
  f_Delta_Res_pMomFD_Corrected->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */

  /*
  myCanvas->cd(1);
  h_Delta_Res_pMomCD_Corrected->Draw();
  TF1 * f_Delta_Res_pMomCD_Corrected = new TF1("Delta_Res_pMomCD",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
  f_Delta_Res_pMomCD_Corrected->SetParameter(0,1000);
  f_Delta_Res_pMomCD_Corrected->SetParameter(1,0);
  f_Delta_Res_pMomCD_Corrected->SetParameter(2,0.05);
  TFitResultPtr p_Delta_Res_pMomCD_Corrected = h_Delta_Res_pMomCD_Corrected->Fit(f_Delta_Res_pMomCD_Corrected,"qeSrn","",-0.2,0.1);
  text.DrawLatex(-0.38,h_Delta_Res_pMomCD_Corrected->GetMaximum()*0.9,Form("#mu = %g %%",p_Delta_Res_pMomCD_Corrected->Parameter(1)*100));
  text.DrawLatex(-0.38,h_Delta_Res_pMomCD_Corrected->GetMaximum()*0.8,Form("#sigma = %g %%",p_Delta_Res_pMomCD_Corrected->Parameter(2)*100));
  */
  /*
  text.DrawLatex(-0.38,h_Delta_Res_pMomCD_Corrected->GetMaximum()*0.9,Form("#mu = %g #pm %g MeV",p_Delta_Res_pMomCD_Corrected->Parameter(1)*1000,p_Delta_Res_pMomCD_Corrected->ParError(1)*1000));
  text.DrawLatex(-0.38,h_Delta_Res_pMomCD_Corrected->GetMaximum()*0.8,Form("#sigma = %g #pm %g MeV",p_Delta_Res_pMomCD_Corrected->Parameter(2)*1000,p_Delta_Res_pMomCD_Corrected->ParError(2)*1000));
  */
  /*
  f_Delta_Res_pMomCD_Corrected->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */

  /*
  for(int j = 0; j < 6; j++){
    myCanvas->Divide(4,4);
    for(int i = 0; i < 14; i++){
      myCanvas->cd(i+1);
      h_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      getGraph(h_phi_corr_binSector_binTheta[j][i],g_phi_corr_binSector_binTheta[j][i]);
      g_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  for(int j = 0; j < 6; j++){
    myCanvas->Divide(4,4);
    for(int i = 0; i < 14; i++){
      myCanvas->cd(i+1);
      h_e_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      getGraph(h_e_phi_corr_binSector_binTheta[j][i],g_e_phi_corr_binSector_binTheta[j][i]);
      g_e_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
    }    
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  for(int j = 0; j < 6; j++){
    myCanvas->Divide(4,4);
    for(int i = 0; i < 14; i++){
      myCanvas->cd(i+1);
      h_p_phi_corr_binSector_binTheta[j][i]->Draw("colz");
      getGraph(h_p_phi_corr_binSector_binTheta[j][i],g_p_phi_corr_binSector_binTheta[j][i]);
      g_p_phi_corr_binSector_binTheta[j][i]->Draw("SAME");
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
  for(int i = 6; i < 12; i++){
    myCanvas->cd(i-5);
    h_phi_corr_binThetaCD[i]->Draw("colz");
    getGraph(h_phi_corr_binThetaCD[i],g_phi_corr_binThetaCD[i]);
    g_phi_corr_binThetaCD[i]->Draw("SAME");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 12; i < 18; i++){
    myCanvas->cd(i-11);
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

  myCanvas->Divide(2,3);
  for(int i = 6; i < 12; i++){
    myCanvas->cd(i-5);
    h_phi_corr_binMomCD[i]->Draw("colz");
    getGraph(h_phi_corr_binMomCD[i],g_phi_corr_binMomCD[i]);
    g_phi_corr_binMomCD[i]->Draw("SAME");
  }    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
