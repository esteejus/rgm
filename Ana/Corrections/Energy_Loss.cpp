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
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "clas12ana.h"
#include "reweighter.h"
#include "TGraphErrors.h"


using namespace std;
using namespace clas12;

const double c = 29.9792458;

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mass_pi = db->GetParticle(-211)->Mass();
double mD = 1.8756;

double beam_E = 5.98;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double func(double x, double a, double b, double c){
  return a + b*x + (c/x); 
}

double funcAB(double x, double a, double b){
  return a + b*x; 
}

double funcC(double x, double a, double b, double c){
  return a - b/(x-c); 
}

double getExp(TLorentzVector balance_ptr, TLorentzVector par){
  double theta_bpar = balance_ptr.Vect().Angle(par.Vect());
  double Eb = balance_ptr.E();
  double Pb = balance_ptr.P();
  double K = ((mass_pi*mass_pi) - balance_ptr.M2() - par.M2()) / 2;
  double a = Pb*Pb*cos(theta_bpar)*cos(theta_bpar) - Eb*Eb;
  double b = -2 * K * Pb * cos(theta_bpar);
  double c = K*K - Eb*Eb*par.M2();
  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
  return x_min;
}

void getGraph(TH2D * h_myhist, TGraphErrors * g_mygraph, TCanvas * myCanvas, char fileName[100]){
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
    if(proj->GetEntries()<15){continue;}

    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.3,0.3,3);
    double mode = proj->GetBinCenter(proj->GetMaximumBin());
    gFit->SetParameter(0,proj->GetMaximum()/G(0,1,0,0.1));
    gFit->SetParameter(1,mode);
    gFit->SetParLimits(1,mode-0.025,mode+0.025);
    gFit->SetParameter(2,0.01);
    gFit->SetParLimits(2,0.001,0.1);
    
    TFitResultPtr gPoint = proj->Fit(gFit,"SrBeqn","",mode-0.05,0.05);

    myCanvas->Divide(1,1);
    myCanvas->cd(1);    
    proj->Draw();
    gFit->Draw("SAME");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  


    if(gPoint == 0){
      g_mygraph->SetPoint(g_mygraph->GetN(),x,gPoint->Parameter(1));
      g_mygraph->SetPointError(g_mygraph->GetN()-1,0,gPoint->Parameter(2));
    }
    proj->Write();
  }
}

void getFunctionMomTFRP(TH2D * h_myhist, TGraphErrors * g_mygraph, TF1 * f_myfunc, TFitResultPtr & p_mypoint, double min, double max,TCanvas * myCanvas, char fileName[100]){
  getGraph(h_myhist,g_mygraph,myCanvas,fileName);

  f_myfunc->SetLineColor(3);
  f_myfunc->SetLineWidth(1);
  f_myfunc->SetParameter(0,0);
  f_myfunc->SetParLimits(0,-0.2,+0.2);
  f_myfunc->SetParameter(1,0);
  f_myfunc->SetParLimits(1,0,+0.01);
  f_myfunc->SetParameter(2,0);
  f_myfunc->SetParLimits(2,-0.5,+0.05);
  p_mypoint = g_mygraph->Fit(f_myfunc,"SrBeqn","",min,max);
      
}

void getABC(TH2D * h_myhist[8], TGraphErrors * g_mygraph[8], TF1 * f_myfunc[8], TFitResultPtr p_mypoint[8], double min, double max, TGraph * g_Pargraph[3], TF1 * f_Parfunc[3], TFitResultPtr p_Parpoint[3], TF1 * f_Combfunc[8], TCanvas * myCanvas, char fileName[100]){

  for(int i=0; i<8; i++){
    double theta = 5.5 + i*5;
    getFunctionMomTFRP(h_myhist[i],g_mygraph[i],f_myfunc[i],p_mypoint[i],min,max,myCanvas,fileName);
    for(int j=0; j<3; j++){
      g_Pargraph[j]->SetPoint(g_Pargraph[j]->GetN(),theta,p_mypoint[i]->Parameter(j));
    }
  }
  for(int j=0; j<2; j++){
    f_Parfunc[j]->SetLineColor(4);
    f_Parfunc[j]->SetParameter(0,0);
    f_Parfunc[j]->SetParLimits(0,-0.01,0.01);
    f_Parfunc[j]->SetParameter(1,0);
    f_Parfunc[j]->SetParLimits(1,-0.01,0.01);
    p_Parpoint[j] = g_Pargraph[j]->Fit(f_Parfunc[j],"SrBeqn","",5,41);
  }  
  f_Parfunc[2]->SetLineColor(4);
  f_Parfunc[2]->SetParameter(0,0);
  f_Parfunc[2]->SetParLimits(0,-0.1,0.1);
  f_Parfunc[2]->SetParameter(1,0.01);
  f_Parfunc[2]->SetParLimits(1,0,10);
  f_Parfunc[2]->SetParameter(2,60);
  f_Parfunc[2]->SetParLimits(2,45,100);
  p_Parpoint[2] = g_Pargraph[2]->Fit(f_Parfunc[2],"SrBeqn","",5,41);

  for(int i=0; i<8; i++){
    double theta = 5.5 + i*5;
    f_Combfunc[i]->SetLineColor(4);
    f_Combfunc[i]->SetLineWidth(1);
    f_Combfunc[i]->SetParameter(0,f_Parfunc[0]->Eval(theta));
    f_Combfunc[i]->SetParameter(1,f_Parfunc[1]->Eval(theta));
    f_Combfunc[i]->SetParameter(2,f_Parfunc[2]->Eval(theta));
  }

  cout<<"("<<p_Parpoint[0]->Parameter(0)<<" + "<<p_Parpoint[0]->Parameter(1)<<"*theta_1p)"<<endl
      <<"+("<<p_Parpoint[1]->Parameter(0)<<" + "<<p_Parpoint[1]->Parameter(1)<<"*theta_1p)*mom_1p"<<endl
      <<"+("<<p_Parpoint[2]->Parameter(0)<<" - "<<p_Parpoint[2]->Parameter(1)<<"/(theta_1p-"<<p_Parpoint[2]->Parameter(2)<<"))/mom_1p"<<endl;

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

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
  double mN = db->GetParticle(2212)->Mass();
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target_ptr(0,0,0,mD);
  //TLorentzVector target_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector proton_ptr_p1(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector proton_ptr_p2(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pim_ptr(0,0,0,db->GetParticle(-211)->Mass());
  reweighter newWeight(beam_E,6,6,kelly,"AV18");

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;

  //////////////////////////////////
  //Electron
  //////////////////////////////////
  TH1D * h_eDeltaP = new TH1D("eDeltaP","eDeltaP",100,-0.2,0.2);
  hist_list.push_back(h_eDeltaP);

  TH2D * h_emom_eDeltaP_int = new TH2D("h_emom_eDeltaP_int","h_emom_eDeltaP_int",50,0,6,50,-0.2,0.2);
  hist_list.push_back(h_emom_eDeltaP_int);
  TGraphErrors * g_emom_eDeltaP_int = new TGraphErrors();
  g_emom_eDeltaP_int->SetName("g_emom_eDeltaP_sector_int");
  g_emom_eDeltaP_int->SetLineColor(2);
  TF1 * f_emom_eDeltaP_int = new TF1("emom_eDeltaP_int",[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);

  TH2D * h_emom_eDeltaP[6];
  TGraphErrors * g_emom_eDeltaP[6];
  TF1 * f_emom_eDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_emom_eDeltaP_sector_%d",j+1);
    h_emom_eDeltaP[j] = new TH2D(temp_name,temp_name,50,0,6,50,-0.2,0.2);
    hist_list.push_back(h_emom_eDeltaP[j]);
    g_emom_eDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_emom_eDeltaP_sector_%d",j+1);
    g_emom_eDeltaP[j]->SetName(temp_name);
    g_emom_eDeltaP[j]->SetLineColor(2);
    sprintf(temp_name,"f_emom_eDeltaP_sector_%d",j+1);
    f_emom_eDeltaP[j] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
  }

  TH2D * h_emom_eDeltaP_phiGroup[6][3];
  TGraphErrors * g_emom_eDeltaP_phiGroup[6][3];
  for(int j=0; j<6; j++){
    for(int i=0; i<3; i++){
      sprintf(temp_name,"h_emom_eDeltaP_sector_%d_phiGroup_%d",j+1,i-1);
      h_emom_eDeltaP_phiGroup[j][i] = new TH2D(temp_name,temp_name,50,0,6,50,-0.2,0.2);
      hist_list.push_back(h_emom_eDeltaP_phiGroup[j][i]);
      g_emom_eDeltaP_phiGroup[j][i] = new TGraphErrors();
      sprintf(temp_name,"g_emom_eDeltaP_sector_%d_phiGroup_%d",j+1,i-1);
      g_emom_eDeltaP_phiGroup[j][i]->SetName(temp_name);
      g_emom_eDeltaP_phiGroup[j][i]->SetLineColor(2);
    }
  }

  TH2D * h_emom_eDeltaP_thetaGroup[6][8];
  TGraphErrors * g_emom_eDeltaP_thetaGroup[6][8];
  TF1 * f_emom_eDeltaP_thetaGroup[6][8];
  for(int j=0; j<6; j++){
    for(int i=0; i<8; i++){
      sprintf(temp_name,"h_emom_eDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      h_emom_eDeltaP_thetaGroup[j][i] = new TH2D(temp_name,temp_name,25,0,6,50,-0.2,0.2);
      hist_list.push_back(h_emom_eDeltaP_thetaGroup[j][i]);
      g_emom_eDeltaP_thetaGroup[j][i] = new TGraphErrors();
      sprintf(temp_name,"g_emom_eDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      g_emom_eDeltaP_thetaGroup[j][i]->SetName(temp_name);
      g_emom_eDeltaP_thetaGroup[j][i]->SetLineColor(2);
      sprintf(temp_name,"f_emom_eDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      f_emom_eDeltaP_thetaGroup[j][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    }
  }

  TH2D * h_emom_eDeltaP_int_thetaGroup[8];
  TGraphErrors * g_emom_eDeltaP_int_thetaGroup[8];
  TF1 * f_emom_eDeltaP_int_thetaGroup[8];
  for(int i=0; i<8; i++){
    sprintf(temp_name,"h_emom_eDeltaP_int_thetaGroup_%d",i+1);
    h_emom_eDeltaP_int_thetaGroup[i] = new TH2D(temp_name,temp_name,25,0,6,100,-0.2,0.2);
    hist_list.push_back(h_emom_eDeltaP_int_thetaGroup[i]);
    g_emom_eDeltaP_int_thetaGroup[i] = new TGraphErrors();
    sprintf(temp_name,"g_emom_eDeltaP_int_thetaGroup_%d",i+1);
    g_emom_eDeltaP_int_thetaGroup[i]->SetName(temp_name);
    g_emom_eDeltaP_int_thetaGroup[i]->SetLineColor(2);
    sprintf(temp_name,"f_emom_eDeltaP_int_thetaGroup_%d",i+1);
    f_emom_eDeltaP_int_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
  }  

  TH2D * h_etheta_eDeltaP[6];
  TGraphErrors * g_etheta_eDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_etheta_eDeltaP_sector_%d",j+1);
    h_etheta_eDeltaP[j] = new TH2D(temp_name,temp_name,50,0,40,50,-0.2,0.2);
    hist_list.push_back(h_etheta_eDeltaP[j]);
    g_etheta_eDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_etheta_eDeltaP_sector_%d",j+1);
    g_etheta_eDeltaP[j]->SetName(temp_name);
    g_etheta_eDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_ephi_eDeltaP[6];
  TGraphErrors * g_ephi_eDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_ephi_eDeltaP_sector_%d",j+1);
    h_ephi_eDeltaP[j] = new TH2D(temp_name,temp_name,50,-40,40,50,-0.2,0.2);
    hist_list.push_back(h_ephi_eDeltaP[j]);
    g_ephi_eDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_ephi_eDeltaP_sector_%d",j+1);
    g_ephi_eDeltaP[j]->SetName(temp_name);
    g_ephi_eDeltaP[j]->SetLineColor(2);
  }


  //////////////////////////////////
  //Proton FD
  //////////////////////////////////
  TH2D * h_pFDmom_pFDDeltaP_int = new TH2D("h_pFDmom_pFDDeltaP_int","h_pFDmom_pFDDeltaP_int",50,0,3.5,50,-0.2,0.2);
  hist_list.push_back(h_pFDmom_pFDDeltaP_int);
  TGraphErrors * g_pFDmom_pFDDeltaP_int = new TGraphErrors();
  g_pFDmom_pFDDeltaP_int->SetName("g_pFDmom_pFDDeltaP_sector_int");
  g_pFDmom_pFDDeltaP_int->SetLineColor(2);
  TF1 * f_pFDmom_pFDDeltaP_int = new TF1("emom_pFDDeltaP_int",[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);

  TH2D * h_pFDmom_pFDDeltaP[6];
  TGraphErrors * g_pFDmom_pFDDeltaP[6];
  TF1 * f_pFDmom_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d",j+1);
    h_pFDmom_pFDDeltaP[j] = new TH2D(temp_name,temp_name,50,0,4.5,50,-0.2,0.2);
    hist_list.push_back(h_pFDmom_pFDDeltaP[j]);
    g_pFDmom_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d",j+1);
    g_pFDmom_pFDDeltaP[j]->SetName(temp_name);
    g_pFDmom_pFDDeltaP[j]->SetLineColor(2);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d",j+1);
    f_pFDmom_pFDDeltaP[j] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
  }

  TH2D * h_pFDmom_pFDDeltaP_phiGroup[6][3];
  TGraphErrors * g_pFDmom_pFDDeltaP_phiGroup[6][3];
  for(int j=0; j<6; j++){
    for(int i=0; i<3; i++){
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_phiGroup_%d",j+1,i-1);
      h_pFDmom_pFDDeltaP_phiGroup[j][i] = new TH2D(temp_name,temp_name,50,0,4.5,50,-0.2,0.2);
      hist_list.push_back(h_pFDmom_pFDDeltaP_phiGroup[j][i]);
      g_pFDmom_pFDDeltaP_phiGroup[j][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_phiGroup_%d",j+1,i-1);
      g_pFDmom_pFDDeltaP_phiGroup[j][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_phiGroup[j][i]->SetLineColor(2);
    }
  }

  TH2D * h_pFDmom_pFDDeltaP_thetaGroup[6][8];
  TGraphErrors * g_pFDmom_pFDDeltaP_thetaGroup[6][8];
  TF1 * f_pFDmom_pFDDeltaP_thetaGroup[6][8];
  for(int j=0; j<6; j++){
    for(int i=0; i<8; i++){
      sprintf(temp_name,"h_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      h_pFDmom_pFDDeltaP_thetaGroup[j][i] = new TH2D(temp_name,temp_name,50,0,4.5,50,-0.2,0.2);
      hist_list.push_back(h_pFDmom_pFDDeltaP_thetaGroup[j][i]);
      g_pFDmom_pFDDeltaP_thetaGroup[j][i] = new TGraphErrors();
      sprintf(temp_name,"g_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      g_pFDmom_pFDDeltaP_thetaGroup[j][i]->SetName(temp_name);
      g_pFDmom_pFDDeltaP_thetaGroup[j][i]->SetLineColor(2);
      sprintf(temp_name,"f_pFDmom_pFDDeltaP_sector_%d_thetaGroup_%d",j+1,i+1);
      f_pFDmom_pFDDeltaP_thetaGroup[j][i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);

    }
  }

  TH2D * h_pFDmom_pFDDeltaP_int_thetaGroup[8];
  TGraphErrors * g_pFDmom_pFDDeltaP_int_thetaGroup[8];
  TF1 * f_pFDmom_pFDDeltaP_int_thetaGroup[8];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_thetaGroup[8];
  TGraph * g_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_Pars[3];
  TFitResultPtr p_pFDmom_pFDDeltaP_int_Pars[3];
  TF1 * f_pFDmom_pFDDeltaP_int_combined_thetaGroup[8];
  TH2D * h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[8];
  for(int i=0; i<8; i++){
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TH2D(temp_name,temp_name,25,0,4.5,200,-0.2,0.2);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_thetaGroup[i]);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TGraphErrors();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->SetLineColor(2);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_combined_thetaGroup_%d",i+1);
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
    sprintf(temp_name,"h_pFDmom_pFDDeltaP_int_corrected_thetaGroup_%d",i+1);
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i] = new TH2D(temp_name,temp_name,25,0,4.5,100,-0.2,0.2);
    hist_list.push_back(h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]);
  }
  for(int j=0; j<3; j++){
    g_pFDmom_pFDDeltaP_int_Pars[j] = new TGraph();
    sprintf(temp_name,"g_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetName(temp_name);
    g_pFDmom_pFDDeltaP_int_Pars[j]->SetLineColor(3);
    sprintf(temp_name,"f_pFDmom_pFDDeltaP_int_Pars_%d",j+1);
    if(j<2){
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcAB(x[0],p[0],p[1]); },5,41,2);
    }
    else{
      f_pFDmom_pFDDeltaP_int_Pars[j] = new TF1(temp_name,[&](double *x, double *p){ return funcC(x[0],p[0],p[1],p[2]); },5,41,3);
    }
  }
  

  TH2D * h_pFDtheta_pFDDeltaP[6];
  TGraphErrors * g_pFDtheta_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDtheta_pFDDeltaP_sector_%d",j+1);
    h_pFDtheta_pFDDeltaP[j] = new TH2D(temp_name,temp_name,50,0,50,50,-0.2,0.2);
    hist_list.push_back(h_pFDtheta_pFDDeltaP[j]);
    g_pFDtheta_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDtheta_pFDDeltaP_sector_%d",j+1);
    g_pFDtheta_pFDDeltaP[j]->SetName(temp_name);
    g_pFDtheta_pFDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pFDphi_pFDDeltaP[6];
  TGraphErrors * g_pFDphi_pFDDeltaP[6];
  for(int j=0; j<6; j++){
    sprintf(temp_name,"h_pFDphi_pFDDeltaP_sector_%d",j+1);
    h_pFDphi_pFDDeltaP[j] = new TH2D(temp_name,temp_name,50,-35,35,50,-0.2,0.2);
    hist_list.push_back(h_pFDphi_pFDDeltaP[j]);
    g_pFDphi_pFDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pFDphi_pFDDeltaP_sector_%d",j+1);
    g_pFDphi_pFDDeltaP[j]->SetName(temp_name);
    g_pFDphi_pFDDeltaP[j]->SetLineColor(2);
  }


  //////////////////////////////////
  //Proton CD
  //////////////////////////////////
  TH2D * h_pCDmom_pCDDeltaP_int = new TH2D("h_pCDmom_pCDDeltaP_int","h_pCDmom_pCDDeltaP_int",50,0,2.5,50,-0.2,0.2);
  hist_list.push_back(h_pCDmom_pCDDeltaP_int);
  TGraphErrors * g_pCDmom_pCDDeltaP_int = new TGraphErrors();
  g_pCDmom_pCDDeltaP_int->SetName("g_pCDmom_pCDDeltaP_sector_int");
  g_pCDmom_pCDDeltaP_int->SetLineColor(2);
  TF1 * f_pCDmom_pCDDeltaP_int = new TF1("emom_pCDDeltaP_int",[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);

  TH2D * h_pCDtheta_pCDDeltaP_int = new TH2D("h_pCDtheta_pCDDeltaP_int","h_pCDtheta_pCDDeltaP_int",50,30,130,50,-0.2,0.2);
  hist_list.push_back(h_pCDtheta_pCDDeltaP_int);
  TGraphErrors * g_pCDtheta_pCDDeltaP_int = new TGraphErrors();
  g_pCDtheta_pCDDeltaP_int->SetName("g_pCDtheta_pCDDeltaP_sector_int");
  g_pCDtheta_pCDDeltaP_int->SetLineColor(2);

  TH2D * h_pCDmom_pCDDeltaP[18];
  TGraphErrors * g_pCDmom_pCDDeltaP[18];
  TF1 * f_pCDmom_pCDDeltaP[18];
  for(int j=0; j<18; j++){
    sprintf(temp_name,"h_pCDmom_pCDDeltaP_bin_%d",j+1);
    h_pCDmom_pCDDeltaP[j] = new TH2D(temp_name,temp_name,50,0,2,50,-0.2,0.2);
    hist_list.push_back(h_pCDmom_pCDDeltaP[j]);
    g_pCDmom_pCDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDmom_pCDDeltaP_bin_%d",j+1);
    g_pCDmom_pCDDeltaP[j]->SetName(temp_name);
    g_pCDmom_pCDDeltaP[j]->SetLineColor(2);
    sprintf(temp_name,"f_pCDmom_pCDDeltaP_sector_%d",j+1);
    f_pCDmom_pCDDeltaP[j] = new TF1(temp_name,[&](double *x, double *p){ return func(x[0],p[0],p[1],p[2]); },0.2,6,3);
  }

  TH2D * h_pCDmomT_pCDDeltaP[18];
  TGraphErrors * g_pCDmomT_pCDDeltaP[18];
  for(int j=0; j<18; j++){
    sprintf(temp_name,"h_pCDmomT_pCDDeltaP_bin_%d",j+1);
    h_pCDmomT_pCDDeltaP[j] = new TH2D(temp_name,temp_name,50,0,2,50,-0.2,0.2);
    hist_list.push_back(h_pCDmomT_pCDDeltaP[j]);
    g_pCDmomT_pCDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDmomT_pCDDeltaP_bin_%d",j+1);
    g_pCDmomT_pCDDeltaP[j]->SetName(temp_name);
    g_pCDmomT_pCDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pCDtheta_pCDDeltaP[18];
  TGraphErrors * g_pCDtheta_pCDDeltaP[18];
  for(int j=0; j<18; j++){
    sprintf(temp_name,"h_pCDtheta_pCDDeltaP_bin_%d",j+1);
    h_pCDtheta_pCDDeltaP[j] = new TH2D(temp_name,temp_name,50,30,130,50,-0.2,0.2);
    hist_list.push_back(h_pCDtheta_pCDDeltaP[j]);
    g_pCDtheta_pCDDeltaP[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDtheta_pCDDeltaP_bin_%d",j+1);
    g_pCDtheta_pCDDeltaP[j]->SetName(temp_name);
    g_pCDtheta_pCDDeltaP[j]->SetLineColor(2);
  }

  TH2D * h_pCDphi_pCDDeltaP_bin[9];
  TGraphErrors * g_pCDphi_pCDDeltaP_bin[9];
  for(int j=0; j<9; j++){
    sprintf(temp_name,"h_pCDphi_pCDDeltaP_bin_%d",j+1);
    h_pCDphi_pCDDeltaP_bin[j] = new TH2D(temp_name,temp_name,50,-180,180,50,-0.2,0.2);
    hist_list.push_back(h_pCDphi_pCDDeltaP_bin[j]);
    g_pCDphi_pCDDeltaP_bin[j] = new TGraphErrors();
    sprintf(temp_name,"g_pCDphi_pCDDeltaP_bin_%d",j+1);
    g_pCDphi_pCDDeltaP_bin[j]->SetName(temp_name);
    g_pCDphi_pCDDeltaP_bin[j]->SetLineColor(2);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next() && counter <100000000000)
    {

      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	wep = original_weight;// * newWeight.get_weight_ep(c12->mcparts());
	wepp = original_weight;// * newWeight.get_weight_epp(c12->mcparts());
      }

      //Display completed  
      counter++;
      if((counter%100000) == 0){
	cerr << "\n" <<counter/100000 <<" hundred thousand completed";
      }    
      if((counter%10000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pims = clasAna.getByPid(-211);
      if(electrons.size() == 1)
	{

	  SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );
	  double WSq = (mN*mN) - Q2 + (2*omega*mN);
	  double W = sqrt(WSq);
	  double vtz_e = electrons[0]->par()->getVz();
	  
	  int sector_e = electrons[0]->getSector()-1;
	  double mom_e = el.P();
	  double theta_e = el.Theta() * 180 / M_PI;
	  double phi_e = el.Phi() * 180 / M_PI;
	  double shift_e = 7.5;
	  shift_e += (sector_e==0)?0:(sector_e==1)?60:(sector_e==2)?120:(sector_e==3 && phi_e>0)?180:(sector_e==3 && phi_e<0)?-180:(sector_e==4)?-120:(sector_e==5)?-60:0;
	  phi_e -= shift_e;


	  if(vtz_e<-5.5){continue;}
	  if(vtz_e>0){continue;}

	  TVector3 pmc_e(c12->mcparts()->getPx(0),c12->mcparts()->getPy(0),c12->mcparts()->getPz(0));
	  TVector3 pmc_p(c12->mcparts()->getPx(1),c12->mcparts()->getPy(1),c12->mcparts()->getPz(1));

	  double eDeltaP = pmc_e.Mag() - el.P();
	  int thetaeGroup = (theta_e - 10.0)/3.0;
	  //Fill histograms
	  //Electron
	  h_eDeltaP->Fill(eDeltaP,wep);	  
	  h_emom_eDeltaP_int->Fill(mom_e,eDeltaP,wep);
	  h_emom_eDeltaP[sector_e]->Fill(mom_e,eDeltaP,wep);
	  h_emom_eDeltaP_phiGroup[sector_e][(phi_e<-5)?0:(phi_e<5)?1:2]->Fill(mom_e,eDeltaP,wep);
	  if((thetaeGroup>=0) && (thetaeGroup<8)){
	    h_emom_eDeltaP_thetaGroup[sector_e][thetaeGroup]->Fill(mom_e,eDeltaP,wep);	      
	    h_emom_eDeltaP_int_thetaGroup[thetaeGroup]->Fill(mom_e,eDeltaP,wep);
	  }
	  h_etheta_eDeltaP[sector_e]->Fill(theta_e,eDeltaP,wep);
	  h_ephi_eDeltaP[sector_e]->Fill(phi_e,eDeltaP,wep);

	  if(protons.size() != 1){continue;}
	  int index_p1 = 0;
	  SetLorentzVector(proton_ptr_p1,protons[index_p1]);
	  double beta_p1 = protons[index_p1]->par()->getBeta();
	  double mom_p1 = proton_ptr_p1.P();
	  double momT_p1 = proton_ptr_p1.Vect().Perp();
	  double theta_p1 = proton_ptr_p1.Theta() * 180 / M_PI;
	  double phi_p1 = proton_ptr_p1.Phi() * 180 / M_PI;
	  double vtz_p1 = protons[index_p1]->par()->getVz();
	  if(fabs(vtz_e-vtz_p1)>1.5){continue;}

	  /*
	  TVector3 pmc_p1;
	  if((c12->mcparts()->getPid(1)==2212) && (c12->mcparts()->getPid(2)==2212)){
	    double theta1L = pmc_pLead.Angle(proton_ptr_p1.Vect());
	    double theta1R = pmc_pRec.Angle(proton_ptr_p1.Vect());
	    if(theta1L<theta1R){
	      pmc_p1 = pmc_pLead;
	    }
	    else{
	      pmc_p1 = pmc_pRec;
	    }
	  }
	  else if((c12->mcparts()->getPid(1)==2212)){
	    pmc_p1 = pmc_pLead;
	  }
	  else if((c12->mcparts()->getPid(2)==2212)){
	    pmc_p1 = pmc_pRec;
	  }
	  else{
	    continue;
	    }*/
	  double p1DeltaP = pmc_p.Mag() - proton_ptr_p1.P();
	  
	  if(protons[index_p1]->getRegion()==FD){
	    int sector_p1 = protons[index_p1]->getSector()-1;
	    double shift_p1 = -22;
	    int thetapFDGroup = (theta_p1 - 3.0)/5.0;
	    double correction = (-0.00345606 + -0.000122651*theta_p1)
	      +(0.00125829 + 4.38197e-05*theta_p1)*mom_p1
	      +(0.000340149 - 0.675201/(theta_p1-59.4993))/mom_p1;
	    
	    shift_p1 += (sector_p1==0)?0:(sector_p1==1)?60:(sector_p1==2)?120:(sector_p1==3 && phi_p1>0)?180:(sector_p1==3 && phi_p1<0)?-180:(sector_p1==4)?-120:(sector_p1==5)?-60:0;
	    phi_p1 -= shift_p1;
	    //Proton FD
	    h_pFDmom_pFDDeltaP_int->Fill(mom_p1,p1DeltaP,wep);
	    h_pFDmom_pFDDeltaP[sector_p1]->Fill(mom_p1,p1DeltaP,wep);
	    h_pFDmom_pFDDeltaP_phiGroup[sector_p1][(phi_p1<-5)?0:(phi_p1<5)?1:2]->Fill(mom_p1,p1DeltaP,wep);
	    if((thetapFDGroup>=0) && (thetapFDGroup<8)){
	      h_pFDmom_pFDDeltaP_thetaGroup[sector_p1][thetapFDGroup]->Fill(mom_p1,p1DeltaP,wep);	      
	      h_pFDmom_pFDDeltaP_int_thetaGroup[thetapFDGroup]->Fill(mom_p1,p1DeltaP,wep);
	      h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[thetapFDGroup]->Fill(mom_p1,p1DeltaP-correction,wep);
	    }

	    h_pFDtheta_pFDDeltaP[sector_p1]->Fill(theta_p1,p1DeltaP,wep);
	    h_pFDphi_pFDDeltaP[sector_p1]->Fill(phi_p1,p1DeltaP,wep);
	    
	  }

	  if(protons[index_p1]->getRegion()==CD){
	    int bin_CD = (phi_p1+180)/20;
	    int bin_mom_CD = (mom_p1<0.35)?0:(mom_p1<0.55)?1:(mom_p1<0.75)?2:(mom_p1<0.95)?3:(mom_p1<1.15)?4:(mom_p1<1.35)?5:(mom_p1<1.55)?6:(mom_p1<1.75)?7:8;
	    //Proton CD
	    h_pCDmom_pCDDeltaP_int->Fill(mom_p1,p1DeltaP,wep);
	    h_pCDtheta_pCDDeltaP_int->Fill(theta_p1,p1DeltaP,wep);
	    h_pCDmom_pCDDeltaP[bin_CD]->Fill(mom_p1,p1DeltaP,wep);
	    h_pCDmomT_pCDDeltaP[bin_CD]->Fill(momT_p1,p1DeltaP,wep);
	    h_pCDtheta_pCDDeltaP[bin_CD]->Fill(theta_p1,p1DeltaP,wep);
	    h_pCDphi_pCDDeltaP_bin[bin_mom_CD]->Fill(phi_p1,p1DeltaP,wep);
	  }	    


	}
    }
  /*
  getFunctionMom(h_emom_eDeltaP_int,g_emom_eDeltaP_int,f_emom_eDeltaP_int,1.0,5.5);
  */
  //getFunctionMom(h_pFDmom_pFDDeltaP_int,g_pFDmom_pFDDeltaP_int,f_pFDmom_pFDDeltaP_int,0.3,3.5);
  /*
  getFunctionMom(h_pCDmom_pCDDeltaP_int,g_pCDmom_pCDDeltaP_int,f_pCDmom_pCDDeltaP_int,0.3,3.5);

  getGraph(h_pCDtheta_pCDDeltaP_int,g_pCDtheta_pCDDeltaP_int);

  for(int j=0; j<6; j++){
    getFunctionMom(h_emom_eDeltaP[j],g_emom_eDeltaP[j],f_emom_eDeltaP[j],1.0,5.5);
    for(int i=0; i<8; i++){
      getFunctionMom(h_emom_eDeltaP_thetaGroup[j][i],g_emom_eDeltaP_thetaGroup[j][i],f_emom_eDeltaP_thetaGroup[j][i],0.3,5.5);
    }
    for(int i=0; i<3; i++){
      getGraph(h_emom_eDeltaP_phiGroup[j][i],g_emom_eDeltaP_phiGroup[j][i]);
    }    
    getGraph(h_etheta_eDeltaP[j],g_etheta_eDeltaP[j]);
    getGraph(h_ephi_eDeltaP[j],g_ephi_eDeltaP[j]);
  }
  for(int i=0; i<8; i++){
    getFunctionMom(h_emom_eDeltaP_int_thetaGroup[i],g_emom_eDeltaP_int_thetaGroup[i],f_emom_eDeltaP_int_thetaGroup[i],0.3,5.5);
  }
  */
  /*
  for(int j=0; j<6; j++){
    getFunctionMom(h_pFDmom_pFDDeltaP[j],g_pFDmom_pFDDeltaP[j],f_pFDmom_pFDDeltaP[j],0.3,3.5);
    for(int i=0; i<8; i++){
      getFunctionMom(h_pFDmom_pFDDeltaP_thetaGroup[j][i],g_pFDmom_pFDDeltaP_thetaGroup[j][i],f_pFDmom_pFDDeltaP_thetaGroup[j][i],0.3,3.5);
    }
    for(int i=0; i<3; i++){
      getGraph(h_pFDmom_pFDDeltaP_phiGroup[j][i],g_pFDmom_pFDDeltaP_phiGroup[j][i]);
    }    
    getGraph(h_pFDtheta_pFDDeltaP[j],g_pFDtheta_pFDDeltaP[j]);
    getGraph(h_pFDphi_pFDDeltaP[j],g_pFDphi_pFDDeltaP[j]);
    }*/
  //for(int i=0; i<8; i++){
  //  getFunctionMomTFRP(h_pFDmom_pFDDeltaP_int_thetaGroup[i],g_pFDmom_pFDDeltaP_int_thetaGroup[i],f_pFDmom_pFDDeltaP_int_thetaGroup[i],p_pFDmom_pFDDeltaP_int_thetaGroup[i],0.3,5.5);
    //}
  /*
  for(int i=0; i<8; i++){
    getFunctionMomTFRP(h_pFDmom_pFDDeltaP_int_thetaGroup[i],g_pFDmom_pFDDeltaP_int_thetaGroup[i],f_pFDmom_pFDDeltaP_int_thetaGroup[i],p_pFDmom_pFDDeltaP_int_thetaGroup[i],0.3,5.5);    
    if(p_pFDmom_pFDDeltaP_int_thetaGroup[i]==0){
      cout<<"It is here "<<p_pFDmom_pFDDeltaP_int_thetaGroup[i]->Parameter(0)<<endl<<endl<<endl<<endl<<endl<<endl<<endl<<endl;
    }
  }
  */
  
  /*  
  for(int j=0; j<18; j++){
    getFunctionMom(h_pCDmom_pCDDeltaP[j],g_pCDmom_pCDDeltaP[j],f_pCDmom_pCDDeltaP[j],0.3,3.5);
    getGraph(h_pCDmomT_pCDDeltaP[j],g_pCDmomT_pCDDeltaP[j]);
    getGraph(h_pCDtheta_pCDDeltaP[j],g_pCDtheta_pCDDeltaP[j]);
  }
  for(int j=0; j<9; j++){
    getGraph(h_pCDphi_pCDDeltaP_bin[j],g_pCDphi_pCDDeltaP_bin[j]);
    }*/
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
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
  
  getABC(h_pFDmom_pFDDeltaP_int_thetaGroup,g_pFDmom_pFDDeltaP_int_thetaGroup,f_pFDmom_pFDDeltaP_int_thetaGroup,p_pFDmom_pFDDeltaP_int_thetaGroup,0.3,5.5,g_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_Pars,p_pFDmom_pFDDeltaP_int_Pars,f_pFDmom_pFDDeltaP_int_combined_thetaGroup,myCanvas,fileName);


  //////////////////////////////////
  //Electron
  //////////////////////////////////
  /*  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_eDeltaP->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_emom_eDeltaP_int->Draw("colz");
  g_emom_eDeltaP_int->Draw("SAME");
  f_emom_eDeltaP_int->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pFDmom_pFDDeltaP_int->Draw("colz");
  g_pFDmom_pFDDeltaP_int->Draw("SAME");
  f_pFDmom_pFDDeltaP_int->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  /*
  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pCDmom_pCDDeltaP_int->Draw("colz");
  g_pCDmom_pCDDeltaP_int->Draw("SAME");
  f_pCDmom_pCDDeltaP_int->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);    
  h_pCDtheta_pCDDeltaP_int->Draw("colz");
  g_pCDtheta_pCDDeltaP_int->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  /*myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_emom_eDeltaP[j]->Draw("colz");
    g_emom_eDeltaP[j]->Draw("SAME");
    f_emom_eDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    g_emom_eDeltaP_phiGroup[j][0]->SetLineColor(1);
    g_emom_eDeltaP_phiGroup[j][0]->Draw();
    for(int i=1; i<3; i++){
      g_emom_eDeltaP_phiGroup[j][i]->SetLineColor(i+1);
      g_emom_eDeltaP_phiGroup[j][i]->Draw("SAME");
    }  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  for(int i=0; i<8; i++){
    myCanvas->Divide(3,2);
    for(int j=0; j<6; j++){
      myCanvas->cd(j+1);    
      h_emom_eDeltaP_thetaGroup[j][i]->Draw("colz");
      g_emom_eDeltaP_thetaGroup[j][i]->Draw("SAME");
      f_emom_eDeltaP_thetaGroup[j][i]->Draw("SAME");
    }  
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_emom_eDeltaP_int_thetaGroup[i]->Draw("colz");
    g_emom_eDeltaP_int_thetaGroup[i]->Draw("SAME");
    f_emom_eDeltaP_int_thetaGroup[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_etheta_eDeltaP[j]->Draw("colz");
    g_etheta_eDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_ephi_eDeltaP[j]->Draw("colz");
    g_ephi_eDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  //////////////////////////////////
  //Proton FD
  //////////////////////////////////
  /*
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pFDmom_pFDDeltaP[j]->Draw("colz");
    g_pFDmom_pFDDeltaP[j]->Draw("SAME");
    f_pFDmom_pFDDeltaP[j]->Draw("SAME");
    }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  */
/*
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    g_pFDmom_pFDDeltaP_phiGroup[j][0]->SetLineColor(1);
    g_pFDmom_pFDDeltaP_phiGroup[j][0]->Draw();
    for(int i=1; i<3; i++){
      g_pFDmom_pFDDeltaP_phiGroup[j][i]->SetLineColor(i+1);
      g_pFDmom_pFDDeltaP_phiGroup[j][i]->Draw("SAME");
    }  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  for(int i=0; i<8; i++){
    myCanvas->Divide(3,2);
    for(int j=0; j<6; j++){
      myCanvas->cd(j+1);    
      h_pFDmom_pFDDeltaP_thetaGroup[j][i]->Draw("colz");
      g_pFDmom_pFDDeltaP_thetaGroup[j][i]->Draw("SAME");
      f_pFDmom_pFDDeltaP_thetaGroup[j][i]->Draw("SAME");
    }  
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }
  */

  myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

    myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("colz");
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  
  myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("colz");
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("SAME");
    f_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("colz");
    g_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("SAME");
    f_pFDmom_pFDDeltaP_int_thetaGroup[i]->Draw("SAME");
    f_pFDmom_pFDDeltaP_int_combined_thetaGroup[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  myCanvas->Divide(2,2);
  for(int j=0; j<3; j++){
    myCanvas->cd(j+1);    
    g_pFDmom_pFDDeltaP_int_Pars[j]->Draw();
    f_pFDmom_pFDDeltaP_int_Pars[j]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int i=0; i<8; i++){
    myCanvas->cd(i+1);    
    h_pFDmom_pFDDeltaP_int_corrected_thetaGroup[i]->Draw("colz");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  


  /*
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pFDtheta_pFDDeltaP[j]->Draw("colz");
    g_pFDtheta_pFDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pFDphi_pFDDeltaP[j]->Draw("colz");
    g_pFDphi_pFDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  //////////////////////////////////
  //Proton CD
  //////////////////////////////////
  /*  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pCDmom_pCDDeltaP[j]->Draw("colz");
    g_pCDmom_pCDDeltaP[j]->Draw("SAME");
    f_pCDmom_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=6; j<12; j++){
    myCanvas->cd(j-5);    
    h_pCDmom_pCDDeltaP[j]->Draw("colz");
    g_pCDmom_pCDDeltaP[j]->Draw("SAME");
    f_pCDmom_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=12; j<18; j++){
    myCanvas->cd(j-11);    
    h_pCDmom_pCDDeltaP[j]->Draw("colz");
    g_pCDmom_pCDDeltaP[j]->Draw("SAME");
    f_pCDmom_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  /*
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pCDmomT_pCDDeltaP[j]->Draw("colz");
    g_pCDmomT_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=6; j<12; j++){
    myCanvas->cd(j-5);    
    h_pCDmomT_pCDDeltaP[j]->Draw("colz");
    g_pCDmomT_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=12; j<18; j++){
    myCanvas->cd(j-11);    
    h_pCDmomT_pCDDeltaP[j]->Draw("colz");
    g_pCDmomT_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  /*
  myCanvas->Divide(3,2);
  for(int j=0; j<6; j++){
    myCanvas->cd(j+1);    
    h_pCDtheta_pCDDeltaP[j]->Draw("colz");
    g_pCDtheta_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=6; j<12; j++){
    myCanvas->cd(j-5);    
    h_pCDtheta_pCDDeltaP[j]->Draw("colz");
    g_pCDtheta_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,2);
  for(int j=12; j<18; j++){
    myCanvas->cd(j-11);    
    h_pCDtheta_pCDDeltaP[j]->Draw("colz");
    g_pCDtheta_pCDDeltaP[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(3,3);
  for(int j=0; j<9; j++){
    myCanvas->cd(j+1);    
    h_pCDphi_pCDDeltaP_bin[j]->Draw("colz");
    g_pCDphi_pCDDeltaP_bin[j]->Draw("SAME");
  }  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  */
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  f->Close();

  return 0;
}

/*
	  TLorentzVector balance_ptr = beam + target_ptr - proton_ptr_p1 - proton_ptr_p2;
	  double theta_be = balance_ptr.Vect().Angle(el.Vect());
	  double Eb = balance_ptr.E();
	  double Pb = balance_ptr.P();
	  double K = (pim_ptr.M2() - balance_ptr.M2() - el.M2()) / 2;
	  double a = Pb*Pb*cos(theta_be)*cos(theta_be) - Eb*Eb;
	  double b = -2 * K * Pb * cos(theta_be);
	  double c = K*K - Eb*Eb*el.M2();
	  double x_min = (-b - sqrt(b*b - 4*a*c))/(2*a);
o	  double x_max = (-b + sqrt(b*b - 4*a*c))/(2*a);
*/
