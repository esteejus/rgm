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
double getPhiDiff(TVector3 A, TVector3 B){
  TVector3 Aperp(A.X(),A.Y(),0);
  TVector3 Bperp(B.X(),B.Y(),0);
  double phidiff = (A.Angle(B))*180/M_PI;
  if(A.Cross(B).Z()<0){phidiff=phidiff*-1;}
  return phidiff;
}

//vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_MomCD = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};//{0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
vector<double> bE_ThetaCD = {35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89};
vector<double> bE_PhiCD = {-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
vector<double> bE_ThetaE = {10,13,16,19,22,25,28,31,34,37};
vector<double> bE_ThetapFD = {19,22,25,28,31,34,37,40,43,46};
double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}

double Squ(double x){
  return x*x;
}


double DoubQuad(double x, double A, double B, double C, double D, double E, double F){
  double X = Quad(x,A,B,C)*Quad(x,A,B,C) - Quad(x,D,E,F)*Quad(x,D,E,F);
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}

double SmearCD[8]={10.049,3.07365,0.461349,1.47116,-0.832429,1.45479,1.67907,3.90999};

double Trig3(double x, double A, double B, double C, double D, double E, double F, double G){
  return A + B*sin((x*1*M_PI/180)+C) + D*sin((x*2*M_PI/180)+E) + F*sin((x*3*M_PI/180)+G); 
}

double DoubTrig3(double x, double A, double B, double C, double D, double E, double F, double G, double H){
  double X = Trig3(x,A,B,C,D,E,F,G)*Trig3(x,A,B,C,D,E,F,G) - H*H;
  if(X>0.01){
    return sqrt(X); 
  }
  return 0.1;
}


int binX(vector<double> XS, double X){
  for(int i = 0; i <= XS.size(); i++){if(X<XS[i]){return i-1;}}
  return -1;
}

auto db=TDatabasePDG::Instance();
double mass_p = db->GetParticle(2212)->Mass();
double mD = 1.8756;
double beam_E = 5.984792;
double beam_E_sigma = 0.00299;
const double me = 0.000511;
const double mU = 0.9314941024;
const double m_4He = 4.00260325415 * mU - 2*me;


double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double SQ(double x){ return x*x;}

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

void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfileData.root inputfileSim.root \n\n\n";

}

int main(int argc, char ** argv)
{

  if(argc < 2)
    {
      Usage();
      return -1;
    }



  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;

  TFile * inFileData = new TFile(argv[2]);
  TFile * inFileSim = new TFile(argv[3]);
  
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

  TH1D * h_Data_e_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"ecorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Data_e_corr_binSector_binTheta[j-1][i] = (TH1D*)inFileData->Get(temp_name);
      hist_list.push_back(h_Data_e_corr_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_Data_p_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pcorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Data_p_corr_binSector_binTheta[j-1][i] = (TH1D*)inFileData->Get(temp_name);
      hist_list.push_back(h_Data_p_corr_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_Data_corr_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"corr_theta_%d",i);
    sprintf(temp_title,"#Delta p (%d< #theta < %d);#Delta p;Counts",min,max);
    h_Data_corr_binThetaCD[i] = (TH1D*)inFileData->Get(temp_name);
    hist_list.push_back(h_Data_corr_binThetaCD[i]);
  }

  TH1D * h_Data_corr_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_Data_corr_binPhiCD[i] = (TH1D*)inFileData->Get(temp_name);
    hist_list.push_back(h_Data_corr_binPhiCD[i]);
  }

  
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  TH1D * h_Sim_e_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"ecorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Sim_e_corr_binSector_binTheta[j-1][i] = (TH1D*)inFileSim->Get(temp_name);
      hist_list.push_back(h_Sim_e_corr_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_Sim_p_corr_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pcorr_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Sim_p_corr_binSector_binTheta[j-1][i] = (TH1D*)inFileSim->Get(temp_name);
      hist_list.push_back(h_Sim_p_corr_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_Sim_corr_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"corr_theta_%d",i);
    sprintf(temp_title,"#Delta p (%d< #theta < %d);#Delta p;Counts",min,max);
    h_Sim_corr_binThetaCD[i] = (TH1D*)inFileSim->Get(temp_name);
    hist_list.push_back(h_Sim_corr_binThetaCD[i]);
  }

  TH1D * h_Sim_corr_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_Sim_corr_binPhiCD[i] = (TH1D*)inFileSim->Get(temp_name);
    hist_list.push_back(h_Sim_corr_binPhiCD[i]);
  }

  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  ////////////////////////////////////////////
  TH1D * h_Sim_e_corr_Smear_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"ecorr_Smear_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Sim_e_corr_Smear_binSector_binTheta[j-1][i] = (TH1D*)inFileSim->Get(temp_name);
      hist_list.push_back(h_Sim_e_corr_Smear_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_Sim_p_corr_Smear_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pcorr_Smear_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"#Delta p Sector %d (%d< #theta < %d);#Delta p [GeV];Counts",j,min,max);
      h_Sim_p_corr_Smear_binSector_binTheta[j-1][i] = (TH1D*)inFileSim->Get(temp_name);
      hist_list.push_back(h_Sim_p_corr_Smear_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_Sim_corr_Smear_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"corr_Smear_theta_%d",i);
    sprintf(temp_title,"#Delta p (%d< #theta < %d);#Delta p;Counts",min,max);
    h_Sim_corr_Smear_binThetaCD[i] = (TH1D*)inFileSim->Get(temp_name);
    hist_list.push_back(h_Sim_corr_Smear_binThetaCD[i]);
  }

  TH1D * h_Sim_corr_Smear_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"corr_Smear_phi_%d",i);
    sprintf(temp_title,"#Delta p (%d< #phi < %d);#Delta p;Counts",min,max);
    h_Sim_corr_Smear_binPhiCD[i] = (TH1D*)inFileSim->Get(temp_name);
    hist_list.push_back(h_Sim_corr_Smear_binPhiCD[i]);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }


  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
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
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  ///////////////////////////////////
  TGraph * g_Data_ep_sigma[6];
  TGraph * g_Data_e_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Data_e_sigma[j] = new TGraph;
    g_Data_ep_sigma[j] = new TGraph;
    //myCanvas->Divide(3,4);
    for(int i = 0; i < 11; i++){
      //myCanvas->cd(i+1);
      //h_Data_e_corr_binSector_binTheta[j][i]->Draw();      
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Data_e_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.006,0.01);
      f_thetabin->SetParameter(2,0.016);
      f_thetabin->SetParLimits(2,0.003,0.025);
      TFitResultPtr point = h_Data_e_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.012,0.01);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Data_e_corr_binSector_binTheta[j][i]->GetEntries()>350)){
	g_Data_e_sigma[j]->SetPoint(g_Data_e_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Data_ep_sigma[j]->SetPoint(g_Data_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }

  
  TGraph * g_Data_p_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Data_p_sigma[j] = new TGraph;
    //myCanvas->Divide(3,3);
    for(int i = 7; i < 14; i++){
      //myCanvas->cd(i-6);
      //h_Data_p_corr_binSector_binTheta[j][i]->Draw();
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Data_p_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.02,0.02);
      f_thetabin->SetParameter(2,0.024);
      f_thetabin->SetParLimits(2,0.003,0.03);
      TFitResultPtr point = h_Data_p_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.03,0.03);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Data_p_corr_binSector_binTheta[j][i]->GetEntries()>190)){
	g_Data_p_sigma[j]->SetPoint(g_Data_p_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Data_ep_sigma[j]->SetPoint(g_Data_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }


  ///////////////////////////////////
  TGraph * g_Sim_ep_sigma[6];
  TGraph * g_Sim_e_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Sim_e_sigma[j] = new TGraph;
    g_Sim_ep_sigma[j] = new TGraph;
    //myCanvas->Divide(3,4);
    for(int i = 0; i < 11; i++){
      //myCanvas->cd(i+1);
      //h_Sim_e_corr_binSector_binTheta[j][i]->Draw();      
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Sim_e_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.006,0.01);
      f_thetabin->SetParameter(2,0.016);
      f_thetabin->SetParLimits(2,0.003,0.025);
      TFitResultPtr point = h_Sim_e_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.012,0.01);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Sim_e_corr_binSector_binTheta[j][i]->GetEntries()>350)){
	g_Sim_e_sigma[j]->SetPoint(g_Sim_e_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Sim_ep_sigma[j]->SetPoint(g_Sim_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }


  TGraph * g_Sim_p_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Sim_p_sigma[j] = new TGraph;
    //myCanvas->Divide(3,3);
    for(int i = 7; i < 14; i++){
      //myCanvas->cd(i-6);
      //h_Sim_p_corr_binSector_binTheta[j][i]->Draw();
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Sim_p_corr_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.02,0.02);
      f_thetabin->SetParameter(2,0.024);
      f_thetabin->SetParLimits(2,0.003,0.03);
      TFitResultPtr point = h_Sim_p_corr_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.03,0.03);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Sim_p_corr_binSector_binTheta[j][i]->GetEntries()>190)){
	g_Sim_p_sigma[j]->SetPoint(g_Sim_p_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Sim_ep_sigma[j]->SetPoint(g_Sim_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }

  ///////////////////////////////////
  TGraph * g_Smear_Sim_ep_sigma[6];
  TGraph * g_Smear_Sim_e_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Smear_Sim_e_sigma[j] = new TGraph;
    g_Smear_Sim_ep_sigma[j] = new TGraph;
    //myCanvas->Divide(3,4);
    for(int i = 0; i < 11; i++){
      //myCanvas->cd(i+1);
      //h_Sim_e_corr_Smear_binSector_binTheta[j][i]->Draw();      
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Sim_e_corr_Smear_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.006,0.01);
      f_thetabin->SetParameter(2,0.016);
      f_thetabin->SetParLimits(2,0.003,0.025);
      TFitResultPtr point = h_Sim_e_corr_Smear_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.012,0.01);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Sim_e_corr_Smear_binSector_binTheta[j][i]->GetEntries()>350)){
	g_Smear_Sim_e_sigma[j]->SetPoint(g_Smear_Sim_e_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Smear_Sim_ep_sigma[j]->SetPoint(g_Smear_Sim_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }

  
  TGraph * g_Smear_Sim_p_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_Smear_Sim_p_sigma[j] = new TGraph;
    //myCanvas->Divide(3,3);
    for(int i = 7; i < 14; i++){
      //myCanvas->cd(i-6);
      //h_Sim_p_corr_Smear_binSector_binTheta[j][i]->Draw();
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.05,0.05,3);
      f_thetabin->SetParameter(0,h_Sim_p_corr_Smear_binSector_binTheta[j][i]->GetMaximum());
      f_thetabin->SetParameter(1,0);
      f_thetabin->SetParLimits(1,-0.02,0.02);
      f_thetabin->SetParameter(2,0.024);
      f_thetabin->SetParLimits(2,0.003,0.03);
      TFitResultPtr point = h_Sim_p_corr_Smear_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-0.03,0.03);
      //f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if((point!=-1)&&(h_Sim_p_corr_Smear_binSector_binTheta[j][i]->GetEntries()>190)){
	g_Smear_Sim_p_sigma[j]->SetPoint(g_Smear_Sim_p_sigma[j]->GetN(),x,100*point->Parameter(2));
	g_Smear_Sim_ep_sigma[j]->SetPoint(g_Smear_Sim_ep_sigma[j]->GetN(),x,100*point->Parameter(2));
      }
    }
    //myCanvas->Print(fileName,"pdf");
    //myCanvas->Clear();
  }
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////

  TF1 * f_Data_FDres[6];
  TFitResultPtr p_Data_FDres[6];
  for(int i = 0; i < 6; i++){
    f_Data_FDres[i] = new TF1("f_Data_res",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },7,45,3);
    f_Data_FDres[i]->SetParameter(0,5.0);
    f_Data_FDres[i]->SetParameter(1,1.0);
    f_Data_FDres[i]->SetParameter(2,1.0);
    p_Data_FDres[i] = g_Data_ep_sigma[i]->Fit(f_Data_FDres[i],"SrBeqn","",7,41);
  }

  TF1 * f_Sim_FDres[6];
  TFitResultPtr p_Sim_FDres[6];
  for(int i = 0; i < 6; i++){
    f_Sim_FDres[i] = new TF1("f_Sim_res",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },7,45,3);
    f_Sim_FDres[i]->SetParameter(0,5.0);
    f_Sim_FDres[i]->SetParameter(1,1.0);
    f_Sim_FDres[i]->SetParameter(2,1.0);
    p_Sim_FDres[i] = g_Sim_ep_sigma[i]->Fit(f_Sim_FDres[i],"SrBeqn","",7,41);
  }

  TF1 * f_Smear_Sim_FDres[6];
  TFitResultPtr p_Smear_Sim_FDres[6];
  for(int i = 0; i < 6; i++){
    f_Smear_Sim_FDres[i] = new TF1("f_Smear_Sim_res",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },7,45,3);
    f_Smear_Sim_FDres[i]->SetParameter(0,5.0);
    f_Smear_Sim_FDres[i]->SetParameter(1,1.0);
    f_Smear_Sim_FDres[i]->SetParameter(2,1.0);
    p_Smear_Sim_FDres[i] = g_Smear_Sim_ep_sigma[i]->Fit(f_Smear_Sim_FDres[i],"SrBeqn","",7,41);
  }


  /////////////
  /////////////

  
  cout<<"{";
  TF1 * f_Diff_FDres[6];
  TF1 * f_Sup_FDres[6];
  for(int i = 0; i < 6; i++){
    f_Diff_FDres[i] = new TF1("f_Diff_res",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },7,45,3);
    double d0 = p_Data_FDres[i]->Parameter(0);
    double d1 = p_Data_FDres[i]->Parameter(1);
    double d2 = p_Data_FDres[i]->Parameter(2);
    double s0 = p_Sim_FDres[i]->Parameter(0);
    double s1 = p_Sim_FDres[i]->Parameter(1);
    double s2 = p_Sim_FDres[i]->Parameter(2);
    f_Diff_FDres[i]->SetParameter(0,d0-s0);
    f_Diff_FDres[i]->SetParameter(1,d1-s1);
    f_Diff_FDres[i]->SetParameter(2,d2-s2);

    f_Sup_FDres[i] = new TF1("f_Sup_res",[&](double *x, double *p){ return DoubQuad(x[0],p[0],p[1],p[2],p[3],p[4],p[5]); },7,45,36);
    f_Sup_FDres[i]->SetParameter(0,d0);
    f_Sup_FDres[i]->SetParameter(1,d1);
    f_Sup_FDres[i]->SetParameter(2,d2);
    f_Sup_FDres[i]->SetParameter(3,s0);
    f_Sup_FDres[i]->SetParameter(4,s1);
    f_Sup_FDres[i]->SetParameter(5,s2);
    cout<<"{"<<d0<<","<<d1<<","<<d2<<","<<s0<<","<<s1<<","<<s2<<"},\n";
  }
  
  double x_ab1[2] = {7,45};
  double y_ab1[2] = {0,3.0};  
  TGraph * r_ab1[6];
  for(int i = 0; i < 6; i++){
    r_ab1[i]= new TGraph(2,x_ab1,y_ab1);
    r_ab1[i]->SetLineColor(0);
    sprintf(temp_title,"#sigma_{#Delta p} vs. #theta (e,e'p) Sector %d;#theta^{#circ};#sigma_{#Delta p} [MeV]",i+1);
    r_ab1[i]->SetTitle(temp_title);
  }
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1[i]->Draw();
    g_Data_e_sigma[i]->SetLineWidth(2);
    g_Data_p_sigma[i]->SetLineWidth(2);
    f_Data_FDres[i]->SetLineWidth(2);
    g_Data_e_sigma[i]->SetLineColor(2);
    g_Data_p_sigma[i]->SetLineColor(2);
    f_Data_FDres[i]->SetLineColor(1);
    g_Data_e_sigma[i]->Draw("SAME");
    g_Data_p_sigma[i]->Draw("SAME");
    f_Data_FDres[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1[i]->Draw();
    g_Sim_e_sigma[i]->SetLineColor(4);
    g_Sim_p_sigma[i]->SetLineColor(4);
    f_Sim_FDres[i]->SetLineColor(1);
    g_Sim_e_sigma[i]->Draw("SAME");
    g_Sim_p_sigma[i]->Draw("SAME");
    f_Sim_FDres[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
    
  TGraph * r_abc1[6];
  for(int i = 0; i < 6; i++){
    r_abc1[i]= new TGraph(2,x_ab1,y_ab1);
    r_abc1[i]->SetLineColor(0);
    sprintf(temp_title,"#sqrt{#sigma_{#Delta p,Data}^{2}-#sigma_{#Delta p,Sim}^{2}} vs. #theta (e,e'p) Sector %d;#theta^{#circ};#sqrt{#sigma_{#Delta p,Data}^{2}-#sigma_{#Delta p,Sim}^{2}} [MeV]",i+1);
    r_abc1[i]->SetTitle(temp_title);
  }
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_abc1[i]->Draw();
    f_Sup_FDres[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1[i]->Draw();
    g_Data_e_sigma[i]->Draw("SAME");
    g_Data_p_sigma[i]->Draw("SAME");
    g_Sim_e_sigma[i]->Draw("SAME");
    g_Sim_p_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1[i]->Draw();
    g_Data_e_sigma[i]->Draw("SAME");
    g_Data_p_sigma[i]->Draw("SAME");
    g_Smear_Sim_e_sigma[i]->SetLineColor(9);
    g_Smear_Sim_p_sigma[i]->SetLineColor(9);
    g_Smear_Sim_e_sigma[i]->Draw("SAME");
    g_Smear_Sim_p_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  //////////////////////////////////////
  
  TGraph * g_Data_sigma_CD = new TGraph;
  //myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    //myCanvas->cd(i+1);
    //h_Data_corr_binThetaCD[i]->Draw();    
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    f_thetabin->SetParameter(0,h_Data_corr_binThetaCD[i]->GetMaximum());
    f_thetabin->SetParameter(1,0);
    f_thetabin->SetParameter(2,0.1);
    TFitResultPtr point = h_Data_corr_binThetaCD[i]->Fit(f_thetabin,"SrBeqn","",-0.15,0.15);
    //f_thetabin->Draw("SAME");
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    if((point!=-1)){
      g_Data_sigma_CD->SetPoint(g_Data_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  TGraph * g_Sim_sigma_CD = new TGraph;
  //myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    //myCanvas->cd(i+1);
    //h_Sim_corr_binThetaCD[i]->Draw();
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.15,0.15,3);
    f_thetabin->SetParameter(0,h_Sim_corr_binThetaCD[i]->GetMaximum());
    f_thetabin->SetParameter(1,0);
    f_thetabin->SetParameter(2,0.1);
    f_thetabin->SetParLimits(2,0,0.06);
    TFitResultPtr point = h_Sim_corr_binThetaCD[i]->Fit(f_thetabin,"SrBeqn","",-0.15,0.15);
    //f_thetabin->Draw("SAME");
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    if((point!=-1)){
      g_Sim_sigma_CD->SetPoint(g_Sim_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  double x_ab2[2] = {35,62};
  double y_ab2[2] = {0,15};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetLineColor(0);
  r_ab2->SetTitle("#sigma_{#Delta p} vs. #theta (e,e'p_{CD});#theta^{#circ};#sigma_{#Delta p} [MeV]");

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab2->Draw();
  g_Data_sigma_CD->SetLineWidth(2);
  g_Data_sigma_CD->SetLineColor(2);
  g_Data_sigma_CD->Draw("SAME");
  g_Sim_sigma_CD->SetLineWidth(1);
  g_Sim_sigma_CD->SetLineColor(4);
  g_Sim_sigma_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  //////////////////////////////////////
  TGraph * g_Phi_Data_sigma_CD = new TGraph;
  //myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    //myCanvas->cd(i+1);
    //h_Data_corr_binPhiCD[i]->Draw();     
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Data_corr_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,0);
    f_phibin->SetParameter(2,0.1);
    f_phibin->SetParLimits(2,0.0,0.15);
    TFitResultPtr point = h_Data_corr_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    if(h_Data_corr_binPhiCD[i]->GetEntries()<30){continue;}
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Data_sigma_CD->SetPoint(g_Phi_Data_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  //myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    //myCanvas->cd(i-8);
    //h_Data_corr_binPhiCD[i]->Draw();     
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Data_corr_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,0);
    f_phibin->SetParameter(2,0.1);
    f_phibin->SetParLimits(2,0.0,0.15);
    TFitResultPtr point = h_Data_corr_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    if(h_Data_corr_binPhiCD[i]->GetEntries()<30){continue;}
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Data_sigma_CD->SetPoint(g_Phi_Data_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  
  TGraph * g_Phi_Sim_sigma_CD = new TGraph;
  //myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    //myCanvas->cd(i+1);
    //h_Sim_corr_binPhiCD[i]->Draw();
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Sim_corr_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,0);
    f_phibin->SetParameter(2,0.1);
    f_phibin->SetParLimits(2,0,0.06);
    if(h_Sim_corr_binPhiCD[i]->GetEntries()<30){continue;}
    TFitResultPtr point = h_Sim_corr_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Sim_sigma_CD->SetPoint(g_Phi_Sim_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  //myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    //myCanvas->cd(i-8);
    //h_Sim_corr_binPhiCD[i]->Draw();
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Sim_corr_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,0);
    f_phibin->SetParameter(2,0.1);
    f_phibin->SetParLimits(2,0,0.06);
    TFitResultPtr point = h_Sim_corr_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    if(h_Sim_corr_binPhiCD[i]->GetEntries()<30){continue;}
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Sim_sigma_CD->SetPoint(g_Phi_Sim_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  TGraph * g_Phi_Smear_Sim_sigma_CD = new TGraph;
  //myCanvas->Divide(3,3);
  for(int i = 0; i < 9; i++){
    //myCanvas->cd(i+1);
    //h_Sim_corr_Smear_binPhiCD[i]->Draw();
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Sim_corr_Smear_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,h_Sim_corr_Smear_binPhiCD[i]->GetMean());
    f_phibin->SetParameter(2,h_Sim_corr_Smear_binPhiCD[i]->GetStdDev());
    f_phibin->SetParLimits(2,h_Sim_corr_Smear_binPhiCD[i]->GetStdDev()-0.02,0.35);
    if(h_Sim_corr_Smear_binPhiCD[i]->GetEntries()<30){continue;}
    TFitResultPtr point = h_Sim_corr_Smear_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Smear_Sim_sigma_CD->SetPoint(g_Phi_Smear_Sim_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  //myCanvas->Divide(3,3);
  for(int i = 9; i < 18; i++){
    //myCanvas->cd(i-8);
    //h_Sim_corr_Smear_binPhiCD[i]->Draw();
    TF1 * f_phibin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.35,0.35,3);
    f_phibin->SetParameter(0,h_Sim_corr_Smear_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(0,h_Sim_corr_Smear_binPhiCD[i]->GetMaximum());
    f_phibin->SetParameter(1,h_Sim_corr_Smear_binPhiCD[i]->GetMean());
    f_phibin->SetParameter(2,h_Sim_corr_Smear_binPhiCD[i]->GetStdDev());
    f_phibin->SetParLimits(2,h_Sim_corr_Smear_binPhiCD[i]->GetStdDev()-0.02,0.35);
    TFitResultPtr point = h_Sim_corr_Smear_binPhiCD[i]->Fit(f_phibin,"SrBeqn","",-0.35,0.35);
    //f_phibin->Draw("SAME");
    if(h_Sim_corr_Smear_binPhiCD[i]->GetEntries()<30){continue;}
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if((point!=-1)){
      g_Phi_Smear_Sim_sigma_CD->SetPoint(g_Phi_Smear_Sim_sigma_CD->GetN(),x,100*point->Parameter(2));
    }
  }
  //myCanvas->Print(fileName,"pdf");
  //myCanvas->Clear();

  
  double x_ab3[2] = {-155,155};
  double y_ab3[2] = {0,15};  
  TGraph * r_ab3 = new TGraph(2,x_ab3,y_ab3);
  r_ab3->SetLineColor(0);
  r_ab3->SetTitle("#sigma_{#Delta p} vs. #phi (e,e'p_{CD});#phi^{#circ};#sigma_{#Delta p} [MeV]");


  TGraph * r_abc3 = new TGraph(2,x_ab3,y_ab3);
  r_abc3->SetLineColor(0);
  r_abc3->SetTitle("#sqrt{#sigma_{#Delta p,Data}^{2}-#sigma_{#Delta p,Sim}^{2}} vs. #phi (e,e'p_{CD});#phi^{#circ};#sqrt{#sigma_{#Delta p,Data}^{2}-#sigma_{#Delta p,Sim}^{2}} [MeV]");

  TF1 * f_Data_CD = new TF1("f_Data_CD",[&](double *x, double *p){ return Trig3(x[0],p[0],p[1],p[2],p[3],p[4],p[5],p[6]); },-180,180,7);
  f_Data_CD->SetParLimits(0,4,16);
  f_Data_CD->SetParLimits(1,0,8);
  f_Data_CD->SetParLimits(3,0,8);
  f_Data_CD->SetParLimits(5,0,8);
  f_Data_CD->SetParLimits(2,-M_PI,M_PI);
  f_Data_CD->SetParLimits(4,-M_PI,M_PI);
  f_Data_CD->SetParLimits(6,-M_PI,M_PI);
  TFitResultPtr p_Data_CD = g_Phi_Data_sigma_CD->Fit(f_Data_CD,"SrBeqn","",-180,180);

    TF1 * f_Sim_CD = new TF1("f_Sim_CD",[&](double *x, double *p){ return p[0]; },-180,180,1);
  f_Data_CD->SetParLimits(0,0,8);
  TFitResultPtr p_Sim_CD = g_Phi_Sim_sigma_CD->Fit(f_Sim_CD,"SrBeqn","",-180,180);

  TF1 * f_Sup_CD = new TF1("f_Sup_CD",[&](double *x, double *p){ return DoubTrig3(x[0],p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]); },-180,180,8);
  for(int i = 0; i < 8; i++){
    f_Sup_CD->SetParameter(i,SmearCD[i]);
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab3->Draw();
  g_Phi_Data_sigma_CD->SetLineWidth(2);
  g_Phi_Data_sigma_CD->SetLineColor(2);
  g_Phi_Data_sigma_CD->Draw("SAME");
  g_Phi_Sim_sigma_CD->SetLineWidth(1);
  g_Phi_Sim_sigma_CD->SetLineColor(4);
  g_Phi_Sim_sigma_CD->Draw("SAME");
  f_Data_CD->SetLineColor(1);
  f_Sim_CD->SetLineColor(1);
  f_Data_CD->Draw("SAME");
  f_Sim_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_abc3->Draw();
  f_Sup_CD->SetLineColor(2);
  f_Sup_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab3->Draw();
  g_Phi_Data_sigma_CD->SetLineWidth(2);
  g_Phi_Data_sigma_CD->SetLineColor(2);
  g_Phi_Data_sigma_CD->Draw("SAME");
  g_Phi_Sim_sigma_CD->SetLineWidth(1);
  g_Phi_Sim_sigma_CD->SetLineColor(4);
  g_Phi_Sim_sigma_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab3->Draw();
  g_Phi_Data_sigma_CD->SetLineWidth(2);
  g_Phi_Data_sigma_CD->SetLineColor(2);
  g_Phi_Data_sigma_CD->Draw("SAME");
  g_Phi_Smear_Sim_sigma_CD->SetLineWidth(1);
  g_Phi_Smear_Sim_sigma_CD->SetLineColor(4);
  g_Phi_Smear_Sim_sigma_CD->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  cout<<"{"<<p_Data_CD->Parameter(0)<<","<<p_Data_CD->Parameter(1)<<","<<p_Data_CD->Parameter(2)<<","<<p_Data_CD->Parameter(3)<<","<<p_Data_CD->Parameter(4)<<","<<p_Data_CD->Parameter(5)<<","<<p_Data_CD->Parameter(6)<<"}\n";
  cout<<p_Sim_CD->Parameter(0)<<"\n";
  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}
