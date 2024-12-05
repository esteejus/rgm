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
#include <TGraph.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TLatex.h"
#include "HipoChain.h"
#include "clas12ana.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
auto db=TDatabasePDG::Instance();
double mproton = db->GetParticle(2212)->Mass();
double mpion = db->GetParticle(211)->Mass();
double mD = 1.8756;
double beam_E = 5.98;

double sq(double x){return x*x;}

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}
double E(double x, double N, double tau){
  return N * exp( x / tau) ; 
}

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void Fit_Delta_Eprime(TH1D * h_Delta_Eprime[9], TF1 * f_Delta_Eprime[9], TFitResultPtr p_Delta_Eprime[9], TF1 * f_G1[9], TF1 * f_G2[9], TF1 * f_E[9], double isMC){

  char temp_name[100];
  char temp_title[100];

  if(!isMC){
    int left_bin[9] = {50,50,50,51,52,53,27,28,29};
    int right_bin[9] = {80,80,80,81,82,83,42,43,44};
    double delta_mu_max[9] = {0.35,0.35,0.35,0.4,0.4,0.4,0.4,0.4,0.4};
    double delta_sigma_max[9] = {0.05,0.05,0.05,0.1,0.1,0.1,0.1,0.1,0.1};
    //double left_sigma[9] = {32,36,42,47,60,70,75,80,85};
    //double right_sigma[9] = {80,90,100,140,160,170,170,170,170};
    //double fit_upper_limit[9] = {0.3,0.3,0.3,0.4,0.3,0.4,0.5,0.5,0.5};
    double fit_tau[9] = {0.5,0.3,0.2,0.15,0.15,0.15,0.1,0.1,0.1};
    for(int i = 0; i < 9; i++){
      sprintf(temp_name,"f_Delta_Eprime_%d",i);
      f_Delta_Eprime[i] = new TF1(temp_name,[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]) + G(x[0],p[3],p[4],p[5]) + E(x[0],p[6],p[7]); },-0.5,0.5,8);

      f_Delta_Eprime[i]->SetParameter(0,h_Delta_Eprime[i]->GetBinContent(left_bin[i]));
      f_Delta_Eprime[i]->SetParameter(1,h_Delta_Eprime[i]->GetBinCenter(left_bin[i]));
      f_Delta_Eprime[i]->SetParameter(2,0.05);

      f_Delta_Eprime[i]->SetParameter(3,h_Delta_Eprime[i]->GetBinContent(right_bin[i]));
      f_Delta_Eprime[i]->SetParameter(4,0.3);
      f_Delta_Eprime[i]->SetParameter(5,0.05);

      f_Delta_Eprime[i]->SetParameter(6,h_Delta_Eprime[i]->GetMaximum());
      f_Delta_Eprime[i]->SetParameter(7,fit_tau[i]);

      /////////////////////////////////////////////////
      f_Delta_Eprime[i]->SetParLimits(0,h_Delta_Eprime[i]->GetBinContent(left_bin[i])*0.05,h_Delta_Eprime[i]->GetBinContent(left_bin[i])*1.1);
      f_Delta_Eprime[i]->SetParLimits(1,-0.1,0.1);
      f_Delta_Eprime[i]->SetParLimits(2,0.01,0.15);

      f_Delta_Eprime[i]->SetParLimits(3,h_Delta_Eprime[i]->GetBinContent(right_bin[i])*0.05,h_Delta_Eprime[i]->GetBinContent(right_bin[i])*1.1);
      f_Delta_Eprime[i]->SetParLimits(4,0.28,delta_mu_max[i]);
      f_Delta_Eprime[i]->SetParLimits(5,0.01,delta_sigma_max[i]);


      f_Delta_Eprime[i]->SetParLimits(6,h_Delta_Eprime[i]->GetMaximum()*0.05,h_Delta_Eprime[i]->GetMaximum()*10);
      f_Delta_Eprime[i]->SetParLimits(7,0.05,1);

      p_Delta_Eprime[i] = h_Delta_Eprime[i]->Fit(f_Delta_Eprime[i],"qesrn","",-0.25,0.4);

      f_G1[i] = new TF1(temp_name,[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
      f_G1[i]->SetParameter(0,p_Delta_Eprime[i]->Parameter(0));
      f_G1[i]->SetParameter(1,p_Delta_Eprime[i]->Parameter(1));
      f_G1[i]->SetParameter(2,p_Delta_Eprime[i]->Parameter(2));
      f_G1[i]->SetLineColor(1);

      f_G2[i] = new TF1(temp_name,[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
      f_G2[i]->SetParameter(0,p_Delta_Eprime[i]->Parameter(3));
      f_G2[i]->SetParameter(1,p_Delta_Eprime[i]->Parameter(4));
      f_G2[i]->SetParameter(2,p_Delta_Eprime[i]->Parameter(5));
      f_G2[i]->SetLineColor(3);

      f_E[i] = new TF1(temp_name,[&](double *x, double *p){ return E(x[0],p[0],p[1]); },-0.5,0.5,2);
      f_E[i]->SetParameter(0,p_Delta_Eprime[i]->Parameter(6));
      f_E[i]->SetParameter(1,p_Delta_Eprime[i]->Parameter(7));
      f_E[i]->SetLineColor(4);
    }
  }
  else{
    for(int i = 0; i < 9; i++){
      sprintf(temp_name,"f_Delta_Eprime_%d",i);
      f_Delta_Eprime[i] = new TF1(temp_name,[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.5,0.5,3);
      f_Delta_Eprime[i]->SetParameter(0,h_Delta_Eprime[i]->GetMaximum());
      f_Delta_Eprime[i]->SetParameter(1,h_Delta_Eprime[i]->GetBinCenter(h_Delta_Eprime[i]->GetMaximumBin()));
      f_Delta_Eprime[i]->SetParameter(2,0.1);

      f_Delta_Eprime[i]->SetParLimits(0,h_Delta_Eprime[i]->GetMaximum()*0.1/G(0,1,0,0.1),h_Delta_Eprime[i]->GetMaximum()*1.1/G(0,1,0,0.1));
      f_Delta_Eprime[i]->SetParLimits(1,h_Delta_Eprime[i]->GetBinCenter(h_Delta_Eprime[i]->GetMaximumBin())-0.05,h_Delta_Eprime[i]->GetBinCenter(h_Delta_Eprime[i]->GetMaximumBin())+0.05);
      f_Delta_Eprime[i]->SetParLimits(2,0.01,0.2);


      if(i!=8){
	p_Delta_Eprime[i] = h_Delta_Eprime[i]->Fit(f_Delta_Eprime[i],"qeSrn","",-0.05,0.05);
      }
    }
  }
}

void Usage()
{
  std::cerr << "Usage: ./code isMC outputfile.pdf inputfiles.root \n\n\n";

}

int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }


  int isMC = atoi(argv[1]);
  char * pdfFile = argv[2];
  cout<<"Ouput PDF file "<< pdfFile <<endl;

 //Creat input histo
  TFile * inFile = new TFile(argv[3]);
  char temp_name[100];
  TH1D * h_Delta_Eprime[9];
  TH1D * h_Delta_Eprime_sectors[6][9];
  for(int i=0; i<9; i++){
    int min = 10+(3*i);
    sprintf(temp_name,"Delta_Eprime_%d",min);
    h_Delta_Eprime[i] = (TH1D*)inFile->Get(temp_name);    

    h_Delta_Eprime[i]->GetXaxis()->CenterTitle();
    h_Delta_Eprime[i]->GetYaxis()->SetTitleOffset(1.0);
    h_Delta_Eprime[i]->GetYaxis()->CenterTitle();  
    h_Delta_Eprime[i]->GetXaxis()->SetTitleSize(0.08);
    h_Delta_Eprime[i]->GetYaxis()->SetTitleSize(0.05);
    h_Delta_Eprime[i]->GetYaxis()->SetTitleOffset(1.0);

    for(int j=1; j<=6; j++){
      sprintf(temp_name,"Delta_Eprime_%d_sector_%d",min,j);
      h_Delta_Eprime_sectors[j-1][i] = (TH1D*)inFile->Get(temp_name);    
      h_Delta_Eprime_sectors[j-1][i]->GetXaxis()->CenterTitle();
      h_Delta_Eprime_sectors[j-1][i]->GetYaxis()->SetTitleOffset(1.0);
      h_Delta_Eprime_sectors[j-1][i]->GetYaxis()->CenterTitle();  
      h_Delta_Eprime_sectors[j-1][i]->GetXaxis()->SetTitleSize(0.08);
      h_Delta_Eprime_sectors[j-1][i]->GetYaxis()->SetTitleSize(0.05);
      h_Delta_Eprime_sectors[j-1][i]->GetYaxis()->SetTitleOffset(1.0);
    }
  }
  /////////////////////////////////////////////////////
  //Fit
  /////////////////////////////////////////////////////

  TF1 * f_Delta_Eprime[9];
  TF1 * f_G1[9];
  TF1 * f_G2[9];
  TF1 * f_E[9];
  TFitResultPtr p_Delta_Eprime[9];
  Fit_Delta_Eprime(h_Delta_Eprime,f_Delta_Eprime,p_Delta_Eprime,f_G1,f_G2,f_E,isMC);


  TF1 * f_Delta_Eprime_sectors[6][9];
  TF1 * f_G1_sectors[6][9];
  TF1 * f_G2_sectors[6][9];
  TF1 * f_E_sectors[6][9];
  TFitResultPtr p_Delta_Eprime_sectors[6][9];
  for(int j=1; j<=6; j++){
    Fit_Delta_Eprime(h_Delta_Eprime_sectors[j-1],f_Delta_Eprime_sectors[j-1],p_Delta_Eprime_sectors[j-1],f_G1_sectors[j-1],f_G2_sectors[j-1],f_E_sectors[j-1],isMC);
  }



  
  TGraph * g_Delta_Eprime_mu;
  TGraph * g_Delta_Eprime_sigma;
  g_Delta_Eprime_mu = new TGraph();
  g_Delta_Eprime_sigma = new TGraph();
  for(int i=0; i<9; i++){
    double x = 11.5 + (3 * (double)i);
    if(!(isMC && (i==8))){
      g_Delta_Eprime_mu->SetPoint(g_Delta_Eprime_mu->GetN(),x,p_Delta_Eprime[i]->Parameter(1)*1000);
      g_Delta_Eprime_sigma->SetPoint(g_Delta_Eprime_sigma->GetN(),x,p_Delta_Eprime[i]->Parameter(2)*1000);
    }
  }


  TGraph * g_Delta_Eprime_mu_sectors[6];
  TGraph * g_Delta_Eprime_sigma_sectors[6];
  for(int j=0; j<6; j++){
    g_Delta_Eprime_mu_sectors[j] = new TGraph();
    g_Delta_Eprime_sigma_sectors[j] = new TGraph();
    for(int i=0; i<9; i++){
      double x = 11.5 + (3 * (double)i);
      if(!(isMC && (i==8))){
	g_Delta_Eprime_mu_sectors[j]->SetPoint(g_Delta_Eprime_mu_sectors[j]->GetN(),x,p_Delta_Eprime_sectors[j][i]->Parameter(1)*1000);
	g_Delta_Eprime_sigma_sectors[j]->SetPoint(g_Delta_Eprime_sigma_sectors[j]->GetN(),x,p_Delta_Eprime_sectors[j][i]->Parameter(2)*1000);
      }
    }
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

  myCanvas->Divide(3,3);
  for(int i = 1; i < 10; i++){
    myCanvas->cd(i);
    h_Delta_Eprime[i-1]->Draw();        
    f_Delta_Eprime[i-1]->Draw("SAME");
    if(!isMC){
      f_G1[i-1]->Draw("SAME");
      f_G2[i-1]->Draw("SAME");
      f_E[i-1]->Draw("SAME");
    }
    if(!(isMC && (i==9))){
      text.DrawLatex(-0.4,h_Delta_Eprime[i-1]->GetMaximum()*0.9,Form("#mu = %g #pm %g MeV",p_Delta_Eprime[i-1]->Parameter(1)*1000,p_Delta_Eprime[i-1]->ParError(1)*1000));
      text.DrawLatex(-0.4,h_Delta_Eprime[i-1]->GetMaximum()*0.8,Form("#sigma = %g #pm %g MeV",p_Delta_Eprime[i-1]->Parameter(2)*1000,p_Delta_Eprime[i-1]->ParError(2)*1000));
    }    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  g_Delta_Eprime_mu->Draw();    
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  g_Delta_Eprime_sigma->Draw();    
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

  for(int j=1; j<=6; j++){
    myCanvas->Divide(3,3);
    for(int i = 1; i < 10; i++){
      myCanvas->cd(i);
      h_Delta_Eprime_sectors[j-1][i-1]->Draw();        
      f_Delta_Eprime_sectors[j-1][i-1]->Draw("SAME");
      if(!isMC){
	f_G1_sectors[j-1][i-1]->Draw("SAME");
	f_G2_sectors[j-1][i-1]->Draw("SAME");
	f_E_sectors[j-1][i-1]->Draw("SAME");
      }
      if(!(isMC && (i==9))){
	text.DrawLatex(-0.4,h_Delta_Eprime_sectors[j-1][i-1]->GetMaximum()*0.9,Form("#mu = %g #pm %g MeV",p_Delta_Eprime_sectors[j-1][i-1]->Parameter(1)*1000,p_Delta_Eprime_sectors[j-1][i-1]->ParError(1)*1000));
	text.DrawLatex(-0.4,h_Delta_Eprime_sectors[j-1][i-1]->GetMaximum()*0.8,Form("#sigma = %g #pm %g MeV",p_Delta_Eprime_sectors[j-1][i-1]->Parameter(2)*1000,p_Delta_Eprime_sectors[j-1][i-1]->ParError(2)*1000));
      }    
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  g_Delta_Eprime_mu_sectors[0]->SetLineColor(1);    
  g_Delta_Eprime_mu_sectors[0]->Draw();    
  for(int j=1; j<6; j++){
    g_Delta_Eprime_mu_sectors[j]->SetLineColor(j+1);    
    g_Delta_Eprime_mu_sectors[j]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  g_Delta_Eprime_sigma_sectors[0]->SetLineColor(1);    
  g_Delta_Eprime_sigma_sectors[0]->Draw();    
  for(int j=1; j<6; j++){
    g_Delta_Eprime_sigma_sectors[j]->SetLineColor(j+1);    
    g_Delta_Eprime_sigma_sectors[j]->Draw("SAME");    
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}
