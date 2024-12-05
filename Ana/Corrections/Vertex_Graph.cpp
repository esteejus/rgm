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
#include "TRandom3.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
//TRandom3 * thisRand = new TRandom3(0);;

double Quad(double x, double A, double B, double C){
  return A + B*x + C*x*x; 
}
double T1(double x, double A, double B){
  return A + B/x; 
}
double T2(double x, double A, double B, double C){
  return A + B/x + C/(180-x); 
}



//vector<double> bE_ThetaCD = {35,40,45,50,55,60,70};
vector<double> bE_MomCD = {0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0};//{0.5,1.0,1.3,1.6,2.0,2.5,3.0};
vector<double> bE_Theta = {8,10,12,14,16,18,20,23,26,29,32,35,38,41,45};
vector<double> bE_ThetaCD = {35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125};
vector<double> bE_PhiCD = {-180,-160,-140,-120,-100,-80,-60,-40,-20,0,20,40,60,80,100,120,140,160,180};
vector<double> bE_Phi = {-35,-15,-5,0,5,10,15,25,35};
vector<double> bE_ThetaE = {10,13,16,19,22,25,28,31,34,37};
vector<double> bE_ThetapFD = {19,22,25,28,31,34,37,40,43,46};

auto db=TDatabasePDG::Instance();

double G(double x, double N, double mu, double sigma){
  return (N/(sigma*sqrt(2*M_PI))) * exp(-0.5 * sq((x-mu)/sigma)) ; 
}

double SQ(double x){ return x*x;}

void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfile.root \n\n\n";

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

  TFile * inFile = new TFile(argv[2]);

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;

  TH1D * h_e_vtz_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"evtz_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Z Vertex Sector %d (%d< #theta < %d);Z Vertex [cm];Counts",j,min,max);
      h_e_vtz_binSector_binTheta[j-1][i] = (TH1D*)inFile->Get(temp_name);
      hist_list.push_back(h_e_vtz_binSector_binTheta[j-1][i]);
    }
  }

  TH1D * h_p_vtz_binSector_binTheta[6][14];
  for(int j=1; j<=6; j++){
    for(int i=0; i<14; i++){
      int min = bE_Theta[i];
      int max = bE_Theta[i+1];
      sprintf(temp_name,"pvtz_sector_%d_theta_%d",j,i);
      sprintf(temp_title,"Z Vertex Sector %d (%d< #theta < %d);Z Vertex [cm];Counts",j,min,max);
      h_p_vtz_binSector_binTheta[j-1][i] = (TH1D*)inFile->Get(temp_name);
      hist_list.push_back(h_p_vtz_binSector_binTheta[j-1][i]);
    }
  }
  
  TH1D * h_vtz_binThetaCD[18];
  for(int i=0; i<18; i++){
    int min = bE_ThetaCD[i];
    int max = bE_ThetaCD[i+1];
    sprintf(temp_name,"vtz_theta_%d",i);
    sprintf(temp_title,"Z Vertex (%d< #theta < %d);Z Vertex [cm];Counts",min,max);
    h_vtz_binThetaCD[i] = (TH1D*)inFile->Get(temp_name);
    h_vtz_binThetaCD[i]->SetTitle(temp_title);
    hist_list.push_back(h_vtz_binThetaCD[i]);
  }

  TH1D * h_vtz_binPhiCD[18];
  for(int i=0; i<18; i++){
    int min = bE_PhiCD[i];
    int max = bE_PhiCD[i+1];
    sprintf(temp_name,"vtz_phi_%d",i);
    sprintf(temp_title,"Z Vertex (%d< #phi < %d);Z Vertex;Counts",min,max);
    h_vtz_binPhiCD[i] = (TH1D*)inFile->Get(temp_name);
    hist_list.push_back(h_vtz_binPhiCD[i]);
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

  TGraph * g_e_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_e_sigma[j] = new TGraph;
    myCanvas->Divide(3,4);
    for(int i = 1; i < 11; i++){
      myCanvas->cd(i);
      myCanvas->GetPad(i)->SetLeftMargin(0.15);
      if(h_e_vtz_binSector_binTheta[j][i]->GetEntries()<100){continue;}
      h_e_vtz_binSector_binTheta[j][i]->Draw();      
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-5,1,3);
      f_thetabin->SetParameter(0,h_e_vtz_binSector_binTheta[j][i]->GetMaximum());
      double m = h_e_vtz_binSector_binTheta[j][i]->GetMean();
      double s = h_e_vtz_binSector_binTheta[j][i]->GetStdDev();
      f_thetabin->SetParameter(1,m);
      f_thetabin->SetParameter(2,s);
      f_thetabin->SetParLimits(2,0.1,2.0);
      TFitResultPtr point = h_e_vtz_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",m-s,m+s);
      f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if(point!=-1){
	g_e_sigma[j]->SetPoint(g_e_sigma[j]->GetN(),x,point->Parameter(2));
      }
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }

  
  TGraph * g_p_sigma[6];
  for(int j = 0; j < 6 ; j++){
    g_p_sigma[j] = new TGraph;
    myCanvas->Divide(4,3);
    for(int i = 2; i < 14; i++){
      myCanvas->cd(i-1);
      if(h_p_vtz_binSector_binTheta[j][i]->GetEntries()<300){h_p_vtz_binSector_binTheta[j][i]->Rebin(2);}
      h_p_vtz_binSector_binTheta[j][i]->Draw();
      TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-7,1,3);
      f_thetabin->SetParameter(0,h_p_vtz_binSector_binTheta[j][i]->GetMaximum());
      double m = h_e_vtz_binSector_binTheta[j][i]->GetMean();
      double s = h_e_vtz_binSector_binTheta[j][i]->GetStdDev();
      f_thetabin->SetParameter(1,m);
      f_thetabin->SetParameter(2,s);
      f_thetabin->SetParLimits(2,0.1,2.0);
      TFitResultPtr point = h_p_vtz_binSector_binTheta[j][i]->Fit(f_thetabin,"SrBeqn","",-6.5,1);
      f_thetabin->Draw("SAME");
      double x = (bE_Theta[i]+bE_Theta[i+1])/2;
      if(point!=-1){
	g_p_sigma[j]->SetPoint(g_p_sigma[j]->GetN(),x,point->Parameter(2));
      }
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
  }
  
  TGraph * g_sigma = new TGraph;
  myCanvas->Divide(4,5);
  for(int i = 0; i < 18; i++){
    myCanvas->cd(i+1);
    h_vtz_binThetaCD[i]->Draw();
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-5,1,3);
    f_thetabin->SetParameter(0,h_vtz_binThetaCD[i]->GetMaximum());
    //f_thetabin->SetParLimits(0,h_vtz_binThetaCD[i]->GetMaximum()*0.5,h_vtz_binThetaCD[i]->GetMaximum()*2);
    double m = h_vtz_binThetaCD[i]->GetMean();
    double s = h_vtz_binThetaCD[i]->GetStdDev();
    f_thetabin->SetParameter(1,m);
    f_thetabin->SetParameter(2,s);
    TFitResultPtr point = h_vtz_binThetaCD[i]->Fit(f_thetabin,"SrBeqn","",m-s,m+s);
    f_thetabin->Draw("SAME");
    double x = (bE_ThetaCD[i]+bE_ThetaCD[i+1])/2;
    if(point!=-1){
      g_sigma->SetPoint(g_sigma->GetN(),x,point->Parameter(2));
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  /*
  TGraph * g_sigma_Phi = new TGraph;
  myCanvas->Divide(4,5);
  for(int i = 0; i < 18; i++){
    myCanvas->cd(i+1);
    h_vtz_binPhiCD[i]->Draw();
    TF1 * f_thetabin = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-5,1,3);
    f_thetabin->SetParameter(0,h_vtz_binPhiCD[i]->GetMaximum());
    double m = h_vtz_binPhiCD[i]->GetMean();
    double s = h_vtz_binPhiCD[i]->GetStdDev();
    f_thetabin->SetParameter(1,m);
    f_thetabin->SetParameter(2,s);
    TFitResultPtr point = h_vtz_binPhiCD[i]->Fit(f_thetabin,"SrBeqn","",m-s,m+s);
    f_thetabin->Draw("SAME");
    double x = (bE_PhiCD[i]+bE_PhiCD[i+1])/2;
    if(point!=-1){
      g_sigma_Phi->SetPoint(g_sigma_Phi->GetN(),x,point->Parameter(2));
    }
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  */
  
  double x_ab1[2] = {12,30};
  double y_ab1[2] = {0,1.0};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetLineColor(0);
  r_ab1->SetTitle("Vertex Z vs. #theta_{e} (e,e'p);#theta_{e}^{#circ};Vertex Z [cm]");

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab1->Draw();
  for(int i = 0; i < 6; i++){
    g_e_sigma[i]->SetLineColor(i+1);
    g_e_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  cout<<"Electrons=\n";
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab1->Draw();
    g_e_sigma[i]->SetLineColor(2);
    g_e_sigma[i]->Draw("SAME");

    TF1 * f_t = new TF1("f_t",[&](double *x, double *p){ return T1(x[0],p[0],p[1]); },0,90,2);
    TFitResultPtr point = g_e_sigma[i]->Fit(f_t,"SrBeqn","",x_ab1[0],x_ab1[1]);
    f_t->SetLineColor(1);
    f_t->Draw("SAME");
    cout<<point->Parameter(0)<<"&"<<point->Parameter(1)<<"\n";
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  cout<<"\n\n\n\n";
  
  double x_ab2[2] = {12,40};
  double y_ab2[2] = {0,2.5};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetLineColor(0);
  r_ab2->SetTitle("Vertex Z vs. #theta_{e} (e,e'p);#theta_{e}^{#circ};Vertex Z [cm]");

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  r_ab2->Draw();
  for(int i = 0; i < 6; i++){
    g_p_sigma[i]->SetLineColor(i+1);
    g_p_sigma[i]->Draw("SAME");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  cout<<"Protons=\n";
  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    r_ab2->Draw();
    g_p_sigma[i]->SetLineColor(2);
    g_p_sigma[i]->Draw("SAME");

    TF1 * f_t = new TF1("f_t",[&](double *x, double *p){ return T1(x[0],p[0],p[1]); },0,90,2);
    TFitResultPtr point = g_p_sigma[i]->Fit(f_t,"SrBeqn","",15,35);
    f_t->SetLineColor(1);
    f_t->Draw("SAME");
    cout<<point->Parameter(0)<<"&"<<point->Parameter(1)<<"\n";
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  cout<<"\n\n\n\n";

  double x_ab3[2] = {40,120};
  double y_ab3[2] = {0,0.35};  
  TGraph * r_ab3 = new TGraph(2,x_ab3,y_ab3);
  r_ab3->SetLineColor(0);
  r_ab3->SetTitle("Vertex Z vs. #theta_{e} (e,e'p);#theta_{e}^{#circ};Vertex Z [cm]");

  myCanvas->Divide(1,1);
  r_ab3->Draw();
  g_sigma->Draw("SAME");

  TF1 * f_t = new TF1("f_t",[&](double *x, double *p){ return Quad(x[0],p[0],p[1],p[2]); },0,180,3);
  TFitResultPtr point = g_sigma->Fit(f_t,"SrBeqn","",0,180);
  f_t->Draw("SAME");
  cout<<"Protons CD \n"<<point->Parameter(0)<<"&"<<point->Parameter(1)<<"&"<<point->Parameter(2)<<"\n";

  
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}
