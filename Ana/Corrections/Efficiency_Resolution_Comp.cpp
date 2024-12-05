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

void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfiledata.root inputfilesim.root \n\n\n";

}

int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  char * pdfFile = argv[1];
  TFile * inFileData = new TFile(argv[2]);
  TFile * inFileSim  = new TFile(argv[3]);

  TH1D * h_thetaFD_eData = (TH1D*)inFileData->Get("thetaFD_e");
  TH1D * h_thetaFD_epData = (TH1D*)inFileData->Get("thetaFD_ep");

  TH1D * h_thetaFD_eSim = (TH1D*)inFileSim->Get("thetaFD_e");
  TH1D * h_thetaFD_epSim = (TH1D*)inFileSim->Get("thetaFD_ep");
  
  
  TH1D * h_Q2FD_eData = (TH1D*)inFileData->Get("Q2FD_e");
  TH1D * h_Q2FD_epData = (TH1D*)inFileData->Get("Q2FD_ep");

  TH1D * h_Q2FD_eSim = (TH1D*)inFileSim->Get("Q2FD_e");
  TH1D * h_Q2FD_epSim = (TH1D*)inFileSim->Get("Q2FD_ep");
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
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_eData->Draw();
  h_thetaFD_eSim->SetLineColor(2);
  h_thetaFD_eSim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_epData->Draw();
  h_thetaFD_epSim->SetLineColor(2);
  h_thetaFD_epSim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_epData->Divide(h_thetaFD_eData);
  h_thetaFD_epSim->Divide(h_thetaFD_eSim);
  h_thetaFD_epData->Draw();
  h_thetaFD_epSim->SetLineColor(2);
  h_thetaFD_epSim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  double x_ab1[2] = {20,45};
  double y_ab1[2] = {0,1.1};  
  TGraph * r_ab1 = new TGraph(2,x_ab1,y_ab1);
  r_ab1->SetLineColor(0);
  r_ab1->SetTitle("#eta_{Data}/#eta_{Sim} vs. #theta_{q}^{#circ};#theta_{q}^{#circ};#eta_{Data}/#eta_{Sim}");

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetaFD_epData->Divide(h_thetaFD_epSim);
  r_ab1->Draw();
  h_thetaFD_epData->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Q2FD_epData->Divide(h_Q2FD_eData);
  h_Q2FD_epSim->Divide(h_Q2FD_eSim);
  h_Q2FD_epData->Draw();
  h_Q2FD_epSim->SetLineColor(2);
  h_Q2FD_epSim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  
  double x_ab2[2] = {1,5};
  double y_ab2[2] = {0,1.1};  
  TGraph * r_ab2 = new TGraph(2,x_ab2,y_ab2);
  r_ab2->SetLineColor(0);
  r_ab2->SetTitle("#eta_{Data}/#eta_{Sim} vs. Q2;Q2;#eta_{Data}/#eta_{Sim}");
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Q2FD_epData->Divide(h_Q2FD_epSim);
  r_ab2->Draw();
  h_Q2FD_epData->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}

