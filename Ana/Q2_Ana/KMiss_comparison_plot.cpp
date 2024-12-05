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
#include "many_plots.h"
#include "reweighter.h"
#include "TGraphErrors.h"
#include "Corrections.h"


using namespace std;
using namespace clas12;


void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputdata.root inputsim.root \n\n\n";
}


int main(int argc, char ** argv)
{

  if(argc < 4)
    {
      Usage();
      return -1;
    }



  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;
  TFile * inDataFile = new TFile(argv[2]);
  TFile * inSimFile = new TFile(argv[3]);

  TH1D * h_kMissZQ_ep_data = (TH1D*)inDataFile->Get("kMissZQ_ep");
  TH1D * h_kMissZQ_ep_sim = (TH1D*)inSimFile->Get("kMissZQ_ep");


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

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->GetPad(1)->SetLeftMargin(0.15);
  h_kMissZQ_ep_data->Draw();
  h_kMissZQ_ep_sim->Scale(h_kMissZQ_ep_data->GetBinContent(45)/h_kMissZQ_ep_sim->GetBinContent(45));
  h_kMissZQ_ep_sim->SetLineColor(2);
  h_kMissZQ_ep_sim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  return 0;
}
