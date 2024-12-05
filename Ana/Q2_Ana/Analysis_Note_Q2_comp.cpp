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

const double c = 29.9792458;


void Usage()
{
  std::cerr << "Usage: ./code outputfile.pdf inputfiledata.root inputfilesim.root \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }



  char * pdfFile = argv[1];
  cout<<"Ouput PDF file "<< pdfFile <<endl;

  TFile * inDataFile = new TFile(argv[2]);
  TFile * inSimFile = new TFile(argv[3]);

  vector<many_plots> hist_Data_list_ep;
  many_plots h_Data_xB("xB","x_{B}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_xB);
  many_plots h_Data_Q2("Q2","Q^{2}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_Q2);
  many_plots h_Data_omega("omega","#omega",inDataFile);
  hist_Data_list_ep.push_back(h_Data_omega);
  many_plots h_Data_thetae("thetae","#theta_{e}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_thetae);
  many_plots h_Data_phie("phie","#phi_{e}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_phie);
  many_plots h_Data_plead("plead","p_{Lead}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_plead);
  many_plots h_Data_thetalead("thetalead","#theta_{Lead}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_thetalead);
  many_plots h_Data_thetalead_FD("thetalead_FD","FD #theta_{Lead}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_thetalead_FD);
  many_plots h_Data_philead("philead","#phi_{Lead}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_philead);
  many_plots h_Data_pmiss("pmiss","p_{miss}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_pmiss);
  many_plots h_Data_mmiss("mmiss","m_{miss}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_mmiss);
  many_plots h_Data_emiss("emiss","E_{miss}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_emiss);
  many_plots h_Data_thetapq("thetapq","#theta_{Lead,q}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_thetapq);
  many_plots h_Data_thetamissq("thetamissq","#theta_{miss,q}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_thetamissq);
  many_plots h_Data_poq("poq","p/q",inDataFile);
  hist_Data_list_ep.push_back(h_Data_poq);
  many_plots h_Data_kmissZQ("kmissZQ","k_{missZQ}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_kmissZQ);
  many_plots h_Data_kmissZB("kmissZB","k_{missZB}",inDataFile);
  hist_Data_list_ep.push_back(h_Data_kmissZB);
 
  vector<many_plots> hist_Data_list_epp;
  many_plots h_Data_precoil("precoil","p_{recoil}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_precoil);
  many_plots h_Data_prel("prel","p_{rel}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_prel);
  many_plots h_Data_thetamissrecoil("thetamissrecoil","#theta_{miss,recoil}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_thetamissrecoil);
  many_plots h_Data_thetacmrel("thetacmrel","#theta_{cm,rel}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_thetacmrel);
  many_plots h_Data_pcm("pcm","p_{cm}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_pcm);
  many_plots h_Data_pcmx("pcmx","p_{X,cm}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_pcmx);
  many_plots h_Data_pcmy("pcmy","p_{||,cm}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_pcmy);
  many_plots h_Data_pcmz("pcmz","p_{miss,cm}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_pcmz);
  many_plots h_Data_E2miss("E2miss","E_{2,miss}",inDataFile);
  hist_Data_list_epp.push_back(h_Data_E2miss);




  vector<many_plots> hist_Sim_list_ep;
  many_plots h_Sim_xB("xB","x_{B}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_xB);
  many_plots h_Sim_Q2("Q2","Q^{2}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_Q2);
  many_plots h_Sim_omega("omega","#omega",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_omega);
  many_plots h_Sim_thetae("thetae","#theta_{e}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_thetae);
  many_plots h_Sim_phie("phie","#phi_{e}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_phie);
  many_plots h_Sim_plead("plead","p_{Lead}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_plead);
  many_plots h_Sim_thetalead("thetalead","#theta_{Lead}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_thetalead);
  many_plots h_Sim_thetalead_FD("thetalead_FD","FD #theta_{Lead}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_thetalead_FD);
  many_plots h_Sim_philead("philead","#phi_{Lead}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_philead);
  many_plots h_Sim_pmiss("pmiss","p_{miss}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_pmiss);
  many_plots h_Sim_mmiss("mmiss","m_{miss}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_mmiss);
  many_plots h_Sim_emiss("emiss","E_{miss}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_emiss);
  many_plots h_Sim_thetapq("thetapq","#theta_{Lead,q}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_thetapq);
  many_plots h_Sim_thetamissq("thetamissq","#theta_{miss,q}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_thetamissq);
  many_plots h_Sim_poq("poq","p/q",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_poq);
  many_plots h_Sim_kmissZQ("kmissZQ","k_{missZQ}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_kmissZQ);
  many_plots h_Sim_kmissZB("kmissZB","k_{missZB}",inSimFile);
  hist_Sim_list_ep.push_back(h_Sim_kmissZB);

  vector<many_plots> hist_Sim_list_epp;
  many_plots h_Sim_precoil("precoil","p_{recoil}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_precoil);
  many_plots h_Sim_prel("prel","p_{rel}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_prel);
  many_plots h_Sim_thetamissrecoil("thetamissrecoil","#theta_{miss,recoil}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_thetamissrecoil);
  many_plots h_Sim_thetacmrel("thetacmrel","#theta_{cm,rel}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_thetacmrel);
  many_plots h_Sim_pcm("pcm","p_{cm}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_pcm);
  many_plots h_Sim_pcmx("pcmx","p_{X,cm}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_pcmx);
  many_plots h_Sim_pcmy("pcmy","p_{||,cm}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_pcmy);
  many_plots h_Sim_pcmz("pcmz","p_{miss,cm}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_pcmz);
  many_plots h_Sim_E2miss("E2miss","E_{2,miss}",inSimFile);
  hist_Sim_list_epp.push_back(h_Sim_E2miss);
  
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
  TH1D * h_xB_ep_data = (TH1D*)inDataFile->Get("xB");
  h_xB_ep_data->Draw();
  TH1D * h_xB_ep_sim = (TH1D*)inSimFile->Get("xB");
  h_xB_ep_sim->SetLineColor(2);
  h_xB_ep_sim->Draw("SAME");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  for(int i=0; i<hist_Data_list_ep.size(); i++){
    myCanvas->Divide(3,4);
    hist_Data_list_ep[i].Draw_hist_set(fileName,myCanvas);
    hist_Sim_list_ep[i].Draw_hist_set_same(fileName,myCanvas,hist_Data_list_ep[i]);  
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }  

  for(int i=0; i<hist_Data_list_ep.size(); i++){
    myCanvas->Divide(3,4);
    hist_Data_list_ep[i].Draw_hist_set_epp(fileName,myCanvas);
    hist_Sim_list_ep[i].Draw_hist_set_epp_same(fileName,myCanvas,hist_Data_list_ep[i]);  
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }  

  for(int i=0; i<hist_Data_list_epp.size(); i++){
    myCanvas->Divide(3,4);
    hist_Data_list_epp[i].Draw_hist_set_epp(fileName,myCanvas);
    hist_Sim_list_epp[i].Draw_hist_set_epp_same(fileName,myCanvas,hist_Data_list_epp[i]);  
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }  

  
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


   
  return 0;
}


