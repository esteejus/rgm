#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "eventcut_old.h"
#include "functions_old.h"

using namespace std;
using namespace clas12;

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/cutfile.txt> <path/to/input.hipo> \n";
}


double binEdges[] = { 0.35, 0.38, 0.41, 0.44, 0.47, 0.5, 0.53, 0.56, 0.59, 0.61, 0.64, 0.67, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };
int binEdgeslength = sizeof(binEdges)/sizeof(binEdges[0]) -1;

int main(int argc, char ** argv)
{

  if(argc < 6)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  
  bool isMC = false;
  if(atoi(argv[1]) == 1){isMC=true;}

  double Ebeam = atof(argv[2]);
  
  TFile * outFile = new TFile(argv[3],"RECREATE");
  eventcut myCut(Ebeam,argv[4]);
  myCut.print_cuts();
  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  
  /////////////////////////////////////
  //Prepare histograms
  /////////////////////////////////////
  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(0.8);

  char temp_name[100];
  char temp_title[100];

  TH1D * h_pmiss_Rec = new TH1D("pmiss_Rec","p_{miss} Rec;p_{miss};Counts",binEdgeslength,binEdges);
  hist_list_1.push_back(h_pmiss_Rec);
  TH1D * h_Q2_Rec = new TH1D("Q2_Rec","Q^{2};Q^{2};Counts",100,1.0,6.0);
  hist_list_1.push_back(h_Q2_Rec);

  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Sumw2();
    hist_list_1[i]->GetXaxis()->CenterTitle();
    hist_list_1[i]->GetYaxis()->CenterTitle();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Sumw2();
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->CenterTitle();
  }


  int counter = 0;

  //Define cut class
  //while((chain.Next()==true) && (counter<100000)){
  while(chain.Next()==true){
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      // get particles by type
      auto allParticles = c12->getDetParticles();
      auto electrons=c12->getByID(11);
      double weight = 1;
      if(isMC){weight=c12->mcevent()->getWeight();}
      TVector3 	p_b(0,0,Ebeam);
      if(!myCut.electroncut(c12)){continue;}      
      TVector3 p_e;
      p_e.SetMagThetaPhi(electrons[0]->getP(),electrons[0]->getTheta(),electrons[0]->getPhi());
      double vtz_e = electrons[0]->par()->getVz();
      TVector3	p_q = p_b - p_e;
      double nu = Ebeam - p_e.Mag();
      double QSq = p_q.Mag2() - (nu*nu);
      double xB = QSq / (2 * mN * nu);
      double WSq = (mN*mN) - QSq + (2*nu*mN);
      double theta_e = p_e.Theta() * 180 / M_PI;
      if(QSq < 1.5){continue;}
      int num_L = 0;
      int index_L = -1;
      bool isFDlead = false;
      bool isCDlead = false;
      for(int j = 0; j < allParticles.size(); j ++){
	if( LeadFDProton_Cut(c12,Ebeam,j) || LeadCDProton_Cut(c12,Ebeam,j)){	  
	  num_L++;
	  index_L=j;
	  isFDlead = LeadFDProton_Cut(c12,Ebeam,j);
	  isCDlead = LeadCDProton_Cut(c12,Ebeam,j);
	}
      }
      if(num_L!=1){continue;}
      TVector3 p_L;
      p_L.SetMagThetaPhi(allParticles[index_L]->getP(),allParticles[index_L]->getTheta(),allParticles[index_L]->getPhi());
      TVector3 p_1 = p_L - p_q;
      TVector3 p_miss = p_1;
      double mmiss = get_mmiss(p_b,p_e,p_L);
      double Loq = p_L.Mag() / p_q.Mag();
      double theta_Lq = p_L.Angle(p_q) * 180 / M_PI;
      double vtz_L = allParticles[index_L]->par()->getVz();
      int numlayers_hit = 0;
      if(xB < 1.2){continue;}
      if(theta_Lq>25){continue;}
      if(Loq<0.6){continue;}
      if(Loq>0.96){continue;}
      if(mmiss>1.1){continue;}
      int num_pos  = 0;
      int index_Rp1 = -1;
      for(int j = 0; j < allParticles.size(); j ++){	
	if(AllProton_Cut(c12,Ebeam,j) && (j!=index_L)){
	  if(num_pos==0){index_Rp1=j;}
	  num_pos++;
	}	
      }
      if(num_pos<1){continue;}
      double mom = allParticles[index_Rp1]->getP();
      TVector3 p_2;
      p_2.SetMagThetaPhi(allParticles[index_Rp1]->getP(),allParticles[index_Rp1]->getTheta(),allParticles[index_Rp1]->getPhi());
      TVector3 p_rel = p_1-p_2;
      p_rel.SetMag(p_rel.Mag()/2);
      TVector3 p_cm = p_1+p_2;
      double vtz_Rp1 = allParticles[index_Rp1]->par()->getVz();
      if(fabs(vtz_L-vtz_Rp1)>1.5){continue;}
      h_pmiss_Rec->Fill(p_miss.Mag(),weight);
      h_Q2_Rec->Fill(QSq,weight);
  }
  cout<<counter<<endl;

  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }
  outFile->Close();
}

