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

using namespace std;
using namespace clas12;

const double c = 29.9792458;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}

void Usage()
{
  std::cerr << "Usage: ./code outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }



  TString outFile = argv[1];
  char * pdfFile = argv[2];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;


  clas12ana clasAna;

  //Read in target parameter files                                                                                                                                                           
  clasAna.readInputParam("../ana.par");
  clasAna.readEcalSFPar("../paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("../paramsPI_40Ca_x2.dat");
  clasAna.printParams();
    

  clas12root::HipoChain chain;
  for(int k = 3; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  //now get reference to (unique)ptr for accessing data in loop
  //this will point to the correct place when file changes
  //  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;

  double beam_E = 5.98;
  const double me = 0.000511;
  const double mU = 0.9314941024;
  const double m_4He = 4.00260325415 * mU - 2*me;
  
  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector nucleus_ptr(0,0,0,m_4He);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  vector<TH1*> hist_list;

  //TH2D * h_phi_theta = new TH2D("phi_theta","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);
  //TH2D * h_poq_thetapq = new TH2D("poq_thetapq","#theta_{pq} vs. p/q ;p/q;#theta_{pq}",100,0.2,1.4,100,0,60);
  TH1D * h_Emiss = new TH1D("Emiss","E_{miss} ",100,-0.1,0.5);
  TH2D * h_omega_Emiss = new TH2D("omega_Emiss","E_{miss} vs. #omega;#omega;E_{miss}",100,0.0,2.5,100,-0.1,0.5);
  TH1D * h_xB = new TH1D("xB","x_{B} ",100,1.2,2);
  TH1D * h_Q2 = new TH1D("Q2","Q^{2} ",100,1,5);
  TH2D * h_omega_Q2 = new TH2D("omega_Q2","Q^{2} vs. #omega;#omega;Q^{2}",100,0.0,2.5,100,1,5);
  TH2D * h_Emiss_Mmiss = new TH2D("Emiss_Mmiss","M_{miss} vs. E_{miss};E_{miss};M_{miss}",100,-0.1,0.5,100,0.8,1.2);
  TH1D * h_pmiss = new TH1D("pmiss","p_{miss}",100,0.3,1);

  //hist_list.push_back(h_poq_thetapq);
  hist_list.push_back(h_Emiss      );
  hist_list.push_back(h_omega_Emiss);
  hist_list.push_back(h_xB         );
  hist_list.push_back(h_Q2         );
  hist_list.push_back(h_omega_Q2   );
  hist_list.push_back(h_Emiss_Mmiss);
  hist_list.push_back(h_pmiss      );

  clasAna.setEcalSFCuts();
  clasAna.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();

  clasAna.setVertexCuts();
  clasAna.setVertexCorrCuts();
  //clasAna.setDCEdgeCuts();
  
  while(chain.Next())
    {

      double weight = c12->mcevent()->getWeight(); //used if MC events have a weight 

      //Display completed  
      counter++;
      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      //auto electrons=c12->getByID(11);
      auto protons = clasAna.getByPid(2212);
      if(electrons.size() == 1 && protons.size() >= 1)
	{
	  SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );

	  clasAna.getLeadRecoilSRC(beam,target,el);
	  auto lead    = clasAna.getLeadSRC();
	  auto recoil  = clasAna.getRecoilSRC();

	  if(lead.size() == 1)
	    {
	      SetLorentzVector(lead_ptr,lead[0]);
	      TLorentzVector miss = q + target - lead_ptr; 
	      TLorentzVector miss_Am1 = q + nucleus_ptr - lead_ptr; 
	      double TB = miss_Am1.E() - miss_Am1.M();
	      double TP = lead_ptr.E() - lead_ptr.M();
	      double Emiss = q.E() - TP - TB;
	      //h_poq_thetapq->Fill(lead_ptr.Rho()/q.Rho(),lead_ptr.Angle(q.Vect())*180/M_PI);
	      h_Emiss->Fill(Emiss);
	      h_omega_Emiss->Fill(omega,Emiss);
	      h_xB->Fill(xB);
	      h_Q2->Fill(Q2);
	      h_omega_Q2->Fill(omega,Q2);
	      h_Emiss_Mmiss->Fill(Emiss,miss.M());
	      h_pmiss->Fill(miss.Rho());
	    }

	}
    }

  //clasAna.WriteDebugPlots();

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
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
  /////////////////////////////////////

  for(int i=0; i<hist_list.size(); i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    hist_list[i]->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

