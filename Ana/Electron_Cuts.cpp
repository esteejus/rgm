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
  clas12ana clasAna2;

  //Read in target parameter files                                                                                                                                                           
  clasAna.readInputParam("ana.par");
  clasAna.readEcalSFPar("paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("paramsPI_40Ca_x2.dat");
  clasAna.printParams();
  
  clasAna2.readInputParam("ana.par");
  clasAna2.readEcalSFPar("paramsSF_40Ca_x2.dat");
  clasAna2.readEcalPPar("paramsPI_40Ca_x2.dat");
  clasAna2.printParams();
  

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

  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  TH1D * h_Q2_bc = new TH1D("Q2_bc","Q^{2} ",1000,0, 5);
  TH1D * h_xB_bc = new TH1D("xB_bc","x_{B} ",1000,0, 2);
  TH2D * h_phi_theta_bc = new TH2D("phi_theta_bc","#phi_{e} vs. #theta_{e} ;#phi_{e};#theta_{e}",100,-180,180,100,5,40);

  TH1D * h_nphe_bc = new TH1D("nphe_bc","#Photo-electrons in HTCC;#Photo-electrons;Counts",40,0,40);

  TH1D * h_PCedep_bc = new TH1D("PCedep_bc","PCal E_{dep};PCal E_{dep} [GeV];Counts",100,0,0.6);

  TH2D * h_mom_SF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_SF_bc_%d",i+1);
    sprintf(temp_title,"p_{e} vs. Sampling Faction Sector=%d;Momentum [GeV];Sampling Fraction",i+1);
    h_mom_SF_bc[i] = new TH2D(temp_name,temp_title,100,0,7,100,0.1,0.35);
  }
  TH2D * h_mom_SF_ac[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"mom_SF_ac_%d",i+1);
    sprintf(temp_title,"p_{e} vs. Sampling Faction Sector=%d;Momentum [GeV];Sampling Fraction",i+1);
    h_mom_SF_ac[i] = new TH2D(temp_name,temp_title,100,0,7,100,0.1,0.35);
  }

  TH2D * h_PCedep_SF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"PCedep_SF_bc_%d",i+1);
    sprintf(temp_title,"PCal E_{dep} vs. Sampling Faction Sector=%d;PCal E_{dep} [GeV];Sampling Fraction",i+1);
    h_PCedep_SF_bc[i] = new TH2D(temp_name,temp_title,100,0,1.25,100,0.1,0.35);
  }
  TH2D * h_PCedep_SF_ac[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"PCedep_SF_ac_%d",i+1);
    sprintf(temp_title,"PCal E_{dep} vs. Sampling Faction Sector=%d;PCal E_{dep} [GeV];Sampling Fraction",i+1);
    h_PCedep_SF_ac[i] = new TH2D(temp_name,temp_title,100,0,1.25,100,0.1,0.35);
  }

  TH2D * h_PCSF_ECINSF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"PCSF_ECINSF_bc_%d",i+1);
    sprintf(temp_title,"PCal Sampling Fraction vs. Inner Cal Sampling Faction Sector=%d;PCal Sampling Fraction; Inner Cal Sampling Fraction",i+1);
    h_PCSF_ECINSF_bc[i] = new TH2D(temp_name,temp_title,100,0.0,0.35,100,0.0,0.35);
  }

  TH1D * h_vtz_e_bc = new TH1D("vtz_e_bc","Electron Z Vertex;Vertex [cm];Counts",100,-10,10);

  TH2D * h_Vcal_SF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"Vcal_SF_bc_%d",i+1);
    sprintf(temp_title,"ECAL V coordinate vs. Sampling Fraction Sector=%d;ECAL V coordinate;Sampling Fraction",i+1);
    h_Vcal_SF_bc[i] = new TH2D(temp_name,temp_title,60,0,30,100,0.1,0.35);
  }

  TH2D * h_Wcal_SF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"Wcal_SF_bc_%d",i+1);
    sprintf(temp_title,"ECAL W coordinate vs. Sampling Fraction Sector=%d;ECAL W coordinate;Sampling Fraction",i+1);
    h_Wcal_SF_bc[i] = new TH2D(temp_name,temp_title,60,0,30,100,0.1,0.35);
  }

  TH1D *h_DCedge_weight_bc[3][6];
  TH1D *h_DCedge_bc[3][6];
  for(int j=0; j<3; j++){
    for(int i=0; i<6; i++){
      sprintf(temp_name,"DCedge_weight_bc_%d_%d",j+1,i+1);
      sprintf(temp_title,"Distance from DC Edge Sector=%d Region=%d;Distance [cm];Average #chi^{2}/DoF",j+1,i+1);
      h_DCedge_weight_bc[j][i] = new TH1D(temp_name,temp_title,50,0,50);
    
      sprintf(temp_name,"DCedge_bc_%d_%d",j+1,i+1);
      h_DCedge_bc[j][i] = new TH1D(temp_name,temp_title,50,0,50);
    }
  }

  TH1D * h_Chi2DoF_bc[6];
  for(int i=0; i<6; i++){
    sprintf(temp_name,"Chi2DoF_bc_%d",i+1);
    sprintf(temp_title,"#chi^{2}/DoF Sector=%d;#chi^{2}/DoF;Counts",i+1);
    h_Chi2DoF_bc[i] = new TH1D(temp_name,temp_title,100,0,100);
  }

  clasAna2.setEcalSFCuts();
  clasAna2.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  //clasAna.setPidCuts();

  clasAna.setVertexCuts();
  //clasAna.setVertexCorrCuts();
  clasAna.setDCEdgeCuts();
  
  //clasAna.setVzcuts(-6,1);
  //clasAna.setVertexCorrCuts(-3,1);

  while(chain.Next())
    {

      double weight = c12->mcevent()->getWeight(); //used if MC events have a weight 

      //Display completed  
      counter++;
      clasAna.Run(c12);
      clasAna2.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);

      auto electrons_check = clasAna2.getByPid(11);

      if(electrons.size() == 1)
	{
	  SetLorentzVector(el,electrons[0]);
	  //	  SetLorentzVector(ptr,protons[0]);

	  TLorentzVector q = beam - el; //photon  4-vector            
          double Q2        = -q.M2(); // Q^2
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) ); //x-borken
	  h_Q2_bc->Fill(Q2);
	  h_xB_bc->Fill(xB);
	  h_phi_theta_bc->Fill(el.Phi()*180/M_PI,el.Theta()*180/M_PI);

	  int nphe = electrons[0]->che(HTCC)->getNphe();
	  h_nphe_bc->Fill(nphe);

	  double PCedep = electrons[0]->cal(PCAL)->getEnergy();
	  double ECINedep = electrons[0]->cal(ECIN)->getEnergy();
	  double ECOUTedep = electrons[0]->cal(ECOUT)->getEnergy();
	  h_PCedep_bc->Fill(PCedep);

	  double SF = (PCedep + ECINedep + ECOUTedep)/el.Rho();
	  int esector = electrons[0]->getSector();
	  h_mom_SF_bc[esector-1]->Fill(el.Rho(),SF);

	  h_PCedep_SF_bc[esector-1]->Fill(PCedep,SF);

	  ///////////////////////////////////////////
	  if(electrons_check.size() != 1){continue;}	

	  h_mom_SF_ac[esector-1]->Fill(el.Rho(),SF);
		  
	  h_PCedep_SF_ac[esector-1]->Fill(PCedep,SF);

	  h_PCSF_ECINSF_bc[esector-1]->Fill(PCedep/el.Rho(),ECINedep/el.Rho());

	  double vtz_e = electrons[0]->par()->getVz();
	  h_vtz_e_bc->Fill(vtz_e);

	  h_Vcal_SF_bc[esector-1]->Fill(electrons[0]->cal(PCAL)->getLv(),SF);

	  h_Wcal_SF_bc[esector-1]->Fill(electrons[0]->cal(PCAL)->getLw(),SF);

	  double DCedge[3];
	  DCedge[0]  = electrons[0]->traj(DC,6 )->getFloat("edge",electrons[0]->traj(DC,6 )->getIndex());
	  DCedge[1]  = electrons[0]->traj(DC,18)->getFloat("edge",electrons[0]->traj(DC,18)->getIndex());
	  DCedge[2]  = electrons[0]->traj(DC,36)->getFloat("edge",electrons[0]->traj(DC,36)->getIndex());
	  double Chi2DoF = electrons[0]->trk(DC)->getChi2()/electrons[0]->trk(DC)->getNDF();

	  for(int k=0; k<3; k++){
	    h_DCedge_bc[k][esector-1]->Fill(DCedge[k]);
	    h_DCedge_weight_bc[k][esector-1]->Fill(DCedge[k],Chi2DoF);
	  }
	  
	  h_Chi2DoF_bc[esector-1]->Fill(Chi2DoF);
	}
    }

  //clasAna.WriteDebugPlots();

  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  h_Q2_bc->Write();
  h_xB_bc->Write();


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

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Q2_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_xB_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_phi_theta_bc->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nphe_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_PCedep_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_SF_bc[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_PCedep_SF_bc[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  ////////////////////////////////////////////

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_mom_SF_ac[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_PCedep_SF_ac[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_PCSF_ECINSF_bc[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_vtz_e_bc->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_Vcal_SF_bc[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();    

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    h_Wcal_SF_bc[i]->Draw("colz");  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();    

  for(int j = 0; j < 3; j++){
    myCanvas->Divide(2,3);
    for(int i = 0; i < 6; i++){
      myCanvas->cd(i+1);
      h_DCedge_weight_bc[j][i]->Divide(h_DCedge_bc[j][i]);
      h_DCedge_weight_bc[j][i]->Draw();  
      h_DCedge_weight_bc[j][i]->SetMaximum(100);
      h_DCedge_weight_bc[j][i]->SetMinimum(0);
    }
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  }

  myCanvas->Divide(2,3);
  for(int i = 0; i < 6; i++){
    myCanvas->cd(i+1);
    myCanvas->SetLogy();  
    h_Chi2DoF_bc[i]->Draw();  
    myCanvas->SetLogy();  
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  
  /////////////////////////////////////
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

