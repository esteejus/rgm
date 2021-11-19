#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"
using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
	      rp->par()->getPz(),p4.M());

}

void e4nu_Monitoring(TString inFile = "", TString outputFile = "", double beamE = 0, double vz_min = -10, double vz_max = 10){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();

  Double_t mD = 1.8756;

  //initialising clas12writer with path to output file
  clas12databases::SetCCDBLocalConnection("ccdb.sqlite");
  clas12databases::SetRCDBRootConnection("rcdb.root");

  clas12root::HipoChain chain;

  chain.Add(inFile);

  chain.SetReaderTags({0});
  auto config_c12=chain.GetC12Reader();

  auto& c12=chain.C12ref();
  auto& rcdbData= config_c12->rcdb()->current();//struct with all relevent rcdb values        
  auto beam_energy = rcdbData.beam_energy/1000;
  cout << "Beam energy is " << beam_energy << endl;

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();

  //some particles
  TLorentzVector beam(0,0,beam_energy, beam_energy);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
   
  gBenchmark->Start("timer");
 
  TH1D *mult_p = new TH1D("mult_p","Proton Multiplicity;Multiplicity p;Counts",5,0,5);
  TH1D *mult_pip = new TH1D("mult_pip","#pi^{+} Multiplicity;Multiplicity #pi^{+};Counts",5,0,5);
  TH1D *mult_pim = new TH1D("mult_pim","#pi^{-} Multiplicity;Multiplicity #pi^{-};Counts",5,0,5);

  TH2D *q2_theta_h = new TH2D("q2_theta_h","Q^{2};Q^{2};Counts",400,0,2,40,0,40);

  TH2D *ecal_h = new TH2D("ecal_h",";Q^{2};Counts",400,0,7,100,0,2);

  TH1D *q2_mc_h = new TH1D("q2_mc_h",";Q^{2};Counts",400,0,2);

  TH1D *q2_h = new TH1D("q2_h",";Q^{2};Counts",400,0,2);
  TH2D *el_theta_mom = new TH2D("el_theta_mom","Electron Angle vs Momentum;Momentum (GeV/c);#theta",400,0.5,2.5,100,0,50);

  TH1D *p_chi2_cd = new TH1D("p_chi2_cd",";Proton Chi2;Counts",100,-4,4);
  TH1D *pip_chi2_cd = new TH1D("pip_chi2_cd",";Pip Chi2;Counts",100,-4,4);
  TH1D *pim_chi2_cd = new TH1D("pim_chi2_cd",";Pim Chi2;Counts",100,-4,4);

  TH1D *p_chi2_fd = new TH1D("p_chi2_fd",";Proton Chi2;Counts",100,-4,4);
  TH1D *pip_chi2_fd = new TH1D("pip_chi2_fd",";Pip Chi2;Counts",100,-4,4);
  TH1D *pim_chi2_fd = new TH1D("pim_chi2_fd",";Pim Chi2;Counts",100,-4,4);

  TH1D *energy_transfer_1 = new TH1D("energy_transfer_1","Energy transfer  5 < #theta < 15; #omega; Counts",1000,4.5,7);
  TH1D *energy_transfer_2 = new TH1D("energy_transfer_2","Energy transfer  15 < #theta < 25; #omega; Counts",1000,4.5,7);
  TH1D *energy_transfer_3 = new TH1D("energy_transfer_3","Energy transfer  25 < #theta < 35; #omega; Counts",1000,4.5,7);

  TH1D *invariant_mass = new TH1D("invariant_mass",";Mass (GeV/c);Counts",1000,0,2);

  TH1D *el_vz_h = new TH1D("el_vz_h",";Z-vertex (cm);Counts",100,-10,10);

  int counter = 0;

  double p_chi2_fd_cut = 1.;
  double p_chi2_cd_cut = 1.;

  double pip_chi2_fd_cut = 1.;
  double pip_chi2_cd_cut = 1.;

  double pim_chi2_fd_cut = 1.;
  double pim_chi2_cd_cut = 1.;



  while(chain.Next())
    {

      beam_energy = rcdbData.beam_energy/1000;
      beam.SetE(beam_energy);
      beam.SetPz(beam_energy);

      if( beamE > 1e-4)
	{
	  if( counter == 0)
	    cout<<"Beam energy set manually: "<<beamE<<endl;
	  beam.SetE(beamE);
	  beam.SetPz(beamE);
	}



      TLorentzVector el_mc;
      c12->mcparts()->setEntry(0);

      el_mc.SetXYZM(c12->mcparts()->getPx(), c12->mcparts()->getPy(),c12->mcparts()->getPz(),db->GetParticle(11)->Mass());


      // get particles by type
      auto electrons=c12->getByID(11);
      auto protons=c12->getByID(2212);
      auto pips=c12->getByID(211);
      auto pims=c12->getByID(-211);
      

      if(electrons.size()==1)
	{
	  
	  double energy =  electrons[0]->cal(PCAL)->getEnergy() +  electrons[0]->cal(ECIN)->getEnergy() +  electrons[0]->cal(ECOUT)->getEnergy();
	  
	  bool ecal = ( (electrons[0]->cal(ECIN)->getLv() >= 14 && electrons[0]->cal(ECIN)->getLw() >= 14) && (energy > 0.18 && energy < 0.28) && electrons[0]->getP() > 0.5 && electrons[0]->getP() < 10 );


	  SetLorentzVector(el,electrons[0]);


	  TLorentzVector q = beam - el;                  //photon  4-vector                  
	  double q2       = -q.M2();
	  double x_b      = q2/(2 * mass_p * (beam.E() - el.E()) ); //x-borken 
	  double vz_e   = electrons[0]->par()->getVz();

	  //Electron quality cuts
	  bool targ_1 = abs( vz_e - (-1.48)) < 1 * 1.25;
	  bool targ_2 = abs( vz_e - (-6.3))  < 1 * 1.25;

	  //	  ecal_h->Fill(electrons[0]->getP(),energy/electrons[0]->getP());
	    //	  if( !(ecal && (targ_1 || targ_2) ) )//empty rgk run

	  if( !(ecal && vz_e < vz_max && vz_e > vz_min))
	    continue;

	  el_vz_h->Fill(vz_e);
	  ecal_h->Fill(electrons[0]->getP(),energy/electrons[0]->getP());

	  //remove FT (seems to be still in the MC; FT is not used in RGM)
	  if( electrons[0]->getRegion() == 1000)
	    continue;

	  //Electron kinematics monitoring
	  //(e,e') cross sections

	  el_theta_mom->Fill(el.P(),el.Theta()*TMath::RadToDeg());
	  q2_h->Fill(q2);

	  if(el.Theta()*TMath::RadToDeg() < 15 && el.Theta()*TMath::RadToDeg() > 5)
	    energy_transfer_1->Fill(q.E());
	  else if(el.Theta()*TMath::RadToDeg() < 25  && el.Theta()*TMath::RadToDeg() >= 15)
	    energy_transfer_2->Fill(q.E());
	  else if(el.Theta()*TMath::RadToDeg() < 35  && el.Theta()*TMath::RadToDeg() >= 25)
	    energy_transfer_3->Fill(q.E());
	  
	  
	  //MULTIPLICITIES 
	  //Monitor the proton, pion multiplicities 
	  int num_p = 0;
	  int num_pip = 0;
	  int num_pim = 0;
	  
	  for(auto &p : protons)
	    {
	      double p_chi2_cd_v = 999999;
	      double p_chi2_fd_v = 999999;
	      
	      if(p->getRegion()==CD)
		{
		  p_chi2_cd->Fill(p->par()->getChi2Pid());
		  p_chi2_cd_v = p->par()->getChi2Pid();
		}
	      if(p->getRegion()==FD)
		{
		  p_chi2_fd->Fill(p->par()->getChi2Pid());
		  p_chi2_fd_v = p->par()->getChi2Pid();
		}
	      
	      if( abs(p_chi2_fd_v ) <  p_chi2_fd_cut || abs(p_chi2_cd_v ) <  p_chi2_cd_cut )
		num_p++;
	      
	    }
	  
	  
	  for(auto &p : pims)
	    {
	      double pim_chi2_cd_v = 999999;
	      double pim_chi2_fd_v = 999999;
	      
	      if(p->getRegion()==CD)
		{
		  pim_chi2_cd->Fill(p->par()->getChi2Pid());
		  pim_chi2_cd_v = p->par()->getChi2Pid();
		}
	      
	      if(p->getRegion()==FD)
		{
		  pim_chi2_fd->Fill(p->par()->getChi2Pid());
		  pim_chi2_fd_v = p->par()->getChi2Pid();
		}
	      
	      if( abs(pim_chi2_fd_v ) <  pim_chi2_fd_cut || abs(pim_chi2_cd_v ) <  pim_chi2_cd_cut )
		num_pim++;
	      
	    }

	  
	  for(auto &p : pips)
	    {
	      double pip_chi2_cd_v = 999999;
	      double pip_chi2_fd_v = 999999;
	      
	      if(p->getRegion()==CD)
		{
		  pip_chi2_cd->Fill(p->par()->getChi2Pid());
		  pip_chi2_cd_v = p->par()->getChi2Pid();
		}           
	      if(p->getRegion()==FD)
		{
		  pip_chi2_fd->Fill(p->par()->getChi2Pid());
		  pip_chi2_fd_v = p->par()->getChi2Pid();
		}
	      
	      if( abs(pip_chi2_fd_v ) < pip_chi2_fd_cut || abs(pip_chi2_cd_v ) <  pip_chi2_cd_cut )
		num_pip++;
	      
	    }
	  
	  
	  mult_p->Fill(num_p);
	  mult_pim->Fill(num_pim);
	  mult_pip->Fill(num_pip);
	  
	  
	  //INVARIANT MASS
	  //Calculating the Invariant mass of the p + pi- system 
	  for(int iPr = 0; iPr < protons.size(); iPr++)
	    {
	      for(int iPim = 0; iPim < pims.size(); iPim++)
		{
		
		  SetLorentzVector(pr,protons[iPr]);
		  SetLorentzVector(pim,pims[iPim]);
		  
		  double vz_e_p = electrons[0]->par()->getVz() - protons[iPr]->par()->getVz();
		  TLorentzVector inv_mass = pr + pim; //missing 4-vector
		  
		  invariant_mass->Fill(inv_mass.M());
		}
	    }

	  
	}

      counter++;
      
    }
  

  TFile *outFile = new TFile(outputFile + ".root","RECREATE");

  TCanvas *c1 = new TCanvas("c1","c1",800,1000);
  TString fileName = outputFile + ".pdf";

  c1->Divide(2,3);
  c1->cd(1);
  q2_h->Draw();

  c1->cd(2);
  el_theta_mom -> Draw("colz");

  c1->Print(fileName+"(");
  c1->Clear();


  c1->Divide(2,3);
  c1->cd(1);
  energy_transfer_1 -> Draw();

  c1->cd(2);
  energy_transfer_2 -> Draw();

  c1->cd(3);
  energy_transfer_3 -> Draw();

  c1->Print(fileName);
  c1->Clear();


  c1->Divide(2,3);
  c1->cd(1);
  c1->cd(1)->SetLogy();
  mult_p->Draw();

  c1->cd(3);
  c1->cd(3)->SetLogy();
  mult_pip->Draw();

  c1->cd(4);
  c1->cd(4)->SetLogy();
  mult_pim->Draw();

  c1->Print(fileName+")");
  c1->Clear();






  mult_p -> Write();
  mult_pip -> Write();
  mult_pim -> Write();

  q2_theta_h-> Write();

  q2_h-> Write();
  q2_mc_h-> Write();

  el_theta_mom -> Write();

  p_chi2_cd -> Write();
  pip_chi2_cd -> Write();
  pim_chi2_cd -> Write();

  p_chi2_fd -> Write();
  pip_chi2_fd -> Write();
  pim_chi2_fd -> Write();

  el_vz_h->Write();

  energy_transfer_1 -> Write();
  energy_transfer_2 -> Write();
  energy_transfer_3 -> Write();
  invariant_mass -> Write();

  ecal_h->Write();
  
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " read events = "<<counter<<endl;
  }
