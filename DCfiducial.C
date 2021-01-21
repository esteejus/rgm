#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TCutG.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"

using namespace clas12;


void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
             rp->par()->getPz(),p4.M());

}


TVector3 rotate(TVector3 vec, int sector) 
{
  double rot_ang = -(sector -1)*60 *TMath::DegToRad();

  vec.RotateZ(rot_ang);
  
  return vec;
}

void DCfiducial(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();


  /////////////////////////////////////
  //ignore this just getting file name!
  TString inputFile;
  TString outputFile;

  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains(".hipo"))){
      inputFile=opt(5,opt.Sizeof());
    }
  }
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////

  //  outputFile = inputFile(inputFile.Index("inc_")+4,inputFile.Index(".hipo"));
  //  outputFile = inputFile(inputFile.Index("qe")+2,inputFile.Index(".hipo"));
  outputFile = inputFile(inputFile.Index("recon_qe")+8,inputFile.Index(".hipo"));
  cout<<"Analysing hipo file "<<inputFile<<endl;

  TFile *out = new TFile("/work/clas12/users/esteejus/output/DCFiducial/MC/hist"+outputFile+".root","RECREATE");
  cout <<"Out file" <<"/work/clas12/users/esteejus/output/DCFiducial/MC/hist"+outputFile+".root"<<endl;

  TTree *tree = new TTree("tree","MC data");

  Int_t sector = -99, region = 0;
  TVector3 mom_in(0,0,0);
  TVector3 mom_out(0,0,0);
  Double_t region1_x, region1_y; 
  Double_t region2_x, region2_y; 
  Double_t region3_x, region3_y; 
  Double_t chi2_ndf = 0;

  tree->Branch("sector",&sector, "sector/I");
  tree->Branch("region1_x",&region1_x, "region1_x/D");
  tree->Branch("region1_y",&region1_y, "region1_y/D");
  tree->Branch("region2_x",&region2_x, "region2_x/D");
  tree->Branch("region2_y",&region2_y, "region2_y/D");
  tree->Branch("region3_x",&region3_x, "region3_x/D");
  tree->Branch("region3_y",&region3_y, "region3_y/D");

  tree->Branch("mom_in","TVector3", &mom_in);
  tree->Branch("mom_out","TVector3", &mom_out);


  TFile *cutf = TFile::Open("cuts.root");
  TCutG *elec_cut = (TCutG *)cutf->Get("elec");
  TCutG *ftof_p = (TCutG *)cutf->Get("ftof_p");
  TCutG *ctof_p = (TCutG *)cutf->Get("ctof_p");
   

  auto pdg_db=TDatabasePDG::Instance();
  TLorentzVector el(0,0,0,pdg_db->GetParticle(11)->Mass());

  TH2D *el_angle = new TH2D("el_angle","Electron angles",180,0,180,360,-180,180);

  TH2D *hit_map[4][7]; //3 regions 6 sectors index starts at 0 for region 1....
  TH2D *hit_map_rot[4][7]; //3 regions 6 sectors index starts at 0 for region 1....
  TH2D *hit_map_angle_rot[4][7]; //3 regions 6 sectors index starts at 0 for region 1....

  TH2D *chi2_ndf_rot[4][7]; //3 regions 6 sectors index starts at 0 for region 1....
  TH2D *chi2_ndf_angle_rot[4][7]; //3 regions 6 sectors index starts at 0 for region 1....

  for(int i = 1; i <= 3; i++)
    {
      for(int j = 1; j <= 6; j++)
	{
	  hit_map[i][j] = new TH2D(Form("hit_map_%d_%d",i,j), Form("Region %d Sector %d",i,j),600,-300,300,600,-300,300);
	  hit_map_rot[i][j] = new TH2D(Form("hit_map_rot_%d_%d",i,j), Form("Region %d Sector %d (Rotated)",i,j),150,0,150,300,-150,150);
	  hit_map_angle_rot[i][j] = new TH2D(Form("hit_map_angle_rot_%d_%d",i,j), Form("Region %d Sector %d Angles(Rotated)",i,j),200,0,100,240,-60,60);

	  chi2_ndf_rot[i][j] = new TH2D(Form("chi2_ndf_rot_%d_%d",i,j), Form("Region %d Sector %d Chi2 (Rotated)",i,j),150,0,150,300,-150,150);
	  chi2_ndf_angle_rot[i][j] = new TH2D(Form("chi2_ndf_angle_rot_%d_%d",i,j), Form("Region %d Sector %d Chi2 Angles(Rotated)",i,j),200,0,100,240,-60,60);

	}
    }

  auto *el_xy_1 = new TH2F("el_xy_1","Region 1 e^{-} xy",400,-200,200,400,-200,200);
  auto *el_xy_2 = new TH2F("el_xy_2","Region 2 e^{-} xy",700,-350,350,400,-350,350);
  auto *el_xy_3 = new TH2F("el_xy_3","Region 3 e^{-} xy",1000,-500,500,400,-500,500);
  
  gBenchmark->Start("timer");
  int counter=0;
 
  clas12reader c12(inputFile.Data());
  //Add some event Pid based selections                                                                                                                                                     
  //  c12.addAtLeastPid(11,1); //at least 1 electron                                                                                                                                      
  //c12.addExactPid(11,1);    //exactly 1 electron                                                                                                                                          
  //c12.addExactPid(211,1);    //exactly 1 pi+                                                                                                                                              
  //c12.addExactPid(-211,1);    //exactly 1 pi-                                                                                                                                             
  //c12.addExactPid(2212,1);    //exactly 1 proton                                                                                                                                          
  //c12.addExactPid(22,2);    //exactly 2 gamma                                                                                                                                             
  //////c12.addZeroOfRestPid();  //nothing else                                                                                                                                             
  //////c12.useFTBased(); //and use the Pids from RECFT   

  TVector3 region_1(0,0,0);      
  TVector3 region_2(0,0,0);      
  TVector3 region_3(0,0,0);      

      while(c12.next()==true)
	{
	  
	  mom_in.SetXYZ(0,0,0);
	  c12.mcparts()->setEntry(0);
	  mom_in.SetXYZ(c12.mcparts()->getPx(), c12.mcparts()->getPy(),c12.mcparts()->getPz());                                                                                

	  for(auto& p : c12.getDetParticles())
	    {
	      //  get predefined selected information
	      p->getTime();
	      p->getDetEnergy();
	      p->getDeltaEnergy();
	      
	      //check trigger bits
	      // if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
	      // else cout<<"NOT"<<endl;
	      
	      // get any detector information (if exists for this particle)
	      // there should be a get function for any entry in the bank
	      //	      int LAYER = 18;
	      p->trk(DC)->getSector();
	      p->trk(DC)->getChi2();
	      /*
	      p->traj(DC,LAYER)->getPindex();
	      p->traj(DC,LAYER)->getDetector();
	      p->traj(DC,LAYER)->getLayer();
	      p->traj(DC,LAYER)->getCx();
	      p->traj(DC,LAYER)->getCy();
	      p->traj(DC,LAYER)->getCz();
	      p->traj(DC,LAYER)->getX();
	      p->traj(DC,LAYER)->getY();
	      p->traj(DC,LAYER)->getZ();
	      p->traj(DC,LAYER)->getPath();  
	      */

	      //MC::Lund                                                                              
	      // For the jth particle in the lund file                                                
	      // 0:e , 1:p ,2:n in this Lund file                                                    /* 

	      //    c12.mcparts()->setEntry(0);                                                            
	      /*
	      c12.mcparts()->getPid();                                                                
	      c12.mcparts()->getPx();                                                                 
	      c12.mcparts()->getPy();                                                                 
	      c12.mcparts()->getPz();                                                                 
	      c12.mcparts()->getVx();                                                                 
	      c12.mcparts()->getVy();                                                                 
	      c12.mcparts()->getVz();                                                                 
	      c12.mcparts()->getMass();                                                               */
	      sector = -99; region = -99;
	      region1_x = -99; region1_y = -99; region2_x = -99; region2_y = -99; region3_x = -99; region3_y = -99;
	      mom_out.SetXYZ(0,0,0);
	      chi2_ndf = -99;

	      // get particles by type
	      auto electrons=c12.getByID(11);
	      //	      auto protons=c12.getByID(2212);


	      if(electrons.size() == 1)
		{

		  int iEl = 0;
		  mom_out.SetXYZ(electrons[iEl]->par()->getPx(),electrons[iEl]->par()->getPy(),electrons[iEl]->par()->getPz());
		  chi2_ndf =  electrons[iEl]->trk(DC)->getChi2()/electrons[iEl]->trk(DC)->getNDF();
		  sector = electrons[iEl]->trk(DC)->getSector();

		  SetLorentzVector(el,electrons[iEl]);
		  el_angle->Fill(el.Theta()*TMath::RadToDeg(), el.Phi()*TMath::RadToDeg());
		  region_1.SetXYZ(electrons[iEl]->traj(DC,6)->getX(),electrons[iEl]->traj(DC,6)->getY(),electrons[iEl]->traj(DC,6)->getZ());
		  region_2.SetXYZ(electrons[iEl]->traj(DC,18)->getX(),electrons[iEl]->traj(DC,18)->getY(),electrons[iEl]->traj(DC,18)->getZ());
		  region_3.SetXYZ(electrons[iEl]->traj(DC,36)->getX(),electrons[iEl]->traj(DC,36)->getY(),electrons[iEl]->traj(DC,36)->getZ());

		  auto region_1_rot = rotate(region_1,sector);
		  auto region_2_rot = rotate(region_2,sector);
		  auto region_3_rot = rotate(region_3,sector);

		  region1_x = region_1_rot.X();
		  region1_y = region_1_rot.Y();

		  region2_x = region_2_rot.X();
		  region2_y = region_2_rot.Y();

		  region3_x = region_3_rot.X();
		  region3_y = region_3_rot.Y();

		  if(electrons[iEl]->getP() > 1)
		    {

		      if( abs(region_1.X()) > 1e-3 && abs(region_1.Y()) > 1e-3 )
			{
			  el_xy_1->Fill(region_1.X(), region_1.Y());
			  hit_map[1][sector]->Fill(region_1.X(),region_1.Y());
			  hit_map_rot[1][sector]->Fill(region_1_rot.X(),region_1_rot.Y());
			  hit_map_angle_rot[1][sector]->Fill(region_1_rot.Theta()*TMath::RadToDeg(),region_1_rot.Phi()*TMath::RadToDeg());

			  chi2_ndf_rot[1][sector]->Fill(region_1_rot.X(),region_1_rot.Y(),chi2_ndf);
			  chi2_ndf_angle_rot[1][sector]->Fill(region_1_rot.Theta()*TMath::RadToDeg(),region_1_rot.Phi()*TMath::RadToDeg(),chi2_ndf);
			}

		      if( abs(region_2.X()) > 1e-3 && abs(region_2.Y()) > 1e-3 )
			{
			  el_xy_2->Fill(region_2.X(), region_2.Y());
			  hit_map[2][sector]->Fill(region_2.X(),region_2.Y());
			  hit_map_rot[2][sector]->Fill(region_2_rot.X(),region_2_rot.Y());
			  hit_map_angle_rot[2][sector]->Fill(region_2_rot.Theta()*TMath::RadToDeg(),region_2_rot.Phi()*TMath::RadToDeg());

			  chi2_ndf_rot[2][sector]->Fill(region_2_rot.X(),region_2_rot.Y(),chi2_ndf);
			  chi2_ndf_angle_rot[2][sector]->Fill(region_2_rot.Theta()*TMath::RadToDeg(),region_2_rot.Phi()*TMath::RadToDeg(),chi2_ndf);
			}

		      if( abs(region_3.X()) > 1e-3 && abs(region_3.Y()) > 1e-3 )
			{
			  el_xy_3->Fill(region_3.X(), region_3.Y());
			  hit_map[3][sector]->Fill(region_3.X(),region_3.Y());
			  hit_map_rot[3][sector]->Fill(region_3_rot.X(),region_3_rot.Y());
			  hit_map_angle_rot[3][sector]->Fill(region_3_rot.Theta()*TMath::RadToDeg(),region_3_rot.Phi()*TMath::RadToDeg());

			  chi2_ndf_rot[3][sector]->Fill(region_3_rot.X(),region_3_rot.Y(),chi2_ndf);
			  chi2_ndf_angle_rot[3][sector]->Fill(region_3_rot.Theta()*TMath::RadToDeg(),region_3_rot.Phi()*TMath::RadToDeg(),chi2_ndf);
			}

		    }

		}// one electron condition
	      
	      counter++;
	    }//particle loop

	  tree->Fill();
	}
      
      gBenchmark->Stop("timer");
      gBenchmark->Print("timer");

      out->cd();
      tree->Write();
      el_xy_1->Write();
      el_xy_2->Write();
      el_xy_3->Write();

      for(int i = 1; i <= 3; i++)
	{
	  for(int j = 1; j <= 6; j++)
	    {
	      chi2_ndf_rot[i][j]->Divide(hit_map_rot[i][j]);
	      chi2_ndf_angle_rot[i][j]->Divide(hit_map_angle_rot[i][j]);

	      hit_map[i][j]->Write();
	      hit_map_rot[i][j]->Write();
	      hit_map_angle_rot[i][j]->Write();	    
	      chi2_ndf_rot[i][j]->Write();
	      chi2_ndf_angle_rot[i][j]->Write();

	    }
	}
      
      el_angle->Write();
      el_angle->Draw("colz");
      out->Close();

      auto finish = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = finish - start;
      std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
      
}
