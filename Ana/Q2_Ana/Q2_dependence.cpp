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

using namespace std;
using namespace clas12;

const double c = 29.9792458;

int binpmiss(double pmiss){
  if( (pmiss>0.4) && (pmiss<0.55)){
    return 0;
  }
  else if( (pmiss>0.55) && (pmiss<0.7)){
    return 1;
  }
  else if( (pmiss>0.7) && (pmiss<0.85)){
    return 2;
  }
  else if( (pmiss>0.85) && (pmiss<1.0)){
    return 3;	 
  }
  return -1;
}

int binQ2(double q2){
  if(q2<1.65){return 0;}
  else if(q2<1.80){return 1;}
  else if(q2<1.95){return 2;}
  else if(q2<2.10){return 3;}
  else if(q2<2.25){return 4;}
  else if(q2<2.40){return 5;}
  else if(q2<2.70){return 6;}
  else if(q2<3.00){return 7;}
  else if(q2<3.50){return 8;}
  else{return 9;}
}

double binEdges_Q2[] = {1.5,1.65,1.80,1.95,2.10,2.25,2.40,2.70,3.00,3.50,5.0};
int binEdgeslength_Q2 = sizeof(binEdges_Q2)/sizeof(binEdges_Q2[0]) -1;

double binEdges[] = { 0.35, 0.38, 0.41, 0.44, 0.47, 0.5, 0.53, 0.56, 0.59, 0.62, 0.65, 0.68, 0.71, 0.76, 0.80, 0.85, 0.90, 0.95, 1.0};

int binEdgeslength = sizeof(binEdges)/sizeof(binEdges[0]) -1;


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
  clasAna.readInputParam("/w/hallb-scshelf2102/clas12/users/awild/RGM/rgm/Ana/ana.par");
  clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas12/users/awild/RGM/rgm/Ana/paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("/w/hallb-scshelf2102/clas12/users/awild/RGM/rgm/Ana/paramsPI_40Ca_x2.dat");
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 3; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

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
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  //////////////////////////////////
  char temp_name[100];
  char temp_title[100];

  
  vector<TH1*> hist_list;

  ////////////////////////////////////////////////
  //ep
  ////////////////////////////////////////////////
  TH1D * h_pmiss_SRC = new TH1D("pmiss_SRC","p_{miss} SRC;p_{miss};Counts",binEdgeslength,binEdges);
  hist_list.push_back(h_pmiss_SRC);
  TH1D * h_Q2_SRC_Q2bin[10];
  TH1D * h_pmiss_SRC_Q2bin[10];
  for(int i=0; i<10; i++){
    sprintf(temp_name,"Q2_SRC_Q2bin_%d",i+1);
    sprintf(temp_title,"Q^{2} Q^{2}bin=%d;Q^{2};Counts",i+1);
    h_Q2_SRC_Q2bin[i] = new TH1D(temp_name,temp_title,90,1.5,6.0);
    hist_list.push_back(h_Q2_SRC_Q2bin[i]);

    sprintf(temp_name,"pmiss_SRC_Q2bin_%d",i+1);
    sprintf(temp_title,"p_{miss} Q^{2}bin=%d;p_{miss};Counts",i+1);
    h_pmiss_SRC_Q2bin[i] = new TH1D(temp_name,temp_title,binEdgeslength,binEdges);
    hist_list.push_back(h_pmiss_SRC_Q2bin[i]);
  }
  TH1D * h_Q2_SRC_pmissbin[4];
  for(int i=0; i<4; i++){
    sprintf(temp_name,"Q2_SRC_pmissbin_%d",i+1);
    sprintf(temp_title,"Q^{2} p_{miss}bin=%d;Q^{2};Counts",i+1);
    h_Q2_SRC_pmissbin[i] = new TH1D(temp_name,temp_title,binEdgeslength_Q2,binEdges_Q2);
    hist_list.push_back(h_Q2_SRC_pmissbin[i]);
  }

  TH1D * h_emiss_SRC_pmissbin[4];
  for(int i=0; i<4; i++){
    sprintf(temp_name,"emiss_SRC_pmissbin_%d",i+1);
    sprintf(temp_title,"E_{miss} p_{miss}bin=%d;E_{miss};Counts",i+1);
    h_emiss_SRC_pmissbin[i] = new TH1D(temp_name,temp_title,100,-0.1,0.5);
    hist_list.push_back(h_emiss_SRC_pmissbin[i]);
  }

  ////////////////////////////////////////////////
  //epp
  ////////////////////////////////////////////////
  TH1D * h_pmiss_Rec = new TH1D("pmiss_Rec","p_{miss} Rec;p_{miss};Counts",binEdgeslength,binEdges);
  hist_list.push_back(h_pmiss_Rec);
  TH1D * h_p_z_cm_Rec = new TH1D("p_z_cm_Rec","p_{z,C.M.};p_{z,C.M.};Counts",50,-1,1);
  hist_list.push_back(h_p_z_cm_Rec);
  TH1D * h_p_y_cm_Rec = new TH1D("p_y_cm_Rec","p_{#perp,C.M.};p_{#perp,C.M.};Counts",50,-1,1);
  hist_list.push_back(h_p_y_cm_Rec);
  TH1D * h_p_x_cm_Rec = new TH1D("p_x_cm_Rec","p_{x,C.M.};p_{x,C.M.};Counts",50,-1,1);
  hist_list.push_back(h_p_x_cm_Rec);

  int bins_Q2bin[10] = {35,35,35,35,35,25,25,25,15,15};
  TH1D * h_Q2_Rec_Q2bin[10];
  TH1D * h_pmiss_Rec_Q2bin[10];
  TH1D * h_p_x_cm_Rec_Q2bin[10];
  TH1D * h_p_y_cm_Rec_Q2bin[10];
  TH1D * h_p_z_cm_Rec_Q2bin[10];
  for(int i=0; i<10; i++){
    sprintf(temp_name,"Q2_Rec_Q2bin_%d",i+1);
    sprintf(temp_title,"Q^{2} Q^{2}bin=%d;Q^{2};Counts",i+1);
    h_Q2_Rec_Q2bin[i] = new TH1D(temp_name,temp_title,90,1.5,6.0);
    hist_list.push_back(h_Q2_Rec_Q2bin[i]);

    sprintf(temp_name,"pmiss_Rec_Q2bin_%d",i+1);
    sprintf(temp_title,"p_{miss} Q^{2}bin=%d;p_{miss};Counts",i+1);
    h_pmiss_Rec_Q2bin[i] = new TH1D(temp_name,temp_title,binEdgeslength,binEdges);
    hist_list.push_back(h_pmiss_Rec_Q2bin[i]);

    sprintf(temp_name,"p_x_cm_Rec_Q2bin_%d",i+1);
    sprintf(temp_title,"p_{x,cm} Q^{2}bin=%d;p_{x,cm};Counts",i+1);
    h_p_x_cm_Rec_Q2bin[i] = new TH1D(temp_name,temp_title,bins_Q2bin[i],-1.0,1.0);
    hist_list.push_back(h_p_x_cm_Rec_Q2bin[i]);

    sprintf(temp_name,"p_y_cm_Rec_Q2bin_%d",i+1);
    sprintf(temp_title,"p_{y,cm} Q^{2}bin=%d;p_{y,cm};Counts",i+1);
    h_p_y_cm_Rec_Q2bin[i] = new TH1D(temp_name,temp_title,bins_Q2bin[i],-1.0,1.0);
    hist_list.push_back(h_p_y_cm_Rec_Q2bin[i]);

    sprintf(temp_name,"p_z_cm_Rec_Q2bin_%d",i+1);
    sprintf(temp_title,"p_{z,cm} Q^{2}bin=%d;p_{z,cm};Counts",i+1);
    h_p_z_cm_Rec_Q2bin[i] = new TH1D(temp_name,temp_title,bins_Q2bin[i],-1.0,1.0);
    hist_list.push_back(h_p_z_cm_Rec_Q2bin[i]);
  }
  TH1D * h_Q2_Rec_pmissbin[4];
  for(int i=0; i<4; i++){
    sprintf(temp_name,"Q2_Rec_pmissbin_%d",i+1);
    sprintf(temp_title,"Q^{2} p_{miss}bin=%d;Q^{2};Counts",i+1);
    h_Q2_Rec_pmissbin[i] = new TH1D(temp_name,temp_title,binEdgeslength_Q2,binEdges_Q2);
    hist_list.push_back(h_Q2_Rec_pmissbin[i]);
  }

  TH1D * h_emiss_Rec_pmissbin[4];
  for(int i=0; i<4; i++){
    sprintf(temp_name,"emiss_Rec_pmissbin_%d",i+1);
    sprintf(temp_title,"E_{miss} p_{miss}bin=%d;E_{miss};Counts",i+1);
    h_emiss_Rec_pmissbin[i] = new TH1D(temp_name,temp_title,100,-0.1,0.5);
    hist_list.push_back(h_emiss_Rec_pmissbin[i]);
  }

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  clasAna.setEcalSFCuts();
  clasAna.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();

  clasAna.setVertexCuts();
  clasAna.setVertexCorrCuts();
  clasAna.setDCEdgeCuts();

  clasAna.setVzcuts(-6,1);
  clasAna.setVertexCorrCuts(-3,1);
  
  clasAna.setPidCuts(false); //clas chi2pid
  clasAna.setProtonPidCuts(true); //tof vs mom pid (proton)

  double num = 0;
  double den = 0;
  
  while(chain.Next())
    {

      double weight = c12->mcevent()->getWeight(); //used if MC events have a weight 
      weight = 1;
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      //Display completed  
      counter++;
      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      if(electrons.size() == 1 && protons.size() >= 1)
	{

	  SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2        = -q.M2();
	  double omega = q.E();
          double xB       = Q2/(2 * mass_p * (beam.E() - el.E()) );

	  clasAna.getLeadRecoilSRC(beam,deut_ptr,el);
	  auto lead    = clasAna.getLeadSRC();
	  auto recoil  = clasAna.getRecoilSRC();

	  if(lead.size() == 1)
	    {
	      den+=1.0;
	      SetLorentzVector(lead_ptr,lead[0]);
	      TLorentzVector miss = q + deut_ptr - lead_ptr;
	      TLorentzVector neg_miss = -miss;
	      bool rec = false;
	      int bp = binpmiss(miss.P());
	      double ei = lead_ptr.E() - omega;
	      double emiss = mass_p - ei;
	      
	      if(recoil.size() == 1){rec = true;}
	      if(miss.P()<0.3){continue;}

	      h_pmiss_SRC->Fill(miss.P(),weight);
	      h_Q2_SRC_Q2bin[binQ2(Q2)]->Fill(Q2,weight);
	      h_pmiss_SRC_Q2bin[binQ2(Q2)]->Fill(miss.P(),weight);	      
	      if(bp!=-1){
		h_Q2_SRC_pmissbin[bp]->Fill(Q2);
		h_emiss_SRC_pmissbin[bp]->Fill(emiss);
	      }

	      if(recoil.size() == 1){

		num+=1.0;
		SetLorentzVector(recoil_ptr,recoil[0]);
		TVector3 v_miss = neg_miss.Vect();
		TVector3 v_rec  = recoil_ptr.Vect();
		TVector3 v_rel  = (v_miss - v_rec) * 0.5;
		TVector3 v_cm   = v_miss + v_rec;

		TVector3 vt = v_miss.Unit();
		TVector3 vy = v_miss.Cross(q.Vect()).Unit();
		TVector3 vx = vt.Cross(vy).Unit();

		h_pmiss_Rec->Fill(miss.P(),weight);
		h_Q2_Rec_Q2bin[binQ2(Q2)]->Fill(Q2,weight);
		h_p_z_cm_Rec->Fill(v_cm.Dot(vt),weight);
		h_p_y_cm_Rec->Fill(v_cm.Dot(vy),weight);
		h_p_x_cm_Rec->Fill(v_cm.Dot(vx),weight);
		if(bp!=-1){
		  h_Q2_Rec_pmissbin[bp]->Fill(Q2);
		  h_emiss_Rec_pmissbin[bp]->Fill(emiss);
		}
		
		h_pmiss_Rec_Q2bin[binQ2(Q2)]->Fill(miss.P(),weight);
		h_p_z_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vt),weight);
		h_p_y_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vy),weight);
		h_p_x_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vx),weight);

	      }
	    }

	}
    }

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  } 

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

  for(int i=0; i<hist_list.size(); i++){
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    hist_list[i]->Draw("colz");
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();  
  } 

  h_pmiss_Rec->Divide(h_pmiss_SRC);
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_Rec->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  

  for(int i = 0; i < 10; i ++){
    h_pmiss_Rec_Q2bin[i]->Divide(h_pmiss_SRC_Q2bin[i]);
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    h_pmiss_Rec_Q2bin[i]->Draw();
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();      
  }

  for(int i = 0; i < 4; i ++){
    h_Q2_Rec_pmissbin[i]->Divide(h_Q2_SRC_pmissbin[i]);
    myCanvas->Divide(1,1);
    myCanvas->cd(1);
    h_Q2_Rec_pmissbin[i]->Draw();
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();      
  }
  

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

