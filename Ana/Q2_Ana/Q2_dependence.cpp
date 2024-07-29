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
#include "TGraph.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "reweighter.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;

int binpmiss(double pmiss){
  if( (pmiss>0.4) && (pmiss<0.5)){
    return 0;
  }
  else if( (pmiss>0.5) && (pmiss<0.6)){
    return 1;
  }
  else if( (pmiss>0.6) && (pmiss<0.7)){
    return 2;
  }
  else if( (pmiss>0.7) && (pmiss<1.0)){
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
  std::cerr << "Usage: ./code isMC outputfile.root outputfile.pdf inputfiles.hipo \n\n\n";

}


int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }


  int isMC = atoi(argv[1]);
  TString outFile = argv[2];
  char * pdfFile = argv[3];

  cout<<"Ouput file "<< outFile <<endl;
  cout<<"Ouput PDF file "<< pdfFile <<endl;
  
  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;
  double beam_E = 5.98;
  const double me = 0.000511;
  const double mU = 0.9314941024;
  const double m_4He = 4.00260325415 * mU - 2*me;
  reweighter newWeight(beam_E,6,6);

  
  clas12ana clasAna;
  clasAna.printParams();
    
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();
  
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




  TH1D * h_pLead_ep = new TH1D("pLead_ep","p_{Lead} SRC;p_{Lead};Counts",100,0,3);
  hist_list.push_back(h_pLead_ep);
  TH1D * h_pMiss_ep = new TH1D("pMiss_ep","p_{Miss} SRC;p_{Miss};Counts",100,0,1);
  hist_list.push_back(h_pMiss_ep);

  ////////////////////////////////////////////////
  //epp
  ////////////////////////////////////////////////
  TH2D * h_vtz_e_l = new TH2D("vtz_e_l","vtz_e_l",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_e_l);
  TH2D * h_vtz_e_r = new TH2D("vtz_e_r","vtz_e_r",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_e_r);
  TH2D * h_vtz_l_r = new TH2D("vtz_l_r","vtz_l_r",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_l_r);
  TH2D * h_vtz_l_r_FD_CD = new TH2D("vtz_l_r_FD_CD","vtz_l_r",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_l_r_FD_CD);
  TH2D * h_vtz_l_r_CD_FD = new TH2D("vtz_l_r_CD_FD","vtz_l_r",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_l_r_CD_FD);
  TH2D * h_vtz_l_r_CD_CD = new TH2D("vtz_l_r_CD_CD","vtz_l_r",100,-8,2,100,-8,2);
  hist_list.push_back(h_vtz_l_r_CD_CD);


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

  TH1D * h_thetamissrec_epp = new TH1D("thetamissrec_epp","#theta_{miss,rec};#theta_{miss,rec};Counts",100,0,180);
  hist_list.push_back(h_thetamissrec_epp);

  TH1D * h_pLead_epp = new TH1D("pLead_epp","p_{Lead} SRC;p_{Lead};Counts",100,0,3);
  hist_list.push_back(h_pLead_epp);
  TH1D * h_pMiss_epp = new TH1D("pMiss_epp","p_{Miss} SRC;p_{Miss};Counts",100,0,1);
  hist_list.push_back(h_pMiss_epp);
  TH1D * h_pRec_epp = new TH1D("pRec_epp","p_{Rec} SRC;p_{Rec};Counts",100,0,1);
  hist_list.push_back(h_pRec_epp);
  TH1D * h_thetaRec_epp = new TH1D("thetaRec_epp","#theta_{Rec} SRC;#theta_{Rec};Counts",100,0,120);
  hist_list.push_back(h_thetaRec_epp);
  TH1D * h_phiRec_epp = new TH1D("phiRec_epp","#phi_{Rec} SRC;#phi_{Rec};Counts",100,-180,180);
  hist_list.push_back(h_phiRec_epp);

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Sumw2();
    hist_list[i]->GetXaxis()->CenterTitle();
    hist_list[i]->GetYaxis()->CenterTitle();
  }

  double num = 0;
  double den = 0;
  
  while(chain.Next())
    {

      double wep = 1;
      double wepp = 1;
      if(isMC==1){
	double original_weight = c12->mcevent()->getWeight(); //used if MC events have a weight
	//wep = original_weight * newWeight.get_weight_ep(c12->mcparts());
	//wepp = original_weight * newWeight.get_weight_epp(c12->mcparts());
	wep = original_weight;
	wepp = original_weight;
      }

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
	  double vtz_e = electrons[0]->par()->getVz();

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
	      double poq = lead_ptr.P()/q.P();
	      double thetapq = lead_ptr.Vect().Angle(q.Vect())*180/M_PI;
	      double thetamissq = q.Vect().Angle(neg_miss.Vect())*180/M_PI;
	      double vtz_l = lead[0]->par()->getVz();

	      if(miss.P()<0.3){continue;}

	      //cout<<"FD="<<electrons[0]->trk(DC)->getStatus() <<endl;
	      //cout<<"CD="<<electrons[0]->trk(CVT)->getStatus()<<endl<<endl;

	      if((lead[0]->getRegion()==CD)){		
		h_vtz_e_l->Fill(vtz_e,vtz_l,wepp);}

	      h_pmiss_SRC->Fill(miss.P(),wep);
	      h_Q2_SRC_Q2bin[binQ2(Q2)]->Fill(Q2,wep);
	      h_pmiss_SRC_Q2bin[binQ2(Q2)]->Fill(miss.P(),wep);	      
	      
	      h_pLead_ep->Fill(lead_ptr.P(),wep);
	      h_pMiss_ep->Fill(miss.P(),wep);

	      if(bp!=-1){
		h_Q2_SRC_pmissbin[bp]->Fill(Q2,wep);
		h_emiss_SRC_pmissbin[bp]->Fill(emiss,wep);
	      }

	      if(recoil.size() == 1){
		//if(c12->mcparts()->getPid(2)==2212){		
		num+=1.0;
		SetLorentzVector(recoil_ptr,recoil[0]);
		TVector3 v_miss = neg_miss.Vect();
		TVector3 v_rec  = recoil_ptr.Vect();
		TVector3 v_rel  = (v_miss - v_rec) * 0.5;
		TVector3 v_cm   = v_miss + v_rec;

		TVector3 vt = v_miss.Unit();
		TVector3 vy = v_miss.Cross(q.Vect()).Unit();
		TVector3 vx = vt.Cross(vy).Unit();
		double vtz_r = recoil[0]->par()->getVz();

		h_vtz_e_r->Fill(vtz_e,vtz_r,wepp);
		h_vtz_l_r->Fill(vtz_l,vtz_r,wepp);
		if((lead[0]->getRegion()==FD) && (recoil[0]->getRegion()==CD)){
		  h_vtz_l_r_FD_CD->Fill(vtz_l,vtz_r,wepp);
		}
		if((lead[0]->getRegion()==CD) && (recoil[0]->getRegion()==FD)){
		  h_vtz_l_r_CD_FD->Fill(vtz_l,vtz_r,wepp);
		}
		if((lead[0]->getRegion()==CD) && (recoil[0]->getRegion()==CD)){
		  h_vtz_l_r_CD_CD->Fill(vtz_l,vtz_r,wepp);
		}

		h_thetamissrec_epp->Fill(v_miss.Angle(v_rec)*180/M_PI,wepp);

		
		h_pmiss_Rec->Fill(miss.P(),wepp);

		h_pLead_epp->Fill(lead_ptr.P(),wepp);
		h_pMiss_epp->Fill(miss.P(),wepp);
		h_pRec_epp->Fill(recoil_ptr.P(),wepp);
		h_thetaRec_epp->Fill(recoil_ptr.Theta()*180/M_PI,wepp);
		h_phiRec_epp->Fill(recoil_ptr.Phi()*180/M_PI,wepp);
		
		h_Q2_Rec_Q2bin[binQ2(Q2)]->Fill(Q2,wepp);
		h_p_z_cm_Rec->Fill(v_cm.Dot(vt),wepp);
		h_p_y_cm_Rec->Fill(v_cm.Dot(vy),wepp);
		h_p_x_cm_Rec->Fill(v_cm.Dot(vx),wepp);
		if(bp!=-1){
		  h_Q2_Rec_pmissbin[bp]->Fill(Q2,wepp);
		  h_emiss_Rec_pmissbin[bp]->Fill(emiss,wepp);
		}
		
		h_pmiss_Rec_Q2bin[binQ2(Q2)]->Fill(miss.P(),wepp);
		h_p_z_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vt),wepp);
		h_p_y_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vy),wepp);
		h_p_x_cm_Rec_Q2bin[binQ2(Q2)]->Fill(v_cm.Dot(vx),wepp);
				
	      }

	    }

	}
    }
  
  /*
  TGraph * g_Q2_p_y_cm = new TGraph();
  g_Q2_p_y_cm->SetName("g_Q2_p_y_cm");
  for(int j = 0; j < 10; j++){
    double x = h_Q2_Rec_Q2bin[j]->GetMean();
    TF1 * gFit = new TF1("GausFit",[&](double *x, double *p){ return G(x[0],p[0],p[1],p[2]); },-0.25,0.25,3);
    gFit->SetParameter(0,h_p_x_cm_Rec_Q2bin[j]->GetMaximum());
    gFit->SetParameter(1,0);
    gFit->SetParameter(2,0.1);
    TFitResultPtr gPoint = h_p_x_cm_Rec_Q2bin[j]->Fit(gFit,"SrBeqn","",-0.25,0.25);
    if(gPoint == 0){
      g_Q2_p_y_cm->SetPoint(g_Q2_p_y_cm->GetN(),x,gPoint->Parameter(2));
    }
  }
  */
  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();
  for(int i=0; i<hist_list.size(); i++){
    hist_list[i]->Write();
  } 
  //g_Q2_p_y_cm->Write();

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

