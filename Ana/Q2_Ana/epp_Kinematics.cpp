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


  //TH1D * h_E2miss = new TH1D("E2miss","E_{2 miss} ",100,-0.1,0.5);
  //TH2D * h_omega_E2miss = new TH2D("omega_E2miss","E_{2 miss} vs. #omega;#omega;E_{2 miss}",100,0.0,2.5,100,-0.1,0.5);

  
  vector<many_plots> hist_list_ep;

  many_plots h_xB("xB","x_{B}",1.1,2);
  hist_list_ep.push_back(h_xB);
  many_plots h_Q2("Q2","Q^{2}",1.0,6.0);
  hist_list_ep.push_back(h_Q2);
  many_plots h_omega("omega","#omega",0.0,3.0);
  hist_list_ep.push_back(h_omega);
  many_plots h_thetae("thetae","#theta_{e}",0.0,40);
  hist_list_ep.push_back(h_thetae);
  many_plots h_phie("phie","#phi_{e}",-180,180);
  hist_list_ep.push_back(h_phie);

  many_plots h_plead("plead","p_{Lead}",0.6,4.0);
  hist_list_ep.push_back(h_plead);
  many_plots h_thetalead("thetalead","#theta_{Lead}",0.0,90);
  hist_list_ep.push_back(h_thetalead);
  many_plots h_thetalead_FD("thetalead_FD","FD #theta_{Lead}",0.0,90);
  hist_list_ep.push_back(h_thetalead_FD);
  many_plots h_thetalead_CD("thetalead_CD","CD #theta_{Lead}",0.0,90);
  hist_list_ep.push_back(h_thetalead_CD);
  many_plots h_philead("philead","#phi_{Lead}",-180,180);
  hist_list_ep.push_back(h_philead);
  many_plots h_pmiss("pmiss","p_{miss}",0.0,1.0);
  hist_list_ep.push_back(h_pmiss);
  many_plots h_mmiss("mmiss","m_{miss}",0.7,1.2);
  hist_list_ep.push_back(h_mmiss);
  many_plots h_emiss("emiss","E_{miss}",-0.1,0.6);
  hist_list_ep.push_back(h_emiss);
  many_plots h_thetapq("thetapq","#theta_{Lead,q}",0,30);
  hist_list_ep.push_back(h_thetapq);
  many_plots h_thetamissq("thetamissq","#theta_{miss,q}",120,180);
  hist_list_ep.push_back(h_thetamissq);
  many_plots h_poq("poq","p/q",0.55,1.0);
  hist_list_ep.push_back(h_poq);

  vector<many_plots> hist_list_epp;

  many_plots h_precoil("precoil","p_{recoil}",0.15,1.0);
  hist_list_epp.push_back(h_precoil);
  many_plots h_prel("prel","p_{rel}",0.15,1.0);
  hist_list_epp.push_back(h_prel);
  many_plots h_thetamissrecoil("thetamissrecoil","#theta_{miss,recoil}",0,180);
  hist_list_epp.push_back(h_thetamissrecoil);
  many_plots h_thetacmrel("thetacmrel","#theta_{cm,rel}",0,180);
  hist_list_epp.push_back(h_thetacmrel);
  many_plots h_pcm("pcm","p_{cm}",0.0,1.0);
  hist_list_epp.push_back(h_pcm);
  many_plots h_pcmx("pcmx","p_{X,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmx);
  many_plots h_pcmy("pcmy","p_{||,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmy);
  many_plots h_pcmz("pcmz","p_{miss,cm}",-0.75,0.75);
  hist_list_epp.push_back(h_pcmz);
  many_plots h_E2miss("E2miss","E_{2,miss}",-0.1,1.0);
  hist_list_epp.push_back(h_E2miss);

  clasAna.setEcalSFCuts();
  clasAna.setEcalPCuts();

  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();

  clasAna.setVertexCuts();
  clasAna.setVertexCorrCuts();
  //clasAna.setDCEdgeCuts();

  double num = 0;
  double den = 0;
  
  while(chain.Next())
    {

      double weight = c12->mcevent()->getWeight(); //used if MC events have a weight 

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
	      bool FD = (lead[0]->sci(clas12::FTOF1A)->getDetector() == 12) || (lead[0]->sci(clas12::FTOF1B)->getDetector() == 12) || (lead[0]->sci(clas12::FTOF2)->getDetector() == 12);
	      bool CD = (lead[0]->sci(clas12::CTOF)->getDetector() == 4);
	      TLorentzVector miss = q + deut_ptr - lead_ptr; 	   
	      TLorentzVector neg_miss = -miss;	      
	      TLorentzVector miss_Am1 = q + nucleus_ptr - lead_ptr; 
	      double TB = miss_Am1.E() - miss_Am1.M();
	      double TP1 = lead_ptr.E() - lead_ptr.M();
	      double Emiss = q.E() - TP1 - TB;
	      
	      bool rec = false;
	      if(recoil.size() == 1){rec = true;}
	      if(miss.P()<0.3){continue;}
	      if(CD && lead_ptr.Theta()*180/M_PI<45){continue;}
	      h_xB.Fill_hist_set(rec,Q2,xB);
	      h_Q2.Fill_hist_set(rec,Q2,Q2);
	      h_omega.Fill_hist_set(rec,Q2,omega);
	      h_thetae.Fill_hist_set(rec,Q2,el.Theta()*180/M_PI);
	      h_phie.Fill_hist_set(rec,Q2,el.Phi()*180/M_PI);

	      h_plead.Fill_hist_set(rec,Q2,lead_ptr.P());
	      h_thetalead.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI);
	      if(FD){
	      h_thetalead_FD.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI);
	      }
	      else if(CD){
		h_thetalead_CD.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI);	      
	      }

	      h_philead.Fill_hist_set(rec,Q2,lead_ptr.Phi()*180/M_PI);

	      h_pmiss.Fill_hist_set(rec,Q2,miss.P());
	      h_mmiss.Fill_hist_set(rec,Q2,miss.M());
	      h_emiss.Fill_hist_set(rec,Q2,Emiss);
	      	 
	      h_thetapq.Fill_hist_set(rec,Q2,lead_ptr.Angle(q.Vect())*180/M_PI);
	      h_thetamissq.Fill_hist_set(rec,Q2,neg_miss.Angle(q.Vect())*180/M_PI);
	      h_poq.Fill_hist_set(rec,Q2,lead_ptr.P()/q.P());
	      
	      if(recoil.size() == 1){

		num+=1.0;
		SetLorentzVector(recoil_ptr,recoil[0]);
		double TP2 = recoil_ptr.E() - recoil_ptr.M();
		TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
		double TB = miss_Am2.E() - miss_Am2.M();
		double E2miss = q.E() - TP1 - TP2 - TB;
		
		TVector3 v_miss = neg_miss.Vect();
		TVector3 v_rec  = recoil_ptr.Vect();
		TVector3 v_rel  = (v_miss - v_rec) * 0.5;
		TVector3 v_cm   = v_miss + v_rec;

		TVector3 vz = v_miss.Unit();
		TVector3 vy = v_miss.Cross(q.Vect()).Unit();
		TVector3 vx = vz.Cross(vy).Unit();

		h_precoil.Fill_hist_set(rec,Q2,recoil_ptr.P());
		h_prel.Fill_hist_set(rec,Q2,v_rel.Mag());
		h_thetamissrecoil.Fill_hist_set(rec,Q2,v_miss.Angle(v_rec)*180/M_PI);
		h_thetacmrel.Fill_hist_set(rec,Q2,v_cm.Angle(v_rel)*180/M_PI);
		h_pcm.Fill_hist_set(rec,Q2,v_cm.Mag());
		h_pcmx.Fill_hist_set(rec,Q2,v_cm.Dot(vx));
		h_pcmy.Fill_hist_set(rec,Q2,v_cm.Dot(vy));
		h_pcmz.Fill_hist_set(rec,Q2,v_cm.Dot(vz));
		h_E2miss.Fill_hist_set(rec,Q2,E2miss);

	      }
	    }

	}
    }

  /////////////////////////////////////////////////////
  //Now create the output PDFs
  /////////////////////////////////////////////////////
  TFile *f = new TFile(outFile,"RECREATE");
  f->cd();

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

  h_pmiss.Write_ratio_set(f,fileName,myCanvas);

  for(int i=0; i<hist_list_ep.size(); i++){
    hist_list_ep[i].Write_hist_set(f,fileName,myCanvas);
  } 

  for(int i=0; i<hist_list_epp.size(); i++){
    hist_list_epp[i].Write_hist_set_epp(f,fileName,myCanvas);
  } 

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

