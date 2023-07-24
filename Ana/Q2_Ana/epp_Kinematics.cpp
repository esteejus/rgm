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

/*
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

struct hist_set{
  TH1D h_ep;
  TH1D h_epp;
  TH1D h_ep_Q2bin[10];
  TH1D h_epp_Q2bin[10];
};


void Init_hist_set(hist_set &my_set, string temp_name, string temp_title, double xmin, double xmax){
  
  my_set.h_ep.SetNameTitle((temp_name+"_ep").c_str(),(temp_title+";"+temp_title+";Counts").c_str());
  my_set.h_ep.SetBins(100,xmin,xmax);

  my_set.h_epp.SetNameTitle((temp_name+"_epp").c_str(),(temp_title+";"+temp_title+";Counts").c_str());
  my_set.h_epp.SetBins(100,xmin,xmax);

  for(int i=0; i<10; i++){
    my_set.h_ep_Q2bin[i].SetNameTitle((temp_name+"_ep_Q2bin"+std::to_string(i)).c_str(),(temp_title+";"+temp_title+";Counts").c_str());
    my_set.h_ep_Q2bin[i].SetBins(100,xmin,xmax);

    my_set.h_epp_Q2bin[i].SetNameTitle((temp_name+"_epp_Q2bin"+std::to_string(i)).c_str(),(temp_title+";"+temp_title+";Counts").c_str());
    my_set.h_epp_Q2bin[i].SetBins(100,xmin,xmax);
  }
}

void Fill_hist_set(hist_set *my_set, bool is_epp, double Q2, double x){
  my_set->h_ep.Fill(x);
  my_set->h_ep_Q2bin[binQ2(Q2)].Fill(x);
  if(is_epp){
    my_set->h_epp.Fill(x);
    my_set->h_epp_Q2bin[binQ2(Q2)].Fill(x);
  }  
}

void Write_hist_set(hist_set &my_set, TFile *f){
  f->cd();
  my_set.h_ep.Write();
  my_set.h_epp.Write();
  for(int i=0; i<10; i++){
    my_set.h_ep_Q2bin[i].Write();
    my_set.h_epp_Q2bin[i].Write();
  }
}
*/

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


  TH1D * h_E2miss = new TH1D("E2miss","E_{2 miss} ",100,-0.1,0.5);
  TH2D * h_omega_E2miss = new TH2D("omega_E2miss","E_{2 miss} vs. #omega;#omega;E_{2 miss}",100,0.0,2.5,100,-0.1,0.5);

  
  vector<many_plots> hist_list;

  many_plots h_xB("xB","x_{B}",1.1,2);
  hist_list.push_back(h_xB);
  many_plots h_Q2("Q2","Q^{2}",1.0,6.0);
  hist_list.push_back(h_Q2);
  many_plots h_omega("omega","#omega",0.0,3.0);
  hist_list.push_back(h_omega);
  many_plots h_thetae("thetae","#theta_{e}",0.0,40);
  hist_list.push_back(h_thetae);
  many_plots h_phie("phie","#phi_{e}",-180,180);
  hist_list.push_back(h_phie);

  many_plots h_plead("plead","p_{Lead}",0.0,4.0);
  hist_list.push_back(h_plead);
  many_plots h_thetalead("thetalead","#theta_{Lead}",0.0,90);
  hist_list.push_back(h_thetalead);
  many_plots h_philead("philead","#phi_{Lead}",-180,180);
  hist_list.push_back(h_philead);

  many_plots h_pmiss("pmiss","p_{miss}",0.0,1.5);
  hist_list.push_back(h_pmiss);
  many_plots h_mmiss("mmiss","m_{miss}",0.7,1.2);
  hist_list.push_back(h_mmiss);
  many_plots h_emiss("emiss","E_{miss}",-0.1,0.6);
  hist_list.push_back(h_emiss);

  many_plots h_thetapq("thetapq","#theta_{Lead,q}",0,30);
  hist_list.push_back(h_thetapq);
  many_plots h_thetamissq("thetamissq","#theta_{miss,q}",0,180);
  hist_list.push_back(h_thetamissq);
  many_plots h_poq("poq","p/q",0.55,1.0);
  hist_list.push_back(h_poq);


  /*
  TH1D * h_prec = new TH1D("prec","p_{rec}",100,0.3,1);
  TH2D * h_prec_pmiss = new TH2D("prec_pmiss","p_{miss} vs. p_{rec};p_{rec};p_{miss}",100,0.3,1.0,100,0.3,1.0);
  TH1D * h_thetamissq = new TH1D("thetamissq","#theta_{miss q}",100,100,180);
  TH1D * h_costhetamissrec = new TH1D("costhetamissrec","cos(#theta_{miss rec})",100,-1,1);
  */
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
	      TLorentzVector miss = q + deut_ptr - lead_ptr; 	   
	      TLorentzVector miss_Am1 = q + nucleus_ptr - lead_ptr; 
	      double TB = miss_Am1.E() - miss_Am1.M();
	      double TP1 = lead_ptr.E() - lead_ptr.M();
	      double Emiss = q.E() - TP1 - TB;
	      
	      bool rec = false;
	      if(recoil.size() == 1){rec = true;}
	      h_xB.Fill_hist_set(rec,Q2,xB);
	      h_Q2.Fill_hist_set(rec,Q2,Q2);
	      h_omega.Fill_hist_set(rec,Q2,omega);
	      h_thetae.Fill_hist_set(rec,Q2,el.Theta()*180/M_PI);
	      h_phie.Fill_hist_set(rec,Q2,el.Phi()*180/M_PI);

	      h_plead.Fill_hist_set(rec,Q2,lead_ptr.P());
	      h_thetalead.Fill_hist_set(rec,Q2,lead_ptr.Theta()*180/M_PI);
	      h_philead.Fill_hist_set(rec,Q2,lead_ptr.Phi()*180/M_PI);

	      h_pmiss.Fill_hist_set(rec,Q2,miss.P());
	      h_mmiss.Fill_hist_set(rec,Q2,miss.M());
	      h_emiss.Fill_hist_set(rec,Q2,Emiss);
	      	 
	      h_thetapq.Fill_hist_set(rec,Q2,lead_ptr.Angle(q.Vect())*180/M_PI);
	      h_thetamissq.Fill_hist_set(rec,Q2,miss.Angle(q.Vect())*180/M_PI);
	      h_poq.Fill_hist_set(rec,Q2,lead_ptr.P()/q.P());
	      
	      if(recoil.size() == 1){

		num+=1.0;
		SetLorentzVector(recoil_ptr,recoil[0]);
		double TP2 = recoil_ptr.E() - recoil_ptr.M();
		TLorentzVector miss_Am2 = q + nucleus_ptr - lead_ptr - recoil_ptr; 
		double TB = miss_Am2.E() - miss_Am2.M();
		double E2miss = q.E() - TP1 - TP2 - TB;
	      
		/*
		h_E2miss->Fill(E2miss);
		h_omega_E2miss->Fill(omega,E2miss);

		h_prec->Fill(recoil_ptr.Rho());
		h_prec_pmiss->Fill(recoil_ptr.Rho(),miss.Rho());
		h_thetamissq->Fill(180-miss.Angle(q.Vect())*180/M_PI);
		h_costhetamissrec->Fill(cos(miss.Angle(recoil_ptr.Vect())));*/
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

  for(int i=0; i<hist_list.size(); i++){
    hist_list[i].Write_hist_set(f,fileName,myCanvas);
  }

  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");


  f->Close();


  return 0;
}

