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
#include "eventcut_old.h" // compare with new in Monitoring!!
#include "functions_old.h" // compare with new in Monitoring!!

using namespace std;
using namespace clas12;

//const double c = 29.9792458;
const double mp = 0.938272;

void printProgress(double percentage);

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/cutfile.txt> <path/to/input.hipo> \n";
}


int main(int argc, char ** argv)
{

  if(argc < 7)
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
  char * pdfFile = argv[4];
  eventcut myCut(Ebeam,argv[5]);
  myCut.print_cuts();
  clas12root::HipoChain chain;
  for(int k = 6; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  chain.db()->turnOffQADB();                 

        
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

  
  /////////////////////////////////////
  //Histos: electrons
  /////////////////////////////////////
  TH1D * h_nphe = new TH1D("nphe","# Photo-electrons in HTCC;# Photo-electrons;Counts",40,0,40);
  hist_list_1.push_back(h_nphe);
  TH1D * h_vtz_e = new TH1D("vtz_e","Electron z Vertex;vertex;Counts",100,-10,10);
  hist_list_1.push_back(h_vtz_e);


  /////////////////////////////////////
  //Histos: lead proton
  /////////////////////////////////////
  TH1D * h_vtzdiff_ep = new TH1D("vtzdiff_ep","Vertex difference z between e and p",100,-6,6);
  hist_list_1.push_back(h_vtzdiff_ep);
  TH1D * h_chi2pid = new TH1D("chi2pid","Proton #chi^{2}_{PID}",100,-5,5);
  hist_list_1.push_back(h_chi2pid);
  TH2D * h_dbetap = new TH2D("dbeta_p","#Delta #beta vs Momentum;Momentum (GeV/c);#beta_{meas} - p/sqrt(p^{2}+m^{2})",100,0,3,100,-0.2,0.2);
  hist_list_2.push_back(h_dbetap);
  TH2D * h_betap = new TH2D("beta_p","#beta vs Momentum;Momentum (GeV/c);#beta",100,0,3,100,0,1.2);
  hist_list_2.push_back(h_betap);

  TH1D * h_pmiss_pL = new TH1D("pmiss_pL","Missing Momentum (e,e'p_{Lead});p_{miss} (GeV/c);Counts",30,0,1);
  hist_list_1.push_back(h_pmiss_pL);
  TH2D * h_pangles = new TH2D("pangles","Proton angles;#phi;#theta",90,-180,180,60,0,60);
  hist_list_2.push_back(h_pangles);
  TH2D * h_pq = new TH2D("pq","#theta_{p,q} vs p/q;p/q;#theta_{pq}",100,0,2,100,0,100);
  hist_list_2.push_back(h_pq);
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  hist_list_1.push_back(h_mmiss);
  TH1D * h_pmisstheta = new TH1D("pmisstheta","Theta Component of p_{miss};#theta_{pmiss};Counts",180,0,180);
  hist_list_1.push_back(h_pmisstheta);

  TH2D * h_q2_xb = new TH2D("q2_xb","Q^{2} vs x_{B};x_{B};Q^{2} (GeV^{2})",100,0,2,100,0,3);
  hist_list_2.push_back(h_q2_xb);


  /////////////////////////////////////
  //Histos: recoil proton
  /////////////////////////////////////
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
  hist_list_1.push_back(h_psize);
  TH1D * h_pcos0 = new TH1D("pcos0","Cosine of #theta_{pmiss,pp}",100,-1.1,1.1);
  hist_list_1.push_back(h_pcos0);
  TH1D * h_prec_p = new TH1D("prec_p","Recoil Proton Momentum",100,0,2);
  hist_list_1.push_back(h_prec_p);
  TH2D * h_prec_angles = new TH2D("prec_angles","Recoil Proton Angular Distribution;phi (deg);theta (deg)",360,-180,180,180,0,180);
  hist_list_2.push_back(h_prec_angles);

  TH2D * h_pp_pmiss = new TH2D("pp_pmiss","Proton Momentum vs Missing Momentum;p_{miss} (GeV/c);Proton Momentum (GeV/c)",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_pp_pmiss);
  TH2D * h_pptheta = new TH2D("pptheta","Proton #theta vs Momentum;#theta_{p};Momentum (GeV/c)",180,0,180,50,0,1.5);
  hist_list_2.push_back(h_pptheta);


  /////////////////////////////////////
  //Histos: recoil neutron
  /////////////////////////////////////
  TH1D * h_nsize = new TH1D("nsize","Number of Neutrons in Event",10,0,10);
  hist_list_1.push_back(h_nsize);
  TH1D * h_ncos0 = new TH1D("ncos0","Cosine of #theta_{pmiss,pn}",100,-1.1,1.1);
  hist_list_1.push_back(h_ncos0);
  TH1D * h_nrec_p = new TH1D("nrec_p","Recoil Neutron Momentum",100,0,2);
  hist_list_1.push_back(h_nrec_p);
  TH2D * h_nrec_angles = new TH2D("nrec_angles","Recoil Neutron Angular Distribution;phi (deg);theta (deg)",48,-180,180,180,0,180);
  hist_list_2.push_back(h_nrec_angles);

  TH2D * h_pn_pmiss = new TH2D("pn_pmiss","Neutron Momentum vs Missing Momentum;p_{miss} (GeV/c);Neutron Momentum (GeV/c)",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_pn_pmiss);
  TH2D * h_nptheta = new TH2D("nptheta","Neutron #theta vs Momentum;#theta_{n};Momentum (GeV/c)",180,0,180,50,0,1.5);
  hist_list_2.push_back(h_nptheta);



  /////////////////////////////////////
  //Histos: Yields
  /////////////////////////////////////
  TH1D * h_e_count = new TH1D("ecount","Electron Yield",20,0.2,1.);
  TH1D * h_pl_count = new TH1D("plcount","Lead Proton Yield",20,0.2,1.);
  TH1D * h_psrc_count = new TH1D("psrccount","SRC Proton Yield",20,0.2,1.);
  TH1D * h_pwrec_count = new TH1D("pwreccount","SRC Proton with Recoil Yield",20,0.2,1.);
  TH1D * h_pn_count = new TH1D("pncount","pn Yield",20,0.2,1.);
  TH1D * h_pp_count = new TH1D("ppcount","pp Yield",20,0.2,1.);





  /////////////////////////////////////
  //Histos: pmiss
  /////////////////////////////////////
  TH1D * h_pmiss_p = new TH1D("pmiss_p","Missing Momentum (e,e'p_{SRC});p_{miss} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pmiss_p);
  TH1D * h_pmiss_p_wrec = new TH1D("pmiss_p_wrec","Missing Momentum (e,e'p_{SRC}N_{rec});p_{miss} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pmiss_p);
  TH1D * h_pmiss_pp = new TH1D("pmiss_pp","Missing Momentum (e,e'p_{SRC}p_{rec});p_{miss} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pmiss_pp); 
  TH1D * h_pmiss_pn = new TH1D("pmiss_pn","Missing Momentum (e,e'p_{SRC}n_{rec});p_{miss} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pmiss_pn);
  TH1D * h_pmiss_pn_corr = new TH1D("pmiss_pn_corr","Missing Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{miss} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pmiss_pn_corr);

  TH1D * h_pn = new TH1D("pn","Momentum (e,e'p_{SRC}n_{rec});p_{n} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pn);
  TH1D * h_pn_corr = new TH1D("pn_corr","Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{n} (GeV/c);Counts",20,0.2,1);
  hist_list_1.push_back(h_pn_corr);



  /////////////////////////////////////
  //Histos: neutron
  /////////////////////////////////////
  



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
  while(chain.Next()==true){
    //int misc = chain.GetNRecords();
    //Display completed  
    counter++;
    if((counter%1000000) == 0){
      cerr << "\n" <<counter/1000000 <<" million completed";
    }    

    if((counter%100000) == 0){
      cerr << ".";
    }    

    // get particles by type
    auto elec=c12->getByID(11);
    auto prot=c12->getByID(2212);
    auto neut=c12->getByID(2112);
    auto allParticles=c12->getDetParticles();
    double weight = 1;
    if(isMC){weight=c12->mcevent()->getWeight();}


    // initial cuts
    if (prot.size()<1) {continue;}
    if (elec.size()<1) {continue;}
    //if (elec.size()!=1) {continue;}

    // ELECTRONS
    if(!myCut.electroncut(c12)){continue;}

    // electron kinematics
    TVector3 pe;
    TVector3 pb(0,0,Ebeam);
    pe.SetMagThetaPhi(elec[0]->getP(),elec[0]->getTheta(),elec[0]->getPhi());
    TVector3 q = pb - pe;
    double vze = elec[0]->par()->getVz();
    double nu = Ebeam - pe.Mag();
    double Q2 = q.Mag2() - (nu*nu);
    double xB = Q2 / (2*mp*nu);
    double W2 = (mp*mp) - Q2 + (2*nu*mp);
    double etheta = elec[0]->getTheta()*180./M_PI;
    int nphe = elec[0]->che(HTCC)->getNphe();
    double vtz_e = elec[0]->par()->getVz();
    //double EoP_e = (elec[0]->cal(PCAL)->getEnergy() + elec[0]->cal(ECIN)->getEnergy() + elec[0]->cal(ECOUT)->getEnergy()) / pe.Mag();
   
    // electron histograms: quality cuts
    h_nphe->Fill(nphe,weight);
    h_vtz_e->Fill(vtz_e,weight);



    //// LEAD PROTON ////
    // LEAD PROTON - FTOF

    for (int k=0; k<prot.size(); k++)
    {
      // proton kinematics
      TVector3 pL_temp;
      pL_temp.SetMagThetaPhi(prot[k]->getP(),prot[k]->getTheta(),prot[k]->getPhi());
      double vzlead = prot[k]->par()->getVz();
      double chi2pid = prot[k]->par()->getChi2Pid();
      double ltheta = prot[k]->getTheta()*180./M_PI;
      double lphi = prot[k]->getPhi()*180./M_PI;
      double theta_pq = pL_temp.Angle(q) * 180./M_PI;
      double mmiss = get_mmiss(pb,pe,pL_temp);
      double dbetap = prot[k]->par()->getBeta() - pL_temp.Mag()/sqrt(pL_temp.Mag2()+mp*mp);
      // lead histos
      h_vtzdiff_ep->Fill(vze-vzlead,weight);
      h_chi2pid->Fill(chi2pid,weight);
      h_dbetap->Fill(pL_temp.Mag(),dbetap,weight);
      h_betap->Fill(pL_temp.Mag(),prot[k]->par()->getBeta(),weight);
    }


    // I need code here to figure out which proton is lead
    // cutfile_ts.txt - proton in FTOF, 0<ltheta<45, 0<theta_pq<25, -3<chipid<3, -3<vtz_diff<3
    int lead = myCut.leadnucleoncut(c12);
    if (lead<0) {continue;}

    TVector3 pL;
    pL.SetMagThetaPhi(prot[lead]->getP(),prot[lead]->getTheta(),prot[lead]->getPhi());
    double vzlead = prot[lead]->par()->getVz();
    double chi2pid = prot[lead]->par()->getChi2Pid();
    double ltheta = prot[lead]->getTheta()*180./M_PI;
    double lphi = prot[lead]->getPhi()*180./M_PI;
    double theta_pq = pL.Angle(q) * 180./M_PI;
    double mmiss = get_mmiss(pb,pe,pL);
    double dbetap = prot[lead]->par()->getBeta() - pL.Mag()/sqrt(pL.Mag2()+mp*mp);



    if ((vze-vzlead)<-3. || (vze-vzlead)>3.) {continue;}
    if (dbetap<-0.05 || dbetap>0.05) {continue;}
    if (chi2pid<-3.0 || chi2pid>3.0) {continue;}


    //// PMISS ////
    TVector3 pmiss = pL - q;
    double pm_theta = pmiss.Theta()*180./M_PI;


// in yield counts, for each channel, we have to account for the extra runs
// with higher prel cutoffs (epja-sim)
//TVector3 prel = pmiss-p
/*if (pmiss.Mag()>0.2 && pmiss.Mag()<0.3)
  {weight_c = weight/1;}
else if (pmiss.Mag()>0.3 && pmiss.Mag()<0.4)
  {weight_c = weight/2;}
else if (pmiss.Mag()>0.4 && pmiss.Mag()<0.5)
  {weight_c = weight/3;}
else if (pmiss.Mag()>0.5 && pmiss.Mag()<0.6)
  {weight_c = weight/4;}
else if (pmiss.Mag()>0.6 && pmiss.Mag()<0.7)
  {weight_c = weight/5;}
else if (pmiss.Mag()>0.7 && pmiss.Mag()<0.8)
  {weight_c = weight/6;}
else if (pmiss.Mag()>0.8)
  {weight_c = weight/7;}*/
//else if (pmiss.Mag()>0.9 && pmiss.Mag()<1.0)
  //{weight_c = weight/8;}


h_pl_count->Fill(pmiss.Mag(),weight);

    h_pmisstheta->Fill(pm_theta,weight);
    h_pmiss_pL->Fill(pmiss.Mag(),weight);

    //if (pm_theta<40.0 || pm_theta>140.0) {continue;} // 2/1/23 - this is only for deuterium!!!!
    //if (pmiss.Mag()<0.2) {continue;}
    h_pangles->Fill(lphi,ltheta,weight);



    h_pq->Fill(pL.Mag()/q.Mag(),theta_pq,weight);

    h_q2_xb->Fill(xB,Q2,weight);



    // lead SRC cuts

    if (pmiss.Mag()<0.3) {continue;}
    if (xB<1.1) {continue;}
    if (Q2<1.5) {continue;}
    if (pL.Mag()/q.Mag()<0.62 || pL.Mag()/q.Mag()>1.1) {continue;}
    if (theta_pq>25) {continue;}

    // fill e'p(SRC) histogram
    h_mmiss->Fill(mmiss,weight);

    if (mmiss>1.1) {continue;}
h_psrc_count->Fill(pmiss.Mag(),weight);

    h_pmiss_p->Fill(pmiss.Mag(),weight);


//// RECOIL P ////
    // look for recoil proton in CTOF
    int rec_p = -1;
    h_psize->Fill(prot.size(),weight);
    TVector3 p_recp;
    TVector3 p_vecX;
    double p_cos0;

    for (int i=0; i<prot.size(); i++)
    {
      // don't accept lead
      if (i==lead) {continue;}
      p_recp.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
      // recoil p must be in CTOF
      bool is_CTOF = prot[i]->sci(CTOF)->getDetector()==4;
      if (!is_CTOF) {continue;}
      // get momenta/angles of recoil candidates
      h_prec_p->Fill(p_recp.Mag(),weight);
      h_prec_angles->Fill(p_recp.Phi()*180./M_PI,p_recp.Theta()*180./M_PI,weight);
      h_pptheta->Fill(p_recp.Theta()*180./M_PI,pmiss.Mag(),weight);
      if (p_recp.Mag()<0.35) {continue;}
      if (p_recp.Theta()*180./M_PI<40 || p_recp.Theta()*180./M_PI>140) {continue;}
      // close in angle to pmiss
      p_vecX.SetXYZ( prot[i]->par()->getPx(), prot[i]->par()->getPy(), prot[i]->par()->getPz() );
      p_cos0 = pmiss.Dot(p_vecX) / (pmiss.Mag() * p_vecX.Mag());
      h_pcos0->Fill(p_cos0,weight);
      //if (p_cos0>-0.8) {continue;} // not needed?
      // fast recoil
      h_pp_pmiss->Fill(pmiss.Mag(),p_recp.Mag(),weight);
      rec_p = i;
    }

    if (rec_p>-1)
    {
      // fill histo
      h_pmiss_p_wrec->Fill(pmiss.Mag(),weight);
      h_pmiss_pp->Fill(pmiss.Mag(),weight);
      h_pwrec_count->Fill(pmiss.Mag(),weight);
      h_pp_count->Fill(pmiss.Mag(),weight);
    }




//// RECOIL N ////
    // look for recoil neutron in CTOF
    //if (rec_p>-1) {continue;} // note: only look for recoil n if no recoil p was found
    int rec_n = -1;
    h_nsize->Fill(neut.size(),weight);
    if (neut.size()<1) {continue;}
    TVector3 p_recn;
    TVector3 n_vecX;
    double n_cos0;

    for (int i=0; i<neut.size(); i++)
    {
      p_recn.SetMagThetaPhi(neut[i]->getP(),neut[i]->getTheta(),neut[i]->getPhi());
      // recoil n must be in CND
      bool is_CND1 = neut[i]->sci(CND1)->getDetector()==3;
      bool is_CND2 = neut[i]->sci(CND2)->getDetector()==3;
      bool is_CND3 = neut[i]->sci(CND3)->getDetector()==3;
      if (!is_CND1 && !is_CND2 && !is_CND3) {continue;}
      // get momenta/angles of recoil candidates
      h_nrec_p->Fill(p_recn.Mag(),weight);
      h_nrec_angles->Fill(p_recn.Phi()*180./M_PI,p_recn.Theta()*180./M_PI,weight);
      h_nptheta->Fill(pm_theta,pmiss.Mag(),weight);
      if (p_recn.Mag()<0.35) {continue;}
      if (p_recn.Theta()*180./M_PI<40 || p_recn.Theta()*180./M_PI>140) {continue;}
      // close in angle to pmiss
      n_vecX.SetXYZ( neut[i]->par()->getPx(), neut[i]->par()->getPy(), neut[i]->par()->getPz() );
      n_cos0 = pmiss.Dot(n_vecX) / (pmiss.Mag() * n_vecX.Mag());
      h_ncos0->Fill(n_cos0,weight);
      //if (n_cos0>-0.8) {continue;} // not needed?
      // compared to pmiss
      h_pn_pmiss->Fill(pmiss.Mag(),p_recn.Mag(),weight);
      rec_n = i;
    }

    if (rec_n>-1)
    {
      // fill histo
      TVector3 pn;
      pn.SetXYZ( neut[rec_n]->par()->getPx(), neut[rec_n]->par()->getPy(), neut[rec_n]->par()->getPz() );
      double neff = 0.1;
      //double neff = 0.16 - 0.1333*pmiss.Mag();
      h_pmiss_p_wrec->Fill(pmiss.Mag(),weight);
      h_pmiss_pn->Fill(pmiss.Mag(),weight);
      h_pmiss_pn_corr->Fill(pmiss.Mag(),weight/neff);
      h_pn->Fill(pn.Mag(),weight);
      h_pn_corr->Fill(pn.Mag(),weight/neff);
      h_pwrec_count->Fill(pmiss.Mag(),weight);
      h_pn_count->Fill(pmiss.Mag(),weight);
    }


  }



  cout<<counter<<endl;

  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
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



  /////////////////////////////////////////////////////
  //Histos: all neutrals
  /////////////////////////////////////////////////////


  myText->cd();
  text.DrawLatex(0.2,0.9,"C-12");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_nphe->Draw();
  myCanvas->cd(2);
  h_vtz_e->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // electrons
  myText->cd();
  text.DrawLatex(0.2,0.9,"(e,e') Cuts:");
  double line = 0.8;
  if(myCut.getDoCut(e_cuts)){
    myCut.print_cut_onPDF(text,e_nphe,line);
    myCut.print_cut_onPDF(text,e_calv,line);
    myCut.print_cut_onPDF(text,e_calw,line);
    myCut.print_cut_onPDF(text,e_SF,line);
    myCut.print_cut_onPDF(text,e_mom,line);
    myCut.print_cut_onPDF(text,e_vtze,line);
  }
  myText->Print(fileName,"pdf");
  myText->Clear();

  // electrons
  // lead proton in FTOF

  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Lead proton in FTOF");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_vtzdiff_ep->Draw();
  myCanvas->cd(2);  h_chi2pid->Draw();
  myCanvas->cd(3);  h_dbetap->Draw("colz");
  myCanvas->cd(4);  h_betap->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myText->cd();
  text.DrawLatex(0.2,0.9,"Lead proton cuts:");
  text.DrawLatex(0.2,0.8,"-3 cm < (vtz_{e}-vtz_{pL}) < 3 cm");
  text.DrawLatex(0.2,0.7,"-0.05 < #Delta #beta < 0.05");
  text.DrawLatex(0.2,0.6,"-3 < #chi^{2}_{PID} < 3");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pL->Draw("colz");
  myCanvas->cd(2);  h_pmisstheta->Draw();
  myCanvas->cd(3);  h_pangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myText->cd();
  text.DrawLatex(0.2,0.9,"Pmiss cuts:");
  text.DrawLatex(0.2,0.8,"40 < #theta_{pmiss} < 140");
  text.DrawLatex(0.2,0.7,"p_{miss} > 0.2 GeV/c");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pq->Draw("colz");
  myCanvas->cd(2);  h_q2_xb->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myText->cd();
  text.DrawLatex(0.2,0.9,"SRC cuts:");
  text.DrawLatex(0.2,0.8,"0.8 GeV/c^{2} < M_{miss} < 1.2 GeV/c^{2}");
  text.DrawLatex(0.2,0.7,"0.62 < p/q < 0.96");
  text.DrawLatex(0.2,0.6,"Q^{2} > 1.5 GeV^{2}");
  text.DrawLatex(0.2,0.5,"x_{B} > 1.1");
  myText->Print(fileName,"pdf");
  myText->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Recoil protons");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_p->Draw();
  myCanvas->cd(2);  h_mmiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_psize->Draw();
  myCanvas->cd(2);  h_pcos0->Draw();
  myCanvas->cd(3);  h_prec_p->Draw();
  myCanvas->cd(4);  h_prec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pp_pmiss->Draw("colz");
  myCanvas->cd(2);  h_pptheta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Recoil neutrons");
  myText->Print(fileName,"pdf");
  myText->Clear();



  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_nsize->Draw();
  myCanvas->cd(2);  h_ncos0->Draw();
  myCanvas->cd(3);  h_nrec_p->Draw();
  myCanvas->cd(4);  h_nrec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pn_pmiss->Draw("colz");
  myCanvas->cd(2);  h_nptheta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myText->cd();
  text.DrawLatex(0.2,0.9,"Recoil cuts:");
  text.DrawLatex(0.2,0.8,"proton in CTOF, neutron in CND");
  text.DrawLatex(0.2,0.7,"cos #theta_{prec,pmiss} < -0.8");
  text.DrawLatex(0.2,0.6,"p_{rec} > 0.2 GeV/c");
  myText->Print(fileName,"pdf");
  myText->Clear();



  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pmiss_p_wrec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pmiss_pn->Draw();
  myCanvas->cd(3);
  TH1D * h_pmiss_pp_p = (TH1D*)h_pmiss_pp->Clone();
  h_pmiss_pp_p->Divide(h_pmiss_p_wrec);
  h_pmiss_pp_p->Draw();
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p");
  myCanvas->cd(4);
  TH1D * h_pmiss_pn_p = (TH1D*)h_pmiss_pn->Clone();
  h_pmiss_pn_p->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_p->Draw();
  h_pmiss_pn_p->GetYaxis()->SetTitle("pn/p");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pmiss_pn_corr->Draw();
  myCanvas->cd(3);  h_pmiss_pp_p->Draw();
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p");
  myCanvas->cd(4);
  TH1D * h_pmiss_pn_p_corr = (TH1D*)h_pmiss_pn_corr->Clone();
  h_pmiss_pn_p_corr->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_p_corr->Draw();
  h_pmiss_pn_p_corr->GetYaxis()->SetTitle("pn/p");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pn->Draw();
  myCanvas->cd(3);
  TH1D * h_pp_p = (TH1D*)h_pmiss_pp->Clone();
  h_pp_p->Divide(h_pmiss_p_wrec);
  h_pp_p->Draw();
  h_pp_p->GetYaxis()->SetTitle("pp/p");
  myCanvas->cd(4);
  TH1D * h_pn_p = (TH1D*)h_pn->Clone();
  h_pn_p->Divide(h_pmiss_p_wrec);
  h_pn_p->Draw();
  h_pn_p->GetYaxis()->SetTitle("pn/p");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pn_corr->Draw();
  myCanvas->cd(3);  h_pmiss_pp_p->Draw();
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p");
  myCanvas->cd(4);
  TH1D * h_pn_p_corr = (TH1D*)h_pn_corr->Clone();
  h_pn_p_corr->Divide(h_pmiss_p_wrec);
  h_pn_p_corr->Draw();
  h_pn_p_corr->GetYaxis()->SetTitle("pn/p");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_pp_pn = (TH1D*)h_pmiss_pp->Clone();
  h_pp_pn->Divide(h_pmiss_pn_corr);
  h_pp_pn->Scale(0.5);
  h_pp_pn->Draw();
  h_pp_pn->GetYaxis()->SetTitle("pp/2pn");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // yields
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_e_count->Draw();
  h_pl_count->SetLineColor(kOrange);
  h_pl_count->Draw("same");
  h_psrc_count->SetLineColor(kGreen);
  h_psrc_count->Draw("same");
  h_pn_count->SetLineColor(kRed);
  h_pn_count->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // write to root file
  h_e_count->Write();
  h_pl_count->Write();
  h_psrc_count->Write();
  h_pwrec_count->Write();
  h_pn_count->Write();
  h_pp_count->Write();
  h_pmiss_p_wrec->Write();



  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
}



void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

