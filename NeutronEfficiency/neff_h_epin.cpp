#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <iomanip>
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
#include <TLegend.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "eventcut.h"
#include "functions.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
const double mp = 0.938272;
const double m_piplus = 0.13957039;

void printProgress(double percentage);
double get_pin_mmiss(TVector3 p_b, TVector3 p_e, TVector3 ppi);
Double_t background_poly(Double_t *x, Double_t *par);
Double_t signal(Double_t *x, Double_t *par); 
Double_t mmiss_signal_gauss(Double_t *x, Double_t *par);
Double_t mmiss_signal_poly(Double_t *x, Double_t *par);
double * hist_projections(TCanvas * can, TH2D * hist2d, int num_hist, char v);

// theta ranges
double t1 = 40;
double t2 = 50;
double t3 = 60;
double t4 = 70;
double t5 = 80;

// misc
int neff_pbins = 16;
int neff_tbins = 16;
int pgrid_x = ceil(sqrt(neff_pbins));
int pgrid_y = ceil((double)neff_pbins/(double)pgrid_x);
int tgrid_x = ceil(sqrt(neff_tbins));
int tgrid_y = ceil((double)neff_tbins/(double)tgrid_x);
double p_lo = 0.2;
double p_hi = 1.1;
double t_lo = 40;
double t_hi = 65;

void Usage()
{
  std::cerr << "Usage: ./code <MC=1,Data=0> <Ebeam(GeV)> <path/to/output.root> <path/to/output.pdf> <path/to/cutfile.txt> <path/to/input.hipo> \n";
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

  /*std::string rootPath(argv[3]);
  rootPath += ".root";
  const char * fullrootpath = rootPath.c_str();*/
  /*std::string neff_t_name(argv[3]);
  neff_t_name[strlen(neff_t_name)-6] = '\0';
  neff_t_filename = neff_t_name.c_str();*/


  
  // output file names
  TFile * outFile = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];
  eventcut myCut(Ebeam,argv[5]);

  char * basename = argv[3];
  basename[strlen(basename)-5] = '\0';
  string theta_name(string(basename) + "_theta.txt");
  theta_name.c_str();
  string p_name(string(basename) + "_p.txt");
  p_name.c_str();


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
  //Histos: basic event selection
  /////////////////////////////////////
  TH1D * h_pid = new TH1D("particle pid","PID (All Particles)",1000,-3000,3000);
  hist_list_1.push_back(h_pid);

  /////////////////////////////////////
  //Histos: pions
  /////////////////////////////////////
  TH1D * h_pivertex = new TH1D("pi vertex","Pion vertex - electron vertex;v_{pi+} - v_{e} (cm);Counts",50,-5,5);
  hist_list_1.push_back(h_pivertex);
  TH2D * h_dbetap = new TH2D("dbeta_p","#Delta #beta vs Momentum;Momentum (GeV/c);#beta_{meas} - p/sqrt(p^{2}+m^{2})",100,0,3,100,-0.2,0.2);
  hist_list_2.push_back(h_dbetap);
  TH1D * h_pitheta = new TH1D("pitheta","Pion #theta;#theta (degrees);Counts",180,0,180);
  hist_list_1.push_back(h_pitheta);


  /////////////////////////////////////
  //Histos: pmiss
  /////////////////////////////////////
  TH1D * h_thetamiss = new TH1D("thetamiss","Missing Momentum Polar Angle;#theta_{miss} (degrees);Counts",90,0,90);
  hist_list_1.push_back(h_thetamiss);


  /////////////////////////////////////
  //Histos: candidates
  /////////////////////////////////////
  TH1D * h_mmiss_cand = new TH1D("mmiss_cand","Missing Mass p(e,e'#pi^{+})n;Missing Mass (GeV/c^{2});Events",100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_cand);
  TH1D * h_neff_pmiss_denom = new TH1D("neff_pmiss_denom","Neutron Candidates;p_{miss} (GeV/c);Counts",16,0,1.6);
  hist_list_1.push_back(h_neff_pmiss_denom);
  // efficiency angular dependence
  TH1D * h_neff_thetamiss_denom = new TH1D("neff_thetamiss_denom","#theta dependence of efficiency;#theta_{miss} (deg);efficiency",30,30,90);
  hist_list_1.push_back(h_neff_thetamiss_denom);
  TH1D * h_neff_phimiss_denom = new TH1D("neff_phimiss_denom","#phi dependence of efficiency;#phi;efficiency",36,-180,180);
  hist_list_1.push_back(h_neff_phimiss_denom);
  TH1D * h_nsize = new TH1D("nsize","Number of Reconstructed Neutrons in Event;Neutron number;Counts",10,0,10);
  hist_list_1.push_back(h_nsize);





  /////////////////////////////////////
  //Histos: mmiss vs pmiss by angle int
  /////////////////////////////////////
  // denominator
  TH2D * h_mmiss_pmiss_allt_denom = new TH2D("mmiss_pmiss_allt_denom","Neutron Candidates (all angles);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_allt_denom);
  TH2D * h_mmiss_pmiss_int1_denom = new TH2D("mmiss_pmiss_int1_denom","Neutron Candidates (int1);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int1_denom);
  TH2D * h_mmiss_pmiss_int2_denom = new TH2D("mmiss_pmiss_int2_denom","Neutron Candidates (int2);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int2_denom);
  TH2D * h_mmiss_pmiss_int3_denom = new TH2D("mmiss_pmiss_int3_denom","Neutron Candidates (int3);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int3_denom);
  // numerator
  TH2D * h_mmiss_pmiss_allt_numer = new TH2D("mmiss_pmiss_allt_numer","Neutron Candidates (all angles);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_allt_numer);
  TH2D * h_mmiss_pmiss_int1_numer = new TH2D("mmiss_pmiss_int1_numer","Neutron Candidates (int1);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int1_numer);
  TH2D * h_mmiss_pmiss_int2_numer = new TH2D("mmiss_pmiss_int2_numer","Neutron Candidates (int2);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int2_numer);
  TH2D * h_mmiss_pmiss_int3_numer = new TH2D("mmiss_pmiss_int3_numer","Neutron Candidates (int3);p_{miss} (GeV/c);M_{miss} (GeV/c^{2})",100,p_lo,p_hi,100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_pmiss_int3_numer);

  /////////////////////////////////////
  //Histos: pmiss by angle int
  /////////////////////////////////////
  // denominator
  TH1D * h_pmiss_allt_denomD = new TH1D("pmiss_allt_denom","Neutron Candidates (all angles);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_allt_denomD);
  TH1D * h_pmiss_int1_denomD = new TH1D("pmiss_int1_denom","Neutron Candidates (int1);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int1_denomD);
  TH1D * h_pmiss_int2_denomD = new TH1D("pmiss_int2_denom","Neutron Candidates (int2);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int2_denomD);
  TH1D * h_pmiss_int3_denomD = new TH1D("pmiss_int3_denom","Neutron Candidates (int3);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int3_denomD);
  // numerator
  TH1D * h_pmiss_allt_numerD = new TH1D("pmiss_allt_numer","Neutrons (all angles);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_allt_numerD);
  TH1D * h_pmiss_int1_numerD = new TH1D("pmiss_int1_numer","Neutrons (int1);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int1_numerD);
  TH1D * h_pmiss_int2_numerD = new TH1D("pmiss_int2_numer","Neutrons (int2);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int2_numerD);
  TH1D * h_pmiss_int3_numerD = new TH1D("pmiss_int3_numer","Neutrons (int3);p_{miss} (GeV/c);Counts",neff_pbins,p_lo,p_hi);
  hist_list_1.push_back(h_pmiss_int3_numerD);



  /////////////////////////////////////
  //Histos: neutron required
  /////////////////////////////////////
  TH1D * h_tof = new TH1D("tof_all","Time of Flight;TOF (ns);Counts",200,-10,20);
  hist_list_1.push_back(h_tof);
   TH2D * h_nangles = new TH2D("n angles","Neutron Angular Distribution;Phi;Theta",48,-180,180,110,35,100);
  hist_list_2.push_back(h_nangles);
  TH2D * h_mmiss_beta = new TH2D("mmiss_beta","Missing Mass vs Beta p(e,e'#pi^{+})n;Beta;Missing Mass (GeV/c^{2})",100,0,0.9,100,0,1.5);
  hist_list_2.push_back(h_mmiss_beta);

  /////////////////////////////////////
  //Histos: good neutrons
  /////////////////////////////////////
  TH2D * h_pmiss_pn = new TH2D("pmiss_pn","Missing Momentum vs Neutron Momentum;p_{neutron} (GeV/c);p_{miss} (GeV/c)",100,0,1.4,100,0,1.4);
  hist_list_2.push_back(h_pmiss_pn);
  TH2D * h_beta_pmiss = new TH2D("beta_pmiss","Beta vs Pmiss;Pmiss (GeV/c);Beta",100,0,1.5,100,0,1);
  hist_list_2.push_back(h_beta_pmiss);
  TH1D * h_cos0 = new TH1D("cos0","cos #theta_{pmiss,pneutron};cos #theta;Counts",50,-1.05,1.05);
  hist_list_1.push_back(h_cos0);
  TH1D * h_dphi = new TH1D("dphi","#Delta #phi = #phi_{pmiss} - #phi_{n};#Delta #phi;Counts",48,-180,180);
  hist_list_1.push_back(h_dphi);


  /////////////////////////////////////
  //Histos: effiency
  /////////////////////////////////////
  TH2D * h_mmiss_theta_denom = new TH2D("mmiss_theta_denom","Missing Mass vs #theta_{miss} (denominator);#theta_{miss} (deg);Missing Mass (GeV/c^{2})",100,t_lo,t_hi,100,0.5,1.5);
  TH2D * h_mmiss_theta_numer = new TH2D("mmiss_theta_numer","Missing Mass vs #theta_{miss} (numerator);#theta_{miss} (deg);Missing Mass (GeV/c^{2})",100,t_lo,t_hi,100,0.5,1.5);



  /////////////////////////////////////
  //Histos: final neutron selection
  /////////////////////////////////////
  TH2D * h_p_theta = new TH2D("p theta","Momentum vs Theta (Neutrons);Theta;Momentum (GeV/c)",180,30,150,100,0,1.5);
  hist_list_2.push_back(h_p_theta);
  TH2D * h_p_phi = new TH2D("p phi","Momentum vs Phi (Neutrons);Phi;Momentum (GeV/c)",48,-180,180,100,0,2);
  hist_list_2.push_back(h_p_phi);
  TH2D * h_pmiss_thetamiss = new TH2D("pmiss_thetamiss","Missing Momentum vs p_{miss} Polar Angle;#theta_{miss} (degrees);p_{miss} (GeV/c)",100,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_thetamiss);
  TH2D * h_pmiss_pn_cut = new TH2D("pmiss_pn","Missing Momentum vs Neutron Momentum;p_{neutron};p_{miss}",100,0,1.4,100,0,1.4);
  hist_list_2.push_back(h_pmiss_pn_cut);
  TH1D * h_mmiss_withn = new TH1D("mmiss_withn","Missing Mass p(e,e'#pi^{+}n);Missing Mass (GeV/c^{2})",100,0.5,1.5);
  hist_list_1.push_back(h_mmiss_withn);
  TH1D * h_neff_pmiss_numer = new TH1D("neff_pmiss_numer","Neutrons;p_{miss} (GeV/c);Counts",16,0,1.6);
  hist_list_1.push_back(h_neff_pmiss_numer);
  // efficiency angular dependence
  TH1D * h_neff_thetamiss_numer = new TH1D("neff_thetamiss_numer","CND Neutron Detection Efficiency;Polar Angle #theta (deg);Detection Efficiency",30,30,90);
  hist_list_1.push_back(h_neff_thetamiss_numer);
  TH1D * h_neff_phimiss_numer = new TH1D("neff_phimiss_numer","#phi dependence of efficiency;#phi;efficiency",36,-180,180);
  hist_list_1.push_back(h_neff_phimiss_numer);




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
    auto phot=c12->getByID(22);
    auto neut=c12->getByID(2112);
    auto pip=c12->getByID(211);
    auto allParticles=c12->getDetParticles();
    double weight = 1;
    if(isMC){weight=c12->mcevent()->getWeight();}

    double ts = c12->event()->getStartTime(); // event start time



    // GENERAL EVENT SELECTION
    if(!myCut.electroncut(c12)){continue;}
    if (prot.size() != 0) {continue;} // before p>0 // hen !=0
    if (pip.size() != 1) {continue;} 
    if (elec.size() != 1) {continue;} // before e>1 // then !=1
    //if ((neut.size() != 0) && (neut.size() != 1)) {continue;} // not there before // TEMP



    // ELECTRONS
    double vze = elec[0]->par()->getVz();
    double e_theta = elec[0]->getTheta()*180./M_PI;
    if (e_theta>35.) {continue;}


    // PIONS +
    TVector3 p_pip;
    p_pip.SetMagThetaPhi(pip[0]->getP(),pip[0]->getTheta(),pip[0]->getPhi());
    double vzpi = pip[0]->par()->getVz();
    double dbeta = pip[0]->par()->getBeta() - p_pip.Mag()/sqrt(p_pip.Mag2()+m_piplus*m_piplus);
    double pitheta = pip[0]->getTheta()*180./M_PI;

    // reject particles with the wrong PID
    bool trash = 0;
    for (int i=0; i<allParticles.size(); i++)
      {
      int pid = allParticles[i]->par()->getPid();
      h_pid->Fill(pid);
      if (pid!=2112 && pid!=11 && pid!=211 && pid!=22 && pid!=0) {trash=1;}
      //if (pid!=2112 && pid!=11 && pid!=211) {trash=1;} // TEMP
      }
    if (trash==1) {continue;}

    // pion histos
    h_pivertex->Fill(vzpi-vze,weight);
    h_dbetap->Fill(p_pip.Mag(),dbeta,weight);
    h_pitheta->Fill(pitheta,weight);


    // pion cuts
    if ((vzpi-vze)<-4. || (vzpi-vze)>2.) {continue;}
    if (dbeta<-0.03 || dbeta>0.03) {continue;}   
    if (p_pip.Mag() < 0.4) {continue;}
    if (p_pip.Mag() > 3.) {continue;}
    if (pitheta>35.) {continue;}

    // Missing Mass and energy deposition
    TVector3 p_e;
    TVector3 p_b(0,0,Ebeam);
    p_e.SetMagThetaPhi(elec[0]->getP(),elec[0]->getTheta(),elec[0]->getPhi());
    double mmiss = get_pin_mmiss(p_b,p_e,p_pip);
    TVector3 pmiss = p_b - p_e - p_pip;
    double thetamiss = pmiss.Theta()*180/M_PI;

    // thetamiss histo
    h_thetamiss->Fill(thetamiss,weight);

    // missing momentum cuts
    if (thetamiss<40.) {continue;}
    if (thetamiss>140.) {continue;}
    if (pmiss.Mag() < 0.2) {continue;} // beta=0.2  // previously 0.1918
    if (pmiss.Mag() > 1.2) {continue;} // beta=0.8  // previoulsy 1.2528

    h_mmiss_theta_denom->Fill(thetamiss,mmiss,weight);




   h_mmiss_cand->Fill(mmiss,weight);   

    if (mmiss>0.85 && mmiss<1.05)
    {
      h_neff_thetamiss_denom->Fill(thetamiss,weight);
      h_neff_phimiss_denom->Fill(pmiss.Phi()*180/M_PI);
      h_neff_pmiss_denom->Fill(pmiss.Mag(),weight);
      h_pmiss_allt_denomD->Fill(pmiss.Mag(),weight);
      if (thetamiss>t1 && thetamiss<t2) {h_pmiss_int1_denomD->Fill(pmiss.Mag(),weight);}
      else if (thetamiss>t2 && thetamiss<t3) {h_pmiss_int2_denomD->Fill(pmiss.Mag(),weight);}
      else if (thetamiss>t3 && thetamiss<t4) {h_pmiss_int3_denomD->Fill(pmiss.Mag(),weight);}

    }


    h_mmiss_pmiss_allt_denom->Fill(pmiss.Mag(),mmiss,weight);
    if (thetamiss>t1 && thetamiss<t2) {h_mmiss_pmiss_int1_denom->Fill(pmiss.Mag(),mmiss,weight);}
    else if (thetamiss>t2 && thetamiss<t3) {h_mmiss_pmiss_int2_denom->Fill(pmiss.Mag(),mmiss,weight);}
    else if (thetamiss>t3 && thetamiss<t4) {h_mmiss_pmiss_int3_denom->Fill(pmiss.Mag(),mmiss,weight);}



    // REQUIRE A NEUTRON HERE


    // NEUTRON NUMBER
    double sz = neut.size();
    h_nsize->Fill(sz);


    int pick = -1;
    if (neut.size() < 1) {continue;}
    else
    {
      double lowest_dphi = 180;
      for (int i=0; i<neut.size(); i++)
      {
        // in CND? if no - skip to next neutron in event
        bool is_CND1 = neut[i]->sci(CND1)->getDetector()==3;
        bool is_CND2 = neut[i]->sci(CND2)->getDetector()==3;
        bool is_CND3 = neut[i]->sci(CND3)->getDetector()==3;
        if (!is_CND1 && !is_CND2 && !is_CND3) {continue;}
  
        // in expected theta range? if no - skip to next neutron in event
        double n_theta = neut[i]->getTheta()*180./M_PI;
        if (n_theta==0) {continue;}
        if (n_theta<40) {continue;}
        if (n_theta>140) {continue;}
  
        // is beta high enough? if no - skip to next neutron in event
        double beta_n = neut[i]->par()->getBeta();
        if (beta_n<0.2) {continue;}
  
        // pick neutron with lowest dphi
        double this_dphi = std::abs( neut[i]->getPhi()*180./M_PI - pmiss.Phi()*180./M_PI );
        if (this_dphi < lowest_dphi)
        {
          pick = i;
          lowest_dphi = this_dphi;
        }
      }
    }

    if (pick==-1) {continue;}


    // CND later
    bool is_CND1 = neut[pick]->sci(CND1)->getDetector()==3;
    bool is_CND2 = neut[pick]->sci(CND2)->getDetector()==3;
    bool is_CND3 = neut[pick]->sci(CND3)->getDetector()==3;
    if (!is_CND1 && !is_CND2 && !is_CND3) {continue;}

    // Neutron kinematics
    double n_theta = neut[pick]->getTheta()*180./M_PI;
    double n_phi = neut[pick]->getPhi()*180./M_PI;
    double beta_n = neut[pick]->par()->getBeta();
    double tof_n = 0;
    double path = 0;
    if (is_CND1)
    {
      tof_n = neut[pick]->sci(CND1)->getTime() - ts;
      path = neut[pick]->sci(CND1)->getPath();
    }
    else if (is_CND2)
    {
      tof_n = neut[pick]->sci(CND2)->getTime() - ts;
      path = neut[pick]->sci(CND2)->getPath();
    }
    else if (is_CND3)
    {
      tof_n = neut[pick]->sci(CND3)->getTime() - ts;
      path = neut[pick]->sci(CND3)->getPath();
    }
    if (tof_n<0) {continue;}
    TVector3 pn;
    pn.SetMagThetaPhi(neut[pick]->getP(),neut[pick]->getTheta(),neut[pick]->getPhi());

    if (pn.Mag()<0.2) {continue;}
    if (pn.Mag()>1.2) {continue;}

    TVector3 vecX( neut[pick]->par()->getPx(), neut[pick]->par()->getPy(), neut[pick]->par()->getPz() );
    double cos0 = pmiss.Dot(vecX) / (pmiss.Mag() * vecX.Mag() );
    double dphi = pmiss.Phi()*180./M_PI - neut[pick]->getPhi()*180./M_PI;

    // neutron histos
    h_tof->Fill(tof_n,weight);
    h_nangles->Fill(n_phi,n_theta,weight);
    h_mmiss_beta->Fill(beta_n,mmiss,weight);

    // more neutron histos
    h_pmiss_pn->Fill(pn.Mag(),pmiss.Mag(),weight);
    h_beta_pmiss->Fill(pmiss.Mag(),beta_n,weight);
    h_cos0->Fill(cos0,weight);
    h_dphi->Fill(dphi,weight);

    // exclusive cuts (requiring info from other final state particles)
    if (cos0 < 0.9) {continue;}
    if (std::abs(dphi)>20.) {continue;}


    // histos after requiring pn along pmiss
    

    h_mmiss_theta_numer->Fill(thetamiss,mmiss,weight);


    h_p_theta->Fill(n_theta,pn.Mag(),weight);
    h_p_phi->Fill(n_phi,pn.Mag(),weight);
    //h_pmiss_theta->Fill(neut[pick]->getTheta()*180./M_PI,pmiss.Mag(),weight);
    h_pmiss_thetamiss->Fill(thetamiss,pmiss.Mag(),weight);
    h_pmiss_pn_cut->Fill(pn.Mag(),pmiss.Mag(),weight);
    h_mmiss_withn->Fill(mmiss,weight);
    if (mmiss>0.85 && mmiss<1.05)
    {
      h_neff_thetamiss_numer->Fill(thetamiss,weight);
      h_neff_phimiss_numer->Fill(pmiss.Phi()*180/M_PI);
      h_neff_pmiss_numer->Fill(pmiss.Mag(),weight);
      h_pmiss_allt_numerD->Fill(pmiss.Mag(),weight);
      if (thetamiss>t1 && thetamiss<t2) {h_pmiss_int1_numerD->Fill(pmiss.Mag(),weight);}
      else if (thetamiss>t2 && thetamiss<t3) {h_pmiss_int2_numerD->Fill(pmiss.Mag(),weight);}
      else if (thetamiss>t3 && thetamiss<t4) {h_pmiss_int3_numerD->Fill(pmiss.Mag(),weight);}

    }


    h_mmiss_pmiss_allt_numer->Fill(pmiss.Mag(),mmiss,weight);
    if (thetamiss>t1 && thetamiss<t2) {h_mmiss_pmiss_int1_numer->Fill(pmiss.Mag(),mmiss,weight);}
    else if (thetamiss>t2 && thetamiss<t3) {h_mmiss_pmiss_int2_numer->Fill(pmiss.Mag(),mmiss,weight);}
    else if (thetamiss>t3 && thetamiss<t4) {h_mmiss_pmiss_int3_numer->Fill(pmiss.Mag(),mmiss,weight);}


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
  //Histos
  /////////////////////////////////////////////////////
 

  // pions
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"final state: 0p, 1pi+, 1e");
  text.DrawLatex(0.2,0.5,"#theta_{e} < 35");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pid->Draw();
  myCanvas->cd(2);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 
 
  // pions
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"allow only PID=0,11,22,211,2112");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pivertex->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
 
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dbetap->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 



  
  // all photons and neutrons
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen run 015017");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"-4 cm < v_{#pi+}-v_{e} < 2 cm");
  text.DrawLatex(0.2,0.5,"-0.03 < #Delta #beta < 0.03");
  text.DrawLatex(0.2,0.4,"p_{#pi+} > 0.4 GeV/c");
  text.DrawLatex(0.2,0.3,"#theta_{#pi} < 35");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetamiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen run 015017");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"pion cuts");
  text.DrawLatex(0.2,0.5,"40 deg < #theta_{miss} < 140 deg");
  text.DrawLatex(0.2,0.4,"0.094 GeV/c < p_{miss} < 1.25 GeV/c");
  myText->Print(fileName,"pdf");
  myText->Clear();


  ///////////////////////////
  // Denom: momentum dependence
  ///////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Get n_{eff} vs p (denominator)");
  myText->Print(fileName,"pdf");
  myText->Clear();

  // all angles
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_allt_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // theta interval 1
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int1_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  double * sigd_int1 = hist_projections(myCanvas,h_mmiss_pmiss_int1_denom,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int1_denom = new TH1D("pmiss_int1_denom","Efficiency denominator (int1);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int1_denom->SetBinContent(i,*(sigd_int1+i));
//std::cout << "int 1, val " << i << " = " << *(sigd_int1+i) << '\n';
    h_pmiss_int1_denom->SetBinError(i,sqrt(*(sigd_int1+i)));
  }

  // theta interval 2
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int2_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  double * sigd_int2 = hist_projections(myCanvas,h_mmiss_pmiss_int2_denom,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int2_denom = new TH1D("pmiss_int2_denom","Efficiency denominator (int2);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int2_denom->SetBinContent(i,*(sigd_int2+i));
//std::cout << "int 2, val " << i << " = " << *(sigd_int2+i) << '\n';
    h_pmiss_int2_denom->SetBinError(i,sqrt(*(sigd_int2+i)));
  }

  // theta interval 3
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int3_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  myCanvas->cd(1);
  double * sigd_int3 = hist_projections(myCanvas,h_mmiss_pmiss_int3_denom,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int3_denom = new TH1D("pmiss_int3_denom","Efficiency denominator (int3);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int3_denom->SetBinContent(i,*(sigd_int3+i));
//std::cout << "int 3, val " << i << " = " << *(sigd_int3+i) << '\n';
    h_pmiss_int3_denom->SetBinError(i,sqrt(*(sigd_int3+i)));
  }

  // denominator - all theta intervals
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_int1_denom->Draw();
  myCanvas->cd(2);
  h_pmiss_int2_denom->Draw();
  myCanvas->cd(3);
  h_pmiss_int3_denom->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  


  ///////////////////////////
  // Denom: theta dependence
  ///////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Get n_{eff} vs #theta (denominator)");
  myText->Print(fileName,"pdf");
  myText->Clear();

  // theta denominator
  myText->cd();
  text.DrawLatex(0.2,0.9,"Theta denominator");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_theta_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(tgrid_x,tgrid_y);
  myCanvas->cd(1);
  double * sigd_t = hist_projections(myCanvas,h_mmiss_theta_denom,neff_tbins, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();







  // histos requiring a neutron start here
  
  myText->cd();
  text.DrawLatex(0.2,0.5,"p(e,e'#pi^{+}n)");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nsize->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen run 015017");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"pion cuts");
  text.DrawLatex(0.2,0.5,"p_{miss}, M_{miss} cuts");
  text.DrawLatex(0.2,0.4,"Require at least 1 neutron in CND");
  text.DrawLatex(0.2,0.3,"Neutron in at least 1 lever of CND");
  text.DrawLatex(0.2,0.2,"Pick neutron closest to p_{miss} in #phi");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nangles->Draw("colz");
  h_nangles->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_beta->Draw("colz");
  h_mmiss_beta->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_tof->Draw("colz");
  h_tof->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();




  // all neutrons
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen run 015017");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"pion cuts");
  text.DrawLatex(0.2,0.5,"p_{miss}, M_{miss} cuts");
  text.DrawLatex(0.2,0.4,"Require 1 neutron in CND");
  text.DrawLatex(0.2,0.3,"exclude #theta_{n}=0, #phi_{n}=0");
  text.DrawLatex(0.2,0.2,"40 deg < #theta_{n} < 140 deg");
  text.DrawLatex(0.2,0.1,"0.1 < #beta_{n} < 0.8");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_pn->Draw("colz");
  TF1 * line = new TF1("line","x",0,2);
  line->Draw("same");
  myCanvas->cd(3);
  h_beta_pmiss->Draw("colz");
  h_beta_pmiss->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,1);
  myCanvas->cd(1);
  h_cos0->Draw();
  h_cos0->SetStats(0);
  myCanvas->cd(2);
  h_dphi->Draw();
  h_dphi->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // neutrons in direction of pmiss
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Hydrogen run 015017");
  text.DrawLatex(0.2,0.7,"p(e,e'#pi^{+})n");
  text.DrawLatex(0.2,0.6,"pion cuts");
  text.DrawLatex(0.2,0.5,"p_{miss}, M_{miss} cuts");
  text.DrawLatex(0.2,0.4,"exclude #theta_{n}=0, #phi_{n}=0");
  text.DrawLatex(0.2,0.3,"40 deg < #theta_{n} < 140 deg");
  text.DrawLatex(0.2,0.2,"0.1 < #beta_{n} < 0.8");
  text.DrawLatex(0.2,0.1,"cos #theta_{pmiss,pn} > 0.9");
  myText->Print(fileName,"pdf");
  myText->Clear();



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_pn_cut->Draw("colz");
  h_pmiss_pn_cut->SetStats(0);
  line->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_theta->Draw("colz");
  h_p_theta->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_phi->Draw("colz");
  h_p_phi->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_thetamiss->Draw("colz");
  h_pmiss_thetamiss->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  



  ///////////////////////////
  // Numer: momentum dependence
  ///////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Get n_{eff} vs p (numerator)");
  myText->Print(fileName,"pdf");
  myText->Clear();

  // all angles
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_allt_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // all angles
  myCanvas->Divide(pgrid_x,pgrid_y);
  double * justhereforthehist = hist_projections(myCanvas,h_mmiss_pmiss_allt_numer,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // theta interval 1
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int1_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  double * sign_int1 = hist_projections(myCanvas,h_mmiss_pmiss_int1_numer,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int1_numer = new TH1D("pmiss_int1_numer","Efficiency numerator (int1);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int1_numer->SetBinContent(i,*(sign_int1+i));
    h_pmiss_int1_numer->SetBinError(i,sqrt(*(sign_int1+i)));
  }

  // theta interval 2
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int2_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  double * sign_int2 = hist_projections(myCanvas,h_mmiss_pmiss_int2_numer,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int2_numer = new TH1D("pmiss_int2_numer","Efficiency numerator (int2);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int2_numer->SetBinContent(i,*(sign_int2+i));
    h_pmiss_int2_numer->SetBinError(i,sqrt(*(sign_int2+i)));
  }

  // theta interval 3
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss_int3_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(pgrid_x,pgrid_y);
  double * sign_int3 = hist_projections(myCanvas,h_mmiss_pmiss_int3_numer,neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  TH1D * h_pmiss_int3_numer = new TH1D("pmiss_int3_numer","Efficiency numerator (int3);p_{miss} (GeV/c)",neff_pbins,p_lo,p_hi);
  for (int i=0; i<neff_pbins; i++)
  {
    h_pmiss_int3_numer->SetBinContent(i,*(sign_int3+i));
    h_pmiss_int3_numer->SetBinError(i,sqrt(*(sign_int3+i)));
  }

  // numerator - all theta intervals
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_int1_numer->Draw();
  myCanvas->cd(2);
  h_pmiss_int2_numer->Draw();
  myCanvas->cd(3);
  h_pmiss_int3_numer->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Theta numerator");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_theta_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(tgrid_x,tgrid_y);
  hist_projections(myCanvas,h_mmiss_theta_numer,neff_tbins, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();






  // efficiency as a function of phi miss
  myCanvas->Divide(2,1);
  myCanvas->cd(1);
  h_neff_phimiss_denom->Draw();
  myCanvas->cd(2);
  h_neff_phimiss_numer->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_neff_phimiss = (TH1D*)h_neff_phimiss_numer->Clone();
  h_neff_phimiss->Divide(h_neff_phimiss_denom);
  h_neff_phimiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // missing mass for candidates and detected
  myText->cd();
  text.DrawLatex(0.2,0.9,"Denominator and numerator missing mass");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_withn->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_cand->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // efficiency results
  myText->cd();
  text.DrawLatex(0.2,0.9,"Efficiency results");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_neff_pmiss_numer->Draw();
  h_neff_pmiss_numer->SetStats(0);
  myCanvas->cd(2);
  h_neff_pmiss_denom->Draw();
  h_neff_pmiss_denom->SetStats(0);
  myCanvas->cd(3);
  TH1D * h_neff_pmiss = (TH1D*)h_neff_pmiss_numer->Clone();
  h_neff_pmiss->Divide(h_neff_pmiss_denom);
  h_neff_pmiss->Draw();
  h_neff_pmiss->SetStats(0);
  h_neff_pmiss->GetYaxis()->SetTitle("efficiency");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TH1D * h_neff_pmiss_int1 = (TH1D*)h_pmiss_int1_numer->Clone();
  h_neff_pmiss_int1->Divide(h_pmiss_int1_denom);
  h_neff_pmiss_int1->SetLineColor(kMagenta);
  h_neff_pmiss_int1->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_pmiss_int1->Draw();
  myCanvas->cd(2);
  TH1D * h_neff_pmiss_int2 = (TH1D*)h_pmiss_int2_numer->Clone();
  h_neff_pmiss_int2->Divide(h_pmiss_int2_denom);
  h_neff_pmiss_int2->SetLineColor(kGreen);
  h_neff_pmiss_int2->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_pmiss_int2->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_pmiss_int3 = (TH1D*)h_pmiss_int3_numer->Clone();
  h_neff_pmiss_int3->Divide(h_pmiss_int3_denom);
  h_neff_pmiss_int3->SetLineColor(kBlue);
  h_neff_pmiss_int3->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_pmiss_int3->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_neff_pmiss_int1->Draw();
  h_neff_pmiss_int2->Draw("same");
  h_neff_pmiss_int3->Draw("same");
  h_neff_pmiss_int3->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_pmiss_int3->SetStats(0);
  TLegend * leg2 = new TLegend(0.65,0.75,0.89,0.89);
  leg2->SetTextFont(72);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h_neff_pmiss_int1,"40<#theta_{miss}<50","l");
  leg2->AddEntry(h_neff_pmiss_int2,"50<#theta_{miss}<60","l");
  leg2->AddEntry(h_neff_pmiss_int3,"60<#theta_{miss}<70","l");
  //leg2->AddEntry(h_neff_pmiss_int4,"70<#theta_{miss}<80","l");
  leg2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();




  // efficiency as a function of theta miss
  myCanvas->Divide(2,1);
  myCanvas->cd(1);
  h_neff_thetamiss_denom->Draw();
  myCanvas->cd(2);
  h_neff_thetamiss_numer->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_neff_thetamiss = (TH1D*)h_neff_thetamiss_numer->Clone();
  h_neff_thetamiss->Divide(h_neff_thetamiss_denom);
  h_neff_thetamiss->Draw();
  h_neff_thetamiss->SetStats(0);
  // print output
  ofstream outtheta(theta_name);
  for (int i=0; i<h_neff_thetamiss->GetNbinsX(); i++) {
    outtheta << h_neff_thetamiss->GetXaxis()->GetBinCenter(i) << ' ' << h_neff_thetamiss->GetBinContent(i) << ' ' << h_neff_thetamiss->GetBinError(i) << '\n';
  }
  outtheta.close();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();




  // original method
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_allt_numerD->Draw();
  myCanvas->cd(2);
  h_pmiss_int1_numerD->Draw();
  myCanvas->cd(3);
  h_pmiss_int2_numerD->Draw();
  myCanvas->cd(4);
  h_pmiss_int3_numerD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // original method
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_allt_denomD->Draw();
  myCanvas->cd(2);
  h_pmiss_int1_denomD->Draw();
  myCanvas->cd(3);
  h_pmiss_int2_denomD->Draw();
  myCanvas->cd(4);
  h_pmiss_int3_denomD->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // original method
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  TH1D * h_neff_allt_D = (TH1D*)h_pmiss_allt_numerD->Clone();
  h_neff_allt_D->Divide(h_pmiss_allt_denomD);
  h_neff_allt_D->Draw();
  // print output
  ofstream outallp(p_name);
  for (int i=0; i<h_neff_allt_D->GetNbinsX(); i++) {
    outallp << h_neff_allt_D->GetXaxis()->GetBinCenter(i) << ' ' << h_neff_allt_D->GetBinContent(i) << ' ' << h_neff_allt_D->GetBinError(i) << '\n';
  }
  outallp.close();


  myCanvas->cd(2);
  TH1D * h_neff_int1_D = (TH1D*)h_pmiss_int1_numerD->Clone();
  h_neff_int1_D->Divide(h_pmiss_int1_denomD);
  h_neff_int1_D->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_int2_D = (TH1D*)h_pmiss_int2_numerD->Clone();
  h_neff_int2_D->Divide(h_pmiss_int2_denomD);
  h_neff_int2_D->Draw();
   myCanvas->cd(4);
  TH1D * h_neff_int3_D = (TH1D*)h_pmiss_int3_numerD->Clone();
  h_neff_int3_D->Divide(h_pmiss_int3_denomD);
  h_neff_int3_D->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // original method
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_neff_int1_D->SetLineColor(kGreen);
  h_neff_int1_D->SetStats(0);
  h_neff_int1_D->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_int1_D->Draw();
  h_neff_int2_D->SetLineColor(kBlue);
  h_neff_int2_D->SetStats(0);
  h_neff_int2_D->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_int2_D->Draw("same");
  h_neff_int3_D->SetLineColor(kRed);
  h_neff_int3_D->SetStats(0);
  h_neff_int3_D->GetYaxis()->SetRangeUser(0,0.16);
  h_neff_int3_D->Draw("same");
  TLegend * leg3 = new TLegend(0.65,0.75,0.89,0.89);
  leg3->SetTextFont(72);
  leg3->SetTextSize(0.04);
  leg3->AddEntry(h_neff_int1_D,"40<#theta_{miss}<50","l");
  leg3->AddEntry(h_neff_int2_D,"50<#theta_{miss}<60","l");
  leg3->AddEntry(h_neff_int3_D,"60<#theta_{miss}<70","l");
  leg3->Draw();
    // print output
  ofstream outp("neff_out/p_hepin_int.txt");
  for (int i=0; i<h_neff_int1_D->GetNbinsX(); i++) {
    outp << h_neff_int1_D->GetXaxis()->GetBinCenter(i) << ' ';
    outp << h_neff_int1_D->GetBinContent(i) << ' ' << h_neff_int1_D->GetBinError(i) << ' ';
    outp << h_neff_int2_D->GetBinContent(i) << ' ' << h_neff_int2_D->GetBinError(i) << ' ';
    outp << h_neff_int3_D->GetBinContent(i) << ' ' << h_neff_int3_D->GetBinError(i) << '\n';
  }
  outp.close();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





std::cout << h_mmiss_withn->Integral(0,200) << '\n';
std::cout << h_mmiss_cand->Integral(0,200) << '\n';
std::cout << "The neutron efficiency is " << h_mmiss_withn->Integral(0,200)/h_mmiss_cand->Integral(0,200) << '\n';



  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
}


Double_t background_poly(Double_t *x, Double_t *par) {
  return ( par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] );
}

Double_t signal(Double_t *x, Double_t *par) {
  return par[0]*exp(-pow((x[0]-par[1]),2.)/(2*pow(par[2],2.))); 
}

Double_t mmiss_signal_gauss(Double_t *x, Double_t *par) {
  return signal(x,par) + signal(x,&par[3]);
}

Double_t mmiss_signal_poly(Double_t *x, Double_t *par) {
  return signal(x,par) + background_poly(x,&par[3]);
}




void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}


double get_pin_mmiss(TVector3 p_b, TVector3 p_e, TVector3 ppi){

  double Ebeam = p_b.Mag();
  double Ee = p_e.Mag();
  double Epi = sqrt(ppi.Mag2() + m_piplus*m_piplus);
  double emiss = Ebeam - Ee + mp - Epi;
  TVector3 pmiss = p_b - p_e - ppi;

  double mmiss = sqrt( (emiss*emiss) - pmiss.Mag2() );

  return mmiss;
}


double * hist_projections(TCanvas * can, TH2D * hist2d, int num_hist, char v)
{
  double p_start_val[num_hist];
  double x_min = hist2d->GetXaxis()->GetXmin();
  double x_max = hist2d->GetXaxis()->GetXmax();
  double dp = (x_max-x_min)/num_hist;
  double * S = new double[num_hist];
  for (int i=0; i<num_hist; i++)
  {
    p_start_val[i] = x_min + i*dp;
    int bin1 = hist2d->GetXaxis()->FindBin(p_start_val[i]);
    int bin2 = hist2d->GetXaxis()->FindBin(p_start_val[i]+dp) - 1;

    // make projection for x interval
    can->cd(i+1);
    TH1D * proj = hist2d->ProjectionY("",bin1,bin2,"d");

    // create name of missing mass histogram for current momentum/theta interval
    std::ostringstream sObj1, sObj2;
    std::string leftTitle = "Missing Mass in ("; std::string midTitle = ",";
    std::string rightTitle;
    if (v=='p')
    {
      rightTitle = ") GeV/c";
      sObj1 << std::fixed << std::setprecision(3) << p_start_val[i];
      sObj2 << std::fixed << std::setprecision(3) << p_start_val[i] + dp;
    }
    else if (v=='t')
    {
      rightTitle = ") deg";
      sObj1 << std::fixed << std::setprecision(0) << p_start_val[i];
      sObj2 << std::fixed << std::setprecision(0) << p_start_val[i] + dp;
    }
    else
    {
      std::cout << "Invalid projection variable for missing mass\n";
    }
    std::string result = leftTitle + sObj1.str() + midTitle + sObj2.str() + rightTitle;
    proj->SetTitle(result.c_str());

    // draw
    proj->Draw();
    S[i] = proj->Integral(proj->GetXaxis()->FindBin(0.85),proj->GetXaxis()->FindBin(1.05));
  }
  return S;
}


