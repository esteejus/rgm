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
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TChain.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TColor.h>
#include <TLine.h>

#include "clas12reader.h"
#include "HipoChain.h"
#include "eventcut/eventcut.h"
#include "eventcut/functions.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

const double c = 29.9792458;
const double mp = 0.938272;
const double m_piplus = 0.13957039;


// efficiency constants
int neff_pbins = 20;
int neff_tbins = 20;
int tbins_ana = 9;
int pgrid_x = ceil(sqrt(neff_pbins));
int pgrid_y = ceil((double)neff_pbins/(double)pgrid_x);
int tgrid_x = ceil(sqrt(neff_tbins));
int tgrid_y = ceil((double)neff_tbins/(double)tgrid_x);
double Mlow = 0.85;
double Mhigh = 1.05;
double theta_lo = 40;  // 40 for CD p
double theta_hi = 140;
double p_lo = 0.25;
double p_hi = 0.8;
double Mdisp_lo = 0.7;
double Mdisp_hi = 1.3;
double xb_cut = 0.4;



void printProgress(double percentage);
double get_pin_mmiss(TVector3 p_b, TVector3 p_e, TVector3 ppi);
Double_t lorentzian(Double_t *x, Double_t *par);
Double_t poly(Double_t *x, Double_t *par);
Double_t signal(Double_t *x, Double_t *par); 
Double_t mmiss_signal_gauss(Double_t *x, Double_t *par);
Double_t mmiss_signal_poly(Double_t *x, Double_t *par);
Double_t mmiss_signal_lorentz(Double_t *x, Double_t *par);
double * hist_projections(TCanvas * can, TH2D * hist2d, int num_hist, char v);
double * hist_projections_backsub(TCanvas * can, TH2D * hist2d, int num_hist, bool subtract_bk, char v);




void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data =0> <bgnd subtraction =1, no =0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/input.hipo> \n";
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

  bool backsub = false;
  if(atoi(argv[2]) == 1){backsub=true;}
 
  double Ebeam = atof(argv[3]);

  // output file names
  TFile * outFile = new TFile(argv[4],"RECREATE");
  char * pdfFile = argv[5];

  char * basename = argv[4];
  basename[strlen(basename)-5] = '\0';
  string theta_name(string(basename) + "_theta.txt");
  theta_name.c_str();
  string p_name(string(basename) + "_p.txt");
  p_name.c_str();

  clas12root::HipoChain chain;

  // clas12ana
  clas12ana clasAna;

  clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsSF_LD2_x2.dat");
  clasAna.readEcalPPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsPI_LD2_x2.dat");

  //clasAna.printParams();

  clasAna.setProtonPidCuts(true);



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
  //Histo: protons
  /////////////////////////////////////
  TH1D * h_pvertex_cd = new TH1D("pvertex_cd","Proton vertex - electron vertex (CD);v_{p} - v_{e} (cm);Counts",50,-4,4);
  hist_list_1.push_back(h_pvertex_cd);
  TH1D * h_pvertex_fd = new TH1D("pvertex_fd","Proton vertex - electron vertex (FD);v_{p} - v_{e} (cm);Counts",50,-4,4);
  hist_list_1.push_back(h_pvertex_fd);
  TH2D * h_dbetap_cd = new TH2D("dbetap_cd","#Delta #beta vs Momentum (CD);Momentum (GeV/c);#beta_{meas} - p/sqrt(p^{2}+m^{2})",100,0,2,100,-0.15,0.15);
  hist_list_2.push_back(h_dbetap_cd);
  TH2D * h_dbetap_fd = new TH2D("dbetap_fd","#Delta #beta vs Momentum (FD);Momentum (GeV/c);#beta_{meas} - p/sqrt(p^{2}+m^{2})",100,0,2,100,-0.15,0.15);
  hist_list_2.push_back(h_dbetap_fd);

  TH1D * h_pchipid = new TH1D("p chipid","Proton #Chi^{2} PID",100,-3,3);
  hist_list_1.push_back(h_pchipid);
  TH2D * h_thetamiss_xb = new TH2D("thetamiss_xb","#theta_{pred} vs xB;xB;#theta_{pred} (deg)",100,0,2,120,30,150);
  hist_list_2.push_back(h_thetamiss_xb);
  TH2D * h_pmiss_xb = new TH2D("pmiss_xb","p_{pred} vs xB;xB;p_{pred} (GeV/c)",100,0,2,100,0,1.5);
  hist_list_2.push_back(h_pmiss_xb);

  TH2D * h_pangles_cd = new TH2D("pangles_cd","Proton Angular Distribution (CD);Phi;Theta",180,-180,180,90,0,180);
  hist_list_2.push_back(h_pangles_cd);
  TH2D * h_pangles_fd = new TH2D("pangles_fd","Proton Angular Distribution (FD);Phi;Theta",180,-180,180,90,0,180);
  hist_list_2.push_back(h_pangles_fd);

  TH2D * h_pmiss_thetamiss_cd = new TH2D("pmiss_thetamiss_cd","Distribution of Expected Neutrons (CD);#theta_{pred} (deg);p_{pred} (GeV/c)",120,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_thetamiss_cd);
  TH2D * h_pmiss_thetamiss_fd = new TH2D("pmiss_thetamiss_fd","Distribution of Expected Neutrons (FD);#theta_{pred} (deg);p_{pred} (GeV/c)",120,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_thetamiss_fd);



  
  /////////////////////////////////////
  //Histo: all neutrals
  /////////////////////////////////////
  TH2D * h_mmiss_xb = new TH2D("mmiss_xb_cand","Missing Mass vs xB;xB;Missing Mass",100,0,2,100,0.2,2);
  hist_list_2.push_back(h_mmiss_xb);
  TH1D * h_mmiss_withn = new TH1D("mmiss_withn","Missing Mass d(e,e'pn);Missing Mass (GeV/c^{2})",100,0.,2.);
  hist_list_1.push_back(h_mmiss_withn);
  TH1D * h_tof = new TH1D("tof_all","TOF (ns)",60,-10,20);
  hist_list_1.push_back(h_tof);


  
  /////////////////////////////////////
  //Histo: all neutrons
  /////////////////////////////////////
  TH1D * h_nsize = new TH1D("nsize","Number of reconstructed neutrons",5,0,5);
  hist_list_1.push_back(h_nsize);
  TH2D * h_nangles = new TH2D("n angles","Neutron Angular Distribution;Phi;Theta",48,-180,180,120,30,150);
  hist_list_2.push_back(h_nangles);
  TH1D * h_nphi = new TH1D("n phi","Neutron Azimuthal Angle;Phi (deg);Counts",40,-180,180);
  hist_list_1.push_back(h_nphi);





  /////////////////////////////////////
  // Histos: missing mass vs momentum
  /////////////////////////////////////
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","Missing Mass vs p_{pred} d(e,e'p)n;p_{pred} (GeV/c);Missing Mass (GeV/c^{2})",100,0.25,1,100,0.2,1.5);
  hist_list_2.push_back(h_mmiss_pmiss);
  TH2D * h_mmiss_pmissCAND = new TH2D("mmiss_pmissCAND","Missing Mass vs p_{pred} d(e,e'p)n;p_{pred} (GeV/c);Missing Mass (GeV/c^{2})",100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_pmissCAND);
  TH2D * h_mmiss_pmissDET = new TH2D("mmiss_pmissDET","Missing Mass vs p_{pred} d(e,e'p)n;p_{pred} (GeV/c);Missing Mass (GeV/c^{2})",100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_pmissDET);
  TH1D * h_thetamiss_before = new TH1D("thetamiss_before","Predicted Momentum Polar Angle;#theta_{pred};Counts",100,40,140);
  hist_list_1.push_back(h_thetamiss_before);

  TH2D * h_pmissthetamiss = new TH2D("pmissthetamiss","Predicted Momentum Magnitude vs Angle;#theta_{pred};p_{pred} (GeV/c)",110,35,145,100,0.2,1.3);
  hist_list_2.push_back(h_pmissthetamiss);
  


  /////////////////////////////////////
  // Histos: missing mass vs momentum
  /////////////////////////////////////
  TH2D * h_mmiss_theta = new TH2D("mmiss_theta","Predicted Momentum vs #theta;#theta_{pred};Missing Mass (GeV/c^{2})",140,0,140,50,0.,2.);
  hist_list_2.push_back(h_mmiss_theta);
  TH2D * h_mmiss_thetaCAND = new TH2D("mmiss_thetaCAND","Predicted Momentum vs #theta;#theta_{pred};Missing Mass (GeV/c^{2})",100,theta_lo,theta_hi,50,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_thetaCAND);
  TH2D * h_mmiss_thetaDET = new TH2D("mmiss_thetaDET","Predicted Momentum vs #theta;#theta_{pred};Missing Mass (GeV/c^{2})",100,theta_lo,theta_hi,50,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_thetaDET);





  /////////////////////////////////////
  //Histo: compare with pmiss
  /////////////////////////////////////
  TH2D * h_pmiss_pn = new TH2D("pmiss_pn","Predicted Momentum vs Neutron Momentum;p_{neutron} (GeV/c);p_{pred} (GeV/c)",100,0.2,1.25,100,0.2,1.25);
  hist_list_2.push_back(h_pmiss_pn);
  TH1D * h_cos0 = new TH1D("cos0","cos #theta_{pred,n}",50,-1.05,1.05);
  hist_list_1.push_back(h_cos0);
  TH2D * h_pmiss_pn_cut = new TH2D("pmiss_pn_cut","Predicted Momentum vs Neutron Momentum;p_{n} (GeV/c);p_{pred} (GeV/c)",100,0.25,1.25,100,0.25,1.25);
  hist_list_2.push_back(h_pmiss_pn_cut);
  TH2D * h_pmiss_theta = new TH2D("pmiss_theta","Predicted Momentum vs #theta;#theta_{n} (deg);Predicted Momentum p_{pred} (GeV/c)",100,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_theta);
  TH2D * h_dtheta_dphi = new TH2D("dtheta_dphi","#Delta #theta vs #Delta #phi;#phi_{pred}-#phi_{n} (deg);#theta_{pred}-#theta_{n} (deg)",100,-180,180,100,-110,110);
  hist_list_2.push_back(h_dtheta_dphi);
  TH2D * h_theta_edep = new TH2D("theta_edep","Energy Deposition vs #theta;#theta_{n} (deg);Energy Deposition (MeVee)",100,0,180,100,0,30);
  hist_list_2.push_back(h_theta_edep);
 


  TH2D * h_dp_theta = new TH2D("dp_theta","Neutron Momentum Resolution vs #theta;#theta_{n} (deg);p_{pred} - p_{n} (GeV/c)",110,35,145,100,-1,1);
  hist_list_2.push_back(h_dp_theta);
  TH2D * h_dp_thetamiss = new TH2D("dp_thetamiss","Neutron Momentum Resolution vs #theta_{pred};#theta_{pred} (deg);p_{pred} - p_{n} (GeV/c)",110,35,145,100,-1,1);
  hist_list_2.push_back(h_dp_thetamiss);
  TH2D * h_dpp_theta = new TH2D("dpp_theta","Neutron Momentum Error vs #theta_{n};#theta_{n} (deg);(p_{pred} - p_{n})/p_{pred}",110,35,145,100,-0.5,0.5);
  hist_list_2.push_back(h_dpp_theta);

  TH2D * h_compare = new TH2D("compare","Comparison between p_{n} and p_{pred};p_{pred}-p_{n};Angle between p_{pred} and p_{n}",100,-1,1,100,0,50);
  hist_list_2.push_back(h_compare);
  TH2D * h_dp_pmiss = new TH2D("dp_pmiss","Neutron Momentum Error vs p_{pred};p_{pred} (GeV/c);p_{pred} - p_{n}",100,0.2,1,100,-0.3,0.3);
  hist_list_2.push_back(h_dp_pmiss);

  TH2D * h_dpp_edep = new TH2D("dpp_edep","Neutron Momentum Error vs Energy;Energy Deposition (MeVee);(p_{pred} - p_{n})/p_{pred}",100,0,50,100,-0.5,0.5);
  hist_list_2.push_back(h_dpp_edep);
  TH2D * h_thetapn_edep = new TH2D("thetapn_edep","Proton-Neutron Angle vs Energy;Energy Deposition (MeVee);Angle Between Proton and Neutron (deg)",100,0,50,90,0,180);
  hist_list_2.push_back(h_thetapn_edep);
  TH2D * h_cos0_edep = new TH2D("cos0_edep","cos #theta_{pred,n};Energy Deposition (MeVee);cos #theta_{pred,n}",100,0,50,100,-1,1);
  hist_list_2.push_back(h_cos0_edep);
  TH2D * h_cos0_pn = new TH2D("cos0_pn","cos #theta_{pred,n};Neutron Momentum p_{n} (GeV/c);cos #theta_{pred,n}",100,0.2,1.2,100,-1,1);
  hist_list_2.push_back(h_cos0_pn);
  TH2D * h_pn_edep = new TH2D("pn_edep","Neutron Momentum vs Energy Deposition;Energy Deposition (MeVee);p_{n} (GeV/c)",100,0,50,100,0.2,1.3);
  hist_list_2.push_back(h_pn_edep);

  TH2D * h_theta_thetamiss = new TH2D("theta_thetamiss","#theta_{n} - #theta_{pred} (deg);#theta_{pred} (deg);#theta_{n}-#theta_{pred} (deg)",100,40,140,100,-30,30);
  hist_list_2.push_back(h_theta_thetamiss);



  /////////////////////////////////////
  //Histo: cutting out background
  /////////////////////////////////////
  TH1D * h_theta_ppmiss = new TH1D("theta_ppmiss","Angle between p_{p} and p_{pred};#theta_{p,pred}",180,0,180);
  hist_list_1.push_back(h_theta_ppmiss);
  TH1D * h_theta_np = new TH1D("theta_np","Angle between neutron and proton",45,0,180);
  hist_list_1.push_back(h_theta_np);





  /////////////////////////////////////
  //Histo: neutron efficiency by angular range
  /////////////////////////////////////
  TH1D * neff_denom_ang[9]; TH1D * neff_numer_ang[9];
  TH2D * mmiss_pmiss_CAND9[9]; TH2D * mmiss_pmiss_DET9[9];
  std::vector<double> ang_range = {50,55,60,65,70,75,80,85,100,120};
  double ang_vals[9] = {52.5,57.5,62.5,67.5,72.5,77.5,82.5,92.5,110};
  double xerr[9] = {2.5,2.5,2.5,2.5,2.5,2.5,2.5,7.5,10};
  for (int i=0; i<9; i++){
    // numerator and denominator vs pmiss
    sprintf(temp_name,"Efficiency (%d-%d deg);Predicted Neutron Momentum (GeV/c);Efficiency",int(ang_range[i]),int(ang_range[i+1]));
    sprintf(temp_title,"Neutron Efficiency (%d-%d deg);Predicted Neutron Momentum (GeV/c);Efficiency",int(ang_range[i]),int(ang_range[i+1]));
    neff_denom_ang[i] = new TH1D(temp_name,temp_title,neff_pbins,0.25,1.0);
    neff_numer_ang[i] = new TH1D(temp_name,temp_title,neff_pbins,0.25,1.0);
    hist_list_1.push_back(neff_denom_ang[i]);
    hist_list_1.push_back(neff_numer_ang[i]);
    // mmiss vs pmiss - for background subtraction
    sprintf(temp_name,"mmiss_pmiss (%d-%d deg);p_{pred} (GeV/c);M_{miss} (GeV/c^{2})",int(ang_range[i]),int(ang_range[i+1]));
    sprintf(temp_title,"M_{miss} vs p_{pred} (%d-%d deg);p_{pred} (GeV/c);M_{miss} (GeV/c^{2})",int(ang_range[i]),int(ang_range[i+1]));
    mmiss_pmiss_CAND9[i] = new TH2D(temp_name,temp_title,100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
    mmiss_pmiss_DET9[i] = new TH2D(temp_name,temp_title,100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
  }


 


  /////////////////////////////////////
  //Histo: neutron efficiency
  /////////////////////////////////////
  TH1D * h_neff_pmiss_numer_ssb = new TH1D("neff_pm_numer_ssb","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_ssb);
  TH1D * h_neff_pmiss_denom_ssb = new TH1D("neff_pm_denom_ssb","Neutron Candidates;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_denom_ssb);

  TH1D * h_neff_pmiss_numer_2 = new TH1D("neff_pm_numer_2","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_2);
  TH1D * h_neff_pmiss_numer_4 = new TH1D("neff_pm_numer_4","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_4);
  TH1D * h_neff_pmiss_numer_6 = new TH1D("neff_pm_numer_6","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_6);
  TH1D * h_neff_pmiss_numer_8 = new TH1D("neff_pm_numer_8","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_8);
  TH1D * h_neff_pmiss_numer_10 = new TH1D("neff_pm_numer_10","Neutrons;p_{pred} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_10);








  /////////////////////////////////////
  //Histo: neutron efficiency by theta
  /////////////////////////////////////
  TH1D * h_neff_thetamiss_denom_ssb = new TH1D("neff_tm_denom_ssb","Neutron Candidates;#theta_{pred} (deg);Counts",neff_tbins,theta_lo,theta_hi);
  hist_list_1.push_back(h_neff_thetamiss_denom_ssb);
  TH1D * h_neff_thetamiss_numer_ssb = new TH1D("neff_tm_numer_ssb","Detected Neutrons;#theta_{pred} (deg);Counts",neff_tbins,theta_lo,theta_hi);
  hist_list_1.push_back(h_neff_thetamiss_numer_ssb);


  /////////////////////////////////////
  //Histos: 2d efficiency stuff
  /////////////////////////////////////
  TH2D * h_cand2d = new TH2D("cand2d","Neutron Candidates;p_{pred} (GeV/c);#theta_{pred} (degrees)",15,0.2,1.2,30,40,140);
  hist_list_2.push_back(h_cand2d);
  h_cand2d->SetStats(0);
  TH2D * h_det2d = new TH2D("det2d","Detected Neutrons;p_{pred} (GeV/c);#theta_{pred} (degrees)",15,0.2,1.2,30,40,140);
  hist_list_2.push_back(h_det2d);
  h_det2d->SetStats(0);





  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Sumw2();
    hist_list_1[i]->GetXaxis()->CenterTitle();
    hist_list_1[i]->GetYaxis()->CenterTitle();
    //hist_list_1[i]->SetStats(0);
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Sumw2();
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->CenterTitle();
    //hist_list_2[i]->SetStats(0);
  }


  int counter = 0;



  // Define cut class
  while(chain.Next()==true){
    // display completed
    counter++;
    if((counter%1000000) == 0){
      cerr << '\n' << counter/1000000 << " million completed";
    }

    if((counter%100000) == 0){
      cerr << ".";
    }

    // get particles by type
    clasAna.Run(c12);
    auto elec = clasAna.getByPid(11);
    auto prot = clasAna.getByPid(2212);
    auto neut = clasAna.getByPid(2112);
    auto phot = clasAna.getByPid(22);
    auto allParticles=c12->getDetParticles();
    double weight = 1;
    if(isMC){weight=c12->mcevent()->getWeight();}

    double ts = c12->event()->getStartTime();



    // GENERAL EVENT SELECTION

    if(prot.size()!=1) {continue;}
    if(elec.size()!=1) {continue;}
    

    // ELECTRONS
    double vze = elec[0]->par()->getVz();
    TVector3 p_e;
    TVector3 p_b(0,0,Ebeam);
    p_e.SetMagThetaPhi(elec[0]->getP(),elec[0]->getTheta(),elec[0]->getPhi());
    TVector3 p_q = p_b - p_e;
    double nu = Ebeam - p_e.Mag();
    double QSq = p_q.Mag2() - (nu*nu);
    double xB = QSq / (2*mN*nu);

    // PROTONS
    TVector3 pp;
    pp.SetMagThetaPhi(prot[0]->getP(),prot[0]->getTheta(),prot[0]->getPhi());
    double vzp = prot[0]->par()->getVz();
    double dbeta = prot[0]->par()->getBeta() - pp.Mag()/sqrt(pp.Mag2()+mp*mp);
    double chipid = prot[0]->par()->getChi2Pid();

    // reject particles with the wrong PID
    bool trash = 0;
    for (int i=0; i<allParticles.size(); i++)
    {
      int pid = allParticles[i]->par()->getPid();
      if (pid!=2112 && pid!=11 && pid!=2212 && pid!=22 && pid!=0) {trash=1;} // 22,0
    }
    if (trash==1) {continue;}


    // Missing mass and energy despotiion
    double mmiss = get_mmiss(p_b,p_e,pp);
    TVector3 pmiss = p_q - pp;
    double thetamiss = pmiss.Theta()*180./M_PI;


    // check out protons in CD and FD separately    
    if (prot[0]->getRegion()==CD)
    {
      h_pvertex_cd->Fill(vzp-vze,weight);
      if (abs(vzp-vze)>2) {continue;} 
      h_dbetap_cd->Fill(pp.Mag(),dbeta,weight);
      if (abs(dbeta)>0.07) {continue;}
      if (pp.Mag()<0.25) {continue;}
      h_pangles_cd->Fill(pp.Phi()*180./M_PI,pp.Theta()*180./M_PI);
      h_pmiss_thetamiss_cd->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag(),weight);

    }
    else if (prot[0]->getRegion()==FD)
    {
      h_pvertex_fd->Fill(vzp-vze,weight);
      if (abs(vzp-vze)>3) {continue;}
      h_dbetap_fd->Fill(pp.Mag(),dbeta,weight);
      if (abs(dbeta)>0.03) {continue;}
      if (pp.Mag()<0.25) {continue;}
      h_pangles_fd->Fill(pp.Phi()*180./M_PI,pp.Theta()*180./M_PI);
      h_pmiss_thetamiss_fd->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag(),weight);

    }

    h_pchipid->Fill(chipid,weight);
    h_thetamiss_xb->Fill(xB,thetamiss,weight);
    h_pmiss_xb->Fill(xB,pmiss.Mag(),weight);




    // cut on theta component of pmiss
    h_thetamiss_before->Fill(thetamiss,weight);
    h_pmissthetamiss->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
    if (thetamiss<40.) {continue;}
    if (thetamiss>140.) {continue;}

    // allow only missing momentum in the correct range
    if (pmiss.Mag() < 0.25) {continue;} // beta = 0.2, p = 0.1918; beta = 0.25, p = 0.243
    if (pmiss.Mag() > 1.25) {continue;} // beta = 0.8, p = 1.2528

    // xb cut
    h_mmiss_xb->Fill(xB,mmiss,weight);
    if (!backsub && xB<xb_cut) {continue;}


    h_theta_ppmiss->Fill(pmiss.Angle(pp)*180./M_PI,weight);
    //if ((pmiss.Angle(pp)*180./M_PI)<40.) {continue;}  // CD ONLY (I THINK)
    



    // fill mmiss vs p/theta
    h_mmiss_pmiss->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_pmissCAND->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_theta->Fill(thetamiss,mmiss,weight);
    h_mmiss_thetaCAND->Fill(thetamiss,mmiss,weight);

    if (backsub) // background subtraction - fill denominator
    {
      h_neff_pmiss_denom_ssb->Fill(pmiss.Mag(),weight);
      h_neff_thetamiss_denom_ssb->Fill(thetamiss,weight);
      h_cand2d->Fill(pmiss.Mag(),thetamiss,weight);
    }
    else // no background subtraction - explicit Mmiss cut
    {
      if (mmiss>Mlow && mmiss<Mhigh)
      {
        h_neff_pmiss_denom_ssb->Fill(pmiss.Mag(),weight);
        h_neff_thetamiss_denom_ssb->Fill(thetamiss,weight);
        h_cand2d->Fill(pmiss.Mag(),thetamiss,weight);
      }
    }

    // fill denominator by angular ranges
    for (int k=0; k<9; k++)
    {
      if (thetamiss>(ang_range[k]) && thetamiss<(ang_range[k+1]))
      {
        if (backsub) {mmiss_pmiss_CAND9[k]->Fill(pmiss.Mag(),mmiss,weight);}
        if (!backsub && mmiss>Mlow && mmiss<Mhigh) {neff_denom_ang[k]->Fill(pmiss.Mag(),weight);}
      }
    }





    // REQUIRE A NEUTRON STARTING HERE


    // NEUTRONS
    double sz = neut.size();
    h_nsize->Fill(sz,weight);

    int pick = -1;

    if (neut.size() < 1){continue;}
    else
    {
      double lowest_dphi = 180;
      for (int i=0; i<neut.size(); i++)
      {
        // in CND or CTOF? if no - skip to next neutron in event
        bool is_CND1 = neut[i]->sci(CND1)->getDetector()==3;
        bool is_CND2 = neut[i]->sci(CND2)->getDetector()==3;
        bool is_CND3 = neut[i]->sci(CND3)->getDetector()==3;
        bool is_CTOF = neut[i]->sci(CTOF)->getDetector()==4;
        if (!is_CND1 && !is_CND2 && !is_CND3 && !is_CTOF) {continue;}

        // in expected theta range? if no - skip to next neutron in event
        double n_theta = neut[i]->getTheta()*180./M_PI;
        if (n_theta==0) {continue;}
        if (n_theta<40) {continue;}
        if (n_theta>140) {continue;}

        // is beta high enough? if no - skip to next neutron in event
        double beta_n = neut[i]->par()->getBeta();
        if (beta_n<0.25) {continue;}
        
        // pick neutron with lowest dphi
        double this_dphi = std::abs( neut[i]->getPhi()*180./M_PI - pmiss.Phi()*180./M_PI );
        if (this_dphi < lowest_dphi)
        {
          pick = i;
          lowest_dphi = this_dphi;
        }
      } // end loop over neutrons
    }

    if (pick==-1) {continue;}


    // get neutron information
    bool is_CND1 = neut[pick]->sci(CND1)->getDetector()==3;
    bool is_CND2 = neut[pick]->sci(CND2)->getDetector()==3;
    bool is_CND3 = neut[pick]->sci(CND3)->getDetector()==3;
    bool is_CTOF = neut[pick]->sci(CTOF)->getDetector()==4;
    double n_theta = neut[pick]->getTheta()*180./M_PI;
    double beta_n = neut[pick]->par()->getBeta();
    double n_phi = neut[pick]->getPhi()*180./M_PI;
    double tof_n = 0;
    double edep = 0;

    if (is_CND1)
    {
      tof_n = neut[pick]->sci(CND1)->getTime() - ts;
      edep = edep + neut[pick]->sci(CND1)->getEnergy();
    }
    else if (is_CND2)
    {
      tof_n = neut[pick]->sci(CND2)->getTime() - ts;
      edep = edep + neut[pick]->sci(CND2)->getEnergy();
    }
    else if (is_CND3)
    {
      tof_n = neut[pick]->sci(CND3)->getTime() - ts;
      edep = edep + neut[pick]->sci(CND3)->getEnergy();
    }
    else if (is_CTOF)
    {
      tof_n = neut[pick]->sci(CTOF)->getTime() - ts;
      edep = edep + neut[pick]->sci(CTOF)->getEnergy();
    }

    h_tof->Fill(tof_n,weight);


    // define neutron kinematics
    TVector3 pn;
    pn.SetMagThetaPhi(neut[pick]->getP(),neut[pick]->getTheta(),neut[pick]->getPhi());
    TVector3 vecX( neut[pick]->par()->getPx(), neut[pick]->par()->getPy(), neut[pick]->par()->getPz() );
    double cos0 = pmiss.Dot(vecX) / (pmiss.Mag() * vecX.Mag() );
    double dphi = pmiss.Phi()*180./M_PI - neut[pick]->getPhi()*180./M_PI;
    double dp = (pmiss.Mag()-pn.Mag());
    double dpp = dp/pmiss.Mag();



    // compare neutrons to predicted neutron momentum
    h_pmiss_pn->Fill(pn.Mag(),pmiss.Mag(),weight);
    h_dtheta_dphi->Fill(dphi,pmiss.Theta()*180./M_PI - neut[pick]->getTheta()*180./M_PI,weight);
    h_theta_edep->Fill(neut[pick]->getTheta()*180./M_PI,edep,weight);
    h_cos0->Fill(cos0,weight);


    // energy deposition cut
    h_pn_edep->Fill(edep,pn.Mag(),weight);
    h_dpp_edep->Fill(edep,(pmiss.Mag()-pn.Mag())/pmiss.Mag(),weight);
    h_thetapn_edep->Fill(edep,pn.Angle(pp)*180./M_PI,weight);
    h_cos0_edep->Fill(edep,cos0,weight);

    if (edep<1.5) {continue;}

    // pn / pmiss magnitude and angular cut
    h_compare->Fill(dp,pn.Angle(pmiss)*180./M_PI);
    h_dp_pmiss->Fill(pmiss.Mag(),dp,weight);
    h_cos0_pn->Fill(pn.Mag(),cos0,weight);

    if (dp<-0.2 || dp>0.2) {continue;}
    if (pn.Angle(pmiss)*180./M_PI>25) {continue;}


    // after requiring pn along pmiss
    h_nangles->Fill(n_phi,n_theta,weight);
    h_nphi->Fill(n_phi,weight);
    h_pmiss_pn_cut->Fill(pn.Mag(),pmiss.Mag(),weight);


    // look at angular
    h_dp_theta->Fill(n_theta,(pmiss.Mag()-pn.Mag()),weight);
    h_dp_thetamiss->Fill(pmiss.Theta()*180./M_PI,(pmiss.Mag()-pn.Mag()),weight);
    h_dpp_theta->Fill(n_theta,dpp,weight);
    h_theta_thetamiss->Fill(pmiss.Theta()*180./M_PI,(n_theta-pmiss.Theta()*180./M_PI),weight);




    h_mmiss_withn->Fill(mmiss,weight);
    h_pmiss_theta->Fill(neut[pick]->getTheta()*180./M_PI,pmiss.Mag(),weight);

    // cut out fake neutron background
    double theta_np = pn.Angle(pp)*180./M_PI;
    h_theta_np->Fill(theta_np,weight);


    // fill mmiss vs p/theta
    h_mmiss_pmissDET->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_thetaDET->Fill(thetamiss,mmiss,weight);

    if (backsub) // background subtraction - fill numerator
    {
      h_neff_pmiss_numer_ssb->Fill(pmiss.Mag(),weight);
      h_neff_thetamiss_numer_ssb->Fill(thetamiss,weight);
      h_det2d->Fill(pmiss.Mag(),thetamiss,weight);
    }
    else // no background subtraction - explicit Mmiss cut
    {
      if (mmiss>Mlow && mmiss<Mhigh)  
      {
        h_neff_pmiss_numer_ssb->Fill(pmiss.Mag(),weight);
        h_neff_thetamiss_numer_ssb->Fill(thetamiss,weight);
        h_det2d->Fill(pmiss.Mag(),thetamiss,weight);
        if (edep>2) {h_neff_pmiss_numer_2->Fill(pmiss.Mag(),weight);}
        if (edep>4) {h_neff_pmiss_numer_4->Fill(pmiss.Mag(),weight);}
        if (edep>6) {h_neff_pmiss_numer_6->Fill(pmiss.Mag(),weight);}
        if (edep>8) {h_neff_pmiss_numer_8->Fill(pmiss.Mag(),weight);}
        if (edep>10) {h_neff_pmiss_numer_10->Fill(pmiss.Mag(),weight);}
      }
    }

    // fill numerator by angular ranges
    for (int k=0; k<9; k++)
    {
      if (thetamiss>(ang_range[k]) && thetamiss<(ang_range[k+1]))
      {
        if (backsub) {mmiss_pmiss_DET9[k]->Fill(pmiss.Mag(),mmiss,weight);}
        if (!backsub && mmiss>Mlow && mmiss<Mhigh) {neff_numer_ang[k]->Fill(pmiss.Mag(),weight);}
      }
    }


  } // end event loop



  cout<<counter<<endl;

  outFile->cd();
  for(int i=0; i<hist_list_1.size(); i++){
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++){
    hist_list_2[i]->Write();
  }

  h_mmiss_pmissDET->Write();






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
  
  /////////////////////
  // PROTONS
  /////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"d(e,e'pn) 6 GeV");
  text.DrawLatex(0.2,0.8,"Protons");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pvertex_cd->Draw();
  TLine * l_vtz_cd1 = new TLine(-2,0,-2,h_pvertex_cd->GetMaximum());
  l_vtz_cd1->SetLineColor(kRed);
  l_vtz_cd1->SetLineWidth(3);
  l_vtz_cd1->Draw("same");
  TLine * l_vtz_cd2 = new TLine(2,0,2,h_pvertex_cd->GetMaximum());
  l_vtz_cd2->SetLineColor(kRed);
  l_vtz_cd2->SetLineWidth(3);
  l_vtz_cd2->Draw("same");
  myCanvas->cd(2);
  h_dbetap_cd->Draw("colz");
  TLine * l_dbcd1 = new TLine(0,-0.07,2,-0.07);
  l_dbcd1->SetLineColor(kRed);
  l_dbcd1->SetLineWidth(3);
  l_dbcd1->Draw("same");
  TLine * l_dbcd2 = new TLine(0,0.07,2,0.07);
  l_dbcd2->SetLineColor(kRed);
  l_dbcd2->SetLineWidth(3);
  l_dbcd2->Draw("same");
  TLine * l_dbcd3 = new TLine(0.25,-0.15,0.25,0.15);
  l_dbcd3->SetLineColor(kRed);
  l_dbcd3->SetLineWidth(3);
  l_dbcd3->Draw("same");
  myCanvas->cd(3);
  h_pangles_cd->Draw("colz");
  myCanvas->cd(4);
  h_pmiss_thetamiss_cd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pvertex_fd->Draw();
  TLine * l_vtz_fd1 = new TLine(-3,0,-3,h_pvertex_fd->GetMaximum());
  l_vtz_fd1->SetLineColor(kRed);
  l_vtz_fd1->SetLineWidth(3);
  l_vtz_fd1->Draw("same");
  TLine * l_vtz_fd2 = new TLine(3,0,3,h_pvertex_fd->GetMaximum());
  l_vtz_fd2->SetLineColor(kRed);
  l_vtz_fd2->SetLineWidth(3);
  l_vtz_fd2->Draw("same");
  myCanvas->cd(2);
  h_dbetap_fd->Draw("colz");
  TLine * l_dbfd1 = new TLine(0,-0.03,2,-0.03);
  l_dbfd1->SetLineColor(kRed);
  l_dbfd1->SetLineWidth(3);
  l_dbfd1->Draw("same");
  TLine * l_dbfd2 = new TLine(0,0.03,2,0.03);
  l_dbfd2->SetLineColor(kRed);
  l_dbfd2->SetLineWidth(3);
  l_dbfd2->Draw("same");
  TLine * l_dbfd3 = new TLine(0.25,-0.15,0.25,0.15);
  l_dbfd3->SetLineColor(kRed);
  l_dbfd3->SetLineWidth(3);
  l_dbfd3->Draw("same");

  myCanvas->cd(3);
  h_pangles_fd->Draw("colz");
  myCanvas->cd(4);
  h_pmiss_thetamiss_fd->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pchipid->Draw();
  myCanvas->cd(3);
  myCanvas->cd(3)->SetLogz();
  h_thetamiss_xb->Draw("colz");
  myCanvas->cd(4);
  myCanvas->cd(4)->SetLogz();
  h_pmiss_xb->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

 


  /////////////////////
  // NEUTRALS
  /////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Expected neutrons / missing momentum");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetamiss_before->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_thetamiss_cd->Draw("colz");
  h_pmiss_thetamiss_cd->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_thetamiss_fd->Draw("colz");
  h_pmiss_thetamiss_fd->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_xb->Draw("colz");
  if (!backsub)
  {
    TLine * l_xb = new TLine(xb_cut,0.2,xb_cut,2.0);
    l_xb->SetLineColor(kRed);
    l_xb->SetLineWidth(3);
    l_xb->Draw("same");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // NEUTRON CANDIDATES
  // vs MOMENTUM
  /////////////////////
  
  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY BY MOMENTUM - CANDIDATES");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmiss->Draw("colz");
  if (!backsub)
  {
    TLine * l_mpmiss1 = new TLine(0.25,0.85,1,0.85);
    l_mpmiss1->SetLineColor(kRed);
    l_mpmiss1->SetLineWidth(3);
    l_mpmiss1->Draw("same");
    TLine * l_mpmiss2 = new TLine(0.25,1.05,1,1.05);
    l_mpmiss2->SetLineColor(kRed);
    l_mpmiss2->SetLineWidth(3);
    l_mpmiss2->Draw("same");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_pmissCAND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
 

  // mmiss in pmiss bins - background subtracted
  myCanvas->Divide(pgrid_x,pgrid_y);
  double *S_pCAND = hist_projections(myCanvas, h_mmiss_pmissCAND, neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // angle between p and pmiss
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_theta_ppmiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY BY THETA - CANDIDATES");
  myText->Print(fileName,"pdf");
  myText->Clear();

  /////////////////////
  // NEUTRON CANDIDATES
  // vs THETA
  /////////////////////
  
  // mmiss by theta bins
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // mmiss by theta bins
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_thetaCAND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  

  // mmiss by theta bins - background subtracted
  myCanvas->Divide(tgrid_x,tgrid_y);
  double *S_tCAND = hist_projections(myCanvas, h_mmiss_thetaCAND, neff_tbins, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // START REQUIRING
  // NEUTRONS
  /////////////////////


  myText->cd();
  text.DrawLatex(0.2,0.9,"Detected neutrons");
  text.DrawLatex(0.2,0.8,"Compare p_{pred} to p_{n}");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nsize->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_tof->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_pmiss_pn->Draw("colz");
  TF1 * line = new TF1("line","x",0,2);
  line->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dtheta_dphi->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_theta_edep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  



  // neutron selection
  myText->cd();
  text.DrawLatex(0.2,0.9,"Detected neutrons");
  text.DrawLatex(0.2,0.8,"Neutron cuts");
  myText->Print(fileName,"pdf");
  myText->Clear();


  // look at energy deposition
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pn_edep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dpp_edep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_thetapn_edep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_cos0_edep->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // momentum resolution
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  //myCanvas->cd(1)->SetLogz();
  h_compare->Draw("colz");
  TLine * l_comp1 = new TLine(-0.2,0,-0.2,50);
  l_comp1->SetLineColor(kRed);
  l_comp1->SetLineWidth(3);
  l_comp1->Draw("same");
  TLine * l_comp2 = new TLine(0.2,0,0.2,50);
  l_comp2->SetLineColor(kRed);
  l_comp2->SetLineWidth(3);
  l_comp2->Draw("same");
  TLine * l_comp3 = new TLine(-1,25,1,25);
  l_comp3->SetLineColor(kRed);
  l_comp3->SetLineWidth(3);
  l_comp3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_dp_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_cos0_pn->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // final neutron selection
  myText->cd();
  text.DrawLatex(0.2,0.9,"Detected neutrons");
  text.DrawLatex(0.2,0.8,"Final neutron selection");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nphi->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_pn_cut->Draw("colz");
  line->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // final neutron selection
  myText->cd();
  text.DrawLatex(0.2,0.9,"Detected neutrons");
  text.DrawLatex(0.2,0.8,"Angular resolution");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dp_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dp_thetamiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dpp_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_theta_thetamiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // neutrons in direction of pmiss
  myText->cd();
  text.DrawLatex(0.2,0.9,"Neutrons in direction of pmiss");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_withn->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pmiss_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_theta_np->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();




  /////////////////////
  // DETECTED NEUTRON
  // vs PMISS
  /////////////////////

  
  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY BY MOMENTUM - DETECTED");
  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmissDET->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // mmiss in pmiss bins - background subtracted
  myCanvas->Divide(pgrid_x,pgrid_y);
  double * S_pDET = hist_projections(myCanvas, h_mmiss_pmissDET, neff_pbins, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // DETECTED NEUTRON
  // vs THETA
  /////////////////////
  
  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY BY THETA - DETECTED");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_thetaDET->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
 


  myCanvas->Divide(tgrid_x,tgrid_y);
  double * S_tDET = hist_projections(myCanvas, h_mmiss_thetaDET, neff_tbins, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // EFFICIENCY RESULTS
  /////////////////////

  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY RESULTS");
  myText->Print(fileName,"pdf");
  myText->Clear();



  // BACKGROUND SUBTRACTION IN 9 ANGULAR RANGES

  // mmiss in pmiss bins - background subtracted
  if (backsub)
  {
  for (int i=0; i<9; i++)
  {
    // background subtraction for candidates
    myCanvas->Divide(pgrid_x,pgrid_y);
    double *Ssub_pCAND00 = hist_projections_backsub(myCanvas, mmiss_pmiss_CAND9[i], neff_pbins, 0, 'p');
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    // fill denominator with background-subtracted values
    for (int j=0; j<neff_pbins; j++)
    {
      neff_denom_ang[i]->SetBinContent(j,*(Ssub_pCAND00+j));
      neff_denom_ang[i]->SetBinError(i,std::sqrt(*(Ssub_pCAND00+j)));
    }

    // background subtraction for detected neutrons
    myCanvas->Divide(pgrid_x,pgrid_y);
    double *Ssub_pDET00 = hist_projections_backsub(myCanvas, mmiss_pmiss_DET9[i], neff_pbins, 0, 'p');
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    // fill numerator with background-subtracted values
    for (int j=0; j<neff_pbins; j++)
    {
      neff_numer_ang[i]->SetBinContent(j,*(Ssub_pDET00+j));
      neff_numer_ang[i]->SetBinError(j,std::sqrt(*(Ssub_pDET00+j)));
    }

  } // end loop over 9 angular ranges
  } // end "if backsub"



  /*if (backsub) // efficiency vs pmiss using background-subtracted numerator and denominator
  {
    myCanvas->Divide(2,2);
    for (int i=0; i<neff_pbins; i++)
    {
      // NUMERATOR
      // background subtraction
      h_neff_pmiss_numer_ssb->SetBinContent(i,*(Ssub_pDET+i));
      h_neff_pmiss_numer_ssb->SetBinError(i,std::sqrt(*(Ssub_pDET+i)));
      // DENOMINATOR
      // background subtraction
      h_neff_pmiss_denom_ssb->SetBinContent(i,*(Ssub_pCAND+i));
      h_neff_pmiss_denom_ssb->SetBinError(i,std::sqrt(*(Ssub_pCAND+i)));
    }
  }*/
  myCanvas->cd(1);
  h_neff_pmiss_numer_ssb->Draw();
  myCanvas->cd(2);
  h_neff_pmiss_denom_ssb->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_pmiss = (TH1D*)h_neff_pmiss_numer_ssb->Clone();
  h_neff_pmiss->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss->Draw();
  h_neff_pmiss->SetTitle("Neutron Efficiency");
  h_neff_pmiss->GetYaxis()->SetTitle("efficiency");
  h_neff_pmiss->GetYaxis()->SetRangeUser(0.,0.16);
  // print output
  ofstream outp(p_name);
  for (int i=0; i<h_neff_pmiss->GetNbinsX(); i++)
  {
    outp << h_neff_pmiss->GetXaxis()->GetBinCenter(i) << ' ';
    outp << h_neff_pmiss->GetBinContent(i) << ' ';
    outp << h_neff_pmiss->GetBinError(i) << '\n';
  }
  outp.close();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // efficiency vs theta using background-subtracted numerator and denominator
  myCanvas->Divide(2,2);
  /*if (backsub)
  {
    for (int i=0; i<neff_tbins; i++)
    {
        // NUMERATOR
        // background subtraction
        h_neff_thetamiss_numer_ssb->SetBinContent(i,*(Ssub_tDET+i));
        h_neff_thetamiss_numer_ssb->SetBinError(i,std::sqrt(*(Ssub_tDET+i)));
        // DENOMINATOR
        // background subtraction
        h_neff_thetamiss_denom_ssb->SetBinContent(i,*(Ssub_tCAND+i));
        h_neff_thetamiss_denom_ssb->SetBinError(i,std::sqrt(*(Ssub_tCAND+i)));
    }
  }*/
  myCanvas->cd(1);
  h_neff_thetamiss_numer_ssb->Draw();
  myCanvas->cd(2);
  h_neff_thetamiss_denom_ssb->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_thetamiss_ssb = (TH1D*)h_neff_thetamiss_numer_ssb->Clone();
  h_neff_thetamiss_ssb->Divide(h_neff_thetamiss_denom_ssb);
  h_neff_thetamiss_ssb->Draw();
  h_neff_thetamiss_ssb->SetTitle("Neutron Efficiency");
  h_neff_thetamiss_ssb->GetYaxis()->SetTitle("efficiency");
  h_neff_thetamiss_ssb->GetYaxis()->SetRangeUser(0.,0.16);
  // print output
  ofstream outtheta(theta_name);
  for (int i=0; i<h_neff_thetamiss_ssb->GetNbinsX(); i++) {
    outtheta << h_neff_thetamiss_ssb->GetXaxis()->GetBinCenter(i) << ' ';
    outtheta << h_neff_thetamiss_ssb->GetBinContent(i) << ' ';
    outtheta << h_neff_thetamiss_ssb->GetBinError(i) << '\n';
  }
  outtheta.close();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  




  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_neff_pmiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_neff_thetamiss_ssb->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // testing Edep dependence
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  // 2 MeVee
  TH1D* h_neff_pmiss_2 = (TH1D*)h_neff_pmiss_numer_2->Clone();
  h_neff_pmiss_2->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss_2->SetLineColor(kRed);
  h_neff_pmiss_2->SetStats(0);
  h_neff_pmiss_2->Draw();
  h_neff_pmiss_2->GetYaxis()->SetRangeUser(0,0.16);
  // 4 MeVee
  TH1D * h_neff_pmiss_4 = (TH1D*)h_neff_pmiss_numer_4->Clone();
  h_neff_pmiss_4->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss_4->SetLineColor(kOrange);
  h_neff_pmiss_4->SetStats(0);
  h_neff_pmiss_4->Draw("same");
  h_neff_pmiss_4->GetYaxis()->SetRangeUser(0,0.16);
  // 6 MeVee
  TH1D * h_neff_pmiss_6 = (TH1D*)h_neff_pmiss_numer_6->Clone();
  h_neff_pmiss_6->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss_6->SetLineColor(kGreen);
  h_neff_pmiss_6->SetStats(0);
  h_neff_pmiss_6->Draw("same");
  h_neff_pmiss_6->GetYaxis()->SetRangeUser(0,0.16);
  // 8 MeVee
  TH1D * h_neff_pmiss_8 = (TH1D*)h_neff_pmiss_numer_8->Clone();
  h_neff_pmiss_8->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss_8->SetLineColor(kBlue);
  h_neff_pmiss_8->SetStats(0);
  h_neff_pmiss_8->Draw("same");
  h_neff_pmiss_8->GetYaxis()->SetRangeUser(0,0.16);
  // 10 MeVee
  TH1D * h_neff_pmiss_10 = (TH1D*)h_neff_pmiss_numer_10->Clone();
  h_neff_pmiss_10->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss_10->SetLineColor(kViolet);
  h_neff_pmiss_10->SetStats(0);
  h_neff_pmiss_10->Draw("same");
  h_neff_pmiss_10->GetYaxis()->SetRangeUser(0,0.16);
  // legend
  TLegend * leg = new TLegend(0.65,0.65,0.89,0.89);
  leg->SetTextFont(72);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_neff_pmiss_2,"2 MeVee");
  leg->AddEntry(h_neff_pmiss_4,"4 MeVee");
  leg->AddEntry(h_neff_pmiss_6,"6 MeVee");
  leg->AddEntry(h_neff_pmiss_8,"8 MeVee");
  leg->AddEntry(h_neff_pmiss_10,"10 MeVee");
  leg->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // 2D EFFICIENCY
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cand2d->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_det2d->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH2D * h_eff2d = (TH2D*)h_det2d->Clone();
  h_eff2d->Divide(h_cand2d);
  h_eff2d->Draw("colz");
  h_eff2d->SetMaximum(0.20);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();








  // draw neff for 9 angular ranges
  TH1D * neff_ang[9];
  double par[3], p0[9], p0_err[9], p1[9], p1_err[9], p2[9], p2_err[9], ang[9];
  myCanvas->Divide(3,3);
  TF1 * f_poly2 = new TF1("poly2","[0]+[1]*x+[2]*x*x",0.3,1.0);
  //TF1 * f_poly2 = new TF1("poly2","[0]+[1]*x",0.3,1.0);
  myCanvas->cd(1);
  for (int i=0; i<9; i++)
  {
    myCanvas->cd(i+1);
    // draw neff for angular range
    neff_ang[i] = (TH1D*)neff_numer_ang[i]->Clone();
    neff_ang[i]->Divide(neff_denom_ang[i]);
    neff_ang[i]->SetMinimum(0.0);
    neff_ang[i]->SetMaximum(0.18);
    neff_ang[i]->SetLineColor(kBlue-i);
    neff_ang[i]->Draw();
    // fit data to function
    if (!backsub) {f_poly2->SetParLimits(2,0,1);} // set limit only with xB cut - low stats
    f_poly2->SetLineColor(kBlue-i);
    neff_ang[i]->Fit("poly2","Q","",0.3,1.0);
    // get parameters
    f_poly2->GetParameters(par);
    p0[i] = par[0];  p0_err[i] = f_poly2->GetParError(0);
    p1[i] = par[1];  p1_err[i] = f_poly2->GetParError(1);
    p2[i] = par[2];  p2_err[i] = f_poly2->GetParError(2);
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // draw fit parameters, fit to degree-2 polynomial
  std::cout << "parameter 0 fit\n";
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TGraphErrors *gr1 = new TGraphErrors(9,ang_vals,p0,xerr,p0_err);
  gr1->Draw("AP");
  gr1->SetTitle("fit parameter a_{0}");
  gr1->GetXaxis()->SetTitle("Neutron Polar Angle (deg)");
  gr1->GetYaxis()->SetTitle("fit parameter 0");
  f_poly2->SetParLimits(0,0,1);  f_poly2->SetParLimits(1,-0.1,0);  f_poly2->SetParLimits(2,0,0.1);
  gr1->Fit("poly2","","",50,140);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  std::cout << "parameter 1 fit\n";
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TGraphErrors *gr2 = new TGraphErrors(9,ang_vals,p1,xerr,p1_err);
  gr2->Draw("AP");
  gr2->SetTitle("fit parameter a_{1}");
  gr2->GetXaxis()->SetTitle("Neutron Polar Angle (deg)");
  gr2->GetYaxis()->SetTitle("fit parameter 1");
  f_poly2->SetParLimits(0,0,0.1);  f_poly2->SetParLimits(1,-0.1,0);  f_poly2->SetParLimits(2,0,0.1);
  gr2->Fit("poly2","","",50,140);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  std::cout << "parameter 2 fit\n";
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TGraphErrors *gr3 = new TGraphErrors(9,ang_vals,p2,xerr,p2_err);
  gr3->Draw("AP");
  gr3->SetTitle("fit parameter a_{2}");
  gr3->GetXaxis()->SetTitle("Neutron Polar Angle (deg)");
  gr3->GetYaxis()->SetTitle("fit parameter 2");
  f_poly2->SetParLimits(0,-1,0);  f_poly2->SetParLimits(1,0,0.1);  f_poly2->SetParLimits(2,-0.1,0);
  gr3->Fit("poly2","","",50,140);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // efficiency by angular range
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  for (int i=0; i<9; i++)
  {
    neff_ang[i]->Draw("same");
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");
  //outFile->Close();
}






Double_t lorentzian(Double_t *x, Double_t *par) { // height, mean, width
  return (0.5*par[0]*par[2]/TMath::Pi()) / TMath::Max( 1.e-10,(x[0]-par[1])*(x[0]-par[1]) + .25*par[2]*par[2] );
}

Double_t poly(Double_t *x, Double_t *par) {
  return ( par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] );
}

Double_t signal(Double_t *x, Double_t *par) { // height, mean, width
  return par[0]*exp(-pow((x[0]-par[1]),2.)/(2*pow(par[2],2.))); 
}

Double_t landau(Double_t *x, Double_t *par) {
  return par[0]*TMath::Landau(x[0],par[1],par[2]); 
}

Double_t mmiss_signal_gauss(Double_t *x, Double_t *par) {
  return signal(x,par) + signal(x,&par[3]);
}

Double_t mmiss_landau_gauss(Double_t *x, Double_t *par) {
  return par[0]*TMath::Landau(x[0],par[1],par[2]) + signal(x,&par[3]);
}

Double_t mmiss_signal_poly(Double_t *x, Double_t *par) {
  return signal(x,par) + poly(x,&par[3]);
}

Double_t mmiss_signal_lorentz(Double_t *x, Double_t *par) {
  return signal(x,par) + lorentzian(x,&par[3]);
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
  double * SSB = new double[num_hist];
  // plot and fit each graph
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
    SSB[i] = proj->Integral(Mlow,Mhigh);
  }
  return SSB;
}



double * hist_projections_backsub(TCanvas * can, TH2D * hist2d, int num_hist, bool subtract_bk, char v)
{
  double p_start_val[num_hist];
  double x_min = hist2d->GetXaxis()->GetXmin();
  double x_max = hist2d->GetXaxis()->GetXmax();
  double dp = (x_max-x_min)/num_hist;
  double * S = new double[num_hist];
  // plot and fit each graph
  for (int i=0; i<num_hist; i++)
  {
    // set momentum/theta interval for this missing mass projection
    p_start_val[i] = x_min + i*dp;
    int bin1 = hist2d->GetXaxis()->FindBin(p_start_val[i]);
    int bin2 = hist2d->GetXaxis()->FindBin(p_start_val[i]+dp) - 1;
    // make projection for x interval
    can->cd(i+1);
    TH1D * proj = hist2d->ProjectionY("",bin1,bin2,"d");

    // create name of missing mass histogram for current momentum/theta interval
    std:ostringstream sObj1, sObj2;
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

    // fit histogram to Gaussian (signal) + Gaussian (background)
    TF1 * cfit = new TF1("cfit",mmiss_signal_gauss,Mdisp_lo,Mdisp_hi,6);
    cfit->SetLineColor(kMagenta);
    cfit->SetParameters(10000,0.94,0.05,1200,1.2,0.1);
    cfit->SetParLimits(0,0,1000000);
    cfit->SetParLimits(1,0.92,1.01);
    cfit->SetParLimits(2,0.02,0.1);
    cfit->SetParLimits(3,0,1000000);
    cfit->SetParLimits(4,0.95,1.8);
    cfit->SetParLimits(5,0.05,0.3);
    proj->Fit("cfit","QN");
    // separate signal/background
    // fit background
    Double_t par[6];
    cfit->GetParameters(par);
    TF1 * bkfit = new TF1("backfit",signal,Mdisp_lo,Mdisp_hi,3);
    bkfit->SetParameters(&par[3]);
    // fit signal
    TF1 * sgfit = new TF1("sgfit",signal,Mdisp_lo,Mdisp_hi,3);
    sgfit->SetParameters(par);
    // background subtraction!
    TH1D * proj_sub = (TH1D*)proj->Clone();
    proj_sub->Add(bkfit,-1);
    // draw stuff
    if (!subtract_bk) // no background subtraction - draw projection and all fits
    {
      proj->Draw();
      sgfit->SetLineColor(kRed);  sgfit->Draw("same");
      cfit->SetLineColor(kGreen);  cfit->Draw("same");
      bkfit->SetLineColor(kBlue);  bkfit->Draw("same");
    }
    else // background subtraction - draw subtracted projection and signal fit
    {
      proj_sub->Draw();
      sgfit->SetLineColor(kRed);  sgfit->Draw("same");
    }
    // find background-subtracted signal
    S[i] = proj_sub->Integral(proj->GetXaxis()->FindBin(Mlow),proj->GetXaxis()->FindBin(Mhigh));
    //std::cout << "integral = " << proj->Integral(Mlow,Mhigh) << '\n';
    //std::cout << "integral and error = " << proj->IntegralAndError(Mlow,Mhigh) << '\n';
  }
  return S;
}
