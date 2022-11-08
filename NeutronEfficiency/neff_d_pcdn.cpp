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
#include <TGraph.h>
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


// efficiency constants
int neff_pbins = 20;
int neff_tbins = 20;
int pgrid_x = ceil(sqrt(neff_pbins));
int pgrid_y = ceil((double)neff_pbins/(double)pgrid_x);
int tgrid_x = ceil(sqrt(neff_tbins));
int tgrid_y = ceil((double)neff_tbins/(double)tgrid_x);
double Mlow = 0.85;
double Mhigh = 1.05;
double theta_lo = 40;  // 40 for CD p
double theta_hi = 120;
double p_lo = 0.25;
double p_hi = 0.8;
double Mdisp_lo = 0.7;
double Mdisp_hi = 1.3;



void printProgress(double percentage);
double get_pin_mmiss(TVector3 p_b, TVector3 p_e, TVector3 ppi);
Double_t lorentzian(Double_t *x, Double_t *par);
Double_t poly(Double_t *x, Double_t *par);
Double_t signal(Double_t *x, Double_t *par); 
Double_t mmiss_signal_gauss(Double_t *x, Double_t *par);
Double_t mmiss_signal_poly(Double_t *x, Double_t *par);
Double_t mmiss_signal_lorentz(Double_t *x, Double_t *par);
double * hist_projections(TCanvas * can, TH2D * hist2d, int num_hist, char v);





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

  // output file names
  TFile * outFile = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];

  char * basename = argv[3];
  basename[strlen(basename)-5] = '\0';
  string theta_name(string(basename) + "_theta.txt");
  theta_name.c_str();
  string p_name(string(basename) + "_p.txt");
  p_name.c_str();

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
  //Histo: protons
  /////////////////////////////////////
  TH1D * h_pvertex = new TH1D("p vertex","Proton vertex - electron vertex;v_{p} - v_{e} (cm);Counts",50,-4,4);
  hist_list_1.push_back(h_pvertex);
  TH2D * h_dbetap = new TH2D("dbeta_p","#Delta #beta vs Momentum;Momentum (GeV/c);#beta_{meas} - p/sqrt(p^{2}+m^{2})",100,0,3,100,-0.2,0.2);
  hist_list_2.push_back(h_dbetap);
  TH1D * h_pchipid = new TH1D("p chipid","Proton #Chi^{2} PID",100,-3,3);
  hist_list_1.push_back(h_pchipid);

  
  /////////////////////////////////////
  //Histo: all neutrals
  /////////////////////////////////////
  TH2D * h_mmiss_xb = new TH2D("mmiss_xb_cand","Missing Mass vs xB;xB;Missing Mass",100,0,2,100,0,2);
  hist_list_2.push_back(h_mmiss_xb);
  TH1D * h_mmiss_withn = new TH1D("mmiss_withn","Missing Mass d(e,e'pn);Missing Mass (GeV/c^{2})",100,0.,2.);
  hist_list_1.push_back(h_mmiss_withn);
  TH1D * h_tof = new TH1D("tof_all","TOF (ns)",50,-10,40);
  hist_list_1.push_back(h_tof);


  
  /////////////////////////////////////
  //Histo: all neutrons
  /////////////////////////////////////
  TH1D * h_nsize = new TH1D("nsize","Number of reconstructed neutrons",10,0,10);
  hist_list_1.push_back(h_nsize);
  TH2D * h_nangles = new TH2D("n angles","Neutron Angular Distribution;Phi;Theta",48,-180,180,140,30,110);
  hist_list_2.push_back(h_nangles);





  /////////////////////////////////////
  // Histos: missing mass vs momentum
  /////////////////////////////////////
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","Missing Mass vs Pmiss d(e,e'p)n;Pmiss (GeV/c);Missing Mass (GeV/c^{2})",100,0,1,100,0.2,1.5);
  hist_list_2.push_back(h_mmiss_pmiss);
  TH2D * h_mmiss_pmissCAND = new TH2D("mmiss_pmissCAND","Missing Mass vs Pmiss d(e,e'p)n;Pmiss (GeV/c);Missing Mass (GeV/c^{2})",100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_pmissCAND);
  TH2D * h_mmiss_pmissDET = new TH2D("mmiss_pmissDET","Missing Mass vs Pmiss d(e,e'p)n;Pmiss (GeV/c);Missing Mass (GeV/c^{2})",100,p_lo,p_hi,30,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_pmissDET);
  TH1D * h_thetamiss_before = new TH1D("thetamiss_before","Missing Momentum Polar Angle;#theta_{miss};Counts",180,0,180);
  hist_list_1.push_back(h_thetamiss_before);
  


  /////////////////////////////////////
  // Histos: missing mass vs momentum
  /////////////////////////////////////
  TH2D * h_mmiss_theta = new TH2D("mmiss_theta","Missing Momentum vs #theta;Theta;Missing Mass (GeV/c^{2})",140,0,140,50,0.,2.);
  hist_list_2.push_back(h_mmiss_theta);
  TH2D * h_mmiss_thetaCAND = new TH2D("mmiss_thetaCAND","Missing Momentum vs #theta;Theta;Missing Mass (GeV/c^{2})",50,theta_lo,theta_hi,50,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_thetaCAND);
  TH2D * h_mmiss_thetaDET = new TH2D("mmiss_thetaDET","Missing Momentum vs #theta;Theta;Missing Mass (GeV/c^{2})",50,theta_lo,theta_hi,50,Mdisp_lo,Mdisp_hi);
  hist_list_2.push_back(h_mmiss_thetaDET);





  /////////////////////////////////////
  //Histo: compare with pmiss
  /////////////////////////////////////
  TH2D * h_pmiss_pn = new TH2D("pmiss_pn","Missing Momentum vs Neutron Momentum;p_{neutron};p_{miss}",100,0,1.25,100,0,1.25);
  hist_list_2.push_back(h_pmiss_pn);
  TH1D * h_cos0 = new TH1D("cos0","cos #theta_{pmiss,pneutron}",50,-1.05,1.05);
  hist_list_1.push_back(h_cos0);
  TH1D * h_dphi = new TH1D("dphi","#Delta #phi = #phi_{pmiss} - #phi_{n}",48,-180,180);
  hist_list_1.push_back(h_dphi);
  TH2D * h_pmiss_pn_cut = new TH2D("pmiss_pn_cut","Missing Momentum vs Neutron Momentum;p_{neutron};p_{miss}",100,0.25,1.25,100,0.25,1.25);
  hist_list_2.push_back(h_pmiss_pn_cut);
  TH2D * h_pmiss_theta = new TH2D("pmiss_theta","Missing Momentum vs #theta;Theta;Missing Momentum (GeV/c)",100,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_theta);
  TH2D * h_dtheta_dphi = new TH2D("dtheta_dphi","#Delta #theta vs #Delta #phi;#phi_{miss}-#phi_{n};#theta_{miss}-#theta_{n}",100,-180,180,100,-180,180);
  hist_list_2.push_back(h_dtheta_dphi);
  TH2D * h_theta_edep = new TH2D("theta_edep","Energy Deposition vs #theta;Theta;Energy Deposition",100,0,180,100,0,30);
  hist_list_2.push_back(h_theta_edep);
 
  TH2D * h_pmiss_pn_sig = new TH2D("pmiss_pn_sig","Missing Momentum vs Neutron Momentum (Signal);p_{neutron};p_{miss}",100,0.25,1.25,100,0,1.25);
  hist_list_2.push_back(h_pmiss_pn_sig);
  TH2D * h_pmiss_pn_bkg = new TH2D("pmiss_pn_bkg","Missing Momentum vs Neutron Momentum (Background);p_{neutron};p_{miss}",100,0.25,1.25,100,0,1.25);
  hist_list_2.push_back(h_pmiss_pn_bkg);
  TH2D * h_dp_p_sig = new TH2D("dp_p_sig","p_{miss} - p_{n} vs p_{miss} (Signal);p_{miss};p_{miss} - p_{n}",100,0.25,1.2,100,-1,1);
  hist_list_2.push_back(h_dp_p_sig);
  TH2D * h_dp_p_bkg = new TH2D("dp_p_bkg","p_{miss} - p_{n} vs p_{miss} (Background);p_{miss};p_{miss} - p_{n}",100,0.25,1.2,100,-1,1);
  hist_list_2.push_back(h_dp_p_bkg);


  /////////////////////////////////////
  //Histo: cutting out background
  /////////////////////////////////////
  TH1D * h_theta_ppmiss = new TH1D("theta_ppmiss","Angle between p_{p} and p_{miss}",180,0,180);
  hist_list_1.push_back(h_theta_ppmiss);
  TH1D * h_theta_np = new TH1D("theta_np","Angle between neutron and proton",45,0,180);
  hist_list_1.push_back(h_theta_np);
  TH2D * h_pmiss_pn_cutbkg = new TH2D("pmiss_pn_cutbkg","p_{miss} vs p_{n} (background cut);p_{neutron};p_{miss}",100,0.25,1.25,100,0.25,1.25); // CTOF p only
  hist_list_2.push_back(h_pmiss_pn_cutbkg);
  TH2D * h_pmiss_theta_cutbkg = new TH2D("pmiss_theta_cutbkg","Missing Momentum vs p_{miss} Polar Angle;#theta_{miss} (degrees);p_{miss} (GeV/c)",100,30,150,100,0,1.5);
  hist_list_2.push_back(h_pmiss_theta_cutbkg);
 
 


  /////////////////////////////////////
  //Histo: neutron efficiency
  /////////////////////////////////////
  TH1D * h_neff_pmiss_numer_ssb = new TH1D("neff_pm_numer_ssb","Neutrons;p_{miss} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_numer_ssb);
  TH1D * h_neff_pmiss_denom_ssb = new TH1D("neff_pm_denom_ssb","Neutron Candidates;p_{miss} (GeV/c);Counts",neff_pbins,0.25,1);
  hist_list_1.push_back(h_neff_pmiss_denom_ssb);


  /////////////////////////////////////
  //Histo: neutron efficiency by theta
  /////////////////////////////////////
  TH1D * h_neff_thetamiss_denom_ssb = new TH1D("neff_tm_denom_ssb","Neutron Candidates;#theta_{miss} (deg);Counts",neff_tbins,theta_lo,theta_hi);
  hist_list_1.push_back(h_neff_thetamiss_denom_ssb);
  TH1D * h_neff_thetamiss_numer_ssb = new TH1D("neff_tm_numer_ssb","Detected Neutrons;#theta_{miss} (deg);Counts",neff_tbins,theta_lo,theta_hi);
  hist_list_1.push_back(h_neff_thetamiss_numer_ssb);


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
    auto elec=c12->getByID(11);
    auto prot=c12->getByID(2212);
    auto phot=c12->getByID(22);
    auto neut=c12->getByID(2112);
    auto allParticles=c12->getDetParticles();
    double weight = 1;
    if(isMC){weight=c12->mcevent()->getWeight();}

    double ts = c12->event()->getStartTime();



    // GENERAL EVENT SELECTION
    if(!myCut.electroncut(c12)) {continue;}
    if(prot.size()!=1) {continue;}
    if(elec.size()!=1) {continue;}
    //if(phot.size()>0) {continue;}
    //if(neut.size()<0 || neut.size()>1) {continue;}
    

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

    // was originally after p cuts, before mmiss
    bool is_CD = (prot[0]->getRegion()==CD);
    if(!is_CD) {continue;}

    // protons histos
    h_pvertex->Fill(vzp-vze,weight);
    h_dbetap->Fill(pp.Mag(),dbeta,weight);
    h_pchipid->Fill(chipid,weight);

    // proton cuts
    if ((vzp-vze)<-4. || (vzp-vze)>4.) {continue;}
    if (dbeta<-0.05 || dbeta>0.05) {continue;}
    if (chipid<-3 || chipid>3) {continue;}
    if (pp.Mag() < 0.2) {continue;}
    if (pp.Mag() > 1.0) {continue;}

    // Missing mass and energy despotiion
    double mmiss = get_mmiss(p_b,p_e,pp);
    TVector3 pmiss = p_b - p_e - pp;
    double thetamiss = pmiss.Theta()*180./M_PI;

    // cut on theta component of pmiss
    h_thetamiss_before->Fill(thetamiss,weight);
    if (thetamiss<40.) {continue;}
    if (thetamiss>140.) {continue;}

    // allow only missing momentum in the correct range
    if (pmiss.Mag() < 0.243) {continue;} // beta = 0.2, p = 0.1918
    if (pmiss.Mag() > 1.25) {continue;} // beta = 0.8, p = 1.2528

    // xb cut
    h_mmiss_xb->Fill(xB,mmiss,weight);
    if (xB<0.6) {continue;}


    h_theta_ppmiss->Fill(pmiss.Angle(pp)*180./M_PI,weight);
    //if ((pmiss.Angle(pp)*180./M_PI)<40.) {continue;}  // CD ONLY (I THINK)
    
    h_mmiss_pmiss->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_pmissCAND->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_theta->Fill(thetamiss,mmiss,weight);
    h_mmiss_thetaCAND->Fill(thetamiss,mmiss,weight);


    if (mmiss>Mlow && mmiss<Mhigh) { h_neff_pmiss_denom_ssb->Fill(pmiss.Mag(),weight);}
    if (mmiss>Mlow && mmiss<Mhigh) { h_neff_thetamiss_denom_ssb->Fill(thetamiss,weight);}


    // REQUIRE A NEUTRON HERE


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
        if (beta_n<0.25) {continue;}  // should be 0.2, not 0.25, right?
        
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




    bool is_CND1 = neut[pick]->sci(CND1)->getDetector()==3;
    bool is_CND2 = neut[pick]->sci(CND2)->getDetector()==3;
    bool is_CND3 = neut[pick]->sci(CND3)->getDetector()==3;
    double n_theta = neut[pick]->getTheta()*180./M_PI;
    double beta_n = neut[pick]->par()->getBeta();
    double n_phi = neut[pick]->getPhi()*180./M_PI;

    double tof_n = 0;
    double path = 0;
    double edep = 0;
    if (is_CND1)
    {
      tof_n = neut[pick]->sci(CND1)->getTime() - ts;
      path = neut[pick]->sci(CND1)->getPath();
      edep = neut[pick]->sci(CND1)->getEnergy();
    }
    else if (is_CND2)
    {
      tof_n = neut[pick]->sci(CND2)->getTime() - ts;
      path = neut[pick]->sci(CND2)->getPath();
      edep = neut[pick]->sci(CND2)->getEnergy();
    }
    else if (is_CND3)
    {
      tof_n = neut[pick]->sci(CND3)->getTime() - ts;
      path = neut[pick]->sci(CND3)->getPath();
      edep = neut[pick]->sci(CND3)->getEnergy();
    }
    if (tof_n<0) {continue;} //{std::cout << "Look out - negative TOF!/n";}




    TVector3 pn;
    pn.SetMagThetaPhi(neut[pick]->getP(),neut[pick]->getTheta(),neut[pick]->getPhi());
    TVector3 vecX( neut[pick]->par()->getPx(), neut[pick]->par()->getPy(), neut[pick]->par()->getPz() );
    double cos0 = pmiss.Dot(vecX) / (pmiss.Mag() * vecX.Mag() );
    double dphi = pmiss.Phi()*180./M_PI - neut[pick]->getPhi()*180./M_PI;


    // fill histos with good neutrons
    h_tof->Fill(tof_n,weight);
    h_nangles->Fill(n_phi,n_theta,weight);
    h_dphi->Fill(dphi,weight);
    h_pmiss_pn->Fill(pn.Mag(),pmiss.Mag(),weight);
    h_dtheta_dphi->Fill(dphi,pmiss.Theta()*180./M_PI - neut[pick]->getTheta()*180./M_PI,weight);
    h_theta_edep->Fill(neut[pick]->getTheta()*180./M_PI,edep,weight);
    //h_mmiss_beta_after->Fill(beta_n,mmiss,weight);
    //h_beta_pmiss->Fill(pmiss.Mag(),beta_n,weight);
    h_cos0->Fill(cos0,weight);



    if (std::abs(dphi)>20.) {continue;}
    if (cos0 < 0.9) {continue;}

    // after requiring pn along pmiss
    h_pmiss_pn_cut->Fill(pn.Mag(),pmiss.Mag(),weight);

    if (mmiss>1.15) {h_pmiss_pn_bkg->Fill(pn.Mag(),pmiss.Mag(),weight);}
    if (mmiss<1.05) {h_pmiss_pn_sig->Fill(pn.Mag(),pmiss.Mag(),weight);}
    if (mmiss>1.15) {h_dp_p_bkg->Fill(pmiss.Mag(),pmiss.Mag()-pn.Mag(),weight);}
    if (mmiss<1.05) {h_dp_p_sig->Fill(pmiss.Mag(),pmiss.Mag()-pn.Mag(),weight);}

    if (std::abs(pmiss.Mag()-pn.Mag())>0.2) {continue;}

    h_mmiss_withn->Fill(mmiss,weight);
    h_pmiss_theta->Fill(neut[pick]->getTheta()*180./M_PI,pmiss.Mag(),weight);

    // cut out fake neutron background
    double theta_np = pn.Angle(pp)*180./M_PI;
    h_theta_np->Fill(theta_np,weight);

    if (is_CD && theta_np<40.) {continue;}
    h_pmiss_pn_cutbkg->Fill(pn.Mag(),pmiss.Mag(),weight); // CTOF p only
    h_pmiss_theta_cutbkg->Fill(n_theta,pmiss.Mag(),weight);

    

    h_mmiss_pmissDET->Fill(pmiss.Mag(),mmiss,weight);
    h_mmiss_thetaDET->Fill(thetamiss,mmiss,weight);



    if (mmiss>Mlow && mmiss<Mhigh)  { h_neff_pmiss_numer_ssb->Fill(pmiss.Mag(),weight); }
    if (mmiss>Mlow && mmiss<Mhigh) { h_neff_thetamiss_numer_ssb->Fill(thetamiss,weight); }

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
  
  /////////////////////
  // PROTONS
  /////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Deuterium run 015053");
  text.DrawLatex(0.2,0.7,"d(e,e'p)n");
  text.DrawLatex(0.2,0.6,"allow only PID=0,11,22,2212,2112");
  text.DrawLatex(0.2,0.5,"exactly 1p, 1e");
  text.DrawLatex(0.2,0.4,"proton in FD");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pvertex->Draw();
  myCanvas->cd(2);
  h_dbetap->Draw("colz");
  myCanvas->cd(3);
  h_pchipid->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 

  
  
  /////////////////////
  // NEUTRALS
  /////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Deuterium run 015053");
  text.DrawLatex(0.2,0.7,"d(e,e'p)n");
  text.DrawLatex(0.2,0.6,"Require proton in FTOF");
  text.DrawLatex(0.2,0.5,"-4 cm < v_{p}-v_{e} < 2 cm");
  text.DrawLatex(0.2,0.4,"-0.02 < #Delta #beta < 0.02");
  text.DrawLatex(0.2,0.3,"p_{p} > 0.3 GeV/c");
  text.DrawLatex(0.2,0.2,"-3 < proton #chi^{2}_{PID} <3");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetamiss_before->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_xb->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Deuterium run 015053");
  text.DrawLatex(0.2,0.7,"d(e,e'pn)");
  text.DrawLatex(0.2,0.6,"proton cuts");
  text.DrawLatex(0.2,0.5,"40 deg < #theta_{miss} < 140 deg");
  text.DrawLatex(0.2,0.4,"0.24 GeV/c < p_{miss} < 1.25 GeV/c");
  myText->Print(fileName,"pdf");
  myText->Clear();


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
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_pmissCAND->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
  

  // mmiss in pmiss bins - landau + gaussian
  myCanvas->Divide(pgrid_x,pgrid_y);
  double * SB_pCAND_ = hist_projections(myCanvas, h_mmiss_pmissCAND, neff_pbins, 'p');
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
  

  // mmiss by theta bins - landau + gaussian
  myCanvas->Divide(tgrid_x,tgrid_y);
  double * SB_tCAND = hist_projections(myCanvas, h_mmiss_thetaCAND, neff_tbins, 't');
  std::cout << "Candidates - as a function of theta\n";
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // START REQUIRING
  // NEUTRONS
  /////////////////////


  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Deuterium run 015053");
  text.DrawLatex(0.2,0.7,"d(e,e'pn)");
  text.DrawLatex(0.2,0.6,"proton cuts");
  text.DrawLatex(0.2,0.5,"p_{miss} cuts");
  text.DrawLatex(0.2,0.4,"Require neutron in CND");
  text.DrawLatex(0.2,0.3,"exclude #theta_{n}=0, #phi_{n}=0");
  text.DrawLatex(0.2,0.2,"40 deg < #theta_{n} < 140 deg");
  text.DrawLatex(0.2,0.1,"0.25 < #beta_{n} < 0.8");
  myText->Print(fileName,"pdf");
  myText->Clear();
 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_nsize->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_pn->Draw("colz");
  TF1 * line = new TF1("line","x",0,2);
  line->Draw("same");
  myCanvas->cd(2);
  h_tof->Draw();
  myCanvas->cd(3);
  h_nangles->Draw("colz");
  myCanvas->cd(4);
  h_dtheta_dphi->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear(); 



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();  



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dphi->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();




  // neutrons in direction of pmiss
  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  text.DrawLatex(0.2,0.8,"Deuterium run 015053");
  text.DrawLatex(0.2,0.7,"d(e,e'pn)");
  text.DrawLatex(0.2,0.6,"proton cuts");
  text.DrawLatex(0.2,0.5,"p_{miss} cuts");
  text.DrawLatex(0.2,0.4,"exclude #theta_{n}=0, #phi_{n}=0");
  text.DrawLatex(0.2,0.3,"40 deg < #theta_{n} < 140 deg");
  text.DrawLatex(0.2,0.2,"0.25 < #beta_{n} < 0.8");
  text.DrawLatex(0.2,0.1,"cos #theta_{pmiss,pn} > 0.9, -20 < #Delta #phi < 20");
  myText->Print(fileName,"pdf");
  myText->Clear();



  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_pn_cut->Draw("colz");
  line->Draw("same");
  myCanvas->cd(2);
  h_mmiss_withn->Draw("colz");
  myCanvas->cd(4);
  h_pmiss_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_pn_sig->Draw("colz");
  line->Draw("same");
  myCanvas->cd(2);
  h_pmiss_pn_bkg->Draw("colz");
  line->Draw("same");
  myCanvas->cd(3);
  h_dp_p_sig->Draw("colz");
  myCanvas->cd(4);
  h_dp_p_bkg->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_theta_np->Draw();
  myCanvas->cd(2);
  h_pmiss_pn_cutbkg->Draw("colz"); // CTOF p only
  myCanvas->cd(3);
  h_pmiss_theta_cutbkg->Draw("colz");
  h_pmiss_theta_cutbkg->SetStats(0);
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


  // mmiss in pmiss bins - landau + gaussian
  myCanvas->Divide(pgrid_x,pgrid_y);
  double * SB_pDET = hist_projections(myCanvas, h_mmiss_pmissDET, neff_pbins, 'p');
  std::cout << "Detected - as a function of pmiss\n";
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
 
  // mmiss in theta bins - landau + gaussian
  myCanvas->Divide(tgrid_x,tgrid_y);
  double * SB_tDET = hist_projections(myCanvas, h_mmiss_thetaDET, neff_tbins, 't');
  std::cout << "Detected - as a function of theta\n";
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  /////////////////////
  // EFFICIENCY RESULTS
  /////////////////////

  myText->cd();
  text.DrawLatex(0.2,0.6,"EFFICIENCY RESULTS");
  myText->Print(fileName,"pdf");
  myText->Clear();


  // efficiency vs pmiss using S/(S+B) filling
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_neff_pmiss_numer_ssb->Draw();
  myCanvas->cd(2);
  h_neff_pmiss_denom_ssb->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_pmiss = (TH1D*)h_neff_pmiss_numer_ssb->Clone();
  h_neff_pmiss->Divide(h_neff_pmiss_denom_ssb);
  h_neff_pmiss->Draw();
  h_neff_pmiss->GetYaxis()->SetTitle("efficiency");
  h_neff_pmiss->GetYaxis()->SetRangeUser(0.,0.16);
  // print output
  ofstream outp(p_name);
  for (int i=0; i<h_neff_pmiss->GetNbinsX(); i++) {
    outp << h_neff_pmiss->GetXaxis()->GetBinCenter(i) << ' ';
    outp << h_neff_pmiss->GetBinContent(i) << ' ';
    outp << h_neff_pmiss->GetBinError(i) << '\n';
  }
  outp.close();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





  // efficiency vs theta using S/(S+B) filling
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_neff_thetamiss_numer_ssb->Draw();
  myCanvas->cd(2);
  h_neff_thetamiss_denom_ssb->Draw();
  myCanvas->cd(3);
  TH1D * h_neff_thetamiss_ssb = (TH1D*)h_neff_thetamiss_numer_ssb->Clone();
  h_neff_thetamiss_ssb->Divide(h_neff_thetamiss_denom_ssb);
  h_neff_thetamiss_ssb->Draw();
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




  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  outFile->Close();
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
  double y_min = hist2d->GetYaxis()->GetXmin();
  double y_max = hist2d->GetYaxis()->GetXmin();
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


