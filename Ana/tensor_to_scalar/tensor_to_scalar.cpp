#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <TDatabasePDG.h>
#include "clas12reader.h"
#include "HipoChain.h"
#include "clas12ana.h"
#include "eventcut/functions.h"
#include "neutron-veto/veto_functions.h"

 
using namespace std;
using namespace clas12;
using namespace TMVA;

//const double c = 29.9792458;
const double mp = 0.938272;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp)
{
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void printProgress(double percentage);

void Usage()
{
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <neff alt?> <mlp alt?> <peff alt?> <path/to/input.hipo> \n";
}


int main(int argc, char ** argv)
{

  if(argc < 6)
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

  // parameters for error analysis
  bool neff_alt = false;
  if(atoi(argv[5]) == 1){neff_alt=true;}
  int mlp_alt = atoi(argv[6]); // 1 for low cutoff, 2 for high cutoff
  bool peff_alt = false;
  if(atoi(argv[7]) == 1){peff_alt=true;}

  // create instance of clas12ana class
  clas12ana clasAna;
  clasAna.printParams();

  clas12root::HipoChain chain;
  for(int k = 8; k < argc; k++){
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
  TH2D * h_pangles = new TH2D("pangles","Proton angles;#phi;#theta",90,-180,180,180,0,180);
  hist_list_2.push_back(h_pangles);
  TH2D * h_p_theta = new TH2D("p_theta","Lead proton p vs #theta;#theta;p (GeV/c)",90,0,180,100,0,2.5);
  hist_list_2.push_back(h_p_theta);
  TH2D * h_pq = new TH2D("pq","#theta_{p,q} vs p/q;p/q;#theta_{pq}",100,0,2,100,0,100);
  hist_list_2.push_back(h_pq);
  TH1D * h_mmiss_all = new TH1D("mmiss_all","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  hist_list_1.push_back(h_mmiss_all);
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass;M_{miss} (GeV/c^{2});Counts",100,0,2);
  hist_list_1.push_back(h_mmiss);
  TH1D * h_pmisstheta = new TH1D("pmisstheta","Theta Component of p_{miss};#theta_{pmiss};Counts",180,0,180);
  hist_list_1.push_back(h_pmisstheta);

  TH2D * h_q2_xb = new TH2D("q2_xb","Q^{2} vs x_{B};x_{B};Q^{2} (GeV^{2})",100,0,2,100,0,3);
  hist_list_2.push_back(h_q2_xb);
  TH1D * h_q2 = new TH1D("q2","Q^{2};Counts;Q^{2} (GeV^{2})",100,0,2);
  hist_list_1.push_back(h_q2);

  TH2D * h_lr_angles = new TH2D("lr_angles","Lead and Recoil Polar Angles;Lead theta;Recoil theta",180,0,180,180,0,180);
  hist_list_2.push_back(h_lr_angles);



  /////////////////////////////////////
  //Histos: recoil proton
  /////////////////////////////////////
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
  hist_list_1.push_back(h_psize);
  TH1D * h_pcos0 = new TH1D("pcos0","Cosine of #theta_{pmiss,pp}",100,-1.1,1.1);
  hist_list_1.push_back(h_pcos0);
  TH1D * h_prec_p = new TH1D("prec_p","Recoil Proton Momentum",100,0,2);
  hist_list_1.push_back(h_prec_p);
  TH2D * h_prec_plead = new TH2D("prec_plead","Lead momentum vs Recoil proton momentum;p_{p} (GeV/c);p_{L} (GeV/c)",50,0,1.5,50,0,2);
  hist_list_2.push_back(h_prec_plead);

  TH2D * h_prec_ptheta = new TH2D("prec_ptheta","Recoil Proton Theta vs Momentum;Momentum (GeV/c);#theta (degrees)",50,0.,1.5,50,0,180);
  hist_list_2.push_back(h_prec_ptheta);
  TH2D * h_prec_angles = new TH2D("prec_angles","Recoil Proton Angular Distribution;phi (deg);theta (deg)",48,-180,180,45,0,180);
  hist_list_2.push_back(h_prec_angles);

  TH1D * h_prec_plead_angle = new TH1D("prec_plead_angle","Angle between lead proton and recoil proton;#theta_{lead,recoil};Counts",45,0,180);
  hist_list_1.push_back(h_prec_plead_angle);
  TH2D * h_lpangle_pmiss = new TH2D("lpangle_pmiss","Angle between lead proton and recoil proton vs p_{miss};p_{miss} (GeV/c);#theta_{lead,recoil}",50,0.2,1,90,0,180);
  hist_list_2.push_back(h_lpangle_pmiss);


  /////////////////////////////////////
  //Histos: recoil neutron
  /////////////////////////////////////
  TH1D * h_nsize = new TH1D("nsize","Number of Neutrons in Event",10,0,10);
  hist_list_1.push_back(h_nsize);
  TH1D * h_ncos0 = new TH1D("ncos0","Cosine of #theta_{pmiss,pn}",100,-1.1,1.1);
  hist_list_1.push_back(h_ncos0);
  TH1D * h_nrec_p = new TH1D("nrec_p","Recoil Neutron Momentum",100,0,2);
  hist_list_1.push_back(h_nrec_p);
  TH2D * h_nrec_plead = new TH2D("nrec_plead","Lead momentum vs Recoil neutron momentum;p_{n} (GeV/c);p_{L} (GeV/c)",50,0,1.5,50,0,2);
  hist_list_2.push_back(h_nrec_plead);
  TH2D * h_nrec_ptheta = new TH2D("nrec_ptheta","Recoil Neutron Theta vs Momentum;Momentum (GeV/c);#theta (degrees)",30,0.,1.5,20,0,180);
  hist_list_2.push_back(h_nrec_ptheta);
  TH2D * h_nrec_angles = new TH2D("nrec_angles","Recoil Neutron Angular Distribution;phi (deg);theta (deg)",48,-180,180,45,0,180);
  hist_list_2.push_back(h_nrec_angles);
  TH2D * h_good_nrec_angles = new TH2D("good_nrec_angles","Recoil Neutron Angular Distribution;phi (deg);theta (deg)",48,-180,180,45,0,180);
  hist_list_2.push_back(h_good_nrec_angles);

  TH2D * h_pn_pmiss = new TH2D("pn_pmiss","Neutron Momentum vs Missing Momentum;p_{miss} (GeV/c);Neutron Momentum (GeV/c)",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_pn_pmiss);
  TH2D * h_nptheta = new TH2D("nptheta","Neutron #theta vs Momentum;#theta_{n};Momentum (GeV/c)",180,0,180,50,0,1.5);
  hist_list_2.push_back(h_nptheta);
  TH1D * h_nrec_plead_angle = new TH1D("nrec_plead_angle","Angle between lead proton and recoil neutron;#theta_{lead,recoil};Counts",45,0,180);
  hist_list_1.push_back(h_nrec_plead_angle);
  TH2D * h_lnangle_pmiss = new TH2D("lnangle_pmiss","Angle between lead proton and recoil neutron vs p_{miss};p_{miss} (GeV/c);#theta_{lead,recoil}",50,0.2,1,90,0,180);
  hist_list_2.push_back(h_lnangle_pmiss);




  /////////////////////////////////////
  //Histos: yields
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
  int pmiss_bins = 7;
  TH1D * h_pmiss_p = new TH1D("pmiss_p","Missing Momentum (e,e'p_{SRC});p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_p);
  TH1D * h_pmiss_p_wrec = new TH1D("pmiss_p_wrec","Missing Momentum (e,e'p_{SRC}N_{rec});p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_p_wrec);
  TH1D * h_pmiss_pp_uncorr = new TH1D("pmiss_pp_uncorr","Missing Momentum (e,e'p_{SRC}p_{rec});p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_pp_uncorr);
  TH1D * h_pmiss_pp_corr = new TH1D("pmiss_pp_corr","Missing Momentum (e,e'p_{SRC}p_{rec}) (efficiency corrected);p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_pp_corr);
  TH1D * h_pmiss_pn_uncorr = new TH1D("pmiss_pn_uncorr","Missing Momentum (e,e'p_{SRC}n_{rec});p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_pn_uncorr);
  TH1D * h_pmiss_pn_corr = new TH1D("pmiss_pn_corr","Missing Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{miss} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pmiss_pn_corr);

  TH1D * h_pn = new TH1D("pn","Momentum (e,e'p_{SRC}n_{rec});p_{n} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pn);
  TH1D * h_pn_corr = new TH1D("pn_corr","Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{n} (GeV/c);Counts",pmiss_bins,0.3,1);
  hist_list_1.push_back(h_pn_corr);



  /////////////////////////////////////
  //Histos: ML features
  /////////////////////////////////////
  TH1D * h_mvaValue_MLP = new TH1D("mvaValue_MLP","MVA Value (MLP);MVA Output Value;Counts",100,0,1);
    hist_list_1.push_back(h_mvaValue_MLP);
  TH1D * h_mvaValue_BDT = new TH1D("mvaValue_BDT","MVA Value (BDT);MVA Output Value;Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_mvaValue_BDT);

  TH1D * h_energy_s = new TH1D("f_energy_s","Neutron Energy",40,0,300);
    hist_list_1.push_back(h_energy_s);
  TH1D * h_layermult_s = new TH1D("f_layermult_s","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_s);
  TH1D * h_size_s = new TH1D("f_size_s","Cluster Size",4,1,4);
    hist_list_1.push_back(h_size_s);
  TH1D * h_cnd_hits_s = new TH1D("f_cnd_hits_s","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_s);
  TH1D * h_cnd_energy_s = new TH1D("f_cnd_energy_s","Nearby CND Energy",50,0,400);
    hist_list_1.push_back(h_cnd_energy_s);
  TH1D * h_ctof_energy_s = new TH1D("f_ctof_energy_s","Nearby CTOF Energy",50,0,200);
    hist_list_1.push_back(h_ctof_energy_s);
  TH1D * h_ctof_hits_s = new TH1D("f_ctof_hits_s","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_s);
  TH1D * h_anglediff_s = new TH1D("f_anglediff_s","CVT Angle Diff",50,0,200);
    hist_list_1.push_back(h_anglediff_s);

  TH1D * h_energy_b = new TH1D("f_energy_b","Neutron Energy",40,0,300);
    hist_list_1.push_back(h_energy_b);
  TH1D * h_layermult_b = new TH1D("f_layermult_b","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_b);
  TH1D * h_size_b = new TH1D("f_size_b","Cluster Size",4,1,4);
    hist_list_1.push_back(h_size_b);
  TH1D * h_cnd_hits_b = new TH1D("f_cnd_hits_b","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_b);
  TH1D * h_cnd_energy_b = new TH1D("f_cnd_energy_b","Nearby CND Energy",50,0,400);
    hist_list_1.push_back(h_cnd_energy_b);
  TH1D * h_ctof_energy_b = new TH1D("f_ctof_energy_b","Nearby CTOF Energy",50,0,200);
    hist_list_1.push_back(h_ctof_energy_b);
  TH1D * h_ctof_hits_b = new TH1D("f_ctof_hits_b","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_b);
  TH1D * h_anglediff_b = new TH1D("f_anglediff_b","CVT Angle Diff",50,0,200);
    hist_list_1.push_back(h_anglediff_b);

  TH1D * h_n_signal = new TH1D("n_signal","ML Signal and Background;Momentum (GeV/c);Counts",pmiss_bins,0.3,1);
    hist_list_1.push_back(h_n_signal);
  TH1D * h_n_background = new TH1D("n_background","ML Signal and Background;Momentum (GeV/c);Counts",pmiss_bins,0.3,1);
    hist_list_1.push_back(h_n_background);
  



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


  // TMVA stuff
  TMVA::Tools::Instance();
  TMVA::Reader * reader = new TMVA::Reader("!Color:!Silent");
  Float_t energy, cnd_energy, ctof_energy, angle_diff, momentum;
  Float_t layermult, size, cnd_hits, ctof_hits;
  reader->AddVariable("energy", &energy);
  reader->AddVariable("layermult", &layermult);
  reader->AddVariable("size", &size);
  reader->AddVariable("cnd_hits", &cnd_hits);
  reader->AddVariable("cnd_energy", &cnd_energy);
  reader->AddVariable("ctof_energy", &ctof_energy);
  reader->AddVariable("ctof_hits", &ctof_hits);
  reader->AddVariable("angle_diff", &angle_diff);

  reader->AddSpectator("momentum", &momentum);

  reader->BookMVA("MLP", "/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm/NeutronVeto/dataset_6gev_pCD/weights/TrainNeutronVeto_TMVA_MLP.weights.xml");

  int n_background = 0;
  int n_signal = 0;


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

    
    clasAna.Run(c12);
    auto elec = clasAna.getByPid(11);
    auto prot = clasAna.getByPid(2212);
    auto neut = clasAna.getByPid(2112);


    auto allParticles=c12->getDetParticles();
    double weight = 1;
    if(isMC){weight=c12->mcevent()->getWeight();}


    // initial cuts
    if (elec.size()!=1) {continue;}

    // ELECTRONS
    // electron kinematics
    TVector3 pe;
    TVector3 pb(0,0,Ebeam);
    pe.SetMagThetaPhi(elec[0]->getP(),elec[0]->getTheta(),elec[0]->getPhi());
    TVector3 q = pb - pe;
    double vze = elec[0]->par()->getVz();
    double nu = Ebeam - pe.Mag();
    double Q2 = q.Mag2() - (nu*nu);
    double xB = Q2 / (2*mp*nu);
    int nphe = elec[0]->che(HTCC)->getNphe();
    double vtz_e = elec[0]->par()->getVz();
   
    // electron histograms: quality cuts
    h_nphe->Fill(nphe,weight);
    h_vtz_e->Fill(vtz_e,weight);



    //// LEAD PROTON ////
    // LEAD PROTON - CTOF
    auto db = TDatabasePDG::Instance();
    double mC = 12.0*0.9315;
    double mD = 1.8756;
    TLorentzVector beam(0,0,Ebeam,Ebeam);
    TLorentzVector target(0,0,0,mD);
    TLorentzVector el(elec[0]->par()->getPx(), elec[0]->par()->getPy(), elec[0]->par()->getPz(), pe.Mag());
    clasAna.getLeadRecoilSRC(beam,target,el);
    auto lead = clasAna.getLeadSRC();
    auto recoil = clasAna.getRecoilSRC();

    if (lead.size()!=1) {continue;}

    TVector3 pL;
    pL.SetMagThetaPhi(lead[0]->getP(),lead[0]->getTheta(),lead[0]->getPhi());
    double vzlead = lead[0]->par()->getVz();
    double chi2pid = lead[0]->par()->getChi2Pid();
    double ltheta = lead[0]->getTheta()*180./M_PI;
    double lphi = lead[0]->getPhi()*180./M_PI;
    double theta_pq = pL.Angle(q) * 180./M_PI;
    double mmiss = get_mmiss(pb,pe,pL);
    double dbetap = lead[0]->par()->getBeta() - pL.Mag()/sqrt(pL.Mag2()+mp*mp);


    // lead histos
    h_vtzdiff_ep->Fill(vze-vzlead,weight);
    h_chi2pid->Fill(chi2pid,weight);
    h_dbetap->Fill(pL.Mag(),dbetap,weight);
    h_betap->Fill(pL.Mag(),lead[0]->par()->getBeta(),weight);


    // look at TOF distributions



    //// PMISS & LEAD SRC CUTS ////
   
    // theta cut, pL cut
    h_pangles->Fill(lphi,ltheta,weight);
    h_p_theta->Fill(ltheta,pL.Mag(),weight);
    if (ltheta<40 || ltheta>140) {continue;}
    if (pL.Mag()<1.0) {continue;}

    // xB cut
    h_q2_xb->Fill(xB,Q2,weight);
    if (xB<1.2) {continue;}
    

    // pmiss cut
    TVector3 pmiss = pL - q;
    double pm_theta = pmiss.Theta()*180./M_PI;
    h_pmisstheta->Fill(pm_theta,weight);
    h_pmiss_pL->Fill(pmiss.Mag(),weight);
    if (pmiss.Mag()<0.3) {continue;}
    if (pmiss.Mag()>1.0) {continue;}

    // Q2 cut
    h_q2->Fill(Q2,weight);
    if (Q2<1.5) {continue;}

    // p/q cut
    h_pq->Fill(pL.Mag()/q.Mag(),theta_pq,weight);
    if (pL.Mag()/q.Mag()>0.96) {continue;}

    // mmiss cut
    h_mmiss->Fill(mmiss,weight);
    if (mmiss>1.1) {continue;}


    h_pmiss_p->Fill(pmiss.Mag(),weight);
    h_pl_count->Fill(pmiss.Mag(),weight);
    h_psrc_count->Fill(pmiss.Mag(),weight);



//// RECOIL P ////
    // look for recoil proton in CTOF

    h_psize->Fill(prot.size(),weight);
    //if (recoil.size()>1) {continue;} // really, if there's one recoil, I should decide which is the best

    TVector3 p_recp;
    double p_cos0;

    if (recoil.size()==1)
    {
      p_recp.SetMagThetaPhi(recoil[0]->getP(),recoil[0]->getTheta(),recoil[0]->getPhi());
      // recoil p must be in CTOF
      bool is_CTOF = recoil[0]->sci(CTOF)->getDetector()==4;
      if (!is_CTOF) {continue;}

      // p and theta cut
      h_prec_ptheta->Fill(p_recp.Mag(),p_recp.Theta()*180./M_PI,weight);
      if (p_recp.Theta()*180./M_PI<45 || p_recp.Theta()*180./M_PI>120) {continue;}
      // limit to central detector acceptance
      if (p_recp.Mag()<0.3) {continue;}
      if (p_recp.Mag()>1.0) {continue;}

      /*// EPJA cuts
      if (p_recp.Theta()*180./M_PI<40 || p_recp.Theta()*180./M_PI>140) {continue;}
      if (p_recp.Mag()<0.3) {continue;}*/


      // CORRECT FOR CVT PROTON EFFICIENCY      
      // quadratic fit
      double peff = -1.14903684*p_recp.Mag2() + 2.78915506*p_recp.Mag() -0.67659258;
      // linear fit
      if (peff_alt==1) {double peff = 1.38158493*p_recp.Mag() -0.30256513;}



      // get momenta/angles of recoil protons
      h_prec_plead->Fill(p_recp.Mag(),pL.Mag(),weight);
      h_prec_p->Fill(p_recp.Mag(),weight);
      h_prec_angles->Fill(p_recp.Phi()*180./M_PI,p_recp.Theta()*180./M_PI,weight);

      h_prec_plead_angle->Fill(p_recp.Angle(pL)*180./M_PI,weight);
      h_lpangle_pmiss->Fill(pmiss.Mag(),p_recp.Angle(pL)*180./M_PI,weight);

      // close in angle to pmiss
      p_cos0 = pmiss.Dot(p_recp) / (pmiss.Mag() * p_recp.Mag());
      h_pcos0->Fill(p_cos0,weight);

      // fill observable histo
      h_pmiss_pp_uncorr->Fill(pmiss.Mag());
      h_pmiss_pp_corr->Fill(pmiss.Mag(),weight/peff); // new
      h_pp_count->Fill(pmiss.Mag(),weight); // fill counts - no peff correction

      // add to "with recoil" p denominator if proton meets recoil conditions
      h_pmiss_p_wrec->Fill(pmiss.Mag(),weight/peff);
    }





//// RECOIL N ////
    // look for recoil neutron in CTOF

    int rec_n = -1;
    h_nsize->Fill(neut.size(),weight);
    TVector3 p_recn;
    double n_cos0;

    for (int i=0; i<neut.size(); i++)
    {
      p_recn.SetMagThetaPhi(neut[i]->getP(),neut[i]->getTheta(),neut[i]->getPhi());
      // recoil n must be in CND
      bool is_CND1 = neut[i]->sci(CND1)->getDetector()==3;
      bool is_CND2 = neut[i]->sci(CND2)->getDetector()==3;
      bool is_CND3 = neut[i]->sci(CND3)->getDetector()==3;
      bool is_CTOF = neut[i]->sci(CTOF)->getDetector()==4;
      if (!is_CND1 && !is_CND2 && !is_CND3 && !is_CTOF) {continue;}


      // limit to central detector acceptance
      h_nrec_ptheta->Fill(p_recn.Mag(),p_recn.Theta()*180./M_PI,weight);
      if (p_recn.Mag()<0.3 || p_recn.Mag()>1.0) {continue;}
      if (p_recn.Theta()*180./M_PI<45 || p_recn.Theta()*180./M_PI>120) {continue;}

      /*// EPJA cuts
      if (p_recn.Theta()*180./M_PI<40 || p_recn.Theta()*180./M_PI>140) {continue;}
      if (p_recn.Mag()<0.3) {continue;}
      if (!is_CND1 && !is_CND2 && !is_CND3) {continue;}*/



      h_nrec_angles->Fill(p_recn.Phi()*180./M_PI,p_recn.Theta()*180./M_PI,weight);
      double n_phi = p_recn.Phi()*180./M_PI;
      if (n_phi>-106 && n_phi<-90) {continue;}
      if (n_phi>14 && n_phi<29) {continue;}
      if (n_phi>134 && n_phi<149) {continue;}


      // calculate features for ML
      Struct ninfo = getFeatures(neut, allParticles, i);
      cnd_hits = ninfo.cnd_hits;
      cnd_energy = ninfo.cnd_energy;
      ctof_hits = ninfo.ctof_hits;
      ctof_energy = ninfo.ctof_energy;
      layermult = ninfo.layermult;
      energy = ninfo.energy;
      size = ninfo.size;
      angle_diff = ninfo.angle_diff;
      // spectator
      momentum = p_recn.Mag();

      // apply ML model
      double mvaValue = reader->EvaluateMVA("MLP");
      h_mvaValue_MLP->Fill(mvaValue,weight);
      h_mvaValue_BDT->Fill(mvaValue,weight);

      // FILL MVA VALUE HISTOGRAM
      double mva_cutoff = 0.45;
      if (mlp_alt==1) {mva_cutoff = 0.4;}
      if (mlp_alt==2) {mva_cutoff = 0.5;}

      // fill counts
      h_pn_count->Fill(pmiss.Mag(),weight);

      

      if (mvaValue>mva_cutoff) // signal
      {
        // ML features
        h_energy_s->Fill(energy,weight);
        h_layermult_s->Fill(layermult,weight);
        h_size_s->Fill(size,weight);
        h_cnd_hits_s->Fill(cnd_hits,weight);
        h_cnd_energy_s->Fill(cnd_energy,weight);
        h_ctof_energy_s->Fill(ctof_energy,weight);
        h_ctof_hits_s->Fill(ctof_hits,weight);
        h_anglediff_s->Fill(angle_diff,weight);
        h_good_nrec_angles->Fill(p_recn.Phi()*180./M_PI,p_recn.Theta()*180./M_PI,weight);

        
        // see if neutron is close in angle to pmiss
        n_cos0 = pmiss.Dot(p_recn) / (pmiss.Mag() * p_recn.Mag());
        h_ncos0->Fill(n_cos0,weight);
        // compared to pmiss
        h_pn_pmiss->Fill(pmiss.Mag(),p_recn.Mag(),weight);


        // get momenta/angles of recoil neutrons (this was originally for candidates)
        h_nrec_plead->Fill(p_recn.Mag(),pL.Mag(),weight);
        h_nrec_p->Fill(p_recn.Mag(),weight);
        h_nptheta->Fill(p_recn.Theta()*180./M_PI,pmiss.Mag(),weight);
        h_nrec_plead_angle->Fill(p_recn.Angle(pL)*180./M_PI,weight);
        h_lnangle_pmiss->Fill(pmiss.Mag(),p_recn.Angle(pL)*180./M_PI,weight);


        rec_n = i; // pick this neutron! ... but what if there's more than 1?

        n_signal = n_signal + 1;
        h_n_signal->Fill(pmiss.Mag(),weight);

      }
      else
      {
        // ML features
        h_energy_b->Fill(energy,weight);
        h_layermult_b->Fill(layermult,weight);
        h_size_b->Fill(size,weight);
        h_cnd_hits_b->Fill(cnd_hits,weight);
        h_cnd_energy_b->Fill(cnd_energy,weight);
        h_ctof_energy_b->Fill(ctof_energy,weight);
        h_ctof_hits_b->Fill(ctof_hits,weight);
        h_anglediff_b->Fill(angle_diff,weight);
//std::cout << rec_n << '\n'; // there's an issue here! sometimes this returns 0

        n_background = n_background + 1;
        h_n_background->Fill(pmiss.Mag(),weight);
      }
    }



    if (rec_n<0) {continue;} // proceed only if we found good neutron
      // fill histo
      TVector3 pn;
      pn.SetXYZ( neut[rec_n]->par()->getPx(), neut[rec_n]->par()->getPy(), neut[rec_n]->par()->getPz() );
      double n_theta = pn.Theta()*180./M_PI;


      // CORRECTION FOR VETO EFFICIENCY
      double e_s = 0.86; double e_b = 0.225;
      if (mlp_alt == 1) {e_s = 0.88; e_b = 0.255;}
      else if (mlp_alt ==2) {e_s = 0.835; e_b = 0.197;}

      double p_to_n[7] = {3.30882, 2.75362, 2.01266, 1.8125, 2.17073, 1.9, 1.55556};
      double this_pton = 0; //double veto_weight = 1;
           if (pmiss.Mag()>=0.3 && pmiss.Mag()<0.4) {this_pton = p_to_n[0];}
      else if (pmiss.Mag()>=0.4 && pmiss.Mag()<0.5) {this_pton = p_to_n[1];}
      else if (pmiss.Mag()>=0.5 && pmiss.Mag()<0.6) {this_pton = p_to_n[2];}
      else if (pmiss.Mag()>=0.6 && pmiss.Mag()<0.7) {this_pton = p_to_n[3];}
      else if (pmiss.Mag()>=0.7 && pmiss.Mag()<0.8) {this_pton = p_to_n[4];}
      else if (pmiss.Mag()>=0.8 && pmiss.Mag()<0.9) {this_pton = p_to_n[5];}
      else if (pmiss.Mag()>=0.9 && pmiss.Mag()<1.0) {this_pton = p_to_n[6];}
      double veto_weight =  ( (1-e_b) - e_b*this_pton ) / (e_s-e_b);
//std::cout << veto_weight << '\n';

      // CORRECTION FOR CND EFFICIENCY
      double p0 = 1; double p1 = 1;

      // regular neff calculation
      double a0[9] = {0.186947, 0.211831, 0.203986, 0.199245, 0.195636, 0.18393, 0.186905, 0.18791, 0.19951 };
      double a1[9] = {-0.117272, -0.144642, -0.144253, -0.145586, -0.144762, -0.132524, -0.142502, -0.150404, -0.1738 };

      // alternative neff calculation
      if (neff_alt)
      {
        double a0_copy[9] = {0.177487, 0.207873, 0.198617, 0.181738, 0.180914, 0.169454, 0.160552, 0.138531, 0.178826 }; // alternative neff
        double a1_copy[9] = {-0.0786727, -0.116044, -0.115987, -0.100439, -0.102689, -0.0958549, -0.0890577, -0.0641776, -0.126627 }; // alternative neff
        memcpy(a0, a0_copy, 9); memcpy(a1, a1_copy, 9);
      }

      if      (n_theta>45 && n_theta<50)  {p0 = a0[0]; p1 = a1[0];}
      else if (n_theta>50 && n_theta<55)  {p0 = a0[1]; p1 = a1[1];}
      else if (n_theta>55 && n_theta<60)  {p0 = a0[2]; p1 = a1[2];}
      else if (n_theta>60 && n_theta<65)  {p0 = a0[3]; p1 = a1[3];}
      else if (n_theta>65 && n_theta<70)  {p0 = a0[4]; p1 = a1[4];}
      else if (n_theta>70 && n_theta<75)  {p0 = a0[5]; p1 = a1[5];}
      else if (n_theta>75 && n_theta<80)  {p0 = a0[6]; p1 = a1[6];}
      else if (n_theta>80 && n_theta<95)  {p0 = a0[7]; p1 = a1[7];}
      else if (n_theta>95 && n_theta<120) {p0 = a0[8]; p1 = a1[8];}

      double neff = p0 + p1*pn.Mag();
      if (neff<=0) {continue;} // just in case

//if (p_recn.Mag()<0.3) {continue;} // THIS IS A HUGE PROBLEM! It eliminates all the events in which the final particle with PID 2112 has momentum of 0, even if they have a good neutron in the event.

      h_pmiss_pn_uncorr->Fill(pmiss.Mag(),weight);
      h_pmiss_pn_corr->Fill(pmiss.Mag(),veto_weight*weight/neff);

      // add to "with recoil" p denominator if neutron meets recoil conditions
      h_pmiss_p_wrec->Fill(pmiss.Mag(),veto_weight*weight/neff);


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

  //h_lr_angles->Write();

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
  // lead proton in FTOF

  myText->cd();
  text.DrawLatex(0.2,0.9,"Basic electron cuts");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_vtzdiff_ep->Draw();
  myCanvas->cd(2);  h_chi2pid->Draw();
  myCanvas->cd(3);  h_dbetap->Draw("colz");
  myCanvas->cd(4);  h_betap->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(2);  h_pmisstheta->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"SRC cuts:");
  text.DrawLatex(0.2,0.8,"p_{miss} > 0.3 GeV/c");
  text.DrawLatex(0.2,0.7,"x_{B} > 1.2");
  text.DrawLatex(0.2,0.6,"0.62 < p/q < 0.96");
  text.DrawLatex(0.2,0.5,"#theta_{pq} > 25 deg");
  text.DrawLatex(0.2,0.4,"M_{miss} < 1.1 GeV/c^{2}");
  //text.DrawLatex(0.2,0.6,"Q^{2} > 1.5 GeV^{2}");

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_p_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_q2_xb->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pmiss_pL->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_q2->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pq->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_mmiss->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  



  myText->Print(fileName,"pdf");
  myText->Clear();

  myText->cd();
  text.DrawLatex(0.2,0.9,"Recoil protons");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_p->Draw();
  myCanvas->cd(3);  h_prec_plead->Draw("colz");
  myCanvas->cd(4);  h_nrec_plead->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_psize->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pcos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_prec_p->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_prec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_prec_ptheta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_prec_plead_angle->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_lpangle_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myText->cd();
  text.DrawLatex(0.2,0.9,"Recoil neutrons");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nsize->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_ncos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nrec_p->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nrec_ptheta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nrec_angles->Draw("colz");
  TLine * l_n1 = new TLine(-106,0,-106,180);
  l_n1->SetLineColor(kRed);
  l_n1->SetLineWidth(3);
  l_n1->Draw("same");
  TLine * l_n2 = new TLine(-90,0,-90,180);
  l_n2->SetLineColor(kRed);
  l_n2->SetLineWidth(3);
  l_n2->Draw("same");
  TLine * l_n3 = new TLine(14,0,14,180);
  l_n3->SetLineColor(kRed);
  l_n3->SetLineWidth(3);
  l_n3->Draw("same");
  TLine * l_n4 = new TLine(29,0,29,180);
  l_n4->SetLineColor(kRed);
  l_n4->SetLineWidth(3);
  l_n4->Draw("same");
  TLine * l_n5 = new TLine(134,0,134,180);
  l_n5->SetLineColor(kRed);
  l_n5->SetLineWidth(3);
  l_n5->Draw("same");
  TLine * l_n6 = new TLine(149,0,149,180);
  l_n6->SetLineColor(kRed);
  l_n6->SetLineWidth(3);
  l_n6->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_good_nrec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_nrec_plead_angle->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // ML output
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_mvaValue_MLP->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_mvaValue_BDT->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // ML features start here
  myCanvas->Divide(3,3);
  myCanvas->cd(1);
  //h_energy_s->SetFillColor(kBlue);
  //h_energy_s->SetFillStyle(1001);
  h_energy_s->Draw();
  h_energy_s->SetStats(0);
  h_energy_b->SetLineColor(kRed);
  //h_energy_b->SetFillColor(kRed);
  //h_energy_b->SetFillStyle(3554);
  h_energy_b->Draw("same");
  h_energy_b->SetStats(0);
  h_energy_b->GetYaxis()->SetRangeUser(0,1.2*max(h_energy_s->GetMaximum(),h_energy_b->GetMaximum()));
  myCanvas->cd(2);
  h_layermult_s->Draw();
  h_layermult_s->SetFillColorAlpha(kBlue,0.35);
  h_layermult_s->SetStats(0);
  h_layermult_s->GetYaxis()->SetRangeUser(0,1.2*max(h_layermult_s->GetMaximum(),h_layermult_b->GetMaximum()));
  h_layermult_b->SetLineColor(kRed);
  h_layermult_b->Draw("same");
  h_layermult_b->SetStats(0);
  h_layermult_b->SetFillColorAlpha(kRed,0.35);
  h_layermult_b->GetYaxis()->SetRangeUser(0,1.2*max(h_layermult_s->GetMaximum(),h_layermult_b->GetMaximum()));
  myCanvas->cd(3);
  h_size_s->Draw();
  h_size_s->SetStats(0);
  h_size_s->GetYaxis()->SetRangeUser(0,1.2*max(h_size_s->GetMaximum(),h_size_b->GetMaximum()));
  h_size_b->SetLineColor(kRed);
  h_size_b->Draw("same");
  h_size_b->SetStats(0);
  h_size_b->GetYaxis()->SetRangeUser(0,1.2*max(h_size_s->GetMaximum(),h_size_b->GetMaximum()));
  myCanvas->cd(4);
  h_cnd_hits_s->Draw();
  h_cnd_hits_s->SetStats(0);
  h_cnd_hits_s->GetYaxis()->SetRangeUser(0,1.2*max(h_cnd_hits_s->GetMaximum(),h_cnd_hits_b->GetMaximum()));
  h_cnd_hits_b->SetLineColor(kRed);
  h_cnd_hits_b->Draw("same");
  h_cnd_hits_b->SetStats(0);
  h_cnd_hits_b->GetYaxis()->SetRangeUser(0,1.2*max(h_cnd_hits_s->GetMaximum(),h_cnd_hits_b->GetMaximum()));
  myCanvas->cd(5);
  h_cnd_energy_s->Draw();
  h_cnd_energy_s->SetStats(0);
  h_cnd_energy_s->GetYaxis()->SetRangeUser(0,1.2*max(h_cnd_energy_s->GetMaximum(),h_cnd_energy_b->GetMaximum()));
  h_cnd_energy_b->SetLineColor(kRed);
  h_cnd_energy_b->Draw("same");
  h_cnd_energy_b->SetStats(0);
  h_cnd_energy_b->GetYaxis()->SetRangeUser(0,1.2*max(h_cnd_energy_s->GetMaximum(),h_cnd_energy_b->GetMaximum()));
  myCanvas->cd(6);
  h_ctof_energy_s->Draw();
  h_ctof_energy_s->SetStats(0);
  h_ctof_energy_s->GetYaxis()->SetRangeUser(0,1.2*max(h_ctof_energy_s->GetMaximum(),h_ctof_energy_b->GetMaximum()));
  h_ctof_energy_b->SetLineColor(kRed);
  h_ctof_energy_b->Draw("same");
  h_ctof_energy_b->SetStats(0);
  h_ctof_energy_b->GetYaxis()->SetRangeUser(0,1.2*max(h_ctof_energy_s->GetMaximum(),h_ctof_energy_b->GetMaximum()));
  myCanvas->cd(7);
  h_ctof_hits_s->Draw();
  h_ctof_hits_s->SetStats(0);
  h_ctof_hits_s->GetYaxis()->SetRangeUser(0,1.2*max(h_ctof_hits_s->GetMaximum(),h_ctof_hits_b->GetMaximum()));
  h_ctof_hits_b->SetLineColor(kRed);
  h_ctof_hits_b->Draw("same");
  h_ctof_hits_b->SetStats(0);
  h_ctof_hits_b->GetYaxis()->SetRangeUser(0,1.2*max(h_ctof_hits_s->GetMaximum(),h_ctof_hits_b->GetMaximum()));
  myCanvas->cd(8);
  h_anglediff_s->Draw();
  h_anglediff_s->SetStats(0);
  h_anglediff_s->GetYaxis()->SetRangeUser(0,1.2*max(h_anglediff_s->GetMaximum(),h_anglediff_b->GetMaximum()));
  h_anglediff_b->SetLineColor(kRed);
  h_anglediff_b->Draw("same");
  h_anglediff_b->SetStats(0);
  h_anglediff_b->GetYaxis()->SetRangeUser(0,1.2*max(h_anglediff_s->GetMaximum(),h_anglediff_b->GetMaximum()));
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();





  // channels
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_p->Draw();
  myCanvas->cd(2);
  h_pmiss_p_wrec->Draw();
  myCanvas->cd(3);
  h_pmiss_pp_uncorr->Draw();
  myCanvas->cd(4);
  h_pmiss_pn_uncorr->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // here we show the charged particle veto correction for pn
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_n_background->SetLineColor(kRed+1);
  h_n_background->Draw();
  h_n_signal->SetLineColor(kCyan-3);
  h_n_signal->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // channels
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_p->Draw();
  myCanvas->cd(2);
  h_pmiss_p_wrec->Draw();
  myCanvas->cd(3);
  h_pmiss_pp_corr->Draw();
  myCanvas->cd(4);
  h_pmiss_pn_corr->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  
  // print yields
  // (e,e'p)
  std::cout << "\n(e,e'p) yields:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_p->GetBinContent(i+1) << ", ";}
  std::cout << "\n(e,e'p) yield errors:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_p->GetBinError(i+1) << ", ";}
  // (e,e'pp)
  std::cout << "\n(e,e'pp) yields:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_pp_corr->GetBinContent(i+1) << ", ";}
  std::cout << "\n(e,e'pp) yield errors:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_pp_corr->GetBinError(i+1) << ", ";}
  // (e,e'pp)
  std::cout << "\n(e,e'pn) yields:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_pn_corr->GetBinContent(i+1) << ", ";}
  std::cout << "\n(e,e'pn) yield errors:\n";
  for (int i=0; i<pmiss_bins; i++) {std::cout << h_pmiss_pn_corr->GetBinError(i+1) << ", ";}
  std::cout << '\n';



  // OBSERVABLE: pp/pN and pn/pN
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  // pn/p (with recoil)
  TH1D * h_pmiss_pn_pwrec = (TH1D*)h_pmiss_pn_corr->Clone();
  h_pmiss_pn_pwrec->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_pwrec->SetLineColor(kRed);
  h_pmiss_pn_pwrec->Draw();
  h_pmiss_pn_pwrec->GetYaxis()->SetTitle("pp/pN and pn/pN");
  h_pmiss_pn_pwrec->SetStats(0);
  h_pmiss_pn_pwrec->GetYaxis()->SetRangeUser(0,1.2);
  // pp/p (with recoil)
  TH1D * h_pmiss_pp_pwrec = (TH1D*)h_pmiss_pp_corr->Clone();
  h_pmiss_pp_pwrec->Divide(h_pmiss_p_wrec);
  h_pmiss_pp_pwrec->SetLineColor(kBlue);
  h_pmiss_pp_pwrec->Draw("same");
  h_pmiss_pp_pwrec->GetYaxis()->SetTitle("pp/pN and pn/pN");
  h_pmiss_pp_pwrec->SetStats(0);
  h_pmiss_pp_pwrec->GetYaxis()->SetRangeUser(0,1.2);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // OBSERVABLE: pp/p and pn/p
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  // pn/p
  TH1D * h_pmiss_pn_p = (TH1D*)h_pmiss_pn_corr->Clone();
  h_pmiss_pn_p->Divide(h_pmiss_p);
  h_pmiss_pn_p->SetLineColor(kRed);
  h_pmiss_pn_p->Draw();
  h_pmiss_pn_p->GetYaxis()->SetTitle("pp/p and pn/p");
  h_pmiss_pn_p->SetStats(0);
  h_pmiss_pn_p->GetYaxis()->SetRangeUser(0,1);
  // pp/p
  TH1D * h_pmiss_pp_p = (TH1D*)h_pmiss_pp_corr->Clone();
  h_pmiss_pp_p->Divide(h_pmiss_p);
  h_pmiss_pp_p->SetLineColor(kBlue);
  h_pmiss_pp_p->Draw("same");
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p and pn/p");
  h_pmiss_pp_p->SetStats(0);
  h_pmiss_pp_p->GetYaxis()->SetRangeUser(0,1);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // OBSERVABLE: pp/2pn
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_pp_pn = (TH1D*)h_pmiss_pp_corr->Clone();
  h_pp_pn->Divide(h_pmiss_pn_corr);
  h_pp_pn->Scale(0.5);
  h_pp_pn->SetLineColor(kMagenta);
  h_pp_pn->Draw();
  h_pp_pn->GetYaxis()->SetTitle("pp/2pn");
  h_pp_pn->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // write yield counts to root file
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

  std::cout << "SIGNAL NEUTRONS:" << '\t' << n_signal << '\n';
  std::cout << "BACKGROUND NEUTRONS:" << '\t' << n_background << '\n';
  std::cout << "NEUTRON BACKGROUND/SIGNAL:" << '\t' << (double)n_background/(double)n_signal << '\n' << '\n';

  for (int j=1; j<pmiss_bins+1; j++) // skip underflow and overflow bins
  {
    double num_b = h_n_background->GetBinContent(j);
    double num_s = h_n_signal->GetBinContent(j);
    std::cout << num_b << '\t' << num_s << '\t' << num_b/num_s << '\n';
  }


  //outFile->Close(); // THIS LINE CAUSES ERRORS!
}



void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

