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
  std::cerr << "Usage: ./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/input.hipo> \n";
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

  // create instance of clas12ana class
  clas12ana clasAna;

  clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsSF_LD2_x2.dat");
  clasAna.readEcalPPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsPI_LD2_x2.dat");

  clasAna.printParams();

  clas12root::HipoChain chain;
  for(int k = 5; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  chain.db()->turnOffQADB();

  // open root file with neutron detection efficiency
  TFile * f_neff = new TFile("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm_build/NeutronEfficiency/neff_2gev_pCDn.root","READ");
  TH2D * h_neff = (TH2D*)f_neff->Get("det2d");

  // open root file with proton detection efficiency
  TFile * f_peff = new TFile("/w/hallb-scshelf2102/clas/clase2/erins/repos/rgm_build/Ana/proton_efficiency/peff_test_101023.root","READ");
  TH1D * h_peff = (TH1D*)f_peff->Get("efficiency_p");

        
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
  TH2D * h_p_theta = new TH2D("p_theta","Lead proton p vs #theta;#theta;p (GeV/c)",90,0,180,50,0,1.5);
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

  TH2D * h_lr_angles = new TH2D("lr_angles","Lead and Recoil Polar Angles;Lead theta;Recoil theta",180,0,180,180,0,180);
  hist_list_2.push_back(h_lr_angles);


  /////////////////////////////////////
  //Histos: CD background
  /////////////////////////////////////
  TH2D * h_beta_p_all = new TH2D("beta_p_all","#beta vs Momentum in the CD;p (GeV/c);#beta",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_beta_p_all);
  TH2D * h_beta_p_pFD = new TH2D("beta_p_pFD","#beta vs Momentum in the CD (assuming p in FD);p (GeV/c);#beta",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_beta_p_pFD);
  TH2D * h_beta_p_pCD = new TH2D("beta_p_pCD","#beta vs Momentum in the CD (assuming p in CD);p (GeV/c);#beta",100,0,1.5,100,0,1.5);
  hist_list_2.push_back(h_beta_p_pCD);
  TH1D * h_numCVT = new TH1D("numCVT","Number of Charged Particles in CVT",10,0,10);
  hist_list_1.push_back(h_numCVT);
  TH1D * h_CVTpid = new TH1D("CVTpid","PID of Charged Particles in CVT",2000,-2500,2500);
  hist_list_1.push_back(h_CVTpid);


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

  TH2D * h_prec_ptheta = new TH2D("prec_ptheta","Recoil Proton Theta vs Momentum;Momentum (GeV/c);#theta (degrees)",30,0.2,1.2,20,40,140);
  hist_list_2.push_back(h_prec_ptheta);
  TH2D * h_prec_angles = new TH2D("prec_angles","Recoil Proton Angular Distribution;phi (deg);theta (deg)",48,-180,180,45,0,180);
  hist_list_2.push_back(h_prec_angles);

  TH2D * h_pptheta = new TH2D("pptheta","Proton #theta vs Momentum;#theta_{p};Momentum (GeV/c)",180,0,180,50,0,1.5);
  hist_list_2.push_back(h_pptheta);
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
  TH2D * h_nrec_ptheta = new TH2D("nrec_ptheta","Recoil Neutron Theta vs Momentum;Momentum (GeV/c);#theta (degrees)",30,0.2,1.2,20,40,140);
  hist_list_2.push_back(h_prec_ptheta);
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
  //Histos: pmiss
  /////////////////////////////////////
  TH1D * h_pmiss_p = new TH1D("pmiss_p","Missing Momentum (e,e'p_{SRC});p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_p);
  TH1D * h_pmiss_p_wrec = new TH1D("pmiss_p_wrec","Missing Momentum (e,e'p_{SRC}N_{rec});p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_p_wrec);
  TH1D * h_pmiss_pp = new TH1D("pmiss_pp","Missing Momentum (e,e'p_{SRC}p_{rec});p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_pp);
  /*TH1D * h_pmiss_pp_corr = new TH1D("pmiss_pp_corr","Missing Momentum (e,e'p_{SRC}p_{rec}) (efficiency corrected);p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_pp_corr);*/
  TH1D * h_pmiss_pn = new TH1D("pmiss_pn","Missing Momentum (e,e'p_{SRC}n_{rec});p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_pn);
  /*TH1D * h_pmiss_pn_corr = new TH1D("pmiss_pn_corr","Missing Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{miss} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pmiss_pn_corr);*/

  TH1D * h_pn = new TH1D("pn","Momentum (e,e'p_{SRC}n_{rec});p_{n} (GeV/c);Counts",10,0.2,1);
  hist_list_1.push_back(h_pn);
  TH1D * h_pn_corr = new TH1D("pn_corr","Momentum (e,e'p_{SRC}n_{rec}) (efficiency corrected);p_{n} (GeV/c);Counts",10,0.2,1);
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


  // set up clas12ana cuts
  //clasAna.setEcalSFCuts();
  //clasAna.setEcalPCuts();
  //clasAna.setEcalEdgeCuts(false); // makes particle PID arrays empty
  //clasAna.setPidCuts(false); // I think I want this to be default?
  //clasAna.setVertexCuts();
  //clasAna.setVertexCorrCuts();
  //clasAna.setDCEdgeCuts(false); // makes particle PID arrays empty
  //clasAna.setCDEdgeCuts(true);
  //  clasAna.setCDRegionCuts();

  clasAna.setVzcuts(-5.3,-0.1);
  //  clasAna.setCDCutRegion(2);  
  //clasAna.setVertexCorrCuts(-3,1);

  // use Andrew's PID for CD protons
  clasAna.setProtonPidCuts(true);

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
    //if (prot.size()<1) {continue;}
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
    double W2 = (mp*mp) - Q2 + (2*nu*mp);
    double etheta = elec[0]->getTheta()*180./M_PI;
    int nphe = elec[0]->che(HTCC)->getNphe();
    double vtz_e = elec[0]->par()->getVz();
   
    // electron histograms: quality cuts
    h_nphe->Fill(nphe,weight);
    h_vtz_e->Fill(vtz_e,weight);
    if (vtz_e<-6 || vtz_e>2) {continue;}


    // CD background - all charged particles in the CVT
    for (int i=0; i<allParticles.size(); i++)
    {
      if (allParticles[i]->par()->getCharge()==0) {continue;} // look at only charged particles
      if (allParticles[i]->trk(CVT)->getDetector()!=5) {continue;} // look at only particles in CVT
      h_beta_p_all->Fill(allParticles[i]->getP(),allParticles[i]->par()->getBeta());
    }

    for (int i=0; i<prot.size(); i++)
    {
      TVector3 pp_temp;
      pp_temp.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
      double mmiss_temp = get_mmiss(pb,pe,pp_temp);
      h_mmiss_all->Fill(mmiss_temp,weight);
    }



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



    // CD background (assuming p in FD)
    for (int i=0; i<allParticles.size(); i++)
    {
      if (allParticles[i]->par()->getCharge()==0) {continue;} // look at only charged particles
      if (allParticles[i]->trk(CVT)->getDetector()!=5) {continue;} // look at only particles in CVT
      if (ltheta<40) {h_beta_p_pFD->Fill(allParticles[i]->getP(),allParticles[i]->par()->getBeta());}
      if (ltheta>40 && ltheta<140) {h_beta_p_pCD->Fill(allParticles[i]->getP(),allParticles[i]->par()->getBeta());}
    }




    // lead histos
    h_vtzdiff_ep->Fill(vze-vzlead,weight);
    if ((vze-vzlead)<-3. || (vze-vzlead)>3.) {continue;}
    h_chi2pid->Fill(chi2pid,weight);
    if (chi2pid<-3.0 || chi2pid>3.0) {continue;}

    h_dbetap->Fill(pL.Mag(),dbetap,weight);
    h_betap->Fill(pL.Mag(),lead[0]->par()->getBeta(),weight);



    // look at TOF distributions



    //// PMISS & LEAD SRC CUTS ////
   
    h_pangles->Fill(lphi,ltheta,weight);
    h_p_theta->Fill(ltheta,pL.Mag(),weight);
    if (lead[0]->par()->getBeta()<0.2) {continue;}
    if (ltheta<40 || ltheta>140) {continue;}
    //if (ltheta>40) {continue;}


    TVector3 pmiss = pL - q;
    double pm_theta = pmiss.Theta()*180./M_PI;
    h_pmisstheta->Fill(pm_theta,weight);

    // src cuts
    h_q2_xb->Fill(xB,Q2,weight);
    if (xB<1.2) {continue;}

    h_pmiss_pL->Fill(pmiss.Mag(),weight);
    if (pmiss.Mag()<0.3) {continue;}

    if (pL.Mag()<1.0) {continue;}

    h_pq->Fill(pL.Mag()/q.Mag(),theta_pq,weight);
    if (pL.Mag()/q.Mag()>0.96) {continue;}

    h_mmiss->Fill(mmiss,weight);
    if (mmiss>1.1) {continue;}




      // get proton efficiency
      //if (pL.Mag() > 1.2) {continue;} // JUST FOR NOW - THERE WILL BE LOTS OF FAST LEADS
      int p_bin = h_peff->GetXaxis()->FindBin(pL.Mag());
      double l_peff = 1;
      //l_peff = h_peff->GetBinContent(p_bin);
      //if (l_peff==0) {l_peff = 1e-6;}

    h_pmiss_p->Fill(pmiss.Mag(),weight/l_peff);


    // count number of recoils in the CVT
    int num_in_CVT = 0;
    for (int i=0; i<allParticles.size(); i++)
    {
      if (allParticles[i]->par()->getCharge()==0) {continue;} // look at only charged particles
      if (allParticles[i]->trk(CVT)->getDetector()!=5) {continue;} // look at only particles in CVT
      num_in_CVT = num_in_CVT + 1;
      h_CVTpid->Fill(allParticles[i]->getPid(),weight);
    }
    h_numCVT->Fill(num_in_CVT);


//// RECOIL P ////
    // look for recoil proton in CTOF

    h_psize->Fill(prot.size(),weight);
    //if (recoil.size()>1) {continue;} // really, if there's one recoil, I should decide which is the best

    TVector3 p_recp;
    TVector3 p_vecX;
    double p_cos0;

    if (recoil.size()==1)
    {
      p_recp.SetMagThetaPhi(recoil[0]->getP(),recoil[0]->getTheta(),recoil[0]->getPhi());
      // recoil p must be in CTOF
      bool is_CTOF = recoil[0]->sci(CTOF)->getDetector()==4;
      if (!is_CTOF) {continue;}

      // limit to central detector acceptance
      if (p_recp.Mag()<0.3) {continue;}
      if (p_recp.Theta()*180./M_PI<50 || p_recp.Theta()*180./M_PI>140) {continue;}

      // get proton efficiency
      int pp_bin = h_peff->GetXaxis()->FindBin(p_recp.Mag());
      double r_peff = 1;
      r_peff = h_peff->GetBinContent(pp_bin); // what if not in range?
      if (r_peff==0) {r_peff = 0.5; std::cout << "peff = " << '\t' << r_peff << '\n';}
      



      // get momenta/angles of recoil protons
      h_prec_plead->Fill(p_recp.Mag(),pL.Mag(),weight);
      h_prec_p->Fill(p_recp.Mag(),weight);
      h_prec_ptheta->Fill(p_recp.Mag(),p_recp.Theta()*180./M_PI,weight);
      h_prec_angles->Fill(p_recp.Phi()*180./M_PI,p_recp.Theta()*180./M_PI,weight);
      h_pptheta->Fill(p_recp.Theta()*180./M_PI,pmiss.Mag(),weight);
      h_prec_plead_angle->Fill(p_recp.Angle(pL)*180./M_PI,weight);
      h_lpangle_pmiss->Fill(pmiss.Mag(),p_recp.Angle(pL)*180./M_PI,weight);

      // close in angle to pmiss
      p_vecX.SetXYZ( recoil[0]->par()->getPx(), recoil[0]->par()->getPy(), recoil[0]->par()->getPz() );
      p_cos0 = pmiss.Dot(p_vecX) / (pmiss.Mag() * p_vecX.Mag());
      h_pcos0->Fill(p_cos0,weight);


      // fill histo
      h_pmiss_pp->Fill(pmiss.Mag(),weight/((l_peff*r_peff)));
      //h_pmiss_pp_corr->Fill(pmiss.Mag(),weight/(peff*ppeff)); // new

      // add to "with recoil" p denominator if proton meets recoil conditions
      h_pmiss_p_wrec->Fill(pmiss.Mag(),weight/(l_peff*r_peff));
    }





//// RECOIL N ////
    // look for recoil neutron in CTOF

    int rec_n = -1;
    h_nsize->Fill(neut.size(),weight);
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

      // limit to central detector acceptance
      if (p_recn.Mag()<0.3) {continue;}
      if (p_recn.Theta()*180./M_PI<50 || p_recn.Theta()*180./M_PI>140) {continue;}


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
      //double mvaErr = reader->GetMVAError(); // returns 0 :(
      //double pSig = reader->GetProba("MLP",0.9); // doesn't exist for MLP
      h_mvaValue_MLP->Fill(mvaValue,weight);
      h_mvaValue_BDT->Fill(mvaValue,weight);

      // FILL MVA VALUE HISTOGRAM
      if (mvaValue>0.45) // signal
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
        n_vecX.SetXYZ( neut[i]->par()->getPx(), neut[i]->par()->getPy(), neut[i]->par()->getPz() );
        n_cos0 = pmiss.Dot(n_vecX) / (pmiss.Mag() * n_vecX.Mag());
        h_ncos0->Fill(n_cos0,weight);
        //if (n_cos0>-0.8) {continue;} // not needed?
        // compared to pmiss
        h_pn_pmiss->Fill(pmiss.Mag(),p_recn.Mag(),weight);
        rec_n = i; // pick this neutron! ... but what if there's more than 1?

        // get momenta/angles of recoil neutrons (this was originally for candidates)
        h_nrec_plead->Fill(p_recn.Mag(),pL.Mag(),weight);
        h_nrec_p->Fill(p_recn.Mag(),weight);
        h_nrec_ptheta->Fill(p_recn.Mag(),p_recn.Theta()*180./M_PI,weight);
        h_nrec_angles->Fill(p_recn.Phi()*180./M_PI,p_recn.Theta()*180./M_PI,weight);
        h_nptheta->Fill(p_recn.Theta()*180./M_PI,pmiss.Mag(),weight);
        h_nrec_plead_angle->Fill(p_recn.Angle(pL)*180./M_PI,weight);
        h_lnangle_pmiss->Fill(pmiss.Mag(),p_recn.Angle(pL)*180./M_PI,weight);


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
      }
    }



    // use the neutron we chose above
    if (rec_n>-1)
    {
      // fill histo
      TVector3 pn;
      pn.SetXYZ( neut[rec_n]->par()->getPx(), neut[rec_n]->par()->getPy(), neut[rec_n]->par()->getPz() );
      int p_bin = h_neff->GetXaxis()->FindBin(pn.Mag());
      int t_bin = h_neff->GetYaxis()->FindBin(pn.Theta()*180./M_PI);
      //double neff = 0;
      //neff = h_neff->GetBinContent(p_bin,t_bin);

      //if (neff==0) neff = 0.10;
      double p0 = 0.21;// 0.172689 + 0.000582989*pn.Theta()*180./M_PI;
      double p1 = -0.23;//0.0638342 - 0.00235126*pn.Theta()*180./M_PI;
      double neff = p0 + p1*pn.Mag();
      if (neff<0) {neff=0.02;}
//std::cout << neff << '\n';

      double veff = 0.85;//0.72/0.90;
//std::cout << neff << '\t' << pn.Mag() << '\t' << pn.Theta()*180./M_PI << '\n';
      h_pmiss_pn->Fill(pmiss.Mag(),weight/(l_peff*neff*veff));
      //h_pmiss_pn_corr->Fill(pmiss.Mag(),weight/(neff*veff));
      //h_pn->Fill(pn.Mag(),weight);
      //h_pn_corr->Fill(pn.Mag(),weight/(neff*veff));

      //h_lr_angles->Fill(ltheta,prot[rec_n]->getTheta()*180./M_PI,weight); // segfault?
      

      // add to "with recoil" p denominator if neutron meets recoil conditions
      h_pmiss_p_wrec->Fill(pmiss.Mag(),weight/(l_peff*neff*veff));
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
  /*myText->cd();
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
  myText->Clear();*/

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

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_beta_p_all->Draw("colz");
  myCanvas->cd(2);  h_beta_p_pFD->Draw("colz");
  myCanvas->cd(3);  h_beta_p_pCD->Draw("colz");
  myCanvas->cd(4);  h_mmiss_all->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pangles->Draw("colz");
  myCanvas->cd(2);  h_p_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pangles->Draw("colz");
  myCanvas->cd(2);  h_pmisstheta->Draw();
  myCanvas->cd(3);  h_q2_xb->Draw("colz");
  myCanvas->cd(4);  h_pmiss_pL->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pq->Draw("colz");
  myCanvas->cd(2);  h_mmiss->Draw();
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

  myText->Print(fileName,"pdf");
  myText->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_numCVT->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_CVTpid->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


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


  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_psize->Draw();
  myCanvas->cd(2);  h_pcos0->Draw();
  myCanvas->cd(3);  h_prec_p->Draw();
  myCanvas->cd(4);  h_prec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);

  myCanvas->cd(2);  h_pptheta->Draw("colz");
  myCanvas->cd(3);  h_prec_ptheta->Draw("colz");
  myCanvas->cd(4);  h_prec_plead_angle->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_lpangle_pmiss->Draw("colz");
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
  
  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_good_nrec_angles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pn_pmiss->Draw("colz");
  myCanvas->cd(2);  h_nptheta->Draw("colz");
  myCanvas->cd(3);  h_nrec_ptheta->Draw("colz");
  myCanvas->cd(4);  h_nrec_plead_angle->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_lnangle_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // ML output
  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_mvaValue_MLP->Draw();
  myCanvas->cd(2);  h_mvaValue_BDT->Draw();
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






/*  myCanvas->Divide(1,1);
  myCanvas->cd(1);  h_pmiss_p_wrec->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/

  // non-efficiency-corrected observable (pp/pN, pn/pN)
/*  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pmiss_pn->Draw();
  myCanvas->cd(3);
  TH1D * h_pmiss_pp_p = (TH1D*)h_pmiss_pp->Clone();
  h_pmiss_pp_p->Divide(h_pmiss_p_wrec);
  h_pmiss_pp_p->Draw();
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p (eff.corrected)");
  myCanvas->cd(4);
  TH1D * h_pmiss_pn_p = (TH1D*)h_pmiss_pn->Clone();
  h_pmiss_pn_p->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_p->Draw();
  h_pmiss_pn_p->GetYaxis()->SetTitle("pn/p (eff. corrected)");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/



/*  // efficiency-corrected observable (pp/pN, pn/pN)
  myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp_corr->Draw();
  myCanvas->cd(2);  h_pmiss_pn_corr->Draw();
  myCanvas->cd(3);
  TH1D * h_pmiss_pp_p_corr = (TH1D*)h_pmiss_pp_corr->Clone();
  h_pmiss_pp_p_corr->Divide(h_pmiss_p_wrec);
  h_pmiss_pp_p_corr->Draw();
  h_pmiss_pp_p_corr->GetYaxis()->SetTitle("pp/p (eff. corrected)");
  myCanvas->cd(4);
  TH1D * h_pmiss_pn_p_corr = (TH1D*)h_pmiss_pn_corr->Clone();
  h_pmiss_pn_p_corr->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_p_corr->Draw();
  h_pmiss_pn_p_corr->GetYaxis()->SetTitle("pn/p (eff. corrected)");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/


  // efficiency-corrected observable (pp/pN, pn/pN)
  // neutron momentum on x-axis instead of missing momentum
  // proton efficiency corrections not added
/*  myCanvas->Divide(2,2);
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
  myCanvas->Clear();*/

  // pp/pN and pn/pN - IDK what this is or why it matters
  /*myCanvas->Divide(2,2);
  myCanvas->cd(1);  h_pmiss_pp->Draw();
  myCanvas->cd(2);  h_pn_corr->Draw();
  myCanvas->cd(3);  h_pmiss_pp_p_corr->Draw();
  h_pmiss_pp_p_corr->GetYaxis()->SetTitle("pp/pN and pn/pN");
  myCanvas->cd(4);  h_pmiss_pn_p_corr->Draw();
  TH1D * h_pn_p_corr = (TH1D*)h_pn_corr->Clone(); // recently changed from pn_p to pmiss_p
  h_pn_p_corr->Divide(h_pmiss_p_wrec);
  h_pn_p_corr->Draw();
  h_pn_p_corr->GetYaxis()->SetTitle("pp/pN and pn/pN");
  h_pmiss_pn_p_corr->GetYaxis()->SetTitle("pp/pN and pn/pN");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/


  // channels
  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pmiss_p->Draw();
  myCanvas->cd(2);
  h_pmiss_p_wrec->Draw();
  myCanvas->cd(3);
  h_pmiss_pp->Draw();
  myCanvas->cd(4);
  h_pmiss_pn->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  // OBSERVABLE: pp/pN and pn/pN
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  // pn/p (with recoil)
  TH1D * h_pmiss_pn_pwrec = (TH1D*)h_pmiss_pn->Clone();
  h_pmiss_pn_pwrec->Divide(h_pmiss_p_wrec);
  h_pmiss_pn_pwrec->SetLineColor(kRed);
  h_pmiss_pn_pwrec->Draw();
  h_pmiss_pn_pwrec->GetYaxis()->SetTitle("pp/pN and pn/pN");
  h_pmiss_pn_pwrec->SetStats(0);
  // pp/p (with recoil)
  TH1D * h_pmiss_pp_pwrec = (TH1D*)h_pmiss_pp->Clone();
  h_pmiss_pp_pwrec->Divide(h_pmiss_p_wrec);
  h_pmiss_pp_pwrec->SetLineColor(kBlue);
  h_pmiss_pp_pwrec->Draw("same");
  h_pmiss_pp_pwrec->GetYaxis()->SetTitle("pp/pN and pn/pN");
  h_pmiss_pp_pwrec->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // OBSERVABLE: pp/p and pn/p
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  // pn/p
  TH1D * h_pmiss_pn_p = (TH1D*)h_pmiss_pn->Clone();
  h_pmiss_pn_p->Divide(h_pmiss_p);
  h_pmiss_pn_p->SetLineColor(kRed);
  h_pmiss_pn_p->Draw();
  h_pmiss_pn_p->GetYaxis()->SetTitle("pp/p and pn/p");
  h_pmiss_pn_p->SetStats(0);
  // pp/p
  TH1D * h_pmiss_pp_p = (TH1D*)h_pmiss_pp->Clone();
  h_pmiss_pp_p->Divide(h_pmiss_p);
  h_pmiss_pp_p->SetLineColor(kBlue);
  h_pmiss_pp_p->Draw("same");
  h_pmiss_pp_p->GetYaxis()->SetTitle("pp/p and pn/p");
  h_pmiss_pp_p->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // OBSERVABLE: pp/2pn
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_pp_pn = (TH1D*)h_pmiss_pp->Clone();
  h_pp_pn->Divide(h_pmiss_pn);
  h_pp_pn->Scale(0.5);
  h_pp_pn->SetLineColor(kMagenta);
  h_pp_pn->Draw();
  h_pp_pn->GetYaxis()->SetTitle("pp/2pn");
  h_pp_pn->SetStats(0);
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  // wrap it up
  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");

  //outFile->Close(); // THIS LINE CAUSES ERRORS!
}



void printProgress(double percentage) {
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

