#include <cstdlib>
#include <iostream>

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
#include <TLine.h>

#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"
#include "efficiency/efficiency.h"
#include "clas12ana.h"


using namespace std;
using namespace clas12;

double func_twolines(double *x, double *par);


void Usage() {
  std::cerr << "Usage: ./ProtonEfficiency Ebeam output-hipo output-root output-pdf input-hipo\n";
}


// constants
int peff_pbins = 16;
double p_min = 0.2;
double p_max = 1.0;
int peff_tbins = 8;
double t_min = 40;
double t_max = 140;


int main(int argc, char ** argv) {

  if(argc<5) {
    std::cerr << "Wrong number of arguments\n";
    Usage();
    return -1;
  }

  // arg 1: beam energy
  double Ebeam = atof(argv[1]);; // 2.07052 or 5.98636

  // arg 2: output hipo
  char * outName = argv[2];
  clas12writer c12writer(outName);

  // args 3-4: output file names
  TFile * f = new TFile(argv[3],"RECREATE");
  char * pdfFile = argv[4];

  // arg 5+: input hipo file
  clas12root::HipoChain chain;
  for (int k=5; k<argc; k++) {
    std::cout << "Input file " << argv[k] << std::endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader(); 
  chain.SetReaderTags({0});
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  chain.db()->turnOffQADB();

  auto currc12=chain.GetC12Reader();

  // clas12ana setup
  clas12ana clasAna;
  clasAna.readEcalSFPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsSF_LD2_x2.dat");
  clasAna.readEcalPPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsPI_LD2_x2.dat");
  clasAna.setProtonPidCuts(true);


  // prepare histograms
  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(0.8);

  char temp_name[100];
  char temp_title[100];

  double momentum;
  int event;

  int counter = 0;


//////////////////////////
/////   HISTOGRAMS   /////
//////////////////////////

  // first proton selection
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
    hist_list_1.push_back(h_psize);
  TH2D * h_pangles1 = new TH2D("pangles1","FD Proton Angles;phi;theta",180,-180,180,90,0,180);
    hist_list_2.push_back(h_pangles1);
  TH1D * h_vtz_ep = new TH1D("vtz_ep","Vertez z difference between proton and electron;Proton vertex - electron vertex;Counts",100,-6,6);
    hist_list_1.push_back(h_vtz_ep);
  TH1D * h_chipid = new TH1D("chipid","Proton #chi_{PID}^{2};#chi_{PID}^{2};Counts",100,-5,5);
    hist_list_1.push_back(h_chipid);
  TH2D * h_dbeta_p = new TH2D("dbeta_p","#Delta #beta vs proton momentum",50,0,3,50,-0.05,0.05);
    hist_list_2.push_back(h_dbeta_p);

  TH2D * h_pt_phi = new TH2D("pt_phi","Transverse Momentum vs #phi;#phi;p_{T} (GeV/c)",180,-180,180,100,0,2);
    hist_list_2.push_back(h_pt_phi);
  TH2D * h_dtof_p = new TH2D("dtof_p","#Delta TOF vs Momentum;Proton momentum (GeV/c);TOF-TOF_{expected} (ns)",100,0,2,100,-1,1);
    hist_list_2.push_back(h_dtof_p);

  // pion selection
  TH2D * h_dbeta_pi = new TH2D("dbeta_pi","#Delta #beta vs proton momentum",50,0,3,50,-0.05,0.05);
    hist_list_2.push_back(h_dbeta_pi);
  TH1D * h_dvz_pi = new TH1D("dvz_pi","Vertex z difference between pi- and e",100,-10,10);
    hist_list_1.push_back(h_dvz_pi);
  TH1D * h_pitheta = new TH1D("pitheta","#pi- Polar Angle;#theta_{#pi-};Counts",90,0,180);
    hist_list_1.push_back(h_pitheta);
  TH2D * h_piangles = new TH2D("piangles","Pi- Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_piangles);

  // denominator stuff (expected proton)
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","Missing Mass vs Missing Momentum;p_{miss} (GeV/c);Missing Mass (GeV/c^{2})",peff_pbins,p_min,p_max,20,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pmiss);
  TH2D * h_mmiss_tmiss = new TH2D("mmiss_tmiss","Missing Mass vs #theta_{miss};#theta_{miss};Missing Mass (GeV/c^{2})",peff_tbins,40,140,30,0.5,1.5);
    hist_list_2.push_back(h_mmiss_tmiss);

  TH2D * h_pt_phi_before = new TH2D("pt_phi_before","p_{pred,T} vs #phi_{pred} (before fiducial cuts);Predicted #phi (deg);Predicted p_{T} (GeV/c)",180,-180,180,50,0,2);
    hist_list_2.push_back(h_pt_phi_before);
  TH2D * h_pt_phi_after = new TH2D("pt_phi_after","p_{pred,T} vs #phi_{pred} (after fiducial cuts);Predicted #phi (deg);Predicted p_{T} (GeV/c)",180,-180,180,50,0,2);
    hist_list_2.push_back(h_pt_phi_after);

  TH2D * h_thetamiss_phimiss = new TH2D("thetamiss_phimiss","#theta_{miss} vs #phi_{miss};#phi_{miss} (degrees);#theta_{miss}",90,-180,180,90,0,180);
    hist_list_2.push_back(h_thetamiss_phimiss);

  TH2D * h_pp_ppi_denom = new TH2D("pp_ppi_denom","Proton and Pion Momenta (Denominator);Proton Momentum (GeV/c);Pion Momentum (GeV/c)",100,0,4,100,0,4);
    hist_list_2.push_back(h_pp_ppi_denom);


  



  // second proton information
  TH2D * h_pangles2 = new TH2D("pangles2","CD Proton Angles;phi;theta",180,-180,180,90,0,180);
    hist_list_2.push_back(h_pangles2);
  TH2D * h_theta_beta = new TH2D("theta_beta","Proton theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    hist_list_2.push_back(h_theta_beta);
  TH2D * h_p_theta = new TH2D("p_theta","Proton Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    hist_list_2.push_back(h_p_theta);
  TH1D * h_tof = new TH1D("tof","Time of Flight",100,-10,20);
    hist_list_1.push_back(h_tof);
  TH2D * h_Edep_beta = new TH2D("Edep_beta","Energy deposition vs #beta;#beta;E_{dep}",50,0,1,50,0,100);
    hist_list_2.push_back(h_Edep_beta);
  TH1D * h_p_all = new TH1D("p_all","Momentum",100,0,1.2);
    hist_list_1.push_back(h_p_all);


  // comparison of second proton to expected proton
  TH1D * h_pxminuspx = new TH1D("pxminuspx","px_{n}-px_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pxminuspx);
  TH1D * h_pyminuspy = new TH1D("pyminuspy","py_{n}-py_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pyminuspy);
  TH1D * h_pzminuspz = new TH1D("pzminuspz","pz_{n}-pz_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pzminuspz);
  TH1D * h_pminusp = new TH1D("pminusp","p_{n}-p_{gen};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pminusp);
  TH2D * h_pvsp = new TH2D("pvsp","Momentum Resolution;p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0.25,1.25,100,0.25,1.25);
    hist_list_2.push_back(h_pvsp);
  TH2D * h_pvsp_mcut = new TH2D("pvsp_mcut","Momentum Resolution (0.85 < M_{miss} < 1.05 GeV/c^{2});p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0.25,1.25,100,0.25,1.25);
    hist_list_2.push_back(h_pvsp_mcut);
  TH2D * h_dpp_p = new TH2D("dpp_p","#Delta p/p vs Momentum;Proton Momentum (GeV/c);#Delta p/p",50,0.3,1.25,50,-1,1);
    hist_list_2.push_back(h_dpp_p);
  TH2D * h_dpp_p_mcut = new TH2D("dpp_p_mcut","#Delta p/p vs Momentum (0.85 < M_{miss} < 1.05 GeV/c^{2});Proton Momentum (GeV/c);#Delta p/p",50,0.3,1.25,50,-1,1);
    hist_list_2.push_back(h_dpp_p_mcut);
  TH1D * h_cos0 = new TH1D("cos0","cos #theta, expected and reconstructed p",50,-1.1,1.1);
    hist_list_1.push_back(h_cos0);
  TH1D * h_cos0_mcut = new TH1D("cos0_mcut","cos #theta, expected and reconstructed p (0.85 < M_{miss} < 1.05 GeV/c^{2})",50,-1.1,1.1);
    hist_list_1.push_back(h_cos0_mcut);
  TH2D * h_pmiss_thetamiss = new TH2D("pmiss_thetamiss","pmiss vs #theta_{pmiss};#theta_{pmiss};pmiss",90,0,180,50,0,1.2);
    hist_list_2.push_back(h_pmiss_thetamiss);
  TH2D * h_selection = new TH2D("selection","(p_{miss}-p_{n})/p_{miss} vs #theta_{n,miss};(p_{miss}-p_{n})/p_{miss};#theta_{n,miss}",100,-3,1,90,0,180);
    hist_list_2.push_back(h_selection);




  // numerator
  TH1D * h_mmiss_cand = new TH1D("mmiss_cand","Missing Mass of Proton Candidates",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss_cand);
  TH1D * h_mmiss_det = new TH1D("mmiss_det","Missing Mass of Detected Protons",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss_det);
  TH2D * h_mmiss_pp = new TH2D("mmiss_pp","Missing Mass vs Measured Proton Momentum",peff_pbins,p_min,p_max,20,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pp);
  TH2D * h_mmiss_pt = new TH2D("mmiss_pt","Missing Mass vs #theta_{miss}",peff_tbins,40,140,30,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pt);
  TH2D * h_pp_ppi_numer = new TH2D("pp_ppi_numer","Proton and Pion Momenta (Numerator);Proton Momentum (GeV/c);Pion Momentum (GeV/c)",100,0,4,100,0,4);
    hist_list_2.push_back(h_pp_ppi_numer);




  // efficiency stuff
  TH1D * h_eff_p_denom = new TH1D("eff_p_denom","Proton Efficiency;Expected Proton Momentum (GeV/c);Counts",peff_pbins,p_min,p_max);
    hist_list_1.push_back(h_eff_p_denom);
  TH1D * h_eff_p_numer = new TH1D("eff_p_numer","Proton Efficiency;Expected Proton Momentum (GeV/c);Counts",peff_pbins,p_min,p_max);
    hist_list_1.push_back(h_eff_p_numer);
  //TH1D * h_eff_p = new TH1D("eff_p","Proton Efficiency;Expected Proton Momentum (GeV/c);Counts",peff_pbins,p_min,p_max);
    //hist_list_1.push_back(h_eff_p);
    
  TH1D * peff_denom_ang[4]; TH1D * peff_numer_ang[4];
  TH2D * mmiss_pmiss_CAND9[4]; TH2D * mmiss_pmiss_DET9[4];
  std::vector<double> ang_range = {45,65,75,90,120};// {45,50,55,60,65,70,75,80,95,120};
  double ang_vals[4] = {55,70,82.5,105};// {47.5,52.5,57.5,62.5,67.5,72.5,77.5,87.5,107.5};
  double xerr[4] = {10,5,7.5,15};// {2.5,2.5,2.5,2.5,2.5,2.5,2.5,7.5,12.5};

  for (int i=0; i<4; i++){
    // numerator and denominator vs pmiss
    sprintf(temp_name,"Proton Efficiency (%d-%d deg);Predicted Proton Momentum (GeV/c);Efficiency",int(ang_range[i]),int(ang_range[i+1]));
    sprintf(temp_title,"Proton Efficiency (%d-%d deg);Predicted Proton Momentum (GeV/c);Efficiency",int(ang_range[i]),int(ang_range[i+1]));
    peff_denom_ang[i] = new TH1D(temp_name,temp_title,peff_pbins,p_min,p_max);
    peff_numer_ang[i] = new TH1D(temp_name,temp_title,peff_pbins,p_min,p_max);
    hist_list_1.push_back(peff_denom_ang[i]);
    hist_list_1.push_back(peff_numer_ang[i]);
    // mmiss vs pmiss - for background subtraction
    sprintf(temp_name,"mmiss_pmiss (%d-%d deg);p_{pred} (GeV/c);M_{miss} (GeV/c^{2})",int(ang_range[i]),int(ang_range[i+1]));
    sprintf(temp_title,"M_{miss} vs p_{pred} (%d-%d deg);p_{pred} (GeV/c);M_{miss} (GeV/c^{2})",int(ang_range[i]),int(ang_range[i+1]));
    mmiss_pmiss_CAND9[i] = new TH2D(temp_name,temp_title,100,p_min,p_max,30,0.5,1.5);
    mmiss_pmiss_DET9[i] = new TH2D(temp_name,temp_title,100,p_min,p_max,30,0.5,1.5);
  }




  TH1D * h_eff_t_denom = new TH1D("eff_t_denom","Proton Efficiency;#theta_{miss};Counts",peff_tbins,t_min,t_max);
    hist_list_1.push_back(h_eff_t_denom);
  TH1D * h_eff_t_numer = new TH1D("eff_t_numer","Proton Efficiency;#theta_{miss};Counts",peff_tbins,t_min,t_max);
    hist_list_1.push_back(h_eff_t_numer);

  TH2D * h_cand2d = new TH2D("cand2d","Proton Candidates;#theta_{p};Momentum (GeV/c)",20,35,145,20,0,1.1);
    hist_list_2.push_back(h_cand2d);
  TH2D * h_det2d = new TH2D("det2d","Detected Protons;#theta_{p};Momentum (GeV/c)",20,35,145,20,0,1.1);
    hist_list_2.push_back(h_det2d);



  // extras
  TH2D * h_mmiss2_w = new TH2D("mmiss2_w","M_{miss} vs W;W;M_{miss}",50,0,6,50,0,3);
    hist_list_2.push_back(h_mmiss2_w);
  TH2D * h_mmiss2_pmiss = new TH2D("mmiss2_pmiss","M_{miss} vs p_{miss};p_{miss};M_{miss}",50,0,1.5,50,0,3);
    hist_list_2.push_back(h_mmiss2_pmiss);
  TH1D * h_mmiss2 = new TH1D("mmiss2","M_{miss} of p #pi-",50,0,4);
    hist_list_1.push_back(h_mmiss2);
  TH2D * h_mmiss2_q2 = new TH2D("mmiss2_q2","M_{miss} vs Q^2;Q^2;M_{miss}",50,0,3,50,0,3);
    hist_list_2.push_back(h_mmiss2_q2);
  TH2D * h_mmiss2_xb = new TH2D("mmiss2_xb","M_{miss} vs x_B;x_B;M_{miss}",50,0,3,50,0,3);
    hist_list_2.push_back(h_mmiss2_xb);



const double mP = 0.93828;
const double mPi = 0.13957;
const double mN = 0.939;
const double mD = 1.8756;
const double c = 29.9792458;



int numevent = 0;
  while(chain.Next())
  {

    // if multiple files in chain
    // we need to update when file changes
    if(currc12!=c12.get()){
      currc12=c12.get();
      // assign a reader to the writer
      c12writer.assignReader(*currc12);
    }

    clasAna.Run(c12);

    // identify particles from REC::Particle
    /*auto elec=c12->getByID(11);
    auto prot = c12->getByID(2212);
    auto piplus = c12->getByID(211);
    auto piminus = c12->getByID(-211);*/
    auto allParticles=c12->getDetParticles();

    auto elec = clasAna.getByPid(11);
    auto prot = clasAna.getByPid(2212);
    auto piplus = clasAna.getByPid(211);
    auto piminus = clasAna.getByPid(-211);

    if (elec.size()!=1) {continue;}
    if (prot.size()<1) {continue;}
    if (piminus.size()<1) {continue;}
    event = c12->runconfig()->getEvent() << '\n';

    // reject particles with the wrong PID
    bool trash = 0;
    for (int i=0; i<allParticles.size(); i++)
    {
      int pid = allParticles[i]->par()->getPid();
      if (pid!=2112 && pid!=11 && pid!=2212 && pid!=0 && pid!=22 && pid!=-211) {trash=1;}
    }
    if (trash==1) {continue;}

    numevent = numevent + 1;
    double starttime = c12->event()->getStartTime();


//////////////////////////
/////    ELECTRONS   /////
//////////////////////////
    TVector3 pe(0.,0.,0.);
    double pe_x = elec[0]->par()->getPx();
    double pe_y = elec[0]->par()->getPy();
    double pe_z = elec[0]->par()->getPz();
    pe.SetXYZ(pe_x,pe_y,pe_z);
    double vze = elec[0]->par()->getVz();
    TVector3 pb(0,0,Ebeam);
    TVector3 pq = pb - pe;
    double nu = Ebeam - pe.Mag();
    double QSq = pq.Mag2() - (nu*nu);
    double xB = QSq / (2*mN*nu);


//////////////////////////
/////     PROTONS    /////
//////////////////////////
    h_psize->Fill(prot.size());
    int p_index = -1;
    TVector3 pp1(0.,0.,0.);
    // technically not optimized - this doesn't address what happens if there are two protons passing cuts
    for (int i=0; i<prot.size(); i++)
    {
      // define quantities
      pp1.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
      double dbeta = prot[i]->par()->getBeta() - pp1.Mag()/sqrt(pp1.Mag2()+mP*mP);
      double vzp = prot[i]->par()->getVz();
      double chipid = prot[i]->par()->getChi2Pid();

      // angle
      double p_theta = pp1.Theta()*180./M_PI;
      h_pangles1->Fill(pp1.Phi()*180./M_PI,p_theta);
      bool p_is_FD = prot[i]->getRegion()==FD;
      if (!p_is_FD) {continue;}

      // momentum and dbeta
      h_dbeta_p->Fill(pp1.Mag(),dbeta);
      if (pp1.Mag() < 0.3) {continue;}
      if (dbeta<-0.015 || dbeta>0.015) {continue;} // make it tighter

      // vertex
      h_vtz_ep->Fill(vzp-vze);
      if ((vzp-vze)<-4. || (vzp-vze)>4.) {continue;}

      /*// chipid
      h_chipid->Fill(chipid);
      if (chipid<-3. || chipid>3.) {continue;}*/

      p_index = i;
    }
    
    if (p_index<0) {continue;}
    pp1.SetMagThetaPhi(prot[p_index]->getP(),prot[p_index]->getTheta(),prot[p_index]->getPhi());



//////////////////////////
//////    PIONS    ///////
//////////////////////////
    TVector3 ppi(0.,0.,0.);
    int pi_index = -1;
    for (int i=0; i<piminus.size(); i++)
    {

      // pion angles
      ppi.SetMagThetaPhi(piminus[i]->getP(),piminus[i]->getTheta(),piminus[i]->getPhi());
      double pi_theta = ppi.Theta()*180./M_PI;
      h_piangles->Fill(ppi.Phi()*180./M_PI,pi_theta);

      // pi theta
      h_pitheta->Fill(pi_theta);
      bool is_FD = piminus[i]->getRegion()==FD;
      if (!is_FD) {continue;} // make consistent with FD proton - theta cut or detector cut?

      // dbeta
      ppi.SetMagThetaPhi(piminus[i]->getP(),piminus[i]->getTheta(),piminus[i]->getPhi());
      double dbeta_pi = piminus[i]->par()->getBeta() - ppi.Mag()/sqrt(ppi.Mag2()+mPi*mPi);
      h_dbeta_pi->Fill(ppi.Mag(),dbeta_pi);
      if (ppi.Mag()<0.3) {continue;}
      if (dbeta_pi<-0.015|| dbeta_pi>0.015) {continue;}

      // vertex
      double vzpi = piminus[i]->par()->getVz();
      h_dvz_pi->Fill(vzpi-vze);
      if ((vzpi-vze)<-4. || (vzpi-vze)>4.) {continue;}

      // pick remaining pion
      pi_index = i;
    }
    if (pi_index<0) {continue;}

    ppi.SetMagThetaPhi(piminus[pi_index]->getP(),piminus[pi_index]->getTheta(),piminus[pi_index]->getPhi());




//////////////////////////
//  MISSING MOMENTUM    //
//////////////////////////

    // missing momentum, energy, mass
    TVector3 pmiss = pq-pp1-ppi;
    momentum = pmiss.Mag();
    double Ep = pow(mP*mP+pp1.Mag2(),0.5);
    double Epi = pow(mPi*mPi+ppi.Mag2(),0.5);
    double mmiss = pow( pow( (nu + mD - Ep - Epi), 2.) - pmiss.Mag2() , 0.5);

    // pmiss should point to CD
    h_thetamiss_phimiss->Fill(pmiss.Phi()*180./M_PI,pmiss.Theta()*180./M_PI);
    if (pmiss.Theta()*180./M_PI<45 || pmiss.Theta()*180./M_PI>120) {continue;}


    // pmiss should have correct magnitude
    if (pmiss.Mag()<p_min || pmiss.Mag()>p_max) {continue;}


    // CD fiducial cuts - doesn't make much of a difference
    double pt = sqrt( pow(pmiss.X(),2) + pow(pmiss.Y(),2) );
    double pmiss_theta = pmiss.Theta()*180./M_PI;
    double phi = pmiss.Phi()*180./M_PI;
    double min_mom_pt = 0.15;
    double cd_edge_cut = 10;
    double fiducial_phi = (-asin(min_mom_pt/pt) - M_PI/2.) * 180./M_PI;
    h_pt_phi_before->Fill(phi,pt);
    if ( (std::abs(phi-fiducial_phi) < cd_edge_cut) || (std::abs(phi-fiducial_phi-120) < cd_edge_cut) || (std::abs(phi-fiducial_phi-240) < cd_edge_cut) ) {continue;}
    h_pt_phi_after->Fill(phi,pt);



    // fill denominator
    h_mmiss_pmiss->Fill(pmiss.Mag(),mmiss);
    h_mmiss_tmiss->Fill(pmiss.Theta()*180./M_PI,mmiss);


    h_mmiss_cand->Fill(mmiss);

    // fill denominator histos
    h_eff_p_denom->Fill(pmiss.Mag());
    // fill denominator by angle
    for (int k=0; k<9; k++)
    {
      if (pmiss_theta>(ang_range[k]) && pmiss_theta<(ang_range[k+1]))
      {
        mmiss_pmiss_CAND9[k]->Fill(pmiss.Mag(),mmiss);
      }
    }
    

    if (mmiss>0.85 && mmiss<1.05) {
      h_cand2d->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
    }


    // EXTRA STUFF
    TVector3 pmiss2 = pp1+ppi;
    double mmiss2 = pow( pow(Ep+Epi,2.) - pmiss2.Mag(), 0.5);
    h_mmiss2_pmiss->Fill(pmiss.Mag(),mmiss2);
    h_mmiss2->Fill(mmiss2);

    double W2 = pow(nu + mD,2.) - pq.Mag2();
    double W = pow(W2,0.5);
    h_mmiss2_w->Fill(W,mmiss2);
    h_mmiss2_q2->Fill(QSq,mmiss);
    h_mmiss2_xb->Fill(xB,mmiss);
    h_pp_ppi_denom->Fill(pp1.Mag(),ppi.Mag());

    // write event to hipo file
    c12writer.writeEvent();


//////////////////////////
////     PROTONS     /////
//////////////////////////

    // LOOP OVER PROTONS  -- MAKE SURE WE'RE NOT LOOKING AT THE SAME PROTON AS BEFORE
    for (int i=0; i<prot.size(); i++) {

      if (i==p_index) {continue;}
    
      // GET PROTON INFORMATION

      // get proton momentum
      double pp2_x = prot[i]->par()->getPx();
      double pp2_y = prot[i]->par()->getPy();
      double pp2_z = prot[i]->par()->getPz();
      TVector3 pp2;
      pp2.SetXYZ(pp2_x,pp2_y,pp2_z);
      
      h_pt_phi->Fill(pp2.Phi()*180./M_PI, sqrt( pow(prot[i]->par()->getPx(),2) + pow(prot[i]->par()->getPy(),2) ));
      double beta_exp = pp2.Mag()/sqrt(pow(pp2.Mag(),2) + pow(mP,2) );
      h_dtof_p->Fill(pp2.Mag(), prot[i]->getPath()/c * (1/(prot[i]->par()->getBeta()) - 1/beta_exp));
       
      // put REC::Scintillator information
      double time; double energy;
      double beta = prot[i]->par()->getBeta();
      
      time =   prot[i]->sci(CTOF)->getTime() - starttime;
      energy = prot[i]->sci(CTOF)->getEnergy();

      double cos0 = pmiss.Dot(pp2) / (pmiss.Mag()*pp2.Mag());
     

      // ESSENTIAL PROTON CUTS
      h_tof->Fill(time);
      if (pp2_x==0 || pp2_y==0 || pp2_z==0) {continue;}

      // proton angular distribution
      double p_theta = pp2.Theta()*180./M_PI;
      h_pangles2->Fill(pp2.Phi()*180./M_PI,p_theta);
      // keep only Central Detector protons
      //bool is_CD = (prot[i]->getRegion()==CD);
      //if (!is_CD) {continue;}
      if (p_theta<45 || p_theta>120) {continue;}

      // proton momentum
      h_p_theta->Fill(p_theta,pp2.Mag());
      if (pp2.Mag()<p_min || pp2.Mag()>p_max) {continue;}



      //if (energy<3) {continue;}


      // histos
      h_Edep_beta->Fill(prot[i]->getBeta(),energy);
    
    h_mmiss2_w->Fill(W,mmiss2);
    h_mmiss2_q2->Fill(QSq,mmiss);
    h_mmiss2_xb->Fill(xB,mmiss);


      // FILL HISTOS FOR PROTON CANDIDATES
      h_cos0->Fill(cos0);
      h_pxminuspx->Fill(pp2_x-pmiss.X());
      h_pyminuspy->Fill(pp2_y-pmiss.Y());
      h_pzminuspz->Fill(pp2_z-pmiss.Z());
      h_pminusp->Fill(pp2.Mag()-pmiss.Mag());
      h_pvsp->Fill(pmiss.Mag(),pp2.Mag());
      h_dpp_p->Fill(pmiss.Mag(),(pmiss.Mag()-pp2.Mag())/pmiss.Mag());
      h_theta_beta->Fill(beta,p_theta);
      h_pmiss_thetamiss->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());

      h_selection->Fill((pmiss.Mag()-pp2.Mag())/pmiss.Mag(),pp2.Angle(pmiss)*180./M_PI);
      h_p_all->Fill(pmiss.Mag());

      if (mmiss>0.85 && mmiss<1.05) {
        h_cos0_mcut->Fill(cos0);
        h_pvsp_mcut->Fill(pmiss.Mag(),pp2.Mag());
        h_dpp_p_mcut->Fill(pmiss.Mag(),(pmiss.Mag()-pp2.Mag())/pmiss.Mag());
      }


      // pick good protons
      //if (pp2.Angle(pmiss)*180./M_PI>20) {continue;}
      //if (cos0<0.8) {continue;}
      //if (abs((pmiss.Mag()-pp2.Mag())/pmiss.Mag())>0.2) {continue;}


      h_mmiss_det->Fill(mmiss);
      //h_mmiss_pp->Fill(pp2.Mag(),mmiss);
      h_mmiss_pp->Fill(pmiss.Mag(),mmiss);
      h_mmiss_pt->Fill(pmiss.Theta()*180./M_PI,mmiss);
      h_pp_ppi_numer->Fill(pp1.Mag(),ppi.Mag());


      // fill numerator histos
      h_eff_p_numer->Fill(pmiss.Mag());
      // fill numerator by angle
      for (int k=0; k<9; k++)
      {
        if (pmiss_theta>(ang_range[k]) && pmiss_theta<(ang_range[k+1]))
        {
          mmiss_pmiss_DET9[k]->Fill(pmiss.Mag(),mmiss);
        }
      }


      if (mmiss>0.85 && mmiss<1.05) {
        h_det2d->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
      }


  }  // closes proton loop

    counter++;

  }  // closes event loop


  c12writer.closeWriter();


  f->cd();
  for(int i=0; i<hist_list_1.size(); i++) {
    hist_list_1[i]->Sumw2();
    hist_list_1[i]->GetXaxis()->CenterTitle();
    hist_list_1[i]->GetYaxis()->CenterTitle();
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++) {
    hist_list_2[i]->Sumw2();
    hist_list_2[i]->GetXaxis()->CenterTitle();
    hist_list_2[i]->GetYaxis()->CenterTitle();
    hist_list_2[i]->SetOption("colz");
    hist_list_2[i]->Write();
  }

  std::cout << '\n' <<counter << " events counted!\n\n";


  ////////////////////////////////////////
  //Create the output PDF
  ////////////////////////////////////////
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


  ////////////////////////////////////////
  //Electron kinematics
  ////////////////////////////////////////
  


  ////////////////////////////////////////
  //FD proton kinematics
  ////////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.8,"FD Proton Selection");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pangles1->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dbeta_p->Draw("colz");
  TLine * l_p_dbeta1 = new TLine(0.3,0.015,3,0.015);
  l_p_dbeta1->SetLineColor(kRed);
  l_p_dbeta1->SetLineWidth(3);
  l_p_dbeta1->Draw("same");
  TLine * l_p_dbeta2 = new TLine(0.3,-0.015,3,-0.015);
  l_p_dbeta2->SetLineColor(kRed);
  l_p_dbeta2->SetLineWidth(3);
  l_p_dbeta2->Draw("same");
  TLine * l_p_dbeta3 = new TLine(0.3,-0.015,0.3,0.015);
  l_p_dbeta3->SetLineColor(kRed);
  l_p_dbeta3->SetLineWidth(3);
  l_p_dbeta3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_vtz_ep->Draw();
  TLine * l_p_vtz1 = new TLine(-4,0,-4,h_vtz_ep->GetMaximum());
  l_p_vtz1->SetLineColor(kRed);
  l_p_vtz1->SetLineWidth(3);
  l_p_vtz1->Draw("same");
  TLine * l_p_vtz2 = new TLine(4,0,4,h_vtz_ep->GetMaximum());
  l_p_vtz2->SetLineColor(kRed);
  l_p_vtz2->SetLineWidth(3);
  l_p_vtz2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_chipid->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  ////////////////////////////////////////
  //Pion kinematics
  ////////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.8,"FD Pion Selection");
  myText->Print(fileName,"pdf");
  myText->Clear(); 

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_piangles->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pitheta->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dbeta_pi->Draw("colz");
  TLine * l_pi_dbeta1 = new TLine(0.3,0.015,3,0.015);
  l_pi_dbeta1->SetLineColor(kRed);
  l_pi_dbeta1->SetLineWidth(3);
  l_pi_dbeta1->Draw("same");
  TLine * l_pi_dbeta2 = new TLine(0.3,-0.015,3,-0.015);
  l_pi_dbeta2->SetLineColor(kRed);
  l_pi_dbeta2->SetLineWidth(3);
  l_pi_dbeta2->Draw("same");
  TLine * l_pi_dbeta3 = new TLine(0.3,-0.015,0.3,0.015);
  l_pi_dbeta3->SetLineColor(kRed);
  l_pi_dbeta3->SetLineWidth(3);
  l_pi_dbeta3->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dvz_pi->Draw();
  TLine * l_pi_vtz1 = new TLine(-4,0,-4,h_vtz_ep->GetMaximum());
  l_pi_vtz1->SetLineColor(kRed);
  l_pi_vtz1->SetLineWidth(3);
  l_pi_vtz1->Draw("same");
  TLine * l_pi_vtz2 = new TLine(4,0,4,h_vtz_ep->GetMaximum());
  l_pi_vtz2->SetLineColor(kRed);
  l_pi_vtz2->SetLineWidth(3);
  l_pi_vtz2->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();






  ////////////////////////////////////////
  //Missing momentum / proton candidates
  ////////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.8,"Expected Proton");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_thetamiss_phimiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pt_phi_before->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pt_phi_after->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pp_ppi_denom->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // background subtraction - p
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_pmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(4,4);
  myCanvas->cd(1);
  double * Scand = hist_projections_backsub(myCanvas, h_mmiss_pmiss, peff_pbins, 1, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // background subtraction - theta
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_tmiss->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(4,4);
  myCanvas->cd(1);
  double * Scand_t = hist_projections_backsub(myCanvas, h_mmiss_tmiss, peff_tbins, 1, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_cand->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cand2d->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  ////////////////////////////////////////
  //CD proton kinematics
  ////////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.8,"Detected CD Proton");
  myText->Print(fileName,"pdf");
  myText->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_tof->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pangles2->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pt_phi->Draw("colz");
  myCanvas->cd(2);
  h_dtof_p->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_p_theta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_mmiss_det->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_Edep_beta->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cos0->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_cos0_mcut->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(2,2);
  myCanvas->cd(1);
  h_pxminuspx->Draw();
  myCanvas->cd(2);
  h_pyminuspy->Draw();
  myCanvas->cd(3);
  h_pzminuspz->Draw();
  myCanvas->cd(4);
  h_pminusp->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pvsp->Draw("colz");
  TLine * l_pvsp = new TLine(0,1.2,0,1.2);
  l_pvsp->SetLineColor(kRed);
  l_pvsp->SetLineWidth(3);
  l_pvsp->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pvsp_mcut->Draw("colz");
  TLine * l_pvsp_mcut = new TLine(0,1.2,0,1.2);
  l_pvsp_mcut->SetLineColor(kRed);
  l_pvsp_mcut->SetLineWidth(3);
  l_pvsp_mcut->Draw("same");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dpp_p->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_dpp_p_mcut->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_selection->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_pp_ppi_numer->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // background subtraction - p
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_pp->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(4,4);
  myCanvas->cd(1);
  double * Sdet = hist_projections_backsub(myCanvas, h_mmiss_pp, peff_pbins, 1, 'p');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  // background subtraction - theta
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  myCanvas->cd(1)->SetLogz();
  h_mmiss_pt->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(4,4);
  myCanvas->cd(1);
  double * Sdet_t = hist_projections_backsub(myCanvas, h_mmiss_pt, peff_tbins, 1, 't');
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_det2d->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



  ////////////////////////////////////////
  //Efficiency
  ////////////////////////////////////////
  myText->cd();
  text.DrawLatex(0.2,0.8,"Proton Efficiency");
  myText->Print(fileName,"pdf");
  myText->Clear();
  
  for (int i=0; i<peff_pbins; i++){
    // numerator - background subtraction
    h_eff_p_numer->SetBinContent(i+1,*(Sdet+i));
    h_eff_p_numer->SetBinError(i+1,std::sqrt(*(Sdet+i)));
    // denominator - background subtraction
    h_eff_p_denom->SetBinContent(i+1,*(Scand+i));
    h_eff_p_denom->SetBinError(i+1,std::sqrt(*(Scand+i)));
  }

  for (int i=0; i<peff_pbins; i++){
    // numerator - background subtraction
    h_eff_t_numer->SetBinContent(i+1,*(Sdet_t+i));
    h_eff_t_numer->SetBinError(i+1,std::sqrt(*(Sdet_t+i)));
    // denominator - background subtraction
    h_eff_t_denom->SetBinContent(i+1,*(Scand_t+i));
    h_eff_t_denom->SetBinError(i+1,std::sqrt(*(Scand_t+i)));
  }

  // efficiency - p
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_eff_p_denom->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_eff_p_numer->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
f->cd();
  TH1D * h_eff_p = (TH1D*)h_eff_p_numer->Clone();
  h_eff_p->SetName("efficiency_p");
  h_eff_p->SetDirectory(0);
  h_eff_p->Divide(h_eff_p_denom);
  h_eff_p->SetStats(0);
  h_eff_p->GetYaxis()->SetTitle("efficiency");
  h_eff_p->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();
hist_list_1.push_back(h_eff_p);
  //h_eff_p->Write();
  

  // fit for proton efficiency: a*sin(bx)+cx
  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TF1 * sinefun = new TF1("sinefun","[0]*sin([1]*x-[2]) + [4]*x + [5]",0.2,1.0);
  TF1 * twogauss = new TF1("twogauss","[0]*exp(-pow((x-[1]),2.)/(2*pow([2],2.))) + [3]*exp(-pow((x-[4]),2.)/(2*pow([5],2.)))",0.2,1.0);
  TF1 * sqrtfun = new TF1("sqrtfun","[0]*sqrt([1]*x-[2])",0.2,1.0);
  h_eff_p->Draw();
  //sinefun->SetParameter(0,0.1);  sinefun->SetParameter(1,12);  sinefun->SetParameter(2,1.3);  sinefun->SetParameter(4,0.9);  sinefun->SetParameter(5,0.1);
  //h_eff_p->Fit("sinefun","","",0.2,1.0);
  //twogauss->SetParameter(0,0.7); twogauss->SetParameter(1,0.45); twogauss->SetParameter(2,0.05);
  //twogauss->SetParameter(3,1);   twogauss->SetParameter(4,0.85); twogauss->SetParameter(5,0.05);
  //h_eff_p->Fit("twogauss","","",0.2,1.0);
  //sqrtfun->SetParameter(0,0.9);  sqrtfun->SetParameter(1,1.2);  sqrtfun->SetParameter(2,0.3);
  //sqrtfun->SetParLimits(0,0.5,1.2); sqrtfun->SetParLimits(1,1,2); sqrtfun->SetParLimits(2,0,1);
  //h_eff_p->Fit("sqrtfun","","",0.2,1.0);

  //TF1 * twolines = new TF1("twolines",func_twolines,0.2,1.0,3);
  //h_eff_p->Fit("twolines","QN");

  //TF1 * linfun = new TF1("poly1",0.2,1.0);
  //h_eff_p->Fit("linfun","","",0.2,1.0);
  //h_eff_p->Fit("pol1");


  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  for (int i=0; i<16; i++) {std::cout << h_eff_p->GetBinContent(i) << ", ";}
  


  // fit for proton efficiency
  




  // efficiency by angular range
  TH1D * h_eff_p_denom_test = new TH1D("test","test",peff_pbins,p_min,p_max);
  TH1D * h_eff_p_numer_test = new TH1D("test","test",peff_pbins,p_min,p_max); 
  for (int i=0; i<4; i++)
  {
    myText->cd();
    text.DrawLatex(0.2,0.8,"next angular range: candidates, then detected");
    myText->Print(fileName,"pdf");
    myText->Clear();
    // background subtraction for candidates
    myCanvas->Divide(4,4);
    double *Ssub_pCAND00 = hist_projections_backsub(myCanvas, mmiss_pmiss_CAND9[i], peff_pbins, 1, 'p');
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();

    // fill denominator with background-subtracted values
    for (int j=0; j<peff_pbins; j++)
    {
      peff_denom_ang[i]->SetBinContent(j,*(Ssub_pCAND00+j));
      peff_denom_ang[i]->SetBinError(j,std::sqrt(*(Ssub_pCAND00+j)));
    }

    // background subtraction for detected neutrons
    myCanvas->Divide(4,4);
    double *Ssub_pDET00 = hist_projections_backsub(myCanvas, mmiss_pmiss_DET9[i], peff_pbins, 1, 'p');
    myCanvas->Print(fileName,"pdf");
    myCanvas->Clear();
    
    // fill numerator with background-subtracted values
    for (int j=0; j<peff_pbins; j++)
    {
      peff_numer_ang[i]->SetBinContent(j,*(Ssub_pDET00+j));
      peff_numer_ang[i]->SetBinError(j,std::sqrt(*(Ssub_pDET00+j)));


h_eff_p_denom_test->Add(peff_denom_ang[i]);
h_eff_p_numer_test->Add(peff_numer_ang[i]);
    }
    
    } // end loop over 9 angular ranges






  // draw peff for 9 angular ranges
  TH1D * peff_ang[4];
  double par[3], p0[9], p0_err[9], p1[9], p1_err[9], p2[9], p2_err[9], ang[9];
  myCanvas->Divide(3,3);
  TF1 * f_poly2 = new TF1("poly2","[0]+[1]*x+[2]*x*x",0.3,1.0);
  myCanvas->cd(1);
  for (int i=0; i<4; i++)
  {
    myCanvas->cd(i+1);
    // draw peff for angular range
    peff_ang[i] = (TH1D*)peff_numer_ang[i]->Clone();
    peff_ang[i]->Divide(peff_denom_ang[i]);
    peff_ang[i]->SetMinimum(0.0);
    peff_ang[i]->SetMaximum(1.0);
    peff_ang[i]->SetLineColor(kBlue-i);
    peff_ang[i]->Draw();
  }
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();



myCanvas->Divide(2,2);
myCanvas->cd(1);
h_eff_p_denom->Draw();
myCanvas->cd(2);
h_eff_p_numer->Draw();
myCanvas->Print(fileName,"pdf");
myCanvas->Clear();

/*myCanvas->Divide(1,1);
myCanvas->cd(1);
h_eff_p_numer_test->Divide(h_eff_p_denom);
h_eff_p_numer_test->Draw();
myCanvas->Print(fileName,"pdf");
myCanvas->Clear();*/





  // efficiency - theta
/*  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_eff_t_denom->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  h_eff_t_numer->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH1D * h_peff_t = (TH1D*)h_eff_t_numer->Clone();
  h_peff_t->Divide(h_eff_t_denom);
  h_peff_t->SetStats(0);
  h_peff_t->Draw();
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();*/



  myCanvas->Divide(1,1);
  myCanvas->cd(1);
  TH2D * h_peff_2d = (TH2D*)h_det2d->Clone();
  h_peff_2d->Divide(h_cand2d);
  h_peff_2d->SetStats(0);
  h_peff_2d->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  myCanvas->Divide(1,1);
  myCanvas->cd(11);
  TH2D * h_mmiss_peff = (TH2D*)h_mmiss_pp->Clone();
  h_mmiss_peff->Divide(h_mmiss_pmiss);
  h_mmiss_peff->SetStats(0);
  h_mmiss_peff->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();

  myCanvas->Divide(1,1);
  myCanvas->cd(11);
  TH2D * h_mmiss_teff = (TH2D*)h_mmiss_pt->Clone();
  h_mmiss_teff->Divide(h_mmiss_tmiss);
  h_mmiss_teff->SetStats(0);
  h_mmiss_teff->Draw("colz");
  myCanvas->Print(fileName,"pdf");
  myCanvas->Clear();


  sprintf(fileName,"%s]",pdfFile);
  myCanvas->Print(fileName,"pdf");



  //f->Close();

  return 0;

}  // closes main function



double func_twolines(double *x, double *par)
{
  double xx = x[0];
  double value = 0;
  if (xx<par[0]) {value = par[1]*xx + par[2];}
  else if (xx>par[0]) {value = par[1]*par[0] + par[2];}

  return value;
}
