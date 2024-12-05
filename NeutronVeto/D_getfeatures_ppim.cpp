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

#include "clas12reader.h"
#include "HipoChain.h"
#include "neutron-veto/veto_functions.h"
#include "clas12ana.h"


using namespace std;
using namespace clas12;

void Usage() {
  std::cerr << "Usage: ./D_getfeatures Ebeam keep_good output-root output-txt input-hipo\n";
}

double getCVTdiff(std::vector<region_part_ptr> &allParticles, TVector3 &pn);

int main(int argc, char ** argv) {

  if(argc<6) {
    std::cerr << "Wrong number of arguments\n";
    Usage();
    return -1;
  }

  // arg 1: beam energy
  double Ebeam = atof(argv[1]);; // 2.07052 or 5.98636

  // arg 2: keep good
  bool keep_good = false;
  if(atoi(argv[2])==1){keep_good=true;}

  // args 3-4: output file names
  TFile * f = new TFile(argv[3],"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");
  std::ofstream outtxt(argv[4]);

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

  // prepare histograms
  vector<TH1*> hist_list_1;
  vector<TH2*> hist_list_2;

  gStyle->SetTitleXSize(0.05);
  gStyle->SetTitleYSize(0.05);

  gStyle->SetTitleXOffset(0.8);
  gStyle->SetTitleYOffset(0.8);

  char temp_name[100];
  char temp_title[100];


  // Set up root tree for TMVA
  Int_t nhits;
  double px, py, pz, momentum;
  Int_t sec[100] = {-1};
  Int_t lay[100] = {-1};
  int event;
  double energy, cnd_energy, ctof_energy, angle_diff;
  int layermult, size, cnd_hits, ctof_hits;
  ntree->Branch("momentum",&momentum,"momentum/D");
  ntree->Branch("energy",&energy,"energy/D");
  ntree->Branch("layermult",&layermult,"layermult/I");
  ntree->Branch("size",&size,"size/I");
  ntree->Branch("cnd_hits",&cnd_hits,"cnd_hits/I");
  ntree->Branch("cnd_energy",&cnd_energy,"cnd_energy/D");
  ntree->Branch("ctof_energy",&ctof_energy,"ctof_energy/D");
  ntree->Branch("ctof_hits",&ctof_hits,"ctof_hits/I");
  ntree->Branch("angle_diff",&angle_diff,"angle_diff/D");

  int counter = 0;

  // set up instance of clas12ana
  clas12ana * clasAna = new clas12ana();

  clasAna->readEcalSFPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsSF_LD2_x2.dat");
  clasAna->readEcalPPar("/w/hallb-scshelf2102/clas12/users/esteejus/rgm/Ana/cutFiles/paramsPI_LD2_x2.dat");

  clasAna->printParams();

  clasAna->setProtonPidCuts(true);


//////////////////////////
/////   HISTOGRAMS   /////
//////////////////////////

  // proton stuff
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
    hist_list_1.push_back(h_psize);
  TH2D * h_pangles = new TH2D("pangles","Proton Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_pangles);
  TH1D * h_vzp = new TH1D("vzp","Vertex difference between proton and electron;z_{p} - z_{e} (cm);Counts",100,-8,8);
    hist_list_1.push_back(h_vzp);
  TH2D * h_dbeta_p = new TH2D("dbeta_p","#Delta #beta_{p} vs proton momentum;p_{p} (GeV/c);#Delta#beta_{p}",50,0,3,50,-0.2,0.2);
    hist_list_2.push_back(h_dbeta_p);

  // pion stuff
  TH2D * h_dbeta_pi = new TH2D("dbeta_pi","#Delta #beta_{#pi} vs pion momentum;p_{#pi} (GeV/c);#Delta#beta_{#pi}",50,0,3,50,-0.2,0.2);
    hist_list_2.push_back(h_dbeta_pi);
  TH1D * h_dvz_pi = new TH1D("dvz_pi","Vertex z difference between #pi^{-} and e;z_{#pi}-z_{e} (cm)",100,-10,10);
    hist_list_1.push_back(h_dvz_pi);
  TH2D * h_piangles = new TH2D("piangles","Pion Angular Distribution;#phi_{#pi} (deg);#theta_{#pi} (deg)",48,-180,180,45,0,180);
    hist_list_2.push_back(h_piangles);


  // neutron stuff
  TH1D * h_nsize = new TH1D("nsize","Number of Neutrons in Event",10,0,10);
    hist_list_1.push_back(h_nsize);


  // reconstructed momentum
  TH2D * h_nangles = new TH2D("nangles","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_nangles);
  TH1D * h_pxminuspx = new TH1D("pxminuspx","px_{n}-px_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pxminuspx);
  TH1D * h_pyminuspy = new TH1D("pyminuspy","py_{n}-py_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pyminuspy);
  TH1D * h_pzminuspz = new TH1D("pzminuspz","pz_{n}-pz_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pzminuspz);
  TH1D * h_pminusp = new TH1D("pminusp","p_{n}-p_{gen};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pminusp);
  TH2D * h_pvsp = new TH2D("pvsp","Momentum Resolution;p_{pred} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    hist_list_2.push_back(h_pvsp);
  TH1D * h_cos0 = new TH1D("cos0","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
    hist_list_1.push_back(h_cos0);
  TH2D * h_dpp = new TH2D("dpp","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    hist_list_2.push_back(h_dpp);
  TH1D * h_energy = new TH1D("energy","Neutron Energy Deposition;Energy (MeV);Counts",100,0,25);
    hist_list_1.push_back(h_energy);
  TH2D * h_sec_phi = new TH2D("sec_phi","Sector vs Phi of CND hits;phi (deg);Sector",90,0,360,25,0,25);
    hist_list_2.push_back(h_sec_phi);
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss);
  TH2D * h_mmiss_pn = new TH2D("mmiss_pn","Missing Mass vs Measured Neutron Momentum;p_{n} (GeV/c);M_{miss} (GeV/c^{2})",50,0.25,1.,50,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pn);
  TH2D * h_mmiss_pmiss = new TH2D("mmiss_pmiss","Missing Mass vs Predicted Neutron Momentum",50,0,1.5,50,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pmiss);

  TH2D * h_theta_beta = new TH2D("theta_beta","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    hist_list_2.push_back(h_theta_beta);
  TH2D * h_p_theta = new TH2D("p_theta","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    hist_list_2.push_back(h_p_theta);
  TH2D * h_pmiss_thetamiss = new TH2D("pmiss_thetamiss","p_{pred} vs #theta_{pred};#theta_{pred};p_{pred} (GeV/c)",90,0,180,50,0,1.3);
    hist_list_2.push_back(h_pmiss_thetamiss);
  TH2D * h_thetapn_pp = new TH2D("thetapn_pp","#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}",40,0,1,40,0,180);
    hist_list_2.push_back(h_thetapn_pp);
  TH1D * h_tof = new TH1D("tof","Time of Flight;TOF (ns)",100,-10,20);
    hist_list_1.push_back(h_tof);
  TH2D * h_compare = new TH2D("compare","(p_{pred}-p_{n})/p_{pred} vs #theta_{n,pred};(p_{pred}-p_{n})/p_{pred};#theta_{n,pred}",100,-3,1,90,0,180);
    hist_list_2.push_back(h_compare);
  TH2D * h_Edep_beta = new TH2D("Edep_beta","Energy deposition vs #beta;#beta;E_{dep}",50,0,1,50,0,100);
    hist_list_2.push_back(h_Edep_beta);
  TH1D * h_p_all = new TH1D("p_all","Momentum",100,0,1.2);
    hist_list_1.push_back(h_p_all);



  // good n / bad n set
  TH2D * h_nangles2 = new TH2D("nangles2","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_nangles2);
  TH1D * h_pxminuspx2 = new TH1D("pxminuspx2","px_{n}-px_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pxminuspx2);
  TH1D * h_pyminuspy2 = new TH1D("pyminuspy2","py_{n}-py_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pyminuspy2);
  TH1D * h_pzminuspz2 = new TH1D("pzminuspz2","pz_{n}-pz_{pred};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pzminuspz2);
  TH1D * h_pminusp2 = new TH1D("pminusp2","p_{n}-p_{gen};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pminusp2);
  TH2D * h_pvsp2 = new TH2D("pvsp2","Momentum Resolution;p_{pred} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    hist_list_2.push_back(h_pvsp2);
  TH1D * h_cos02 = new TH1D("cos02","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
    hist_list_1.push_back(h_cos02);
  TH2D * h_dpp2 = new TH2D("dpp2","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    hist_list_2.push_back(h_dpp2);
  TH1D * h_mmiss2 = new TH1D("mmiss2","Missing Mass",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss2);
  TH2D * h_mmiss_pn2 = new TH2D("mmiss_pn2","Missing Mass vs Neutron Momentum",50,0,1.5,50,0.5,1.5);
    hist_list_2.push_back(h_mmiss_pn2);
  TH1D * h_energy2 = new TH1D("energy2","Neutron Energy Deposition;Energy (MeV);Counts",100,0,25);
    hist_list_1.push_back(h_energy2);
  TH2D * h_theta_beta2 = new TH2D("theta_beta2","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    hist_list_2.push_back(h_theta_beta2);
  TH2D * h_p_theta2 = new TH2D("p_theta2","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    hist_list_2.push_back(h_p_theta2);
  TH2D * h_pmiss_thetamiss2 = new TH2D("pmiss_thetamiss2","p_{pred} vs #theta_{pred};#theta_{pred};p_{pred} (GeV/c)",90,0,180,50,0,1.2);
    hist_list_2.push_back(h_pmiss_thetamiss2);
  TH2D * h_thetapn_pp2 = new TH2D("thetapn_pp2","#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}",40,0,1,40,0,180);
    hist_list_2.push_back(h_thetapn_pp2);
  TH1D * h_tof2 = new TH1D("tof2","Time of Flight;TOF (ns)",100,-10,20);
    hist_list_1.push_back(h_tof2);
  TH2D * h_compare2 = new TH2D("compare2","(p_{pred}-p_{n})/p_{pred} vs #theta_{n,pred};(p_{pred}-p_{n})/p_{pred};#theta_{n,pred}",100,-3,1,90,0,180);
    hist_list_2.push_back(h_compare2);
  TH2D * h_Edep_beta2 = new TH2D("Edep_beta2","Energy deposition vs #beta;#beta;E_{dep}",50,0,1,50,0,100);
    hist_list_2.push_back(h_Edep_beta2);
  TH1D * h_p_cut = new TH1D("p_cut","Momentum",100,0,1.2);
    hist_list_1.push_back(h_p_cut);

  TH2D * h_ptheta_pred = new TH2D("ptheta_pred","Predicted Momentum vs Angle of Final Background Sample;#theta_{pred} (deg);p_{pred} (GeV/c)",110,35,145,100,0.2,1.3);
    hist_list_2.push_back(h_ptheta_pred);
  TH2D * h_ptheta = new TH2D("ptheta","Measured Momentum vs Angle of Final Background Sample;#theta_{p} (deg);p_{p} (GeV/c)",110,35,145,100,0.2,1.3);
    hist_list_2.push_back(h_ptheta);


  // ML features
  TH1D * h_energy_1 = new TH1D("f_energy_1","Neutron Energy",100,0,100);
    hist_list_1.push_back(h_energy_1);
  TH1D * h_layermult_1 = new TH1D("f_layermult_1","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_1);
  TH1D * h_size_1 = new TH1D("f_size_1","Cluster Size",5,0,5);
    hist_list_1.push_back(h_size_1);
  TH1D * h_cnd_hits_1 = new TH1D("f_cnd_hits_1","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_1);
  TH1D * h_cnd_energy_1 = new TH1D("f_cnd_energy_1","Nearby CND Energy",100,0,100);
    hist_list_1.push_back(h_cnd_energy_1);
  TH1D * h_ctof_energy_1 = new TH1D("f_ctof_energy_1","Nearby CTOF Energy",100,0,100);
    hist_list_1.push_back(h_ctof_energy_1);
  TH1D * h_ctof_hits_1 = new TH1D("f_ctof_hits_1","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_1);
  TH1D * h_anglediff_1 = new TH1D("f_anglediff_1","CVT Angle Diff",200,0,200);
    hist_list_1.push_back(h_anglediff_1);

  TH1D * h_energy_2 = new TH1D("f_energy_2","Neutron Energy",100,0,100);
    hist_list_1.push_back(h_energy_2);
  TH1D * h_layermult_2 = new TH1D("f_layermult_2","CND Layer Mult",4,0,4);
    hist_list_1.push_back(h_layermult_2);
  TH1D * h_size_2 = new TH1D("f_size_2","Cluster Size",5,0,5);
    hist_list_1.push_back(h_size_2);
  TH1D * h_cnd_hits_2 = new TH1D("f_cnd_hits_2","Nearby CND Hits",10,0,10);
    hist_list_1.push_back(h_cnd_hits_2);
  TH1D * h_cnd_energy_2 = new TH1D("f_cnd_energy_2","Nearby CND Energy",100,0,100);
    hist_list_1.push_back(h_cnd_energy_2);
  TH1D * h_ctof_energy_2 = new TH1D("f_ctof_energy_2","Nearby CTOF Energy",100,0,100);
    hist_list_1.push_back(h_ctof_energy_2);
  TH1D * h_ctof_hits_2 = new TH1D("f_ctof_hits_2","Nearby CTOF Hits",10,0,10);
    hist_list_1.push_back(h_ctof_hits_2);
  TH1D * h_anglediff_2 = new TH1D("f_anglediff_2","CVT Angle Diff",200,0,200);
    hist_list_1.push_back(h_anglediff_2);




const double mP = 0.93828;
const double mPi = 0.13957;
const double mN = 0.939;
const double mD = 1.8756;



int numevent = 0;
  while(chain.Next())
  {

    // initialize features
    energy = 0; cnd_energy = 0; ctof_energy = 0; angle_diff = 180;
    layermult = 0; size = 0; cnd_hits = 0; ctof_hits = 0;

    // particle pid
    clasAna->Run(c12);
    auto elec = clasAna->getByPid(11);
    auto prot = clasAna->getByPid(2212);
    auto neut = clasAna->getByPid(2112);
    auto piplus = clasAna->getByPid(211);
    auto piminus = clasAna->getByPid(-211);
    auto allParticles=c12->getDetParticles();
    if (elec.size()!=1) {continue;}
    if (prot.size()<1) {continue;}
    if (neut.size()<1) {continue;}
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
    TVector3 pp(0.,0.,0.);
    // technically not optimized - this doesn't address what happens if there are two protons passing cuts
    for (int i=0; i<prot.size(); i++) // technically not optimized - 
    {
      // define quantities
      pp.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
      double dbeta = prot[i]->par()->getBeta() - pp.Mag()/sqrt(pp.Mag2()+mP*mP);
      double vzp = prot[i]->par()->getVz();
      double chipid = prot[i]->par()->getChi2Pid();


      // make cuts
      h_dbeta_p->Fill(pp.Mag(),dbeta);
      if (pp.Mag()<0.5) {continue;}
      if (pp.Mag()>3.0) {continue;}
      if (abs(dbeta)>0.05) {continue;}

      h_vzp->Fill(vzp-vze);
      if (abs(vzp-vze)>5.) {continue;}
      //if (chipid<-3. || chipid>3.) {continue;}
      
      double p_theta = pp.Theta()*180./M_PI;
      h_pangles->Fill(pp.Phi()*180./M_PI,p_theta);
      //if (p_theta>40) {continue;}

      p_index = i;
    }
    
    if (p_index<0) {continue;}
    pp.SetMagThetaPhi(prot[p_index]->getP(),prot[p_index]->getTheta(),prot[p_index]->getPhi());
    double p_theta = pp.Theta()*180./M_PI;



//////////////////////////
//////    PIONS    ///////
//////////////////////////
    TVector3 ppi(0.,0.,0.);
    int pi_index = -1;
    for (int i=0; i<piminus.size(); i++)
    {
      ppi.SetMagThetaPhi(piminus[i]->getP(),piminus[i]->getTheta(),piminus[i]->getPhi());

      double vzpi = piminus[i]->par()->getVz();
      h_dvz_pi->Fill(vzpi-vze);
      if ((vzpi-vze)<-4. || (vzpi-vze)>4.) {continue;}

      double dbeta_pi = piminus[i]->par()->getBeta() - ppi.Mag()/sqrt(ppi.Mag2()+mPi*mPi);
      h_dbeta_pi->Fill(ppi.Mag(),dbeta_pi);
      if (ppi.Mag()<0.5) {continue;}
      if (ppi.Mag()>3.0) {continue;}
      if (abs(dbeta_pi)>0.05) {continue;}

      double pi_theta = ppi.Theta()*180./M_PI;
      h_piangles->Fill(ppi.Phi()*180./M_PI,pi_theta);
      //if (pi_theta>40) {continue;}

      pi_index = i;
    }
    if (pi_index<0) {continue;}

    ppi.SetMagThetaPhi(piminus[pi_index]->getP(),piminus[pi_index]->getTheta(),piminus[pi_index]->getPhi());
    double pi_theta = ppi.Theta()*180./M_PI;


    //if (p_theta<40 && pi_theta<40) {continue;}



//////////////////////////
//  MISSING MOMENTUM    //
//////////////////////////

    // missing momentum, energy, mass
    TVector3 pmiss = pq-pp-ppi;
    momentum = pmiss.Mag();
    double Ep = pow(mP*mP+pp.Mag2(),0.5);
    double Epi = pow(mPi*mPi+ppi.Mag2(),0.5);
    double mmiss = pow( pow( (nu + mD - Ep - Epi), 2.) - pmiss.Mag2() , 0.5);


//////////////////////////
////     NEUTRONS    /////
//////////////////////////

    // LOOP OVER NEUTRONS
    h_nsize->Fill(neut.size());
    for (int i=0; i<neut.size(); i++) {
    
      // GET NEUTRON INFORMATION

      // get neutron momentum
      double pn_x = neut[i]->par()->getPx();
      double pn_y = neut[i]->par()->getPy();
      double pn_z = neut[i]->par()->getPz();
      TVector3 pn;
      pn.SetXYZ(pn_x,pn_y,pn_z);
      
      // figure out what layer the hit is in
      bool is_CND1 = (neut[i]->sci(CND1)->getLayer()==1);
      bool is_CND2 = (neut[i]->sci(CND2)->getLayer()==2);
      bool is_CND3 = (neut[i]->sci(CND3)->getLayer()==3);
      bool is_CTOF = (neut[i]->sci(CTOF)->getDetector()==4); 
      if (!is_CND1 && !is_CND2 && !is_CND3 && !is_CTOF) {continue;}
       
      // put REC::Scintillator information
      int sector; int status = 0;
      double time;
      double beta = neut[i]->par()->getBeta();
      
      if (is_CND1)
      {
        sector = neut[i]->sci(CND1)->getSector();
        time =   neut[i]->sci(CND1)->getTime() - starttime;
        status = status + neut[i]->sci(CND1)->getStatus();
      }
      
      if (is_CND3)
      {
        sector = neut[i]->sci(CND3)->getSector();
        time =   neut[i]->sci(CND3)->getTime() - starttime;
        status = status + neut[i]->sci(CND3)->getStatus();
      }
      
      if (is_CND2)
      {
        sector = neut[i]->sci(CND2)->getSector();
        time =   neut[i]->sci(CND2)->getTime() - starttime;
        status = status + neut[i]->sci(CND2)->getStatus();
      }
      // PROBLEM: this gives preference to 2nd-layer hits
      if (is_CTOF)
      {
        sector = (neut[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
        time =   neut[i]->sci(CTOF)->getTime() - starttime;
      }
      if (status!=0) {continue;}

      // GET ML FEATURES FOR THIS NEUTRON
      Struct ninfo = getFeatures(neut, allParticles, i);
      cnd_hits = ninfo.cnd_hits;
      ctof_hits = ninfo.ctof_hits;
      cnd_energy = ninfo.cnd_energy;
      ctof_energy = ninfo.ctof_energy;
      layermult = ninfo.layermult;
      energy = ninfo.energy;
      size = ninfo.size;
      angle_diff = ninfo.angle_diff;


      double cos0 = pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag());
     
      // ESSENTIAL NEUTRONS CUTS
      h_tof->Fill(time);
      if (pn_x==0 || pn_y==0 || pn_z==0) {continue;}
      if (time>10) {continue;}


      h_pvsp->Fill(pmiss.Mag(),pn.Mag());

      // select particle in momentum and angle accepted by CND
      double n_theta = pn.Theta()*180./M_PI;
      if (pn.Mag()<0.25 || pn.Mag()>1) {continue;}
      if (n_theta<45 || n_theta>140) {continue;}


      h_mmiss_pn->Fill(pn.Mag(),mmiss);
      h_mmiss_pmiss->Fill(pmiss.Mag(),mmiss);
      h_mmiss->Fill(mmiss);

      if (mmiss>1.) {continue;}



      h_pmiss_thetamiss->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
      h_thetapn_pp->Fill(pp.Mag(),pp.Angle(pn)*180./M_PI);

      if (cnd_energy>1000) {continue;}
      if (pmiss.Mag()<0.3 || pmiss.Mag()>1.) {continue;}
      if (pmiss.Theta()*180./M_PI<45 || pmiss.Theta()*180./M_PI>140) {continue;}




      if (energy<5) {continue;}



      h_compare->Fill((pmiss.Mag()-pn.Mag())/pmiss.Mag(),pn.Angle(pmiss)*180./M_PI);




      h_tof2->Fill(time);


      // histos
      h_nangles->Fill(pn.Phi()*180./M_PI,n_theta);
      h_energy->Fill(energy);
      h_Edep_beta->Fill(neut[i]->getBeta(),energy);
    
   



      // FILL HISTOS FOR NEUTRON CANDIDATES
      h_cos0->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
      h_pxminuspx->Fill(pn_x-pmiss.X());
      h_pyminuspy->Fill(pn_y-pmiss.Y());
      h_pzminuspz->Fill(pn_z-pmiss.Z());
      h_pminusp->Fill(pn.Mag()-pmiss.Mag());
      h_dpp->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
      h_theta_beta->Fill(beta,n_theta);
      h_p_theta->Fill(n_theta,pn.Mag());

      h_p_all->Fill(pmiss.Mag());
      // ML features
      h_energy_1->Fill(energy);
      h_layermult_1->Fill(layermult);
      h_size_1->Fill(size);
      h_cnd_hits_1->Fill(cnd_hits);
      h_cnd_energy_1->Fill(cnd_energy);
      h_ctof_energy_1->Fill(ctof_energy);
      h_ctof_hits_1->Fill(ctof_hits);
      h_anglediff_1->Fill(angle_diff);

    


  // physics cuts - not being used
  //if (angle_diff<30) {continue;}
  //if (cnd_hits>2) {continue;}
  //if (size>1) {continue;}
  //if (ctof_hits>2) {continue;}


//////////////////////////
/////     SORT       /////
////   GOOD / BAD    /////
////    NEUTRONS     /////
//////////////////////////

  // Determine whether to write to "good (signal) neutron" or "bad (background) neutron" file

  // good_N not used!
  bool good_N = pn.Angle(pmiss)*180./M_PI<20 && abs((pmiss.Mag()-pn.Mag())/pmiss.Mag())<0.2 && cnd_energy<1000; // this never gets used!! which is good - there are no neutrons in this channel.

  bool bad_N = (pn.Angle(pmiss)*180./M_PI<50 && abs((pmiss.Mag()-pn.Mag())/pmiss.Mag())<0.5) && cnd_energy<1000 && (pmiss.Mag()>0.25 && pmiss.Mag()<1.) && (pmiss.Theta()*180./M_PI>45 && pmiss.Theta()*180./M_PI<140);

  bool keep_this_one = keep_good ? good_N : bad_N;

  if (keep_this_one)
  {
    // all neutrons - print features
    outtxt << pmiss.Mag() << ' ';
    outtxt << energy << ' ';
    outtxt << layermult << ' ';
    outtxt << size << ' ';
    outtxt << cnd_hits << ' ';
    outtxt << cnd_energy << ' ';
    outtxt << ctof_energy << ' ';
    outtxt << ctof_hits << ' ';
    outtxt << angle_diff << ' ';
    outtxt << '\n';

    ntree->Fill();

  // FILL HISTOS FOR SIGNAL/BACKGROUND EVENTS
  h_nangles2->Fill(pn.Phi()*180./M_PI,n_theta);
  h_cos02->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
  h_pxminuspx2->Fill(pn_x-pmiss.X());
  h_pyminuspy2->Fill(pn_y-pmiss.Y());
  h_pzminuspz2->Fill(pn_z-pmiss.Z());
  h_pminusp2->Fill(pn.Mag()-pmiss.Mag());
  h_pvsp2->Fill(pmiss.Mag(),pn.Mag());
  h_dpp2->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
  h_mmiss2->Fill(mmiss);
  h_mmiss_pn->Fill(pn.Mag(),mmiss);
  h_energy2->Fill(energy);
  h_theta_beta2->Fill(beta,n_theta);
  h_p_theta2->Fill(n_theta,pn.Mag());
  h_pmiss_thetamiss2->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
  h_thetapn_pp2->Fill(pp.Mag(),pp.Angle(pn)*180./M_PI);

  h_compare2->Fill((pmiss.Mag()-pn.Mag())/pmiss.Mag(),pn.Angle(pmiss)*180./M_PI);
  h_Edep_beta2->Fill(neut[i]->getBeta(),energy);
  h_p_cut->Fill(pmiss.Mag());

  h_ptheta_pred->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
  h_ptheta->Fill(pn.Theta()*180/M_PI,pn.Mag());

  // ML features
  h_energy_2->Fill(energy);
  h_layermult_2->Fill(layermult);
  h_size_2->Fill(size);
  h_cnd_hits_2->Fill(cnd_hits);
  h_cnd_energy_2->Fill(cnd_energy);
  h_ctof_energy_2->Fill(ctof_energy);
  h_ctof_hits_2->Fill(ctof_hits);
  h_anglediff_2->Fill(angle_diff);

  }


  }  // closes neutron loop

    //chain.WriteEvent();

    counter++;

  }  // closes event loop


  f->cd();
  for(int i=0; i<hist_list_1.size(); i++) {
    hist_list_1[i]->Write();
  }
  for(int i=0; i<hist_list_2.size(); i++) {
    hist_list_2[i]->SetOption("colz");
    hist_list_2[i]->Write();
  }

  std::cout << '\n' <<counter << " events counted!\n\n";


  h_vzp->Write();

  outtxt.close();
  ntree->Write();
  f->Close();

  return 0;

}  // closes main function
