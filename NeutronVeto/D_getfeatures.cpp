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


using namespace std;
using namespace clas12;

void Usage() {
  std::cerr << "Usage: ./D_getfeatures Ebeam proton-detector(F/D) output-root output-txt keep_good input-hipo\n";
}

int main(int argc, char ** argv) {

  if(argc<7) {
    std::cerr << "Wrong number of arguments\n";
    Usage();
    return -1;
  }

  // arg 1: beam energy
  double Ebeam = atof(argv[1]);; // 2.07052 or 5.98636

  // arg 2: keep good
  bool keep_good = false;
  if(atoi(argv[2])==1){keep_good=true;}

  // arg 3: proton in Forward Detector or Central Detector
  char pDet = argv[3][0];

  // args 4-5: output file names
  TFile * f = new TFile(argv[4],"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");
  std::ofstream outtxt(argv[5]);

  // arg 6+: input hipo file
  clas12root::HipoChain chain;
  for (int k=6; k<argc; k++) {
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


  // charge is 0 for neutrons and 1 for protons
  /*TString hipo_name = "2gev_hipo/" + outName + ".hipo"; 
  clas12root::HipoChainWriter chain(hipo_name.Data()); // output hipo
  chain.Add(inName.Data()); // input hipo file


  */


  // create output root file and tree
  /*TString root_name = "2gev_root/" + outName + ".root";
  TFile * f = new TFile(root_name.Data(),"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");

  // create output txt file
  TString txt_name = "2gev_txt/" + outName + ".txt";
  ofstream outtxt(txt_name.Data());*/


  Int_t nhits;
  double px, py, pz, energy, time, path;
  Int_t sec[100] = {-1};
  Int_t lay[100] = {-1};
  int event;
  ntree->Branch("px",&px,"momentum x/D");
  ntree->Branch("py",&py,"momentum y/D");
  ntree->Branch("pz",&pz,"momentum z/D");
  ntree->Branch("nhits",&nhits,"number of hits/I");
  ntree->Branch("sec",sec,"sec[10]/I");
  ntree->Branch("lay",lay,"lay[10]/I");
  ntree->Branch("energy",&energy,"energy/D");
  ntree->Branch("event",&event,"event/I");



  // REC::Scintillator
  auto rec_scint = config_c12->addBank("REC::Scintillator");
  auto scint_detector = config_c12->getBankOrder(rec_scint,"detector");
  auto scint_sector = config_c12->getBankOrder(rec_scint,"sector");
  auto scint_layer = config_c12->getBankOrder(rec_scint,"layer");
  auto scint_component = config_c12->getBankOrder(rec_scint,"component");
  auto scint_energy = config_c12->getBankOrder(rec_scint,"energy");

  auto rec_scintx = config_c12->addBank("REC::ScintExtras");
  auto scint_size = config_c12->getBankOrder(rec_scintx,"size");
  

  // other banks
  auto rec_part = config_c12->addBank("REC::Particle");

 
  int counter = 0;
  //auto& c12=chain.C12ref();



  // histos


  // proton stuff
  TH1D * h_psize = new TH1D("psize","Number of Protons in Event",10,0,10);
    hist_list_1.push_back(h_psize);
  TH2D * h_dbeta_p = new TH2D("dbeta_p","#Delta #beta vs proton momentum",50,0,3,50,-0.2,0.2);
    hist_list_2.push_back(h_dbeta_p);


  // neutron stuff
  TH1D * h_nsize = new TH1D("nsize","Number of Neutrons in Event",10,0,10);
    hist_list_1.push_back(h_nsize);


  // reconstructed momentum
  TH2D * h_pangles = new TH2D("pangles","Proton Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_pangles);
  TH2D * h_nangles = new TH2D("nangles","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_nangles);
  TH1D * h_pxminuspx = new TH1D("pxminuspx","px_{n}-px_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pxminuspx);
  TH1D * h_pyminuspy = new TH1D("pyminuspy","py_{n}-py_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pyminuspy);
  TH1D * h_pzminuspz = new TH1D("pzminuspz","pz_{n}-pz_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pzminuspz);
  TH1D * h_pminusp = new TH1D("pminusp","p_{n}-p_{gen};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pminusp);
  TH2D * h_pvsp = new TH2D("pvsp","Momentum Resolution;p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    hist_list_2.push_back(h_pvsp);
  TH1D * h_cos0 = new TH1D("cos0","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
    hist_list_1.push_back(h_cos0);
  TH2D * h_dpp = new TH2D("dpp","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    hist_list_2.push_back(h_dpp);
  TH1D * h_energy = new TH1D("energy","Neutron Energy Deposition;Energy (MeV);Counts",100,0,50);
    hist_list_1.push_back(h_energy);
  TH2D * h_sec_phi = new TH2D("sec_phi","Sector vs Phi of CND hits;phi (deg);Sector",90,0,360,25,0,25);
    hist_list_2.push_back(h_sec_phi);
  TH1D * h_mmiss = new TH1D("mmiss","Missing Mass",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss);
  TH2D * h_mmiss_xb = new TH2D("mmiss_xb","Missing Mass vs x_{B}",50,0,3,50,0.5,1.5);
    hist_list_2.push_back(h_mmiss_xb);
  TH2D * h_theta_beta = new TH2D("theta_beta","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    hist_list_2.push_back(h_theta_beta);
  TH2D * h_p_theta = new TH2D("p_theta","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    hist_list_2.push_back(h_p_theta);
  TH2D * h_pmiss_thetamiss = new TH2D("pmiss_thetamiss","pmiss vs #theta_{pmiss};#theta_{pmiss};pmiss",90,0,180,50,0,1.2);
    hist_list_2.push_back(h_pmiss_thetamiss);
  TH2D * h_thetapn_pp = new TH2D("thetapn_pp","#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}",40,0,1,40,0,180);
    hist_list_2.push_back(h_thetapn_pp);
  TH2D * h_radiusz = new TH2D("radius_z","Radius vs z;z;Radius (cm)",80,-70,70,100,25,40);
    hist_list_2.push_back(h_radiusz);


  // good n / bad n set
  TH2D * h_nangles2 = new TH2D("nangles2","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    hist_list_2.push_back(h_nangles2);
  TH1D * h_pxminuspx2 = new TH1D("pxminuspx2","px_{n}-px_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pxminuspx2);
  TH1D * h_pyminuspy2 = new TH1D("pyminuspy2","py_{n}-py_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pyminuspy2);
  TH1D * h_pzminuspz2 = new TH1D("pzminuspz2","pz_{n}-pz_{miss};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pzminuspz2);
  TH1D * h_pminusp2 = new TH1D("pminusp2","p_{n}-p_{gen};Counts",100,-0.5,0.5);
    hist_list_1.push_back(h_pminusp2);
  TH2D * h_pvsp2 = new TH2D("pvsp2","Momentum Resolution;p_{miss} (GeV/c);p_{measured} (GeV/c)",100,0,1,100,0,1);
    hist_list_2.push_back(h_pvsp2);
  TH1D * h_cos02 = new TH1D("cos02","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
    hist_list_1.push_back(h_cos02);
  TH2D * h_dpp2 = new TH2D("dpp2","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    hist_list_2.push_back(h_dpp2);
  TH1D * h_mmiss2 = new TH1D("mmiss2","Missing Mass",50,0.5,1.5);
    hist_list_1.push_back(h_mmiss2);
  TH2D * h_mmiss_xb2 = new TH2D("mmiss_xb2","Missing Mass vs x_{B}",50,0,3,50,0.5,1.5);
    hist_list_2.push_back(h_mmiss_xb2);
  TH1D * h_energy2 = new TH1D("energy2","Neutron Energy Deposition;Energy (MeV);Counts",100,0,50);
    hist_list_1.push_back(h_energy2);
  TH2D * h_theta_beta2 = new TH2D("theta_beta2","Neutron theta vs beta;#beta;#theta",50,-0.1,1.1,55,35,145);
    hist_list_2.push_back(h_theta_beta2);
  TH2D * h_p_theta2 = new TH2D("p_theta2","Neutron Momentum vs Theta;#theta;p (GeV/c)",55,35,145,50,0,1.2);
    hist_list_2.push_back(h_p_theta2);
  TH2D * h_pmiss_thetamiss2 = new TH2D("pmiss_thetamiss2","pmiss vs #theta_{pmiss};#theta_{pmiss};pmiss",90,0,180,50,0,1.2);
    hist_list_2.push_back(h_pmiss_thetamiss2);
  TH2D * h_thetapn_pp2 = new TH2D("thetapn_pp2","#theta_{pn} vs p_{p};p_{p} (GeV/c);#theta_{pn}",40,0,1,40,0,180);
    hist_list_2.push_back(h_thetapn_pp2);
  TH2D * h_radiusz2 = new TH2D("radius_z2","Radius vs z;z;Radius (cm)",80,-70,70,100,25,40);
    hist_list_2.push_back(h_radiusz2);


const double mP = 0.93828;
const double mN = 0.939;
const double mD = 1.8756;



int numevent = 0;
  //while(chain.Next() && numevent<200)
  while(chain.Next())
  {

    // identify particles from REC::Particle
    //if (!myCut.electroncut(c12)) {continue;}
    auto elec=c12->getByID(11);
    auto prot = c12->getByID(2212);
    auto neut = c12->getByID(2112);
    auto piplus = c12->getByID(211);
    auto piminus = c12->getByID(-211);
    auto allParticles=c12->getDetParticles();
    if (elec.size()!=1) {continue;}
    if (prot.size()!=1) {continue;}
    if (neut.size()<1) {continue;}
    event = c12->runconfig()->getEvent() << '\n';

    // reject particles with the wrong PID
    bool trash = 0;
    for (int i=0; i<allParticles.size(); i++)
    {
      int pid = allParticles[i]->par()->getPid();
      if (pid!=2112 && pid!=11 && pid!=2212 && pid!=0 && pid!=22 && pid!=211) {trash=1;}
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
    double QSq = pq.Mag() - (nu*nu);
    double xB = QSq / (2*mN*nu);


//////////////////////////
/////     PROTONS    /////
//////////////////////////
    h_psize->Fill(prot.size());
    int p_index = -1;
    TVector3 pp(0.,0.,0.);
    for (int i=0; i<prot.size(); i++)
    {
      pp.SetMagThetaPhi(prot[i]->getP(),prot[i]->getTheta(),prot[i]->getPhi());
      double dbeta = prot[i]->par()->getBeta() - pp.Mag()/sqrt(pp.Mag2()+mP*mP);
      h_dbeta_p->Fill(pp.Mag(),dbeta);
      double vzp = prot[i]->par()->getVz();
      double chipid = prot[i]->par()->getChi2Pid();
      if ((vzp-vze)<-4. || (vzp-vze)>4.) {continue;}
      if (chipid<-3. || chipid>3.) {continue;}
      if (dbeta<-0.05 || dbeta>0.05) {continue;}
      // require proton to be in correct angle and momentum range for the requested etector
      double p_theta = pp.Theta()*180./M_PI;
      if (pDet=='F' && ((p_theta>40)                || (pp.Mag()<0.5 || pp.Mag()>3.0))) {continue;}
      if (pDet=='C' && ((p_theta<40 || p_theta>140) || (pp.Mag()<0.2 || pp.Mag()>1.2))) {continue;}
      p_index = i;
    }
    // NOT YET OPTIMIZED - what do I do if there are two protons?
    
    if (p_index<0) {continue;}
    pp.SetMagThetaPhi(prot[p_index]->getP(),prot[p_index]->getTheta(),prot[p_index]->getPhi());
    
    double p_theta = pp.Theta()*180./M_PI;
    h_pangles->Fill(pp.Phi()*180./M_PI,p_theta);


//////////////////////////
//  MISSING MOMENTUM    //
//////////////////////////

    // missing momentum
    TVector3 pmiss = pq-pp;


//////////////////////////
////     NEUTRONS    /////
//////////////////////////

    // LOOP OVER NEUTRONS
    h_nsize->Fill(neut.size());
    for (int i=0; i<neut.size(); i++) {
    
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
      
      int num_hits_inc = 0;
      if (is_CND1) {num_hits_inc = num_hits_inc + 1;}
      if (is_CND2) {num_hits_inc = num_hits_inc + 1;}
      if (is_CND3) {num_hits_inc = num_hits_inc + 1;}
       
      // put REC::Scintillator information
      int pindex, sector, layer, component, clusterid;
      int layermult = 0; int size = -1; int status = -1;
      double time, energy, path, x, y, z, dedx;
      double beta = neut[i]->par()->getBeta();
      
      if (is_CND1)
      {
        pindex = neut[i]->sci(CND1)->getPindex();
        sector = neut[i]->sci(CND1)->getSector();
        layer =  neut[i]->sci(CND1)->getLayer();
        component =  neut[i]->sci(CND1)->getComponent();
        time =   neut[i]->sci(CND1)->getTime() - starttime;
        energy = neut[i]->sci(CND1)->getEnergy();
        path =   neut[i]->sci(CND1)->getPath();
        status = neut[i]->sci(CND1)->getStatus();
        x =      neut[i]->sci(CND1)->getX();
        y =      neut[i]->sci(CND1)->getY();
        z =      neut[i]->sci(CND1)->getZ();
        dedx =   neut[i]->sci(CND1)->getDedx();
        size =   neut[i]->sci(CND1)->getSize();
        layermult = neut[i]->sci(CND1)->getLayermulti();
      }
      
      if (is_CND3)
      {
        pindex = neut[i]->sci(CND3)->getPindex();
        sector = neut[i]->sci(CND3)->getSector();
        layer =  neut[i]->sci(CND3)->getLayer();
        component =  neut[i]->sci(CND3)->getComponent();
        time =   neut[i]->sci(CND3)->getTime() - starttime;
        energy = neut[i]->sci(CND3)->getEnergy();
        path =   neut[i]->sci(CND3)->getPath();
        status = neut[i]->sci(CND3)->getStatus();
        x =      neut[i]->sci(CND3)->getX();
        y =      neut[i]->sci(CND3)->getY();
        z =      neut[i]->sci(CND3)->getZ();
        dedx =   neut[i]->sci(CND3)->getDedx();
        size =   neut[i]->sci(CND3)->getSize();
        layermult = neut[i]->sci(CND3)->getLayermulti();
      }
      
      if (is_CND2)
      {
        pindex = neut[i]->sci(CND2)->getPindex();
        sector = neut[i]->sci(CND2)->getSector();
        layer =  neut[i]->sci(CND2)->getLayer();
        component =  neut[i]->sci(CND2)->getComponent();
        time =   neut[i]->sci(CND2)->getTime() - starttime;
        energy = neut[i]->sci(CND2)->getEnergy();
        path =   neut[i]->sci(CND2)->getPath();
        status = neut[i]->sci(CND2)->getStatus();
        x =      neut[i]->sci(CND2)->getX();
        y =      neut[i]->sci(CND2)->getY();
        z =      neut[i]->sci(CND2)->getZ();
        dedx =   neut[i]->sci(CND2)->getDedx();
        size =   neut[i]->sci(CND2)->getSize();
        layermult = neut[i]->sci(CND2)->getLayermulti();
      }
      // PROBLEM: this gives preference to 2nd-layer hits
      if (!is_CND1 && !is_CND2 && !is_CND3)
      {
        sector = (neut[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
        layer = 0;
        component = 1; // value doesn't matter - not used
        time =   neut[i]->sci(CTOF)->getTime() - starttime;
        energy = neut[i]->sci(CTOF)->getEnergy();
        path =   neut[i]->sci(CTOF)->getPath();
        //status = neut[i]->sci(CTOF)->getStatus();
        x =      neut[i]->sci(CTOF)->getX();
        y =      neut[i]->sci(CTOF)->getY();
        z =      neut[i]->sci(CTOF)->getZ();
        dedx =   energy/3; // getDedx() is 0 (true for sim - is it also true for data?)
        size = neut[i]->sci(CTOF)->getSize();
        layermult = layermult + 1;
        layermult = neut[i]->sci(CTOF)->getLayermulti();
      }

      double n_phi = pn.Phi()*180./M_PI;
      double Ep = sqrt(mN*mN + pp.Mag2());
      double Emiss = Ebeam + mD - pe.Mag() - Ep;
      double mmiss = sqrt((Emiss*Emiss) - pmiss.Mag2());
    
      // calculate layer multiplicity by hand
      layermult=0; // default to 0 for CTOF
      if (is_CND1) {layermult = layermult+1;}
      if (is_CND2) {layermult = layermult+1;}
      if (is_CND3) {layermult = layermult+1;}
    
    
    
      // BASIC NEUTRONS CUTS
      double n_theta = pn.Theta()*180./M_PI;
      if (n_theta<40 || n_theta>140) {continue;}
      if (pn_x==0 || pn_y==0 || pn_z==0) {continue;}
      if (pn.Mag()<0.2) {continue;}
      if (energy<3) {continue;}
      //if (pmiss.Theta()*180./M_PI<40 || pmiss.Theta()*180./M_PI>140) {continue;}
      if (pmiss.Mag()<0.3 || pmiss.Mag()>1.2) {continue;}
    
      //if (xB<0.6) {continue;}
      double cos0 = pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag());
    
    
      // fill histos for neutron PID before good/bad selection
      h_nangles->Fill(pn.Phi()*180./M_PI,n_theta);
      h_cos0->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
      h_pxminuspx->Fill(pn_x-pmiss.X());
      h_pyminuspy->Fill(pn_y-pmiss.Y());
      h_pzminuspz->Fill(pn_z-pmiss.Z());
      h_pminusp->Fill(pn.Mag()-pmiss.Mag());
      h_pvsp->Fill(pmiss.Mag(),pn.Mag());
      h_energy->Fill(energy);
      h_dpp->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
      h_mmiss->Fill(mmiss);
      h_mmiss_xb->Fill(xB,mmiss);
      h_theta_beta->Fill(beta,n_theta);
      h_p_theta->Fill(n_theta,pn.Mag());
      h_pmiss_thetamiss->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
      h_thetapn_pp->Fill(pp.Mag(),pp.Angle(pn)*180./M_PI);
      h_radiusz->Fill(z,pow(x*x+y*y,0.5));
   //std::cout << x << '\t' <<  pow(x*x+y*y,0.5) << endl;
    
  // CND & CTOF NEARBY HITS
  int hits_nearby7 = 0; int ctof_nearby7 = 0;
  double cluster_energy7 = 0; double ctof_energy7 = 0;

  for (int j=0; j<c12->getBank(rec_scint)->getRows(); j++)
  {
    int rec_detector = c12->getBank(rec_scint)->getInt(scint_detector,j);
    if (rec_detector!=3 && rec_detector!=4) {continue;}

    int rec_sector = c12->getBank(rec_scint)->getInt(scint_sector,j);
    int rec_layer = c12->getBank(rec_scint)->getInt(scint_layer,j);
    int rec_component = c12->getBank(rec_scint)->getInt(scint_component,j);
    double rec_energy = c12->getBank(rec_scint)->getFloat(scint_energy,j);

    if (rec_detector==3 && (abs(rec_sector-sector)<3)) // hits in CND
    {
      hits_nearby7 = hits_nearby7 + c12->getBank(rec_scintx)->getInt(scint_size,j) ;
      cluster_energy7 = cluster_energy7 + rec_energy;
    }
    else if (rec_detector==3 && (abs(rec_sector-sector)>21)) // hits in CND, boundary
    {
      hits_nearby7 = hits_nearby7 + c12->getBank(rec_scintx)->getInt(scint_size,j);
      cluster_energy7 = cluster_energy7 + rec_energy;
    }
    else if (rec_detector==4 && (abs(rec_component-2*sector)<3)) // hits in CTOF //technically asymmetric
    {
      ctof_nearby7 = ctof_nearby7 + c12->getBank(rec_scintx)->getInt(scint_size,j) ;
      ctof_energy7 = ctof_energy7 + rec_energy;
    }
    else if (rec_detector==4 && abs(rec_component-2*sector)>44) // hits in CTOF, boundary //technically asymmetric
    {
      ctof_nearby7 = ctof_nearby7 + c12->getBank(rec_scintx)->getInt(scint_size,j) ;
      ctof_energy7 = ctof_energy7 + rec_energy;
    }
  }




  // CVT TRACKS
  double hit12_phi = 180;
  double angle_diff = 360;

  for (int j=0; j<allParticles.size(); j++)
  {
    // want k=1,3,5,7,12
    TVector3 traj1( allParticles[j]->traj(CVT,1)->getX(), allParticles[j]->traj(CVT,1)->getY(), allParticles[j]->traj(CVT,1)->getZ() ); 
    TVector3 traj3( allParticles[j]->traj(CVT,3)->getX(), allParticles[j]->traj(CVT,3)->getY(), allParticles[j]->traj(CVT,3)->getZ() );
    TVector3 traj5( allParticles[j]->traj(CVT,5)->getX(), allParticles[j]->traj(CVT,5)->getY(), allParticles[j]->traj(CVT,5)->getZ() );
    TVector3 traj7( allParticles[j]->traj(CVT,7)->getX(), allParticles[j]->traj(CVT,7)->getY(), allParticles[j]->traj(CVT,7)->getZ() );
    TVector3 traj12( allParticles[j]->traj(CVT,12)->getX(), allParticles[j]->traj(CVT,12)->getY(), allParticles[j]->traj(CVT,12)->getZ() );


    if (traj12.X()==0 || traj12.Y()==0 || traj12.Z()==0) {continue;}

    // take the track that is closest in angle to the neutron hit
    if ( (pn.Angle(traj12)*180./M_PI) < angle_diff )
    {
      hit12_phi = pn.Angle(traj12)*180./M_PI;
      angle_diff = hit12_phi;
    }


    double drdz = 0;
    double r_t1 = traj1.X()*traj1.X() + traj1.Y()*traj1.Y();
    double r_t3 = traj3.X()*traj3.X() + traj3.Y()*traj3.Y();
    double r_t5 = traj5.X()*traj5.X() + traj5.Y()*traj5.Y();
    double r_t7 = traj7.X()*traj7.X() + traj7.Y()*traj7.Y();
    double r_t12 = traj12.X()*traj12.X() + traj12.Y()*traj12.Y();
    /*drdz = drdz + (r_t1 - r_t3) / (traj1.Z() - traj3.Z());
    drdz = drdz + (r_t3 - r_t5) / (traj3.Z() - traj5.Z());
    drdz = drdz + (r_t5 - r_t7) / (traj5.Z() - traj7.Z());
    drdz = drdz + (r_t7 - r_t12) / (traj7.Z() - traj12.Z());*/

  }






//////////////////////////
/////     SORT       /////
////   GOOD / BAD    /////
////    NEUTRONS     /////
//////////////////////////

  // Determine whether to write to "good neutron" or "bad neutron" file

  bool good_N = (cos0>0.9 && abs(pmiss.Mag()-pn.Mag())<0.1 && mmiss< 1.05);

  bool bad_N = (cos0<0.8 && abs(pmiss.Mag()-pn.Mag())>0.2); // shown in paris
  //bool bad_N = (mmiss>0.8 && mmiss<1.05 && (pmiss.Theta()*180./M_PI<40 || pmiss.Theta()*180./M_PI>140));  // Justin's idea
  //bool bad_N = (pmiss.Theta()*180./M_PI<40 || pmiss.Theta()*180./M_PI>140);
  //bool bad_N = !good_N;
  //bool bad_N = (prot.size()==2 && piminus.size()==1); // no stats

  bool keep_this_one = keep_good ? good_N : bad_N;

  if (keep_this_one)
  {
    // all neutrons - print features
    outtxt << energy << ' ';
    outtxt << layermult << ' '; //outtxt << z << ' ';
    outtxt << size << ' '; //outtxt << beta << ' ';
    outtxt << hits_nearby7 << ' ';
    outtxt << cluster_energy7 << ' ';
    outtxt << ctof_energy7 << ' ';
    outtxt << ctof_nearby7 << ' ';
    //outtxt << abs(hit12_phi-n_phi) << ' ';
    outtxt << angle_diff << ' ';
    outtxt << '\n';


  h_nangles2->Fill(pn.Phi()*180./M_PI,n_theta);
  h_cos02->Fill(pmiss.Dot(pn) / (pmiss.Mag()*pn.Mag()));
  h_pxminuspx2->Fill(pn_x-pmiss.X());
  h_pyminuspy2->Fill(pn_y-pmiss.Y());
  h_pzminuspz2->Fill(pn_z-pmiss.Z());
  h_pminusp2->Fill(pn.Mag()-pmiss.Mag());
  h_pvsp2->Fill(pmiss.Mag(),pn.Mag());
  h_dpp2->Fill(pmiss.Mag(),(pmiss.Mag()-pn.Mag())/pmiss.Mag());
  h_mmiss2->Fill(mmiss);
  h_mmiss_xb2->Fill(xB,mmiss);
  h_energy2->Fill(energy);
  h_theta_beta2->Fill(beta,n_theta);
  h_p_theta2->Fill(n_theta,pn.Mag());
  h_pmiss_thetamiss2->Fill(pmiss.Theta()*180./M_PI,pmiss.Mag());
  h_thetapn_pp2->Fill(pp.Mag(),pp.Angle(pn)*180./M_PI);
  h_radiusz2->Fill(z,pow(x*x+y*y,0.5));
  }

  }  // closes neutron loop



    //chain.WriteEvent();
    ntree->Fill();

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


/*
  // write histograms
  h_psize->Write();
  h_pxminuspx->Write();
  h_pyminuspy->Write();
  h_pzminuspz->Write();
  h_pminusp->Write();
  h_pvsp->Write();
  h_pangles->Write();
  h_nangles->Write();
  h_cos0->Write();
  h_energy->Write();
  h_dpp->Write();
  h_mmiss->Write();
  h_mmiss_xb->Write();
  h_dbeta_p->Write();
  h_theta_beta->Write();
  h_p_theta->Write();
  h_pmiss_thetamiss->Write();
  h_thetapn_pp->Write();



  h_nangles2->Write();
  h_cos02->Write();
  h_pxminuspx2->Write();
  h_pyminuspy2->Write();
  h_pzminuspz2->Write();
  h_pminusp2->Write();
  h_pvsp2->Write();
  h_dpp2->Write();
  h_mmiss2->Write();
  h_mmiss_xb2->Write();
  h_energy2->Write();
  h_theta_beta2->Write();
  h_p_theta2->Write();
  h_pmiss_thetamiss2->Write();
  h_thetapn_pp2->Write();
*/

  outtxt.close();
  ntree->Write();
  f->Close();

}  // closes main function