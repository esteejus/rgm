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
  std::cerr << "Usage: ./N_getfeatures charge output-root output-txt input-hipo\n";
}

double getCVTdiff(std::vector<region_part_ptr> &allParticles, TVector3 &pn);

int main(int argc, char ** argv)
{

  if(argc<5) {
    std::cerr << "Wrong number of arguments\n";
    Usage();
    return -1;
  }

  // argument 1: particle charge (1 for proton, 0 for neutron)
  int charge = atoi(argv[1]);

  // argument 2-3: output file names
  TFile * f = new TFile(argv[2],"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");
  std::ofstream outtxt(argv[3]);

  // argument 4+: input hipo files
  clas12root::HipoChain chain;
  for (int k=4; k<argc; k++) {
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


/*  // charge is 0 for neutrons and 1 for protons
  TString hipo_name = "bknd_hipo/" + outName + ".hipo"; 
  clas12root::HipoChainWriter chain(hipo_name.Data()); // output hipo
  chain.Add(inName.Data()); // input hipo file
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();

  auto config_c12=chain.GetC12Reader();


  // create output root file and tree
  TString root_name = "bknd_root/" + outName + ".root";
  TFile * f = new TFile(root_name.Data(),"RECREATE");
  TTree * ntree = new TTree("T","NeutronTree");

  // create output txt file
  TString txt_name = "bknd_txt/" + outName + ".txt";
  ofstream outtxt(txt_name.Data());*/


  Int_t nhits;
  double px, py, pz; //double energy;
  //Float_t energy[100] = {-1};
  Int_t sec[100] = {-1};
  Int_t lay[100] = {-1};
  int event;
  double energy, cnd_energy, ctof_energy, angle_diff;
  int layermult, size, cnd_hits, ctof_hits;
  //int sec, lay, event;
  /*ntree->Branch("px",&px,"momentum x/D");
  ntree->Branch("py",&py,"momentum y/D");
  ntree->Branch("pz",&pz,"momentum z/D");
  ntree->Branch("nhits",&nhits,"number of hits/I");
  ntree->Branch("sec",sec,"sec[10]/I");
  ntree->Branch("lay",lay,"lay[10]/I");
  ntree->Branch("energy",&energy,"energy/D");
  ntree->Branch("event",&event,"event/I");*/
  ntree->Branch("energy",&energy,"energy/D");
  ntree->Branch("layermult",&layermult,"layermult/I");
  ntree->Branch("size",&size,"size/I");
  ntree->Branch("cnd_hits",&cnd_hits,"cnd_hits/I");
  ntree->Branch("cnd_energy",&cnd_energy,"cnd_energy/D");
  ntree->Branch("ctof_energy",&ctof_energy,"ctof_energy/D");
  ntree->Branch("ctof_hits",&ctof_hits,"ctof_hits/I");
  ntree->Branch("angle_diff",&angle_diff,"angle_diff/D");
  
  

  // MC banks
  auto mc_p = config_c12->addBank("MC::Particle");
  auto mc_pid = config_c12->getBankOrder(mc_p,"pid");
  auto mc_px = config_c12->getBankOrder(mc_p,"px");
  auto mc_py = config_c12->getBankOrder(mc_p,"py");
  auto mc_pz = config_c12->getBankOrder(mc_p,"pz");


/*  // CND hits
  auto cnd_hits = config_c12->addBank("CND::hits");
  auto cnd_id = config_c12->getBankOrder(cnd_hits,"id");
  auto cnd_status = config_c12->getBankOrder(cnd_hits,"status");
  auto cnd_trkID = config_c12->getBankOrder(cnd_hits,"trkID");
  auto cnd_sector = config_c12->getBankOrder(cnd_hits,"sector");
  auto cnd_layer = config_c12->getBankOrder(cnd_hits,"layer");
  auto cnd_component = config_c12->getBankOrder(cnd_hits,"component");
  auto cnd_energy = config_c12->getBankOrder(cnd_hits,"energy");
  auto cnd_time = config_c12->getBankOrder(cnd_hits,"time");
  auto cnd_energy_unc = config_c12->getBankOrder(cnd_hits,"energy_unc");
  auto cnd_time_unc = config_c12->getBankOrder(cnd_hits,"time_unc");
  auto cnd_x = config_c12->getBankOrder(cnd_hits,"x");
  auto cnd_y = config_c12->getBankOrder(cnd_hits,"y");
  auto cnd_z = config_c12->getBankOrder(cnd_hits,"z");
  auto cnd_x_unc = config_c12->getBankOrder(cnd_hits,"x_unc");
  auto cnd_y_unc = config_c12->getBankOrder(cnd_hits,"y_unc");
  auto cnd_z_unc = config_c12->getBankOrder(cnd_hits,"z_unc");
  auto cnd_tx = config_c12->getBankOrder(cnd_hits,"tx");
  auto cnd_ty = config_c12->getBankOrder(cnd_hits,"ty");
  auto cnd_tz = config_c12->getBankOrder(cnd_hits,"tz");
  auto cnd_tlength = config_c12->getBankOrder(cnd_hits,"tlength");
  auto cnd_pathlength = config_c12->getBankOrder(cnd_hits,"pathlength");
  auto cnd_indexLadc = config_c12->getBankOrder(cnd_hits,"indexLadc");
  auto cnd_indexRadc = config_c12->getBankOrder(cnd_hits,"indexRadc");
  auto cnd_indexLtdc = config_c12->getBankOrder(cnd_hits,"indexLtdc");
  auto cnd_indexRtdc = config_c12->getBankOrder(cnd_hits,"indexRtdc");

  // CND adc/tdc
  auto cnd_adc = config_c12->addBank("CND::adc");
  auto cnd_tdc = config_c12->addBank("CND::tdc");

  // CND clusters
  auto cnd_clusters = config_c12->addBank("CND::clusters");
  auto clust_id = config_c12->getBankOrder(cnd_clusters,"id");
  auto clust_sector = config_c12->getBankOrder(cnd_clusters,"sector");
  auto clust_layer = config_c12->getBankOrder(cnd_clusters,"layer");
  auto clust_component = config_c12->getBankOrder(cnd_clusters,"component");
  auto clust_nhits = config_c12->getBankOrder(cnd_clusters,"nhits");
  auto clust_energy = config_c12->getBankOrder(cnd_clusters,"energy");
  auto clust_x = config_c12->getBankOrder(cnd_clusters,"x");
  auto clust_y = config_c12->getBankOrder(cnd_clusters,"y");
  auto clust_z = config_c12->getBankOrder(cnd_clusters,"z");
  auto clust_time = config_c12->getBankOrder(cnd_clusters,"time");
  auto clust_status = config_c12->getBankOrder(cnd_clusters,"status");
  auto clust_size = config_c12->getBankOrder(cnd_clusters,"size");

  // CTOF hits
  auto ctof_hits = config_c12->addBank("CTOF::hits");
  auto ctof_id = config_c12->getBankOrder(ctof_hits,"id");
  auto ctof_layer = config_c12->getBankOrder(ctof_hits,"layer");
  auto ctof_sector = config_c12->getBankOrder(ctof_hits,"sector");
  auto ctof_component = config_c12->getBankOrder(ctof_hits,"component");
  auto ctof_energy = config_c12->getBankOrder(ctof_hits,"energy");
  auto ctof_x = config_c12->getBankOrder(ctof_hits,"x");
  auto ctof_y = config_c12->getBankOrder(ctof_hits,"y");
  auto ctof_z = config_c12->getBankOrder(ctof_hits,"z");

  // CTOF adc/tdc
  auto ctof_adc = config_c12->addBank("CTOF::adc");
  auto ctof_tdc = config_c12->addBank("CTOF::tdc");

  // CTOF clusters
  auto ctof_clusters = config_c12->addBank("CTOF::clusters");
  auto ctof_clus_size = config_c12->getBankOrder(ctof_clusters,"size");
  auto ctof_clus_sector = config_c12->getBankOrder(ctof_clusters,"sector");
  auto ctof_clus_layer = config_c12->getBankOrder(ctof_clusters,"layer");
  auto ctof_clus_component = config_c12->getBankOrder(ctof_clusters,"component");
  auto ctof_clus_energy = config_c12->getBankOrder(ctof_clusters,"energy");
  auto ctof_clus_time = config_c12->getBankOrder(ctof_clusters,"time"); // try this one!!!
  auto ctof_clus_status = config_c12->getBankOrder(ctof_clusters,"status");*/

  // REC::Scintillator
  auto rec_scint = config_c12->addBank("REC::Scintillator");
  auto scint_detector = config_c12->getBankOrder(rec_scint,"detector");
  auto scint_sector = config_c12->getBankOrder(rec_scint,"sector");
  auto scint_layer = config_c12->getBankOrder(rec_scint,"layer");
  auto scint_component = config_c12->getBankOrder(rec_scint,"component");
  auto scint_energy = config_c12->getBankOrder(rec_scint,"energy");

  // ScintExtras
  auto rec_scintx = config_c12->addBank("REC::ScintExtras");
  auto scint_dedx = config_c12->getBankOrder(rec_scintx,"dedx");
  auto scint_size = config_c12->getBankOrder(rec_scintx,"size");
  auto scint_layermult = config_c12->getBankOrder(rec_scintx,"layermult");
  



 
  int counter = 0;



  // histos

  // generated momentum
  TH2D * h_px = new TH2D("px","px;generated px;reconstructed px",100,0,1.5,100,0,1.5);
    h_px->SetOption("colz");
  TH2D * h_py = new TH2D("py","py;generated py;reconstructed py",100,0,1.5,100,0,1.5);
    h_py->SetOption("colz");
  TH2D * h_pz = new TH2D("pz","pz;generated pz;reconstructed pz",100,0,1.5,100,0,1.5);
    h_pz->SetOption("colz");
  TH2D * h_p = new TH2D("p","p;generated p;reconstructed p",100,0,1.5,100,0,1.5);
    h_p->SetOption("colz");
  TH1D * h_pg_theta = new TH1D("pg_theta","Generated Theta",180,0,180);

  // reconstructed momentum
  TH2D * h_nangles = new TH2D("nangles","Neutron Angles;phi;theta",48,-180,180,45,0,180);
    h_nangles->SetOption("colz");
  TH1D * h_pxminuspx = new TH1D("pxminuspx","(px_{n}-px_{gen})/px_{gen};Counts",100,-0.5,0.5);
  TH1D * h_pyminuspy = new TH1D("pyminuspy","(py_{n}-py_{gen})/py_{gen};Counts",100,-0.5,0.5);
  TH1D * h_pzminuspz = new TH1D("pzminuspz","(pz_{n}-pz_{gen})/pz_{gen};Counts",100,-0.5,0.5);
  TH1D * h_pminusp = new TH1D("pminusp","p_{n}-p_{gen};Counts",100,-0.5,0.5);
  TH2D * h_pvsp = new TH2D("pvsp","Momentum Resolution;p_{generated} (GeV/c);g_{measured} (GeV/c)",100,0,1,100,0,1);
    h_pvsp->SetOption("colz");
  TH2D * h_dpp = new TH2D("dpp","Momentum Resolution;p_{generated} (GeV/c);#Delta p/p",100,0,1,100,-0.4,0.4);
    h_dpp->SetOption("colz");
  TH1D * h_cos0 = new TH1D("cos0","Cosine of angle between generated and reconstructed p",50,-1.1,1.1);
  TH1D * h_hitsec = new TH1D("hitsec","CND hit sector",25,0,25);
  TH1D * h_hitlay = new TH1D("hitlay","CND hit layer",10,0,4);
  TH1D * h_cos1 = new TH1D("cos1","Cosine of angle between generated p and cluster hit",50,-1.1,1.1);

  TH1D * h_energy = new TH1D("energy","Neutron energy deposition;Energy (MeV);Counts",1000,0,1000);
  TH2D * h_sec_phi = new TH2D("sec_phi","Sector vs Phi of CND hits;phi (deg);Sector",90,0,360,25,0,25);
    h_sec_phi->SetOption("colz");




int numevent = 0;
  //while(chain.Next() && numevent<20)
  while(chain.Next())
  {

    // initialize features
    energy = 0; cnd_energy = 0; ctof_energy = 0; angle_diff = 180;
    layermult = -1; size = 0; cnd_hits = 0; ctof_hits = 0;


    // define particles
    TVector3 p_g(0.,0.,0.);
    TVector3 p(0.,0.,0.);
    TVector3 pe(0.,0.,0.);


    // identify particles from REC::Particle
    auto elec=c12->getByID(11);
    auto nucl= (charge==0) ? c12->getByID(2112) : c12->getByID(2212);
    auto allParticles = c12->getDetParticles();
    double weight = c12->mcevent()->getWeight();
    if (elec.size()!=1) {continue;}
    if (nucl.size()<1) {continue;}
    event = c12->runconfig()->getEvent() << '\n';



    // electron momentum
    double pe_x = elec[0]->par()->getPx();
    double pe_y = elec[0]->par()->getPy();
    double pe_z = elec[0]->par()->getPz();
    pe.SetXYZ(pe_x,pe_y,pe_z);


    // read Monte Carlo nucleon PID and momentum
    double px_g, py_g, pz_g;
    px_g = c12->getBank(mc_p)->getFloat(mc_px,1);
    py_g = c12->getBank(mc_p)->getFloat(mc_py,1);
    pz_g = c12->getBank(mc_p)->getFloat(mc_pz,1);
    p_g.SetXYZ(px_g,py_g,pz_g);
    h_pg_theta->Fill(p_g.Theta()*180./M_PI);

    // read reconstructed nucleon PID and momentum
    int n0 = -1;
    double max_cos0 = 0.5;


  //std::cout << "number of neutrons " << nucl.size() << endl;
  numevent = numevent + 1;
  double starttime = c12->event()->getStartTime();


  // PRINT BANK INFO //
  // LOOP OVER NEUTRONS
  for (int i=0; i<nucl.size(); i++)
  {

    // get neutron momentum
    px = nucl[i]->par()->getPx();
    py = nucl[i]->par()->getPy();
    pz = nucl[i]->par()->getPz();
    p.SetXYZ(px,py,pz);
    double n_theta = p.Theta()*180./M_PI;
    

    // reject neutrons that we have no interest in
    if (px==0 || py==0 || pz==0) {continue;}
    if (p.Mag()<0.2) {continue;}
    if (p_g.Mag()<0.2) {continue;}
    if (p_g.Theta()*180./M_PI<40 || p_g.Theta()*180./M_PI>135) {continue;}
    //if (n_theta<40 || n_theta>135) {continue;}
    bool is_CD = nucl[i]->getRegion()==CD;
    if (!is_CD) {continue;}

    // figure out what layer the hit is in - check 0 if not found
    bool is_CND1 = (nucl[i]->sci(CND1)->getLayer()==1);
    bool is_CND2 = (nucl[i]->sci(CND2)->getLayer()==2);
    bool is_CND3 = (nucl[i]->sci(CND3)->getLayer()==3);

    int num_hits_inc = 0;
    if (is_CND1) {num_hits_inc = num_hits_inc + 1;}
    if (is_CND2) {num_hits_inc = num_hits_inc + 1;}
    if (is_CND3) {num_hits_inc = num_hits_inc + 1;}  

    // put REC::Scintillator information
    int status = -1;
    layermult = -1;
    double time, path, x, y, z;
    //int sector = -1; int layer = -1; int component = -1;
    int sector = 0; int layer = 0; int component = 0;
    double beta = nucl[i]->par()->getBeta();


    // same as cluster information
    if (is_CND1)
    {
      sector = nucl[i]->sci(CND1)->getSector();
      layer =  nucl[i]->sci(CND1)->getLayer();
      component =  nucl[i]->sci(CND1)->getComponent();
      time =   nucl[i]->sci(CND1)->getTime() - starttime;
      energy = nucl[i]->sci(CND1)->getEnergy();
      path =   nucl[i]->sci(CND1)->getPath();
      status = nucl[i]->sci(CND1)->getStatus();
      x =      nucl[i]->sci(CND1)->getX();
      y =      nucl[i]->sci(CND1)->getY();
      z =      nucl[i]->sci(CND1)->getZ();
      size =   nucl[i]->sci(CND1)->getSize();
      //layermult = nucl[i]->sci(CND1)->getLayermulti();
    }

    if (is_CND3)
    {
      sector = nucl[i]->sci(CND3)->getSector();
      layer =  nucl[i]->sci(CND3)->getLayer();
      component =  nucl[i]->sci(CND3)->getComponent();
      time =   nucl[i]->sci(CND3)->getTime() - starttime;
      energy = nucl[i]->sci(CND3)->getEnergy();
      path =   nucl[i]->sci(CND3)->getPath();
      status = nucl[i]->sci(CND3)->getStatus();
      x =      nucl[i]->sci(CND3)->getX();
      y =      nucl[i]->sci(CND3)->getY();
      z =      nucl[i]->sci(CND3)->getZ();
      size =   nucl[i]->sci(CND3)->getSize();
      //layermult = nucl[i]->sci(CND3)->getLayermulti();
    }

    if (is_CND2)
    {
      sector = nucl[i]->sci(CND2)->getSector();
      layer =  nucl[i]->sci(CND2)->getLayer();
      component =  nucl[i]->sci(CND2)->getComponent();
      time =   nucl[i]->sci(CND2)->getTime() - starttime;
      energy = nucl[i]->sci(CND2)->getEnergy();
      path =   nucl[i]->sci(CND2)->getPath();
      status = nucl[i]->sci(CND2)->getStatus();
      x =      nucl[i]->sci(CND2)->getX();
      y =      nucl[i]->sci(CND2)->getY();
      z =      nucl[i]->sci(CND2)->getZ();
      size =   nucl[i]->sci(CND2)->getSize();
      //layermult = nucl[i]->sci(CND2)->getLayermulti();
    }
    // PROBLEM: this gives preference to 2nd-layer hits
    if (!is_CND1 && !is_CND2 && !is_CND3)
    {
      sector = (nucl[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
      layer =  0;
      component = 1; // value doesn't matter - not used
      time =   nucl[i]->sci(CTOF)->getTime() - starttime;
      energy = nucl[i]->sci(CTOF)->getEnergy();
      path =   nucl[i]->sci(CTOF)->getPath();
      //status = nucl[i]->sci(CTOF)->getStatus();
      x =      nucl[i]->sci(CTOF)->getX();
      y =      nucl[i]->sci(CTOF)->getY();
      z =      nucl[i]->sci(CTOF)->getZ();
      size =   nucl[i]->sci(CTOF)->getSize();
      //layermult = nucl[i]->sci(CTOF)->getLayermulti();
    }


    // calculate layer multiplicity by hand
    if (is_CND1) {layermult = layermult+1;}
    if (is_CND2) {layermult = layermult+1;}
    if (is_CND3) {layermult = layermult+1;}

    // fill histos
    h_nangles->Fill(p.Phi()*180./M_PI,n_theta,weight);
    h_cos0->Fill(p_g.Dot(p) / (p_g.Mag()*p.Mag()));
    h_pxminuspx->Fill((px-px_g)/px_g,weight);
    h_pyminuspy->Fill((py-py_g)/py_g,weight);
    h_pzminuspz->Fill((pz-pz_g)/pz_g,weight);
    h_pminusp->Fill(p.Mag()-p_g.Mag(),weight);
    h_pvsp->Fill(p_g.Mag(),p.Mag(),weight);
    h_dpp->Fill(p_g.Mag(),(p_g.Mag()-p.Mag())/p_g.Mag(),weight);
    h_energy->Fill(energy,weight);


if (energy<3) {continue;}




/*    // CTOF BANK INFO
    double ctof_event_energy = 0;
    bool ctof_hit_match = 0;
    int ctof_hit_component = -1;
    int ctof_nearby = 0; // CHECK THAT THIS IS CALCULATED CORRECTLY!!!
    //double ctof_neutron_angle = 0; // HOW TO CHECK FOR SMALLEST?
    for (int j=0; j<c12->getBank(ctof_hits)->getRows(); j++)
    {
      if (c12->getBank(ctof_hits)->getRows()==0) {continue;}
      ctof_hit_component = c12->getBank(ctof_hits)->getInt(ctof_component,j);
      //check nonzero, just in expected range
      if (abs(sector*2-ctof_hit_component)<3) // ADD BOUNDARY CONDITION
      {
        ctof_hit_match = 1;
        ctof_event_energy = ctof_event_energy + c12->getBank(ctof_hits)->getFloat(ctof_energy,j);
        ctof_nearby = ctof_nearby + 1;
      }
      else if (abs(sector*2-ctof_hit_component)>44)
      {
        ctof_hit_match = 1;
        ctof_event_energy = ctof_event_energy + c12->getBank(ctof_hits)->getFloat(ctof_energy,j);
        ctof_nearby = ctof_nearby + 1;
      }
      // check if there is a ctof hit in a sector corresponding to the cnd hit
    }


    // CTOF CLUSTERS INFO
    double ctofclus_event_energy = 0;
    int ctofclus_component = -1;
    int ctofclus_nearby = 0;
    double ctof_time = 0;

    for (int j=0; j<c12->getBank(ctof_clusters)->getRows(); j++)
    {
      if (c12->getBank(ctof_clusters)->getRows()==0) {continue;}

      ctofclus_component = c12->getBank(ctof_clusters)->getInt(ctof_clus_component,j);
      if (abs(sector*2-ctofclus_component)<3)
      {
        ctofclus_event_energy = ctofclus_event_energy + c12->getBank(ctof_clusters)->getFloat(ctof_clus_energy,j);
        ctofclus_nearby = ctofclus_nearby + 1;
        ctof_time = ctof_time + c12->getBank(ctof_clusters)->getFloat(ctof_clus_time,j);
      }
      else if (abs(sector*2-ctofclus_component)>44)
      {
        ctofclus_event_energy = ctofclus_event_energy + c12->getBank(ctof_clusters)->getFloat(ctof_clus_energy,j);
        ctofclus_nearby = ctofclus_nearby + 1;
        ctof_time = ctof_time + c12->getBank(ctof_clusters)->getFloat(ctof_clus_time,j);
      }
    }
    if (ctofclus_nearby>0) {ctof_time = ctof_time/ctofclus_nearby;}



    // CND CLUSTERS BANK INFO
    int hits_nearby7 = 0;
    double cluster_energy7 = 0;
    int layer_width = 0;
    int sector_depth = 0;

    for (int j=0; j<c12->getBank(cnd_clusters)->getRows(); j++)
    {
      int cluster_sector = c12->getBank(cnd_clusters)->getInt(clust_sector,j);
      int cluster_layer = c12->getBank(cnd_clusters)->getInt(clust_layer,j);
      // count nearby hit and nearby energy
      if (abs(cluster_sector-sector)<4)
      {
        hits_nearby7 = hits_nearby7 + 1;
        cluster_energy7 = cluster_energy7 + c12->getBank(cnd_clusters)->getFloat(clust_energy,j);
      }
      else if (abs(cluster_sector-sector)>20)
      {
        hits_nearby7 = hits_nearby7 + 1;
        cluster_energy7 = cluster_energy7 + c12->getBank(cnd_clusters)->getFloat(clust_energy,j);
      }

      // measure width of response in neutron layer
      if (abs(cluster_sector-sector)<3 && cluster_layer==layer) { layer_width = layer_width + 1; }
      else if (abs(cluster_sector-sector)==22 && cluster_layer==layer) { layer_width = layer_width + 1; }
      else if (abs(cluster_sector-sector)==23 && cluster_layer==layer) { layer_width = layer_width + 1; }

      // count how many layers in sector have response
      if (cluster_sector==sector)  { sector_depth = sector_depth + 1; }
    }




    // CND hits (not only clusters) near neutron
    int cndhits_nearby7 = 0;
    double cndhits_energy7 = 0;
    for (int j=0; j<c12->getBank(cnd_hits)->getRows(); j++)
    {
      int cndhits_sector = c12->getBank(cnd_hits)->getInt(cnd_sector,j);
      int cndhits_layer = c12->getBank(cnd_hits)->getInt(cnd_layer,j);
      // count nearby hit and nearby energy
      if (abs(cndhits_sector-sector)<3)  // ADD BOUNDARY CONDITION!!
      {
        cndhits_nearby7 = cndhits_nearby7 + 1;
        cndhits_energy7 = cndhits_energy7 + c12->getBank(cnd_hits)->getFloat(cnd_energy,j);
      }
      if (abs(cndhits_sector-sector)>21)
      {
        cndhits_nearby7 = cndhits_nearby7 + 1;
        cndhits_energy7 = cndhits_energy7 + c12->getBank(cnd_hits)->getFloat(cnd_energy,j);
      }
    }*/


    // CND & CTOF HEARBY HITS
    for (int j=0; j<c12->getBank(rec_scint)->getRows(); j++)
    {
      int rec_detector = c12->getBank(rec_scint)->getInt(scint_detector,j);
      if (rec_detector!=3 && rec_detector!=4) {continue;}
      
      int rec_sector = c12->getBank(rec_scint)->getInt(scint_sector,j);
      int rec_layer = c12->getBank(rec_scint)->getInt(scint_layer,j);
      int rec_component = c12->getBank(rec_scint)->getInt(scint_component,j);
      double rec_energy = c12->getBank(rec_scint)->getFloat(scint_energy,j);
      // note - rec_energy sometimes comes out insanely big
      
      if (rec_detector==3 && (abs(rec_sector-sector)<3)) // hits in CND
      {
        cnd_hits = cnd_hits + c12->getBank(rec_scintx)->getInt(scint_size,j);
        cnd_energy = cnd_energy + rec_energy;
      }
      else if (rec_detector==3 && (abs(rec_sector-sector)>21)) // hits in CND, boundary
      {
        cnd_hits = cnd_hits + c12->getBank(rec_scintx)->getInt(scint_size,j);
        cnd_energy = cnd_energy + rec_energy;
      }
      else if (rec_detector==4 && (abs(rec_component-2*sector)<3)) // hits in CTOF // technically asymmetric
      {
        ctof_hits = ctof_hits + c12->getBank(rec_scintx)->getInt(scint_size,j);
        ctof_energy = ctof_energy + rec_energy;
      }
      else if (rec_detector==4 && (abs(rec_component-2*sector)>44)) // hits in CTOF, boundary // technically asymmetric
      {
        ctof_hits = ctof_hits + c12->getBank(rec_scintx)->getInt(scint_size,j);
        ctof_energy = ctof_energy + rec_energy;
      }

    }



    // CVT Tracks - does basically nothing, almost no tracks are formed
    double hit12_phi = 180;
    for (int j=0; j<allParticles.size(); j++)
    {
      TVector3 traj1( allParticles[j]->traj(CVT,1)->getX(), allParticles[j]->traj(CVT,1)->getY(), allParticles[j]->traj(CVT,1)->getZ() );
      TVector3 traj3( allParticles[j]->traj(CVT,3)->getX(), allParticles[j]->traj(CVT,3)->getY(), allParticles[j]->traj(CVT,3)->getZ() );
      TVector3 traj5( allParticles[j]->traj(CVT,5)->getX(), allParticles[j]->traj(CVT,5)->getY(), allParticles[j]->traj(CVT,5)->getZ() );
      TVector3 traj7( allParticles[j]->traj(CVT,7)->getX(), allParticles[j]->traj(CVT,7)->getY(), allParticles[j]->traj(CVT,7)->getZ() );
      TVector3 traj12( allParticles[j]->traj(CVT,12)->getX(), allParticles[j]->traj(CVT,12)->getY(), allParticles[j]->traj(CVT,12)->getZ() );

      if (traj12.X()==0 || traj12.Y()==0 || traj12.Z()==0) {continue;}

      // take the track that is closest to the neutron hit
      /*if (abs(traj12.Phi()*180./M_PI-n_phi)<angle_diff)
      {
        hit12_phi = traj12.Phi()*180./M_PI;
        angle_diff = abs(hit12_phi-n_phi);
      }*/

      // take the track that is closest in angle to the neutron hit
      if ( (p.Angle(traj12)*180./M_PI) < angle_diff)
      {
        hit12_phi = p.Angle(traj12)*180./M_PI;
        angle_diff = hit12_phi;
      }

    }



/*    // CND & CTOF NEARBY HITS (FROM REC::SCINTILLATOR)
    int hits_nearby7 = 0; int ctof_nearby7 = 0;
    double cluster_energy7 = 0; double ctof_energy7 = 0;
    for (int j=0; j<c12->getBank(rec_scint)->getRows(); j++)
    {
      int rec_detector = c12->getBank(rec_scint)->getInt(scint_detector,j);
      if (rec_detector!=3 && rec_detector!=4) {continue;}

      int rec_sector = c12->getBank(rec_scint)->getInt(scint_sector,j);
      int rec_layer = c12->getBank(rec_scint)->getInt(scint_layer,j);
      int rec_component = c12->getBank(rec_scint)->getFloat(scint_component,j);
      double rec_energy = c12->getBank(rec_scint)->getFloat(scint_energy,j);

      if (rec_detector==3 && (abs(rec_sector-sector)<3)) // hits in CND
      {
        hits_nearby7 = hits_nearby7 + 1;
        cluster_energy7 = cluster_energy7 + rec_energy;
      }
      else if (rec_detector==3 && (abs(rec_sector-sector)>21)) // hits in CND, boundary
      {
        hits_nearby7 = hits_nearby7 + 1;
        cluster_energy7 = cluster_energy7 + rec_energy;
      }
      else if (rec_detector==4 && (abs(rec_component-2*sector)<3)) // hits in CTOF //technically asymmetric
      {
        ctof_nearby7 = ctof_nearby7 + 1;
        ctof_energy7 = ctof_energy7 + rec_energy;
      }
      else if (rec_detector==4 && (abs(rec_component-2*sector)>44)) // hits in CTOF, boundary //technically asymmetric
      {
        ctof_nearby7 = ctof_nearby7 + 1;
        ctof_energy7 = ctof_energy7 + rec_energy;
      }
    }*/



/*    // CVT TRACKS
    //std::cout << "new event\n";
    for (int j=0; j<allParticles.size(); j++)
    {
      //TVector3 traj1( allParticles[j]->traj(CVT,1)->getX(), allParticles[j]->traj(CVT,1)->getY(), allParticles[j]->traj(CVT,1) );
      for (int k=1; k<13; k++)
      {
        //std::cout << "k = " << k << '\t';
        if (allParticles[j]->traj(CVT,k)->getLayer()!=0)
        {
          //if (allParticles[j]->traj(CVT,k)->getX()<1 || allParticles[j]->traj(CVT,k)->getX()>30) {continue;}
          //std::cout << "particle here\n";
          //std::cout << pow( pow(allParticles[j]->traj(CVT,k)->getX(),2) + pow(allParticles[j]->traj(CVT,k)->getY(),2), 0.5) << '\n';
          //if (allParticles[j]->trk(CVT)->getDetector()!=5) {continue;}
          //std::cout << allParticles[j]->trk(CVT)->getSector() << '\n';
          //std::cout << allParticles[j]->traj(CVT,k)->getDetector() << '\n';
          std::cout << allParticles[j]->traj(CVT,k)->getCx() << '\t';
          std::cout << allParticles[j]->traj(CVT,k)->getCy() << '\t';
          std::cout << allParticles[j]->traj(CVT,k)->getCz() << '\t';
          std::cout << pow( pow(allParticles[j]->traj(CVT,k)->getCx(),2) + pow(allParticles[j]->traj(CVT,k)->getCy(),2) + pow(allParticles[j]->traj(CVT,k)->getCz(),2)  ,0.5) << '\n';
        }
      }
        //std::cout << '\n';
    }
//std::cout << '\n';*/




    // Determine whether to write to "good nucleon" or "bad nucleon" file
    double cos0 = p_g.Dot(p) / (p_g.Mag()*p.Mag());
    bool good_N = (cos0>0.9 && p.Mag()>0.2 && abs(px-px_g)/px_g<0.2 && abs(py-py_g)/py_g<0.2 && abs(pz-pz_g)/pz_g<0.2 && abs(p.Mag()-p_g.Mag())/p_g.Mag()<0.1);
    bool bad_N = cos0>0.7 && (p.Mag()>0.2) && abs(p.Mag()-p_g.Mag())/p_g.Mag()<0.2;
    //bool bad_N = (cos0<0.8 || abs(px-px_g)>0.2 || abs(py-py_g)>0.2 || abs(pz-pz_g)>0.2 || abs(p.Mag()-p_g.Mag())>0.2);


    bool keep_this_one = (charge==0) ? good_N : bad_N;

    if (keep_this_one)
    {
      // all nucleons - print features
      outtxt << p_g.Mag() << ' ';
      outtxt << energy << ' ';
      outtxt << layermult << ' '; //////outtxt << z << ' ';
      outtxt << size << ' ';  /////outtxt << beta << ' ';
      outtxt << cnd_hits << ' ';
      outtxt << cnd_energy << ' ';
      outtxt << ctof_energy << ' ';
      outtxt << ctof_hits << ' ';
      outtxt << angle_diff << ' ';
      outtxt << '\n';
    }



  } // end loop over nucleons

    //chain.WriteEvent();
    // THIS WRITES ALL n AND p, NOT JUST GOOD N AND BAD N
    // set features to default values, then cut if not good n/bad n
    ntree->Fill();

    counter++;

  } // end loop over events

  std::cout << '\n' <<counter << " events counted!\n\n";



  // write histograms
  h_p->Write();
  h_pxminuspx->Write();
  h_pyminuspy->Write();
  h_pzminuspz->Write();
  h_pminusp->Write();
  h_pvsp->Write();
  h_dpp->Write();
  h_nangles->Write();
  h_px->Write();
  h_py->Write();
  h_pz->Write();
  h_cos0->Write();
  h_sec_phi->Write();
  h_cos1->Write();
  h_energy->Write();
  h_pg_theta->Write();


  outtxt.close();
  ntree->Write();
  f->Close();


  return 0;

} // closes main function



double getCVTdiff(std::vector<region_part_ptr> &allParticles, TVector3 &pn)
{
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

  /*double drdz = 0;
  double r_t1 = traj1.X()*traj1.X() + traj1.Y()*traj1.Y();
  double r_t3 = traj3.X()*traj3.X() + traj3.Y()*traj3.Y();
  double r_t5 = traj5.X()*traj5.X() + traj5.Y()*traj5.Y();
  double r_t7 = traj7.X()*traj7.X() + traj7.Y()*traj7.Y();
  double r_t12 = traj12.X()*traj12.X() + traj12.Y()*traj12.Y();
  drdz = drdz + (r_t1 - r_t3) / (traj1.Z() - traj3.Z());
  drdz = drdz + (r_t3 - r_t5) / (traj3.Z() - traj5.Z());
  drdz = drdz + (r_t5 - r_t7) / (traj5.Z() - traj7.Z());
  drdz = drdz + (r_t7 - r_t12) / (traj7.Z() - traj12.Z());*/
  }

   return angle_diff;
}
