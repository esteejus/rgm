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

using namespace std;
using namespace clas12;

void Usage() {
  std::cerr << "Usage: ./N_getfeatures charge output-root output-txt input-hipo\n";
}


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



  Int_t nhits;
  double px, py, pz, momentum; //double energy;
  //Float_t energy[100] = {-1};
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
  
  

  // MC banks
  auto mc_p = config_c12->addBank("MC::Particle");
  auto mc_px = config_c12->getBankOrder(mc_p,"px");
  auto mc_py = config_c12->getBankOrder(mc_p,"py");
  auto mc_pz = config_c12->getBankOrder(mc_p,"pz");

  // REC::Scintillator
  auto rec_scint = config_c12->addBank("REC::Scintillator");
  auto scint_detector = config_c12->getBankOrder(rec_scint,"detector");
  auto scint_sector = config_c12->getBankOrder(rec_scint,"sector");
  auto scint_layer = config_c12->getBankOrder(rec_scint,"layer");
  auto scint_component = config_c12->getBankOrder(rec_scint,"component");
  auto scint_energy = config_c12->getBankOrder(rec_scint,"energy");

  // ScintExtras
  auto rec_scintx = config_c12->addBank("REC::ScintExtras");
  //auto scint_dedx = config_c12->getBankOrder(rec_scintx,"dedx");
  auto scint_size = config_c12->getBankOrder(rec_scintx,"size");
  //auto scint_layermult = config_c12->getBankOrder(rec_scintx,"layermult");
  
 
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




  while(chain.Next())
  {

    // initialize features
    energy = 0; cnd_energy = 0; ctof_energy = 0; angle_diff = 180;
    layermult = 0; size = 0; cnd_hits = 0; ctof_hits = 0;


    // define particles
    TVector3 p_g(0.,0.,0.);
    TVector3 p(0.,0.,0.);
    TVector3 pe(0.,0.,0.);


    // identify particles from REC::Particle
    auto elec=c12->getByID(11);
    auto nucl = c12->getByID(2112); // looking for neutrons in e'n and e'p simulations
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
    momentum = p_g.Mag();



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
    bool is_CD = nucl[i]->getRegion()==CD;
    if (!is_CD) {continue;}

    // figure out what layer the hit is in - check 0 if not found
    bool is_CND1 = (nucl[i]->sci(CND1)->getLayer()==1);
    bool is_CND2 = (nucl[i]->sci(CND2)->getLayer()==2);
    bool is_CND3 = (nucl[i]->sci(CND3)->getLayer()==3);


    // put REC::Scintillator information
    double time;
    int sector;
    double beta = nucl[i]->par()->getBeta();


    // same as cluster information
    if (is_CND1)
    {
      sector = nucl[i]->sci(CND1)->getSector();
      //time =   nucl[i]->sci(CND1)->getTime() - starttime;
      energy = nucl[i]->sci(CND1)->getEnergy();
    }

    if (is_CND3)
    {
      sector = nucl[i]->sci(CND3)->getSector();
      //time =   nucl[i]->sci(CND3)->getTime() - starttime;
      energy = nucl[i]->sci(CND3)->getEnergy();
    }

    if (is_CND2)
    {
      sector = nucl[i]->sci(CND2)->getSector();
      //time =   nucl[i]->sci(CND2)->getTime() - starttime;
      energy = nucl[i]->sci(CND2)->getEnergy();
    }
    // PROBLEM: this gives preference to 2nd-layer hits
    if (!is_CND1 && !is_CND2 && !is_CND3)
    {
      sector = (nucl[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
      //time =   nucl[i]->sci(CTOF)->getTime() - starttime;
      energy = nucl[i]->sci(CTOF)->getEnergy();
    }


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


    // function for CND & CTOF nearby hits and energy
    Struct ninfo = getFeatures(nucl, allParticles, i);
    cnd_hits = ninfo.cnd_hits;
    ctof_hits = ninfo.ctof_hits;
    cnd_energy = ninfo.cnd_energy;
    ctof_energy = ninfo.ctof_energy;
    layermult = ninfo.layermult;
    energy = ninfo.energy;
    size = ninfo.size;
    angle_diff = ninfo.angle_diff;


    // Determine whether to write to "good nucleon" or "bad nucleon" file
    double cos0 = p_g.Dot(p) / (p_g.Mag()*p.Mag());
    bool good_N = (cos0>0.9 && p.Mag()>0.2 && abs(px-px_g)/px_g<0.2 && abs(py-py_g)/py_g<0.2 && abs(pz-pz_g)/pz_g<0.2 && abs(p.Mag()-p_g.Mag())/p_g.Mag()<0.1);
    bool bad_N = cos0>0.7 && (p.Mag()>0.2) && abs(p.Mag()-p_g.Mag())/p_g.Mag()<0.2;


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


      ntree->Fill();
    }



  } // end loop over nucleons

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
