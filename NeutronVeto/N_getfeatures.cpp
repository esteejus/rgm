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
    layermult = -1; size = 0; cnd_hits = 0; ctof_hits = 0;


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
    int sector = 0; int component = 0;
    double beta = nucl[i]->par()->getBeta();


    // same as cluster information
    if (is_CND1)
    {
      sector = nucl[i]->sci(CND1)->getSector();
      //time =   nucl[i]->sci(CND1)->getTime() - starttime;
      energy = nucl[i]->sci(CND1)->getEnergy();
      size =   nucl[i]->sci(CND1)->getSize();
    }

    if (is_CND3)
    {
      sector = nucl[i]->sci(CND3)->getSector();
      //time =   nucl[i]->sci(CND3)->getTime() - starttime;
      energy = nucl[i]->sci(CND3)->getEnergy();
      size =   nucl[i]->sci(CND3)->getSize();
    }

    if (is_CND2)
    {
      sector = nucl[i]->sci(CND2)->getSector();
      //time =   nucl[i]->sci(CND2)->getTime() - starttime;
      energy = nucl[i]->sci(CND2)->getEnergy();
      size =   nucl[i]->sci(CND2)->getSize();
    }
    // PROBLEM: this gives preference to 2nd-layer hits
    if (!is_CND1 && !is_CND2 && !is_CND3)
    {
      sector = (nucl[i]->sci(CTOF)->getComponent())/2; // rounded down, ctof component mapped onto cnd sector
      //time =   nucl[i]->sci(CTOF)->getTime() - starttime;
      energy = nucl[i]->sci(CTOF)->getEnergy();
      size =   nucl[i]->sci(CTOF)->getSize();
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
  }

   return angle_diff;
}
