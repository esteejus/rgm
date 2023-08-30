#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>
#include <sstream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include "HipoChain.h"
#include "HipoChainWriter.h"
#include "clas12ana.h"

using namespace std;
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());

}

void Usage()
{
  std::cerr << "Usage: ./testAna inputfiles.hipo outputfile.root \n\n\n";

}



int main(int argc, char ** argv)
{

  if(argc < 2)
    {
      Usage();
      return -1;
    }



  TString inFile  = argv[1];
  TString outFile = argv[2];

  cout<<"Ouput file "<< outFile <<endl;


  clas12ana clasAna;

  //Read in target parameter files                                                                                                                                                           
  clasAna.readInputParam("ana.par");
  clasAna.readEcalSFPar("paramsSF_40Ca_x2.dat");
  clasAna.readEcalPPar("paramsPI_40Ca_x2.dat");

  clasAna.printParams();


  clas12root::HipoChainWriter chain(outFile);
  chain.Add(inFile);

  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();

  //now get reference to (unique)ptr for accessing data in loop
  //this will point to the correct place when file changes
  //  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();

  int counter = 0;
  int cutcounter = 0;

  auto &c12=chain.C12ref();

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;

  double beam_E = 5.98;

  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  TH1D *q2_h = new TH1D("q2_h","Q^2 ",1000,0, 4);
  TH1D *xb_h = new TH1D("xb_h","x_B ",1000,0, 4);
  TH1D *px_com = new TH1D("px_com","Px COM",1000,-500,500);
  TH1D *py_com = new TH1D("py_com","Py COM",1000,-500,500);
  TH1D *pz_com = new TH1D("pz_com","Pz COM",1000,-500,500);

  TH1D *epp_h = new TH1D("epp_h","(e,e'pp)",100,0,2);
  TH1D *ep_h  = new TH1D("ep_h","(e,e'p)",100,0,2);


  clasAna.setEcalSFCuts();
  clasAna.setEcalEdgeCuts();
  clasAna.setPidCuts();
  clasAna.setVertexCuts();
  clasAna.setVertexCorrCuts();
  //  clasAna.setDCEdgeCuts();
  
  clasAna.setVzcuts(-6,1);
  clasAna.setVertexCorrCuts(-3,1);




  while(chain.Next())
    {
      //Display completed  
      counter++;

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);


      if(electrons.size() == 1)
	{
          SetLorentzVector(el,electrons[0]);

	  clasAna.getLeadRecoilSRC(beam,target,el);
	  auto lead    = clasAna.getLeadSRC();
	  auto recoil  = clasAna.getRecoilSRC();

	  if(lead.size() == 1)
	    chain.WriteEvent();

	}

    }


  return 0;
}

