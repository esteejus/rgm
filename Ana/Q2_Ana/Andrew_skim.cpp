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
  std::cerr << "Usage: ./code outputfile.hipo inputfiles.hipo  \n\n\n";

}



int main(int argc, char ** argv)
{

  if(argc < 3)
    {
      Usage();
      return -1;
    }


  //make c12writer for the output hipo file
  char * outName = argv[1];
  cout<<"Ouput file "<< outName <<endl;
  clas12root::HipoChainWriter chain(outName);
  for(int k = 2; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  //chain.GetWriter().writeSpecialBanks(true);
  chain.SetReaderTags({0});
  chain.db()->turnOffQADB();
  auto config_c12=chain.GetC12Reader();
  auto &c12=chain.C12ref();

  ////////////////////////////////////////////////

  int counter = 0;
  int cutcounter = 0;

  auto db=TDatabasePDG::Instance();
  double mass_p = db->GetParticle(2212)->Mass();
  double mD = 1.8756;

  double beam_E = 5.98636;

  //some particles
  TLorentzVector beam(0,0,beam_E,beam_E);
  TLorentzVector target(0,0,0,mD);
  TLorentzVector deut_ptr(0,0,0,mD);
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector lead_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector recoil_ptr(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector ntr(0,0,0,db->GetParticle(2112)->Mass());

  clas12ana clasAna;
  clasAna.printParams();
  
  while(chain.Next())
    {
      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cerr << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cerr << ".";
      }    

      clasAna.Run(c12);
      auto electrons = clasAna.getByPid(11);
      auto protons = clasAna.getByPid(2212);
      auto pip = clasAna.getByPid(211);
      auto pim = clasAna.getByPid(-211);

      if(electrons.size() == 1)
	{
          SetLorentzVector(el,electrons[0]);
	  TLorentzVector q = beam - el;
          double Q2 = -q.M2();
          double xB = Q2/(2 * mass_p * (beam.E() - el.E()));
	  if(xB<1){continue;}
	  if(Q2<1){continue;}
	  int lead_ctr = 0;
	  for(auto p = protons.begin(); p != protons.end();++p){	    
	    SetLorentzVector(lead_ptr,(*p));
	    TLorentzVector miss = q + deut_ptr - lead_ptr;
	    if(lead_ptr.P()<1){continue;}
	    //if(miss.P()<0.3){continue;}
	    lead_ctr++;
	  }
	  
	  if(lead_ctr==0){continue;}
	  chain.WriteEvent();	  
	}
    }
  return 0;
}

