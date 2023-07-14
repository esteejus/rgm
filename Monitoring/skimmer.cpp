#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <vector>
#include <typeinfo>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "clas12reader.h"
#include "clas12writer.h"
#include "HipoChain.h"
#include "eventcut/eventcut.h"
//#include "functions.h"

using namespace std;
using namespace clas12;

void Usage()
{
  std::cerr << "Usage: ./code <Ebeam(GeV)> <path/to/cutfile.txt> <path/to/output.hipo> <path/to/input.hipo> \n";
}


int main(int argc, char ** argv)
{

  if(argc < 5)
    {
      std::cerr<<"Wrong number of arguments.\n";
      Usage();
      return -1;
    }

  /////////////////////////////////////
  //Set cut object with Ebeam and cutfile
  eventcut myCut(atof(argv[1]),argv[2]);
  myCut.print_cuts();

  //make c12writer for the output hipo file
  char * outName = argv[3];
  cout<<"Ouput file "<< outName <<endl;
  clas12writer c12writer(outName);

  //make hipochain that contains input files
  clas12root::HipoChain chain;
  for(int k = 4; k < argc; k++){
    cout<<"Input file "<<argv[k]<<endl;
    chain.Add(argv[k]);
  }
  auto config_c12=chain.GetC12Reader();
  chain.SetReaderTags({0});

  
  //Additional setup of reader and writer
  //get pointer to check for file change
  auto currc12=chain.GetC12Reader();

  //now get reference to (unique)ptr for accessing data in loop
  //this will point to the correct place when file changes
  const std::unique_ptr<clas12::clas12reader>& c12=chain.C12ref();
  chain.db()->turnOffQADB();
  int counter = 0;
  int cutcounter = 0;



/*// CND hits
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
auto ctof_clus_suze = config_c12->getBankOrder(ctof_clusters,"size");
auto ctof_clus_sector = config_c12->getBankOrder(ctof_clusters,"sector");
auto ctof_clus_layer = config_c12->getBankOrder(ctof_clusters,"layer");
auto ctof_clus_component = config_c12->getBankOrder(ctof_clusters,"component");
auto ctof_clus_energy = config_c12->getBankOrder(ctof_clusters,"energy");
auto ctof_clus_time = config_c12->getBankOrder(ctof_clusters,"time");


// ScintExtras
auto scintextras = config_c12->addBank("RECHB::ScintExtras");
auto scint_dedx = config_c12->getBankOrder(scintextras,"dedx");
auto scint_size = config_c12->getBankOrder(scintextras,"size");
auto scint_layermult = config_c12->getBankOrder(scintextras,"layermult");*/


  

  while(chain.Next()==true){

      //Display completed  
      counter++;
      if((counter%1000000) == 0){
	cout << "\n" <<counter/1000000 <<" million completed";
      }    
      if((counter%100000) == 0){
	cout << ".";
      }    

      //if multiple files in chain
      //we need to update when file changes
      if(currc12!=c12.get()){
	currc12=c12.get();
	//assign a reader to the writer
	c12writer.assignReader(*currc12);
      }

  /////////////////////////////////////
  //Electron fiducials and Pid
  //Lead Proton Checks
  //Lead SRC Proton Checks
  //Recoil Proton Checks
  /////////////////////////////////////      
      if(!myCut.electroncut(c12)){continue;}      
      int index_L = myCut.leadnucleoncut(c12);
      if(index_L < 0){ continue; }
      if(!myCut.leadSRCnucleoncut(c12,index_L)){continue;}      
      int index_R = myCut.recoilSRCnucleoncut(c12,index_L);
      if(index_R < 0){ continue; }
      
      cutcounter++;
      c12writer.writeEvent(); 
  }

  c12writer.closeWriter();
  cout<<cutcounter<<" events written to:\n" << outName <<endl;
}

