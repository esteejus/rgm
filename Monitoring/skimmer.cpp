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
#include "eventcut.h"
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

