#ifndef REWEIGHTER_HH
#define REWEIGHTER_HH

#include <vector>
#include <TGraphErrors.h>
#include <TRandom3.h>
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "clas12reader.h"
#include "HipoChain.h"

#include "gcfSRC.hh"
#include "eNCrossSection.hh"

#define REWEIGHTER_DIR _REWEIGHTER_DIR

using namespace std;
using namespace clas12;


class reweighter
{
public:
  reweighter(double E, int Z, int N, ffModel thisMod, char * input_uType_fin);
  ~reweighter();
  
  double get_weight_noT(clas12::mcparticle* mcInfo);
  double get_weight_ep(clas12::mcparticle* mcInfo);
  double get_weight_epp(clas12::mcparticle* mcInfo);
  double Gauss(double x, double mu, double sigma);
    
private:
  
  int Z_nuc;
  int N_nuc;
  double Ebeam;
  
  char * uType_init = new char;
  gcfSRC * gcf_config_init;
  double sigma_cm_init;
  eNCrossSection * CS_config_init;
  
  char * uType_fin = new char;
  gcfSRC * gcf_config_fin;
  double sigma_cm_fin;
  eNCrossSection * CS_config_fin;
  double P[4][2];
  double TN;
  double TNN;  
  
  
};

#endif

