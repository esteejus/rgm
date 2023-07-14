#ifndef VETO_FUNCTIONS_H
#define VETO_FUNCTIONS_H

#include "clas12reader.h"
#include "TVector3.h"

using namespace clas12;

double getCVTdiff(std::vector<region_part_ptr> neutron_list, std::vector<region_part_ptr>& allParticles_list, int i);

struct neutronInfo{
  double cnd_hits;
  double ctof_hits;
  double layermult;
  double size;
  double cnd_energy;
  double ctof_energy;
  double energy;
  double angle_diff;
};

typedef struct neutronInfo Struct;

Struct getFeatures(std::vector<region_part_ptr> neutron_list, std::vector<region_part_ptr>& allParticles_list, int i);



#endif
