#include "veto_functions.h"

//using namespace std;
//using namespace clas12;



double getCVTdiff(std::vector<region_part_ptr> neutron_list, std::vector<region_part_ptr> &allParticles_list, int i)
{
  double hit12_phi = 180;
  double angle_diff = 360;

  // get neutron momentum
  TVector3 pn;
  pn.SetXYZ( neutron_list[i]->par()->getPx(), neutron_list[i]->par()->getPy(), neutron_list[i]->par()->getPz() );

  for (int j=0; j<allParticles_list.size(); j++)
  {
    // want k=1,3,5,7,12
    TVector3 traj1( allParticles_list[j]->traj(CVT,1)->getX(), allParticles_list[j]->traj(CVT,1)->getY(), allParticles_list[j]->traj(CVT,1)->getZ() ); 
    TVector3 traj3( allParticles_list[j]->traj(CVT,3)->getX(), allParticles_list[j]->traj(CVT,3)->getY(), allParticles_list[j]->traj(CVT,3)->getZ() );
    TVector3 traj5( allParticles_list[j]->traj(CVT,5)->getX(), allParticles_list[j]->traj(CVT,5)->getY(), allParticles_list[j]->traj(CVT,5)->getZ() );
    TVector3 traj7( allParticles_list[j]->traj(CVT,7)->getX(), allParticles_list[j]->traj(CVT,7)->getY(), allParticles_list[j]->traj(CVT,7)->getZ() );
    TVector3 traj12( allParticles_list[j]->traj(CVT,12)->getX(), allParticles_list[j]->traj(CVT,12)->getY(), allParticles_list[j]->traj(CVT,12)->getZ() );


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




Struct getFeatures(std::vector<region_part_ptr> neutron_list, std::vector<region_part_ptr>& allParticles_list, int i)
{

  // initialize variables to return
  Struct info;
  info.cnd_hits = 0;
  info.ctof_hits = 0;
  info.cnd_energy = 0;
  info.ctof_energy = 0;
  info.energy = 0;
  info.layermult = 0;
  info.size = 0;
  info.angle_diff = getCVTdiff(neutron_list, allParticles_list, i);


  // determine which CND layer(s) neutron is in
  bool n_isCND1 = (neutron_list[i]->sci(CND1)->getLayer()==1);
  bool n_isCND2 = (neutron_list[i]->sci(CND2)->getLayer()==2);
  bool n_isCND3 = (neutron_list[i]->sci(CND3)->getLayer()==3);
  bool n_isCND = (n_isCND1 || n_isCND2 || n_isCND3);
  bool n_isCTOF = (neutron_list[i]->sci(CTOF)->getDetector()==4);


  // initialize subdetector/layer-dependent quantities
  double n_phi = -360; // neutron phi (range -176.25 to 176.25 degrees) -- they occur at intervals of exactly 7.5 degrees :D
  int sector = -1; // sector in which neutron is localized
  // define subdetector/layer-dependent quantities
  if (n_isCND1)
  {
    n_phi = atan2(neutron_list[i]->sci(CND1)->getY(),neutron_list[i]->sci(CND1)->getX())*180/M_PI;
    sector = neutron_list[i]->sci(CND1)->getSector();
    info.energy = neutron_list[i]->sci(CND1)->getEnergy();
    info.size = neutron_list[i]->sci(CND1)->getSize();
    info.layermult = info.layermult + 1;
  }
  if (n_isCND2)
  {
    n_phi = atan2(neutron_list[i]->sci(CND2)->getY(),neutron_list[i]->sci(CND2)->getX())*180/M_PI;
    sector = neutron_list[i]->sci(CND2)->getSector();
    info.energy = neutron_list[i]->sci(CND2)->getEnergy();
    info.size = neutron_list[i]->sci(CND2)->getSize();
    info.layermult = info.layermult + 1;
  }
  if (n_isCND3)
  {
    n_phi = atan2(neutron_list[i]->sci(CND3)->getY(),neutron_list[i]->sci(CND3)->getX())*180/M_PI;
    sector = neutron_list[i]->sci(CND3)->getSector();
    info.energy = neutron_list[i]->sci(CND3)->getEnergy();
    info.size = neutron_list[i]->sci(CND3)->getSize();
    info.layermult = info.layermult + 1;
  }
  if (n_isCTOF)
  {
    n_phi = atan2(neutron_list[i]->sci(CTOF)->getY(),neutron_list[i]->sci(CTOF)->getX())*180/M_PI;
    sector = neutron_list[i]->sci(CTOF)->getComponent();
    info.energy = neutron_list[i]->sci(CTOF)->getEnergy();
    info.size = neutron_list[i]->sci(CTOF)->getSize();
  }


  // for all particles, look for hits near neutron
  for (int j=0; j<allParticles_list.size(); j++)
  {
    // skip particles that are not in CND or CTOF
    bool part_isCND1 = (allParticles_list[j]->sci(CND1)->getLayer()==1);
    bool part_isCND2 = (allParticles_list[j]->sci(CND2)->getLayer()==2);
    bool part_isCND3 = (allParticles_list[j]->sci(CND3)->getLayer()==3);
    bool part_isCND = (part_isCND1 || part_isCND2 || part_isCND3);
    bool part_isCTOF = (allParticles_list[j]->sci(CTOF)->getDetector()==4);
    if ( !part_isCND && !part_isCTOF) {continue;}

    double part_phi = -360;
    if (part_isCND1) {part_phi = atan2(allParticles_list[j]->sci(CND1)->getY(),allParticles_list[j]->sci(CND1)->getX())*180/M_PI;}
    if (part_isCND2) {part_phi = atan2(allParticles_list[j]->sci(CND2)->getY(),allParticles_list[j]->sci(CND2)->getX())*180/M_PI;}
    if (part_isCND3) {part_phi = atan2(allParticles_list[j]->sci(CND3)->getY(),allParticles_list[j]->sci(CND3)->getX())*180/M_PI;}
    if (part_isCTOF) {part_phi = atan2(allParticles_list[j]->sci(CTOF)->getY(),allParticles_list[j]->sci(CTOF)->getX())*180/M_PI;}


    // look for nearby CND and CTOF hits
    double phi_diff = abs(part_phi-n_phi);
    double tolerance = 30+1; // angular range (degrees) within which to look for hits
    if (part_isCND1 && (phi_diff<tolerance || phi_diff>(360-tolerance)) )
    {
      info.cnd_hits = info.cnd_hits + allParticles_list[j]->sci(CND1)->getSize();
      info.cnd_energy = info.cnd_energy + allParticles_list[j]->sci(CND1)->getEnergy();
    }
    if (part_isCND2 && (phi_diff<tolerance || phi_diff>(360-tolerance)) )
    {
      info.cnd_hits = info.cnd_hits + allParticles_list[j]->sci(CND2)->getSize();
      info.cnd_energy = info.cnd_energy + allParticles_list[j]->sci(CND2)->getEnergy();
    }
    if (part_isCND3 && (phi_diff<tolerance || phi_diff>(360-tolerance)) )
    {
      info.cnd_hits = info.cnd_hits + allParticles_list[j]->sci(CND3)->getSize();
      info.cnd_energy = info.cnd_energy + allParticles_list[j]->sci(CND3)->getEnergy();
    }
    if (part_isCTOF && (phi_diff<tolerance || phi_diff>(360-tolerance)) )
    {
      info.ctof_hits = info.ctof_hits + allParticles_list[j]->sci(CTOF)->getSize();
      info.ctof_energy = info.ctof_energy + allParticles_list[j]->sci(CTOF)->getEnergy();
    }
  } // end loop over all particles

  return info;

} // end function
