#include "clas12ana.h"

struct cutpar{
  std::string id;
  std::vector<double> par = {}; //pi- parameters
};

//helper function for DC fiducials                                                           
TVector3 rotate(TVector3 vec, int sector)
{
  double rot_ang = -(sector -1)*60 *TMath::DegToRad();

  vec.RotateZ(rot_ang);

  return vec;
}


int clas12ana::getCDRegion(const region_part_ptr &p)
{
  int region = -1;

 //Only defined for protons right now!!!!
 if(p->getRegion() == CD) 
    {
      auto px = p -> par() -> getPx();
      auto py = p -> par() -> getPy();
      double pt = sqrt( pow(px,2) + pow(py,2) );
      double fiducial_phi = (-asin(min_mom_pt/pt) - (pi/2)) * 180/pi;
      double phi = p->getPhi() * 180/pi;
      
      if( phi > fiducial_phi &&  phi < (fiducial_phi+120))
	region = 1;
      else if( (phi > (fiducial_phi+120)) &&  (phi < (fiducial_phi+240)) )
	region = 2;
      else
	region = 3;
    }

  return region;
}



void clas12ana::Clear()
 {
   electrons.clear();
   protons.clear();
   deuterons.clear();
   neutrals.clear();
   neutrons.clear();
   piplus.clear();
   piminus.clear();
   kplus.clear();
   kminus.clear();
   otherpart.clear();

   lead_proton.clear();
   recoil_proton.clear();

   current_run = -1;
   beam_energy = 0;
   event_mult = 0;
 }


void clas12ana::Run(const std::unique_ptr<clas12::clas12reader>& c12)
{
  Clear();
  current_run = c12->runconfig()->getRun();
  checkCutParameters(); //check run number has the right cuts 



  auto particles = c12->getDetParticles(); //particles is now a std::vector of particles for this event
  auto electrons_det = c12->getByID(11);

  //DEBUG plots
  if(debug_plots)
    {
      for(auto el : electrons_det)
	debug_c.fillBeforeEl(el);
    }
  
  std::for_each(electrons_det.begin(),electrons_det.end(),[this](auto el)
		{
		  if(!((el->che(HTCC)->getNphe() <= 2)           || //Photo electron min cut
		       (!checkEcalSFCuts(el)   && f_ecalSFCuts)  || //ECAL SF cuts
		       (!checkEcalPCuts(el)    && f_ecalPCuts)   || //ECAL SF cuts
		       (!checkEcalDiagCuts(el) && f_ecalDiagCuts)|| //ECAL Diagonoal SF cuts
		       (!EcalEdgeCuts(el) && f_ecalEdgeCuts)     || //ECAL edge cuts
		       (!checkVertex(el)  && f_vertexCuts)       || //Vertex cut
		       (!DCEdgeCuts(el)   && f_DCEdgeCuts)       || //DC edge cut
		       (el->par()->getP() < 0.8)) ) // minium 800 MeV/c cut for electrons in class     
		    setByPid(el);
		});
  
  
  if(debug_plots)
    {
      for(auto el : electrons)
	debug_c.fillAfterEl(el);
    }
  
  
  if(electrons.size() == 1) //good trigger electron
    {
      
      if(debug_plots)
	{
	  for(auto p : particles)
	    if(p->par()->getPid() == 2212 || p->par()->getPid() == -211 || p->par()->getPid() == 211)
	      debug_c.fillBeforePart(p);
	}
      
      
      /*
	This may be a strange way to check the cuts, maybe there is a better way
	We need to ensure that the flag f_cuts is on otherwise we don't want to apply any cut
	The below logic will return particles that did not pass any cut (for only cut flags that are on
	Then I will invert the logic using !(logic) to return when the particle do pass all cuts
	
	(!checkPidCut(p) && f_pidCuts)    ||	         //PID cuts
	(!checkVertex(p) && f_vertexCuts) ||  //Vertex cut
	(!CDEdgeCuts(p)  && f_CDEdgeCuts) ||  //CD edge cut
	(!CDRegionCuts(p)  && f_CDRegionCuts) ||  //CD edge cut
	(!DCEdgeCuts(p) && f_DCEdgeCuts)  || //DC edge cut
	(!checkVertexCorrelation(electrons_det[0],p) && f_corr_vertexCuts) //Vertex correlation cut between electron
      */
      
      std::for_each(particles.begin(),particles.end(),[this,electrons_det](auto p)
		    {
		      //neutrals and electrons don't follow cuts below, skip them 
		      if(p->par()->getCharge() == 0 && p->par()->getPid() != 11 )
			{
			  setByPid(p);
			  return;
			}
		      else if(p->par()->getPid() != 11  && electrons.size() > 0)
			{
			  ++event_mult;	//charge particles


			  bool check_pid_cuts = ((f_protonpidCuts && checkProtonPidCut(p)) || //check if in proton PID cuts or chi2pid cuts
						 (f_pidCuts && checkPidCut(p))             || // if proton pid cuts if off but pid cuts on just use chi2pid
						 (!f_protonpidCuts && !f_pidCuts));           // if no pid cuts are specified let all particles pass pid 
			  
			  if( !( (!check_pid_cuts)                    || //PID cuts
				 (!checkVertex(p)  && f_vertexCuts)   || //Vertex cut
				 (!CDEdgeCuts(p)   && f_CDEdgeCuts)   || //CD edge cut
				 (!CDRegionCuts(p) && f_CDRegionCuts) || //CD edge cut
				 (!DCEdgeCuts(p)   && f_DCEdgeCuts)   || //DC edge cut
				 (!checkVertexCorrelation(electrons_det[0],p) && f_corr_vertexCuts)) ) //Vertex correlation cut between electron
			    setByPid(p);
			}
		    });
      
      
      if(debug_plots)
	{
	  
	  for(auto p : protons)
	    debug_c.fillAfterPart(p);
	  for(auto p : piplus)
	    debug_c.fillAfterPart(p);
	  for(auto p : piminus)
	    debug_c.fillAfterPart(p);
	  
	  for(auto el : electrons)
	    debug_c.fillAfterEl(el);
	  
	}

    }//good electron loop
  
  
}



void clas12ana::InitSFEcalCuts()
{
  //  cout<<"PARAMETERS for SF vs Ecal cuts"<<endl;
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  if(j < 4)
	    ecal_sf_mean_fcn[i]->SetParameter(j,ecal_sf_fcn_par[i][j]); //only first 3 parameters involve mean

	  //	  cout<<"sector "<<i <<" j "<<j<<" par "<<ecal_sf_fcn_par[i][j]<<endl;
	  ecal_sf_fcn[0][i]->SetParameter(j,ecal_sf_fcn_par[i][j]);
	  ecal_sf_fcn[1][i]->SetParameter(j,ecal_sf_fcn_par[i][j]);
	}

      ecal_sf_fcn[0][i]->SetParameter(6,sigma_cut);
      ecal_sf_fcn[1][i]->SetParameter(6,sigma_cut);
    }
}



void clas12ana::WriteSFEcalCuts()
{
  TFile file_ecal("ecalSFcuts.root","RECREATE");
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  ecal_sf_fcn[0][i]->Write();
	  ecal_sf_fcn[1][i]->Write();
	  ecal_p_fcn[0][i]->Write();
	  ecal_p_fcn[1][i]->Write();

	}

      ecal_sf_mean_fcn[i]->Write();
      ecal_p_mean_fcn[i]->Write();
    }

   file_ecal.Close();
}




void clas12ana::InitSFPCuts()
{
  //   cout<<"PARAMETERS for SF vs P cuts"<<endl;
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  if(j < 4)
	    ecal_p_mean_fcn[i]->SetParameter(j,ecal_p_fcn_par[i][j]); //only first 3 parameters for mean value

	  //	  cout<<"sector "<<i <<" j "<<j<<" par "<<ecal_p_fcn_par[i][j]<<endl;
	  ecal_p_fcn[0][i]->SetParameter(j,ecal_p_fcn_par[i][j]);
	  ecal_p_fcn[1][i]->SetParameter(j,ecal_p_fcn_par[i][j]);
	}

      ecal_p_fcn[0][i]->SetParameter(6,sigma_cut);
      ecal_p_fcn[1][i]->SetParameter(6,sigma_cut);
    }


}

void clas12ana::checkCutParameters()
{
  //This next section sets the Sampling Fraction cuts; pass1 data has two SF regions
  //< 15542 and >=15542 based on the SF timelines
  //Defualt SF cuts are set in Init() function and is the < 15542 events
  if(current_run >= 15542)
    {
      if(previous_run < 15542) //default; do nothing already been set to >=15542
	{
	  //Set new SF cuts, new run range
	  this -> readEcalSFPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsSF_40Ca_x2.dat").c_str() );
	  this -> readEcalPPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsPI_40Ca_x2.dat").c_str());
	  std::cerr << "WARNING:: Run number changed to " << current_run <<". The SF cuts are changed to reflect this new run range" << std::endl;
	}
    }

  //note run ranges cover all possible ranges for a given target
  // the production runs is some subet of the full range
  //Here we check the previous run and the current run, if the current run differs
  //then we need to set a new set of par files
  bool h_runrange    = (current_run >= 15016 && current_run <= 15042);
  bool d_runrange    = ((current_run >= 15043 && current_run <= 15106) || (current_run >= 15433 && current_run <= 15456));
  bool he_runrange   = ((current_run >= 15108 && current_run <= 15164) || (current_run >= 15458 && current_run <= 15490));
  bool cx4_runrange  = (current_run >= 15178 && current_run <= 15317);
  bool snx4_runrange = (current_run >= 15318 && current_run <= 15328);
  bool ca40_runrange = (current_run >= 15355 && current_run <= 15432);
  bool ca48_runrange = (current_run >= 15829 && current_run <= 15884);

  bool h_runrange_prev    = (previous_run >= 15016 && previous_run <= 15042);
  bool d_runrange_prev    = ((previous_run >= 15043 && previous_run <= 15106) || (previous_run >= 15433 && previous_run <= 15456));
  bool he_runrange_prev   = ((previous_run >= 15108 && previous_run <= 15164) || (previous_run >= 15458 && previous_run <= 15490));
  bool cx4_runrange_prev  = (previous_run >= 15178 && previous_run <= 15317);
  bool snx4_runrange_prev = (previous_run >= 15318 && previous_run <= 15328);
  bool ca40_runrange_prev = (previous_run >= 15355 && previous_run <= 15432);
  bool ca48_runrange_prev = (previous_run >= 15829 && previous_run <= 15884);


  //for some reason first 3-4 events are always run 0
  //they seem to always have 0 particles so maybe they are headers; do nothing here
  if(current_run == 0)
    return;

  //MC runs are always run 11 per CLAS default
  //We assume the USER must supply the parameter file inline in the analysis code
  //Since the MC run# is always 11 there is no automatic feature
  
  if(he_runrange && !he_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_he4.par file." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_he4.par").c_str() );
    }
  else if(ca40_runrange && !ca40_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_ca40.par file." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_ca40.par").c_str() );
    }
  else if(ca48_runrange && !ca48_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_ca48.par file." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_ca48.par").c_str() );
    }
  else if(cx4_runrange && !cx4_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_cx4.par file." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_cx4.par").c_str() );
    }
  else if(snx4_runrange && !snx4_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_cx4.par file which is the same as the snx4 cell." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_cx4.par").c_str() );
    }
  else if(d_runrange && !d_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_he4.par file. Which is std. liquid cell default." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_he4.par").c_str() );
    }
  else if(h_runrange && !h_runrange_prev)
    {
      std::cerr << "WARNING:: Run range changed for run " << current_run << ". Setting ana_he4.par file. Which is std. liquid cell default." << std::endl;
      this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_he4.par").c_str() );
    }

  previous_run = current_run;

}

void clas12ana::Init()
{
  if(debug_plots)
    debug_c.InitDebugPlots();

  proton_pid_mean->SetParameters(0.0152222,0.816844,-0.0950375,0.255628);
  proton_pid_sigma->SetParameters(0.0760525,0.240862,-0.000276433,0.229085);

  for(int i = 0; i < 7; i++)
    {
      ecal_sf_mean_fcn[i] = new TF1(Form("ecal_sf_mean_fcn_%d",i),"[0] + [1]/x + [2]/pow(x,2)",0,1.5);
      ecal_sf_fcn[0][i] = new TF1(Form("ecal_sf_fcn_0_%d",i),"[0] + [1]/x + [2]/pow(x,2) - [6]*( [3] + [4]/x + [5]/pow(x,2))",0,1.5);
      ecal_sf_fcn[1][i] = new TF1(Form("ecal_sf_fcn_1_%d",i),"[0] + [1]/x + [2]/pow(x,2) + [6]*( [3] + [4]/x + [5]/pow(x,2))",0,1.5);
      
      ecal_p_mean_fcn[i] = new TF1(Form("ecal_p_mean_fcn_%d",i),"[0] + [1]/x + [2]/pow(x,2)",0,10);
      ecal_p_fcn[0][i] = new TF1(Form("ecal_p_fcn_0_%d",i),"[0] + [1]/x + [2]/pow(x,2) - [6]*( [3] + [4]/x + [5]/pow(x,2))",0,10);
      ecal_p_fcn[1][i] = new TF1(Form("ecal_p_fcn_1_%d",i),"[0] + [1]/x + [2]/pow(x,2) + [6]*( [3] + [4]/x + [5]/pow(x,2))",0,10);
    }
  
  for(int i = 0; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  if(j == 3)
	    {
	      ecal_sf_fcn_par[i][j] = 9999; 
	      ecal_p_fcn_par[i][j] = 9999; 
	    }

	  ecal_sf_fcn_par[i][j] = 0;
	  ecal_p_fcn_par[i][j] = 0;
	}
    }


  //As defualt load 4He analysis cuts and the SF cuts fit on liquid deuterium which apply to runs < 15542
  previous_run = 15108; //set to a defualt helium run
  this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana_he4.par").c_str() );
  this -> readEcalSFPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsSF_LD2_x2.dat").c_str() );
  this -> readEcalPPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsPI_LD2_x2.dat").c_str());

  this -> readInputSRCParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/src_cuts.par").c_str() );
  //  this -> printParams();

}


bool clas12ana::DCEdgeCuts(const region_part_ptr &p)
{
  //true if inside cut
  //cut all charged particles
  if(p->par()->getCharge() != 0)
    {
      auto traj_index_1 = p->traj(DC,6)->getIndex();  //layer 1 
      auto traj_index_2 = p->traj(DC,18)->getIndex(); //layer 2
      auto traj_index_3 = p->traj(DC,36)->getIndex(); //layer 3

      auto traj_edge_1  = p->traj(DC,6)->getFloat("edge",traj_index_1);
      auto traj_edge_2  = p->traj(DC,18)->getFloat("edge",traj_index_2);
      auto traj_edge_3  = p->traj(DC,36)->getFloat("edge",traj_index_3);

      //PUT DC EDGE CUTS IN PARAMETER FILE

      //electron DC cuts
      if(p->par()->getCharge() < 0 && (dc_edge_cut_el.size() == 3  && traj_edge_1 > dc_edge_cut_el[0] && traj_edge_2 > dc_edge_cut_el[1] && traj_edge_3 > dc_edge_cut_el[2]) )
	return true;
      //proton DC cuts
      else if(p->par()->getCharge() > 0 && (dc_edge_cut_ptr.size() == 3  && traj_edge_1 > dc_edge_cut_ptr[0] && traj_edge_2 > dc_edge_cut_ptr[1] && traj_edge_3 > dc_edge_cut_ptr[2]) )
	return true;
      else
	return false;
    }
  else
    return true;
}



bool clas12ana::EcalEdgeCuts(const region_part_ptr &p)
{
  //true if inside cut
  double sampling_frac = getSF(p);

  if(p->par()->getPid() == 11)
    {
      if(p->cal(clas12::PCAL)->getLv() > ecal_edge_cut && p->cal(clas12::PCAL)->getLw() > ecal_edge_cut)
	return true;
      else
	return false;
    }
  
  else
    return true;
}


bool clas12ana::checkGhostTrackCD(const region_part_ptr &p)
{
  /*
    Function returns true if track is a suspected ghost track

    There are two sources of ghost tracks. A ghost track is a particle which gets reconstructed twice
    one of the tracks should only be kept. Two cases have been identified:
    1. CD protons have the exact same CTOF hit component (remove by removing if a tracks shares same CTOF component in CD
    2. Proton reconstructed by CD and FD separately, they have < 5deg angle between tracks. 
  */

  //check ghost tracks only apply to charge particles. Investigated for protons. Need to check pions
  if(p->par()->getCharge() == 0)
    return false;

  for(auto &p2 : protons)
    {
      //case 1
      if(p->getRegion() == clas12::CD && p2->getRegion() == clas12::CD 
	 && p2->sci(CTOF)->getComponent() == p->sci(CTOF)->getComponent())
	return true;
      
      //case 2
      else if( ((p->getRegion() == clas12::FD && p2->getRegion() == clas12::CD) || (p->getRegion() == clas12::CD && p2->getRegion() == clas12::FD)) &&
	       abs(p2->getTheta() - p->getTheta())*TMath::RadToDeg() < ghost_track_cut )
	return true;
    }


  return false;
}

bool clas12ana::checkEcalDiagCuts(const region_part_ptr &p)
{
  double mom = p->par()->getP();
  //true if inside cut
  if(p->par()->getPid() == 11)
    {
      if( (p->cal(clas12::PCAL)->getEnergy() + p->cal(clas12::ECIN)->getEnergy())/mom > ecal_diag_cut && mom > 4.5)
	return true;
      else if(mom <= 4.5)
	return true;
      else
	return false;
    }
  
  else
    return true;
}


bool clas12ana::checkEcalSFCuts(const region_part_ptr &p)
{
  //true if inside cut

  if(p->par()->getPid() == 11)
    {
      double sampling_frac = getSF(p);
      double energy =  p->cal(PCAL)->getEnergy();

      if(energy < pcal_energy_cut)
	return false;

      int sector = p->getSector();
      
      //Turn on for functional form 
      double sf_max_cut = ecal_sf_fcn[1][sector]->Eval(energy);
      double sf_min_cut = ecal_sf_fcn[0][sector]->Eval(energy);
      
      if(sampling_frac < sf_max_cut && sampling_frac > sf_min_cut)
	return true;

      else
	return false;
    }
  
  else
    return false;
}


bool clas12ana::checkEcalPCuts(const region_part_ptr &p)
{
  //true if inside cut

  if(p->par()->getPid() == 11)
    {
      double sampling_frac = getSF(p);
      int sector = p->getSector();
      
      //Turn on for functional form 
      double sf_max_cut = ecal_p_fcn[1][sector]->Eval(p->par()->getP() );
      double sf_min_cut = ecal_p_fcn[0][sector]->Eval(p->par()->getP() );
      
      if(sampling_frac < sf_max_cut && sampling_frac > sf_min_cut)
	return true;

      else
	return false;
    }
  
  else
    return false;
}



double clas12ana::getSF(const region_part_ptr &p)
{
  if(p->par()->getPid() == 11)
    return (p->cal(clas12::PCAL)->getEnergy() +  p->cal(clas12::ECIN)->getEnergy() +  p->cal(clas12::ECOUT)->getEnergy()) / p->par()->getP();
  else
    return -9999.;
}


bool clas12ana::CDRegionCuts(const region_part_ptr &p)
{
  //true if inside cut
  //cut all charged particles
  if(p->par()->getCharge() != 0 &&  p->getRegion() == CD) //neutral particles don't follow cuts
    {
      auto px = p -> par() -> getPx();
      auto py = p -> par() -> getPy();
      double pt = sqrt( pow(px,2) + pow(py,2) );
      double fiducial_phi = (-asin(min_mom_pt/pt) - (pi/2)) * 180/pi;
      double phi = p->getPhi() * 180/pi;
      
      int region = -1;

      if( phi > fiducial_phi &&  phi < (fiducial_phi+120))
	region = 1;
      else if( (phi > (fiducial_phi+120)) &&  (phi < (fiducial_phi+240)) )
	region = 2;
      else
	region = 3;

      if(region == region_cut)
	return false; //inside bad region
      else
	return true; //inside good region CD
    }
  else
    return true; //neutrals dont apply
}





bool clas12ana::CDEdgeCuts(const region_part_ptr &p)
{
  //true if inside cut
  //cut all charged particles
  if(p->par()->getCharge() != 0 &&  p->getRegion() == CD) //neutral particles don't follow cuts
    {
      double edge_first = p->traj(CVT,7)->getEdge();
      TVector3 hit_first(p->traj(CVT,7)->getX(),p->traj(CVT,7)->getY(),p->traj(CVT,7)->getZ());
      double hp_first = hit_first.Phi()*180/M_PI;
      int hit_reg_first = hp_first<-90?1:hp_first<30?2:hp_first<150?3:1;
      
      double edge_last = p->traj(CVT,12)->getEdge();
      TVector3 hit_last(p->traj(CVT,12)->getX(),p->traj(CVT,12)->getY(),p->traj(CVT,12)->getZ());
      double hp_last = hit_last.Phi()*180/M_PI;
      int hit_reg_last = hp_last<-90?1:hp_last<30?2:hp_last<150?3:1;

      if(!((edge_first>cd_edge_cut) && (edge_last>cd_edge_cut) && (hit_reg_first == hit_reg_last))){
	return false;} //inside bad region
      else{
	return true;} //inside good region CD
    }
  else
    return true; //neutrals dont apply
}

bool clas12ana::checkVertexCorrelation(const region_part_ptr &el, const region_part_ptr &p)
{
  //true if inside cut
  if(p->getRegion() == clas12::CD)
    return ( (el->par()->getVz() - p->par()->getVz()) > vertex_corr_cuts_cd.at(0) &&  (el->par()->getVz() - p->par()->getVz()) < vertex_corr_cuts_cd.at(1) );
  else if(p->getRegion() == clas12::FD)
    return ( (el->par()->getVz() - p->par()->getVz()) > vertex_corr_cuts_fd.at(0) &&  (el->par()->getVz() - p->par()->getVz()) < vertex_corr_cuts_fd.at(1) );
  else
    return true;
}

bool clas12ana::checkVertex(const region_part_ptr &p)
{
  //function returns true if inside vertex cuts
  int pid = p->par()->getPid();

  bool in_vxvy = (p->par()->getVx() > vertex_x_cuts.at(0) && p->par()->getVx() < vertex_x_cuts.at(1))
    && (p->par()->getVy() > vertex_y_cuts.at(0) && p->par()->getVy() < vertex_y_cuts.at(1));

  if(!in_vxvy)
    return false;    

  //need to change the PID for protons which are idientified with TOF cuts if turned on
  //otherwise vertex cuts will not properly be done

  if(checkProtonPidCut(p) && f_protonpidCuts)
    pid = 2212;

  if(p->getRegion() == FD) //forward detector cuts
    {
      auto itter = vertex_z_cuts_fd.find(pid);
      if(itter != vertex_z_cuts_fd.end())
	return p->par()->getVz() > itter->second.at(0) && p->par()->getVz() < itter->second.at(1);
      else
	return true;
    }
  else if(p->getRegion() == CD) //central detector cuts
    {
      auto itter = vertex_z_cuts_cd.find(pid);
      if(itter != vertex_z_cuts_cd.end())
	return p->par()->getVz() > itter->second.at(0) && p->par()->getVz() < itter->second.at(1);
      else
	return true;
    }

  return true;
}


bool clas12ana::checkPidCut(const region_part_ptr &p)
{
  //function returns true if inside PID cuts

  //electron pid is handled by ECal sampling fractions cuts NOT here
  if(p->par()->getPid() == 11)
    return true;

  if(p->getRegion() == FD) //forward detector cuts
    {
      auto itter = pid_cuts_fd.find(p->par()->getPid());
      if(itter != pid_cuts_fd.end())
	return (abs(p->par()->getChi2Pid() - itter->second.at(0)) < itter->second.at(1));
      else
	return false;
    }
  else if(p->getRegion() == CD) //central detector cuts
    {
      auto itter = pid_cuts_cd.find(p->par()->getPid());
      if(itter != pid_cuts_cd.end())
	return (abs(p->par()->getChi2Pid() - itter->second.at(0)) < itter->second.at(1));
      else
	return false;
    }

  return true;
}


bool clas12ana::checkProtonPidCut(const region_part_ptr &p)
{
  //returns true if inside cut
  //this is the pid done "by hand" where we cut on TOF vs momentum for protons
  //we only apply to CD where PID default from CLAS getByID is not as good

  //PID cuts only apply to CD
  if(p->getRegion() == clas12::FD) 
    return false;

  //only applies to + charged particles
  if(p->par()->getCharge() <= 0) 
    return false;

  //positive particles are reconstructed even without CTOF hit
  //Assigned a beta = -9999, here we throw away any tracks that could have no beta 
  //these tracks also have a path = 0 which would mae the tof_diff below 0. So we throw away all beta <= 0 tracks
  if(p->par()->getBeta() <= 1e-3)
    return false;

  //get the # of sigma away in PID cut from par file
  auto itter = pid_cuts_cd.find(2212);
  if(itter != pid_cuts_cd.end())
    {
      double mom = p->par()->getP();
      //proton_sigma=2;
      double exp_beta  = mom/sqrt(pow(mom,2) + pow(mass_proton,2)); //expected beta of particle assuming proton mass
      double tof_diff = (p->getPath()/c)*(1/p->par()->getBeta() - 1/exp_beta); //TOF difference measured - expected
      double up_lim  =  proton_pid_mean->Eval(mom) + proton_sigma*proton_pid_sigma->Eval(mom); 
      double low_lim =  proton_pid_mean->Eval(mom) - proton_sigma*proton_pid_sigma->Eval(mom);

      if( !(tof_diff < up_lim && tof_diff > low_lim))
	return false;
    }

  return true;
}



void clas12ana::readEcalSFPar(const char* filename)
{
  int num_par = 6; 
  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                   
      for(int i = 0; i < 2; ++i)
	getline(infile, tp);

      for(int i = 1; i < 7; ++i)
	{
	  getline(infile, tp);  //read data from file object and put it into string.       
          stringstream ss(tp);
          double parameter;
	  //get parameters for a given sector
	  for(int j = 0; j < num_par; j++)
	    {	 
	      ss >> parameter;
	      ecal_sf_fcn_par[i][j] = parameter;
	    }
	}


      InitSFEcalCuts();
    }
  else
    std::cout<<"ECal parameter files does not exist!!!"<<endl;


}


void clas12ana::readEcalPPar(const char* filename)
{
  int num_par = 6; 
  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                   
      for(int i = 0; i < 2; ++i)
	getline(infile, tp);

      for(int i = 1; i < 7; ++i)
	{
	  getline(infile, tp);  //read data from file object and put it into string.       
          stringstream ss(tp);
          double parameter;
	  //get parameters for a given sector
	  for(int j = 0; j < num_par; j++)
	    {	 
	      ss >> parameter;
	      ecal_p_fcn_par[i][j] = parameter;
	    }
	}


      InitSFPCuts();
    }
  else
    std::cout<<"ECal parameter files does not exist!!!"<<endl;

}


void clas12ana::clearInputParam()
{
  pid_cuts_cd.clear();
  pid_cuts_fd.clear();
  vertex_z_cuts_cd.clear();
  vertex_z_cuts_fd.clear();
}

void clas12ana::readInputParam(const char* filename)
{
  clearInputParam();

  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                                                                                              
      for(int i = 0; i < 3; ++i)
        getline(infile, tp);

      while(getline(infile, tp))  //read data from file object and put it into string.                                                                                                       
        {
          stringstream ss(tp);
          string parameter,parameter2;
          double value;
	  //get cut identifier 
	  ss >> parameter;
          if(parameter == "pid_cuts")
            {
	      //get cut values
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string pid_v;
	      string detector;
	      int count = 0; //parameter number
	      int pid = -99;
	      vector<double> par;
	      
	      while(getline(ss2, pid_v, ':'))
		{
		  if(count == 0)
		    pid = stoi(pid_v);
		  else if(count < 3)
		    par.push_back(atof(pid_v.c_str() ));
		  else if(count == 3)
		    detector = pid_v;
		  
		  count++;
		}
	      if(pid != -99) //if pid cut exists in file
		{
		  if(detector == "FD")
		    pid_cuts_fd.insert(pair<int, vector<double> >(pid, par));
		  else if(detector == "CD")
		    pid_cuts_cd.insert(pair<int, vector<double> >(pid, par));
		}
            }//end PID cuts section

          else if(parameter == "vertex_z_cut")
            {
	      //get cut values
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string pid_v;
	      string detector;
	      int count = 0; //parameter number
	      int pid = -99;
	      vector<double> par;
	      
	      while(getline(ss2, pid_v, ':'))
		{
		  if(count == 0)
		    pid = stoi(pid_v);
		  else if(count < 3)
		    par.push_back(atof(pid_v.c_str() ));
		  else if(count == 3)
		    detector = pid_v;
		  
		  count++;
		}
	      if(pid != -99) //if pid cut exists in file
		{
		  if(detector == "FD")
		    vertex_z_cuts_fd.insert(pair<int, vector<double> >(pid, par));
		  else if(detector == "CD")
		    vertex_z_cuts_cd.insert(pair<int, vector<double> >(pid, par));
		}
            }//end vertex z cuts
        }
    }
  else
    cout<<"Parameter file didn't read in "<<endl;

  this->printParams();

  return;
}


void clas12ana::readInputSRCParam(const char* filename)
{
  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                                                                                              
      for(int i = 0; i < 3; ++i)
        getline(infile, tp);

      while(getline(infile, tp))  //read data from file object and put it into string.     
        {
          stringstream ss(tp);
          string parameter,parameter2;
          double value;
	  //get cut identifier 
	  ss >> parameter;

          if(parameter == "q2")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;
		  
		  ++count;
		}
	      q2_cut = {min,max};
            }
          else if(parameter == "xb")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      xb_cut = {min,max};
            }
          else if(parameter == "pmiss")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      pmiss_cut = {min,max};
            }
          else if(parameter == "recoil_mom")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      recoil_mom_cut = {min,max};
            }
          else if(parameter == "miss_mass")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      mmiss_cut = {min,max};
            }
          else if(parameter == "p/q")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      pq_cut = {min,max};
            }
          else if(parameter == "theta_pq")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      theta_pq_cut = {min,max};
            }
          else if(parameter == "lead_mom")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string value;
	      double min,max;
	      int count = 0;
	      while(getline(ss2, value, ':'))
		{
		  stringstream ss3(value);
		  if (count == 0)
		    ss3 >> min;
		  else
		    ss3 >> max;

		  ++count;
		}
	      mom_lead_cut = {min,max};
            }

        }
    }
  else
    cout<<"SRC parameter file didn't read in "<<endl;

  return;
}


void clas12ana::printParams()
{
  cout<<endl;
  cout<<"Target Parameters:"<<endl;
  
  cout<<"Central Detector PID cuts:"<<endl;
  for (auto itr = pid_cuts_cd.begin(); itr != pid_cuts_cd.end(); ++itr) 
    {
      cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
      for(auto a : itr->second)
	cout  << '\t' << a ;
      cout<< '\n';
    }
  
  cout<<"Forward Detector PID cuts:"<<endl;
  for (auto itr = pid_cuts_fd.begin(); itr != pid_cuts_fd.end(); ++itr) 
    {
      cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
      for(auto a : itr->second)
	cout  << '\t' << a ;
      cout<< '\n';
    }

  cout<<"Central Detector Vz cuts:"<<endl;
  for (auto itr = vertex_z_cuts_cd.begin(); itr != vertex_z_cuts_cd.end(); ++itr) 
    {
      cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
      for(auto a : itr->second)
	cout  << '\t' << a ;
      cout<< '\n';
    }

  cout<<"Forward Detector Vz cuts:"<<endl;
  for (auto itr = vertex_z_cuts_fd.begin(); itr != vertex_z_cuts_fd.end(); ++itr)
    {
      cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
      for(auto a : itr->second)
	cout  << '\t' << a ;
      cout<< '\n';
    }

  cout <<"SRC lead and recoil cuts:"<<endl;
  cout << "Q2 {max,min}: "    << q2_cut[0]<<","<<q2_cut[1]<<endl;
  cout << "xB {max,min}: "    << xb_cut[0]<<","<<xb_cut[1]<<endl;
  cout << "Pmiss {max,min}: " << pmiss_cut[0]<<","<<pmiss_cut[1]<<endl;
  cout << "Recoil mom {max,min}: " << recoil_mom_cut[0]<<","<<recoil_mom_cut[1]<<endl;
  cout << "Misisng mass {max,min}: " << mmiss_cut[0]<<","<<mmiss_cut[1]<<endl;
  cout << "|p|/|q| {max,min}: " << pq_cut[0]<<","<<pq_cut[1]<<endl;
  cout << "Theta_pq {max,min}: " << theta_pq_cut[0]<<","<<theta_pq_cut[1]<<endl;
  cout << "Mom. lead {max,min}: " << mom_lead_cut[0]<<","<<mom_lead_cut[1]<<endl;

}

TVector3 clas12ana::getCOM(TLorentzVector lead, TLorentzVector recoil, TLorentzVector q)
{
  TVector3 com = (lead + recoil - q).Vect(); //center of momentum vector
  TVector3 pmiss = (lead - q).Vect();

  //vz along pmiss, vx,vy transverse
  TVector3 vz = pmiss.Unit(); //also called vt
  TVector3 vy = pmiss.Cross(q.Vect()).Unit();
  TVector3 vx = vz.Cross(vy).Unit();

  com*= 1000; //GeV to MeV

  return TVector3(com.Dot(vx),com.Dot(vy),com.Dot(vz));
}



void clas12ana::getLeadRecoilSRC(TLorentzVector beam, TLorentzVector target, TLorentzVector el)
{

  lead_proton.clear();
  recoil_proton.clear();

  TLorentzVector ptr(0,0,0,0);
  TLorentzVector q = beam - el;                  //photon  4-vector	
  double q2        = -q.M2();
  double xb        = q2/(2 * mass_proton * (beam.E() - el.E()) ); //x-borken       
  
  if( !(q2 > q2_cut[0] && xb > xb_cut[0]) )
    return; 
  
  int lead_idx   = -1;
  int lead_mult  = 0;


  for(int idx_ptr = 0; idx_ptr != protons.size(); ++idx_ptr)
    {

      ptr.SetXYZM(protons.at(idx_ptr)->par()->getPx(),protons.at(idx_ptr)->par()->getPy(),protons.at(idx_ptr)->par()->getPz(),mass_proton);

      TLorentzVector miss = beam + target - el - ptr; //missing 4-vector                   
      double pmiss    = miss.P();
      double mmiss    = miss.M();
      double theta_pq = ptr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); //angle between vectors p_miss and q                                                                               
      double p_q      = ptr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|                               
      if( ptr.P() > mom_lead_cut[0] &&  ptr.P() < mom_lead_cut[1]   && 
	  pmiss > pmiss_cut[0] && pmiss < pmiss_cut[1]              && 
	  mmiss > mmiss_cut[0] && mmiss < mmiss_cut[1]              && 
	  theta_pq > theta_pq_cut[0]  && theta_pq < theta_pq_cut[1] &&
	  p_q > pq_cut[0] && p_q < pq_cut[1])
	{
	  lead_idx = idx_ptr;
	  lead_mult++; //check for double lead
	}
    }

  if(lead_idx == -1 || lead_mult != 1)
    return;
  
  lead_proton.push_back(protons.at(lead_idx));


  int recoil_idx = -1;

  for(int idx_ptr = 0; idx_ptr != protons.size(); ++idx_ptr)
    {
      if(idx_ptr == lead_idx)
	continue;

      if(protons[idx_ptr]->par()->getP() > recoil_mom_cut[0] && protons[idx_ptr]->par()->getP() < recoil_mom_cut[1])
	recoil_proton.push_back(protons.at(idx_ptr));

    }

  //need to check if recoil momentum < lead momentum? 


  return;  

}
