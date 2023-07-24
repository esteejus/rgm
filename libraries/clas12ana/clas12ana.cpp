#include "clas12ana.h"

struct cutpar{
  string id;
  vector<double> par = {}; //pi- parameters                                                                                                                                                  
};

//helper function for DC fiducials                                                                                                                                                            
TVector3 rotate(TVector3 vec, int sector)
{
  double rot_ang = -(sector -1)*60 *TMath::DegToRad();

  vec.RotateZ(rot_ang);

  return vec;
}


void clas12ana::Clear()
 {
   //  particles.clear();
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

   event_mult = 0;
 }

void clas12ana::Run(const std::unique_ptr<clas12::clas12reader>& c12)
{
  Clear();
  
  auto particles = c12->getDetParticles(); //particles is now a std::vector of particles for this event
  auto electrons_det = c12->getByID(11);
  
      //DEBUG plots
      if(debug_plots)
	{
	  for (auto el = electrons_det.begin(); el != electrons_det.end();el++)
	    {
	      
	      int sector = (*el)->getSector();
	      double el_mom = (*el)->getP();
	      double el_sf = getSF(*el);
	      double el_pcal_energy = (*el)->cal(PCAL)->getEnergy();
	      
	      double ecin_v = (*el)->cal(ECIN)->getLv();
	      double ecin_w = (*el)->cal(ECIN)->getLw();
	      double ecout_v = (*el)->cal(ECOUT)->getLv();
	      double ecout_w = (*el)->cal(ECOUT)->getLw();
	      double pcal_v = (*el)->cal(PCAL)->getLv();
	      double pcal_w = (*el)->cal(PCAL)->getLw();

	      sf_v_ecalIN_debug->Fill(ecin_v,el_sf);
	      sf_w_ecalIN_debug->Fill(ecin_w,el_sf);
	      sf_v_ecalOUT_debug->Fill(ecout_v,el_sf);
	      sf_w_ecalOUT_debug->Fill(ecout_w,el_sf);
	      sf_v_pcal_debug->Fill(pcal_v,el_sf);
	      sf_w_pcal_debug->Fill(pcal_w,el_sf);
	      
	      if(sector <= 6 && sector >= 1)
		{
		  sf_e_debug_b[sector]->Fill(el_pcal_energy,el_sf); 
		  sf_p_debug_b[sector]->Fill(el_mom,el_sf); 
		}
	      
	      fillDCdebug(*el,dc_hit_map_b);
	    }
	}
  
       //Cut out electrons that do not pass cuts
      electrons_det.erase(std::remove_if(electrons_det.begin(),electrons_det.end(),[this](auto &el){
	    if( el->che(HTCC)->getNphe() <= 2)           //photo electron min cut
	      return true;
	    else if(!checkEcalSFCuts(el) && f_ecalSFCuts)//ECAL SF cuts
	      return true;
	    else if(!checkEcalPCuts(el) && f_ecalPCuts)  //ECAL SF cuts
	      return true;
	    else if(!EcalEdgeCuts(el) && f_ecalEdgeCuts) //ECAL edge cuts
	      return true;
	    else if(!checkVertex(el)  && f_vertexCuts)   //Vertex cut
	      return true;
	    else if(!DCEdgeCuts(el)   && f_DCEdgeCuts)   //DC edge cut
	      return true;
	    else
	      return false;
	  }), electrons_det.end());
	       
	      
      //DEBUG plots
      if(debug_plots)
	{
	  for (auto el = electrons_det.begin(); el != electrons_det.end();el++)
	    {
	      int sector = (*el)->getSector();
	      double el_mom = (*el)->getP();
	      double el_sf = getSF(*el);
	      double el_pcal_energy = (*el)->cal(PCAL)->getEnergy();
	      
	      double ecin_v = (*el)->cal(ECIN)->getLv();
	      double ecin_w = (*el)->cal(ECIN)->getLw();
	      double ecout_v = (*el)->cal(ECOUT)->getLv();
	      double ecout_w = (*el)->cal(ECOUT)->getLw();
	      double pcal_v = (*el)->cal(PCAL)->getLv();
	      double pcal_w = (*el)->cal(PCAL)->getLw();
	      
	      
	      //DEBUG plots
	      if(debug_plots && sector <= 6 && sector >= 1)
		{
		  sf_e_debug_a[sector]->Fill((*el)->cal(PCAL)->getEnergy(),el_sf);
		  sf_p_debug_a[sector]->Fill(el_mom,el_sf);
		}
	      
	      el_vz_debug->Fill( (*el)->par()->getVz());
	      
	      sf_v_ecalIN_a_debug->Fill(ecin_v,el_sf);
	      sf_w_ecalIN_a_debug->Fill(ecin_w,el_sf);
	      sf_v_ecalOUT_a_debug->Fill(ecout_v,el_sf);
	      sf_w_ecalOUT_a_debug->Fill(ecout_w,el_sf);
	      sf_v_pcal_a_debug->Fill(pcal_v,el_sf);
	      sf_w_pcal_a_debug->Fill(pcal_w,el_sf);
	      
	    }
	}
  
   if(electrons_det.size() == 1) //good trigger electron
     {
       
       //       setByPid(electrons_det[0]); //set good trigger electron

       if(debug_plots)
	 {
	   fillDCdebug(electrons_det[0],dc_hit_map_a); //electron DC hit debug maps

	   for (auto p = particles.begin(); p != particles.end();++p)
	     {
	       
	       if( (*p)->par()->getPid() == 2212)
		 fillDCdebug(*p,dc_hit_map_b_proton);
	       if( (*p)->par()->getPid() == 211)
		 fillDCdebug(*p,dc_hit_map_b_pion);
	   
	       double par_mom  = (*p)->par()->getP();
	       double par_beta = (*p)->par()->getBeta();
	       
	       bool is_cd = ( (*p)->getRegion()==CD);
	       bool is_fd = ( (*p)->getRegion()==FD);
	       
	       //DEBUG plots
	       if(debug_plots && ( (*p)->par()->getCharge() >= 1) && ( (*p)->par()->getPid() != 11) )
		 {
		   if(is_cd)
		     pid_cd_debug->Fill(par_mom,par_beta);
		   if(is_fd)
		     pid_fd_debug->Fill(par_mom,par_beta);

		   if((*p)->getRegion() == CD)
		     cd_particles_b->Fill((*p)->getPhi()*180/pi,sqrt( pow((*p)->par()->getPx(),2) + pow((*p)->par()->getPy(),2)));
		 }
	     }
	 }


       //Cut out particles that do not pass cuts
       particles.erase(std::remove_if(particles.begin(),particles.end(),[this,electrons_det](auto &p){
	     //neutrals and electrons don't follow cuts below, skip them 
	     if((p)->par()->getCharge() == 0 || (p)->par()->getPid() == 11 )
	       {
		 setByPid(p);
		 return true;
	       }
	     else
	       {
		 if((p)->par()->getPid() != 11)
		   event_mult++;
	       }
	     
	     
	     if(!checkPidCut(p) && f_pidCuts)	         //PID cuts
	       return true;
	     else if(!checkVertex(p) && f_vertexCuts)  //Vertex cut
	       return true;
	     else if(!CDEdgeCuts(p)  && f_CDEdgeCuts)  //CD edge cut
	       return true;
	     else if(!CDRegionCuts(p)  && f_CDRegionCuts)  //CD edge cut
	       return true;
	     else if(!DCEdgeCuts(p) && f_DCEdgeCuts)   //DC edge cut
	       return true;
	     else if(!checkVertexCorrelation(electrons_det[0],p) && f_corr_vertexCuts) //Vertex correlation cut between electron
	       return true;
	     else{
	       setByPid(p);
	       return false;
	     }
	   }), particles.end());
       

	   if(debug_plots)
	     {
	       for (auto p = particles.begin(); p != particles.end();++p)
		 {
		   if((*p)->par()->getCharge() != 0 && (*p)->par()->getPid() != 11)
		     el_vz_p_debug->Fill( (*p)->par()->getVz() - electrons_det[0]->par()->getVz() );

		   if((*p)->par()->getCharge() >= 1 && (*p)->par()->getPid() != 11 && (*p)->par()->getPid() == 2212)
		     {
		     debugByPid(*p);
		     if( (*p)->par()->getPid() == 2212)
		       fillDCdebug(*p,dc_hit_map_a_proton);
		     if( (*p)->par()->getPid() == 211)
		       fillDCdebug(*p,dc_hit_map_a_pion);

		     if( (*p)->getRegion() == CD)
		       cd_particles_a->Fill((*p)->getPhi()*180/pi,sqrt( pow((*p)->par()->getPx(),2) + pow((*p)->par()->getPy(),2)));

		     }
		 }
	     }
	   
	   
     }//good electron loop
       
   
}


void clas12ana::fillDCdebug(region_part_ptr p,TH2D **h)
{
  h[1]->Fill(p->traj(DC,6)->getX(),p->traj(DC,6)->getY());
  h[2]->Fill(p->traj(DC,18)->getX(),p->traj(DC,18)->getY());
  h[3]->Fill(p->traj(DC,36)->getX(),p->traj(DC,36)->getY());
}


void clas12ana::plotDebug()
{

  TCanvas *c1 = new TCanvas("c1","c1",1000,2000);
  c1->Divide(2,6);

  for(int i = 1; i <=6; i++)
    {   
      c1->cd(i);
      sf_p_debug_b[i]->Draw("colz");
    }
  for(int i = 7; i <=12; i++)
    {   
      c1->cd(i);
      sf_p_debug_a[i]->Draw("colz");
    }


}


void clas12ana::InitSFEcalCuts()
{
  cout<<"PARAMETERS for SF vs Ecal cuts"<<endl;
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  cout<<"sector "<<i <<" j "<<j<<" par "<<ecal_sf_fcn_par[i][j]<<endl;
	  ecal_sf_fcn[0][i]->SetParameter(j,ecal_sf_fcn_par[i][j]);
	  ecal_sf_fcn[1][i]->SetParameter(j,ecal_sf_fcn_par[i][j]);
	}

      ecal_sf_fcn[0][i]->SetParameter(6,sigma_cut);
      ecal_sf_fcn[1][i]->SetParameter(6,sigma_cut);
    }
}



void clas12ana::InitSFPCuts()
{
   cout<<"PARAMETERS for SF vs P cuts"<<endl;
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  cout<<"sector "<<i <<" j "<<j<<" par "<<ecal_p_fcn_par[i][j]<<endl;
	  ecal_p_fcn[0][i]->SetParameter(j,ecal_p_fcn_par[i][j]);
	  ecal_p_fcn[1][i]->SetParameter(j,ecal_p_fcn_par[i][j]);
	}

      ecal_p_fcn[0][i]->SetParameter(6,sigma_cut);
      ecal_p_fcn[1][i]->SetParameter(6,sigma_cut);
    }


}



void clas12ana::Init()
{
      for(int i = 0; i < 7; i++)
	{
	  ecal_sf_fcn[0][i] = new TF1(Form("ecal_sf_fcn_0_%d",i),"[0] + [1]/x + [2]/pow(x,2) - [6]*( [3] + [4]/x + [5]/pow(x,2))",0,1.5);
	  ecal_sf_fcn[1][i] = new TF1(Form("ecal_sf_fcn_1_%d",i),"[0] + [1]/x + [2]/pow(x,2) + [6]*( [3] + [4]/x + [5]/pow(x,2))",0,1.5);

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

  if(debug_plots)
    InitDebugPlots();


  this -> readInputParam( (std::string(CLAS12ANA_DIR) + "/Ana/cutFiles/ana.par").c_str() );
  this ->  readEcalSFPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsSF_40Ca_x2.dat").c_str() );
  this ->   readEcalPPar( (std::string(CLAS12ANA_DIR) +"/Ana/cutFiles/paramsPI_40Ca_x2.dat").c_str());
  this -> printParams();

}


bool clas12ana::DCEdgeCuts(region_part_ptr p)
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

      if(traj_edge_1 > dc_edge_cut && traj_edge_2 > dc_edge_cut && traj_edge_3 > dc_edge_cut)
	  return true;
      else
	return false;
    }
  else
    return true;
}



bool clas12ana::EcalEdgeCuts(region_part_ptr p)
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


bool clas12ana::checkEcalSFCuts(region_part_ptr p)
{
  //true if inside cut

  if(p->par()->getPid() == 11)
    {
      double sampling_frac = getSF(p);
      double energy =  p->cal(PCAL)->getEnergy();

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


bool clas12ana::checkEcalPCuts(region_part_ptr p)
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



double clas12ana::getSF(region_part_ptr p)
{
  if(p->par()->getPid() == 11)
    return (p->cal(clas12::PCAL)->getEnergy() +  p->cal(clas12::ECIN)->getEnergy() +  p->cal(clas12::ECOUT)->getEnergy()) / p->par()->getP();
  else
    return -9999.;
}


bool clas12ana::CDRegionCuts(region_part_ptr p)
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





bool clas12ana::CDEdgeCuts(region_part_ptr p)
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
      
      if( (std::abs(phi-fiducial_phi) < cd_edge_cut) || (std::abs(phi-fiducial_phi-120) < cd_edge_cut) || (std::abs(phi-fiducial_phi-240) < cd_edge_cut) )
	return false; //inside bad region
      else
	return true; //inside good region CD
    }
  else
    return true; //neutrals dont apply
}


bool clas12ana::checkVertex(region_part_ptr p)
{
  //true if inside cut
  return ( (p->par()->getVx() > vertex_x_cuts.at(0) && p->par()->getVx() < vertex_x_cuts.at(1)) 
	&& (p->par()->getVy() > vertex_y_cuts.at(0) && p->par()->getVy() < vertex_y_cuts.at(1))
        && (p->par()->getVz() > vertex_z_cuts.at(0) && p->par()->getVz() < vertex_z_cuts.at(1)) );

}

bool clas12ana::checkVertexCorrelation(region_part_ptr el, region_part_ptr p)
{
  //true if inside cut
  return ( (p->par()->getVz() - el->par()->getVz()) > vertex_corr_cuts.at(0) &&  (p->par()->getVz() - el->par()->getVz()) < vertex_corr_cuts.at(1) );
}

bool clas12ana::checkPidCut(region_part_ptr p)
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


void clas12ana::readEcalSFPar(const char* filename)
{
  int num_par = 6; 
  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                   
      for(int i = 0; i < 2; i++)
	getline(infile, tp);
      cout<<tp<<endl;


      for(int i = 1; i < 7; i++)
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
      for(int i = 0; i < 2; i++)
	getline(infile, tp);
      cout<<tp<<endl;


      for(int i = 1; i < 7; i++)
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


void clas12ana::readInputParam(const char* filename)
{
  ifstream infile;
  infile.open(filename);

  if (infile.is_open())
    {  
      string tp;

      //remove 3 lines of header                                                                                                                                              
      for(int i = 0; i < 3; i++)
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

          else if(parameter == "vertex_cut")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string pid_v;
	      int count = 0;
	      string pid = "";
	      vector<double> par;

	      while(getline(ss2, pid_v, ':'))
		{
		  if(count == 0)
		    pid = pid_v;
		  else
		    par.push_back(atof(pid_v.c_str() ));

		  count++;
		}

	      if(pid != "")
		vertex_cuts.insert(pair<string, vector<double> >(pid, par));
            }




	  /*


          else if(parameter == "cell_pos")
            {
	      ss >> parameter2;
              stringstream ss2(parameter2);
              string cell_v;
              while(getline(ss2, cell_v, ':'))
                cell.push_back(atof(cell_v.c_str()));
            }

          ss >> value;

	  if(parameter == "raster_x")
            raster_x = value;
          else if(parameter == "raster_y")
            raster_y = value;
          else if(parameter == "cell_limit")
            {
              cell.push_back(value); //first limit                                                                                                                                            
              ss >> value; //get second limit                                                                                                                                                 
              cell.push_back(value); //second limit                                                                                                                                           
            }
          else if(pannrameter == "target_A")
            target_A = value;
          else if(parameter == "target_Z")
            target_Z = value;
          else if(parameter == "beamspot_x")
            beamspot_x = value;
          else if(parameter == "beamspot_y")
            beamspot_y = value;

	  */
        }
    }
  else
    cout<<"Parameter file didn't read in "<<endl;

  return;
}


void clas12ana::printParams()
{
  cout<<endl;
  cout<<"Target Parameters:"<<endl;

  cout<<"Central Detector PID cuts:"<<endl;
  for (auto itr = pid_cuts_cd.begin(); itr != pid_cuts_cd.end(); ++itr) {
    cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
	  for(auto a : itr->second)
	    cout  << '\t' << a ;
	  cout<< '\n';
  }

  cout<<"Forward Detector PID cuts:"<<endl;
  for (auto itr = pid_cuts_fd.begin(); itr != pid_cuts_fd.end(); ++itr) {
    cout << '\t' << "Particle type: " << itr->first << '\t' <<"{mean,sigma}: ";
	  for(auto a : itr->second)
	    cout  << '\t' << a ;
	  cout<< '\n';
  }


  cout<<"Vertex cuts:"<<endl;
  for (auto itr = vertex_cuts.begin(); itr != vertex_cuts.end(); ++itr) {
    cout << '\t' << "Particle type: " << itr->first << '\t' <<"{min,max}: ";
	  for(auto a : itr->second)
	    cout  << '\t' << a ;
	  cout<< '\n';
  }

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

  //  if(abs(com.Dot(vy)) < 20)
  //  cout<<com.Mag()<<" "<<com.Dot(vx)<<" "<<com.Dot(vy)<<" "<<com.Dot(vz)<<" "<<lead.P()<<" "<<recoil.P()<<endl;

  return TVector3(com.Dot(vx),com.Dot(vy),com.Dot(vz));
}



void clas12ana::getLeadRecoilSRC(TLorentzVector beam, TLorentzVector target, TLorentzVector el)
{
  lead_proton.clear();
  recoil_proton.clear();

  TLorentzVector ptr;
  TLorentzVector q = beam - el;                  //photon  4-vector	
  double q2        = -q.M2();
  double xb        = q2/(2 * mass_proton * (beam.E() - el.E()) ); //x-borken       
  
  if( !(q2 > q2_cut && xb > xb_cut) )
    return; 
  
  int lead_idx   = -1;
  int lead_mult  = 0;

  for(int idx_ptr = 0; idx_ptr < protons.size(); idx_ptr++)
    {


      ptr.SetXYZM(protons.at(idx_ptr)->par()->getPx(),protons.at(idx_ptr)->par()->getPy(),protons.at(idx_ptr)->par()->getPz(),mass_proton);

      TLorentzVector miss = beam + target - el - ptr; //missing 4-vector                   
      double pmiss = miss.P();
      double mmiss   = miss.M2();
      double theta_pq = ptr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); //angle between vectors p_miss and q                                                                               
      double p_q       = ptr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|                               
      if( pmiss > pmiss_cut && mmiss > mmiss_cut[0] && mmiss < mmiss_cut[1] && theta_pq < theta_pq_cut && (p_q < pq_cut[1] && p_q > pq_cut[0]) )
	{
	  lead_idx = idx_ptr;
	  lead_mult++; //check for double lead
	}
    }

  ptr.SetXYZM(0,0,0,0);
  
  if(lead_idx == -1 || lead_mult != 1)
    return;
  
  lead_proton.push_back(protons.at(lead_idx));

  int recoil_idx = -1;

  for(int idx_ptr = 0; idx_ptr < protons.size(); idx_ptr++)
    {
      if(idx_ptr == lead_idx)
	continue;

      if(protons[idx_ptr]->par()->getP() > recoil_mom_cut)
	recoil_proton.push_back(protons.at(idx_ptr));

    }


  return;  
}



/*

void clas12ana::getLeadRecoilSRC(TLorentzVector beam, TLorentzVector target, TLorentzVector el)
{
  lead_proton.clear();
  recoil_proton.clear();

  TLorentzVector ptr;
  TLorentzVector q = beam - el;                  //photon  4-vector	
  double q2        = -q.M2();
  double xb        = q2/(2 * mass_proton * (beam.E() - el.E()) ); //x-borken       

  if( !(q2 > q2_cut && xb > xb_cut) )
    return; 
  
  int lead_idx   = -1;
  int lead_mult  = 0;
  int recoil_idx = -1;

  for(int idx_ptr = 0; idx_ptr < protons.size(); idx_ptr++)
    {
      ptr.SetXYZM(protons.at(idx_ptr)->par()->getPx(),protons.at(idx_ptr)->par()->getPy(),protons.at(idx_ptr)->par()->getPz(),mass_proton);
      double theta_pq = ptr.Vect().Angle(q.Vect()) * TMath::RadToDeg(); //angle between vectors p_miss and q                                                                               
      double p_q       = ptr.Vect().Mag()/q.Vect().Mag(); // |p|/|q|                               
      if(theta_pq < theta_pq_cut && (p_q < pq_cut[1] && p_q > pq_cut[0]) )
	{
	  lead_idx = idx_ptr;
	  lead_mult++; //check for double lead
	}
    }

  ptr.SetXYZM(0,0,0,0);
  
  if(lead_idx == -1)
    return;
  
  
  ptr.SetXYZM(protons.at(lead_idx)->par()->getPx(),protons.at(lead_idx)->par()->getPy(),protons.at(lead_idx)->par()->getPz(),mass_proton);
  
  TLorentzVector miss = beam + target - el - ptr; //missing 4-vector                   
  double pmiss = miss.P();
  double mmiss   = miss.M();
  
  if( pmiss > pmiss_cut && mmiss > mmiss_cut[0] && mmiss < mmiss_cut[1] )
      lead_proton.push_back(protons.at(lead_idx));



  for(int idx_ptr = 0; idx_ptr < protons.size(); idx_ptr++)
    {
      if(idx_ptr == lead_idx)
	continue;

      if(protons[idx_ptr]->par()->getP() > recoil_mom_cut)
	recoil_proton.push_back(protons.at(idx_ptr));

    }


  return;  
}
*/




void clas12ana::debugByPid(region_part_ptr p)
     {
       int pid = p->par()->getPid();
       double par_mom  = p->par()->getP();
       double par_beta = p->par()->getBeta();

       bool is_cd = ( p->getRegion()==CD);
       bool is_fd = ( p->getRegion()==FD);

       if(is_fd)
	 {
	   if( pid == 2212)
	     pid_proton_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == 45)
	     pid_deuteron_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == 211)
	     pid_piplus_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == -211)
	     pid_piminus_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == 321)
	     pid_kplus_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == -321)
	     pid_kminus_fd_debug->Fill(par_mom,par_beta);
	   else if( pid == 0 || pid == 2112)
	     pid_neutrals_fd_debug->Fill(par_mom,par_beta);
	 }
       else if(is_cd)
	 {
	   if( pid == 2212)
	     pid_proton_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == 45)
	     pid_deuteron_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == 211)
	     pid_piplus_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == -211)
	     pid_piminus_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == 321)
	     pid_kplus_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == -321)
	     pid_kminus_cd_debug->Fill(par_mom,par_beta);
	   else if( pid == 0 || pid == 2112)
	     pid_neutrals_cd_debug->Fill(par_mom,par_beta);
	 }

       
     }


 void clas12ana::InitDebugPlots()
 {

   for(int i = 1; i <= 6; i++)
     {
       sf_p_debug_b[i] = new TH2D(Form("sf_p_debug_b_sector_%d",i),Form("Sampling Fraction Before Cuts Sector_%d",i),100,0,6,100,0,.4);
       sf_p_debug_a[i] = new TH2D(Form("sf_p_debug_a_sector_%d",i),Form("Sampling Fraction  After Cuts Sector_%d",i),100,0,6,100,0,.4);

       sf_e_debug_b[i] = new TH2D(Form("sf_e_debug_b_sector_%d",i),Form("Sampling Fraction Before Cuts Sector_%d",i),100,0,1.5,100,0,.4);
       sf_e_debug_a[i] = new TH2D(Form("sf_e_debug_a_sector_%d",i),Form("Sampling Fraction  After Cuts Sector_%d",i),100,0,1.5,100,0,.4);
     }


   //DC hit maps
   for(int i = 1; i <=3 ; i++)
     {
       dc_hit_map_b[i] = new TH2D(Form("dc_hitmap_before_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_a[i] = new TH2D(Form("dc_hitmap_after_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_a_proton[i] = new TH2D(Form("dc_hitmap_after_proton_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_b_proton[i] = new TH2D(Form("dc_hitmap_before_proton_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_a_pion[i] = new TH2D(Form("dc_hitmap_after_pion_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_b_pion[i] = new TH2D(Form("dc_hitmap_before_pion_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);

     }


 }


 void clas12ana::WriteDebugPlots()
 {
   TFile *f_debugOut = new TFile(debug_fileName,"RECREATE");

   for(int i = 1; i <= 6; i++)
     {
       sf_p_debug_b[i]->Write();
       sf_p_debug_a[i]->Write();
       sf_e_debug_b[i]->Write();
       sf_e_debug_a[i]->Write();
     }

   for(int i = 1; i <= 6; i++)
     {
       ecal_sf_fcn[0][i]->Write();
       ecal_sf_fcn[1][i]->Write();
       ecal_p_fcn[0][i]->Write();
       ecal_p_fcn[1][i]->Write();
     }

   for(int i = 1; i <= 3; i++)
       dc_hit_map_b[i]->Write();

   for(int i = 1; i <= 3; i++)
       dc_hit_map_a[i]->Write();

   for(int i = 1; i <= 3; i++)
       dc_hit_map_b_proton[i]->Write();

   for(int i = 1; i <= 3; i++)
       dc_hit_map_a_proton[i]->Write();

   for(int i = 1; i <= 3; i++)
       dc_hit_map_b_pion[i]->Write();

   for(int i = 1; i <= 3; i++)
       dc_hit_map_a_pion[i]->Write();


   sf_v_ecalIN_debug->Write();
   sf_w_ecalIN_debug->Write();
   sf_v_ecalOUT_debug->Write();
   sf_w_ecalOUT_debug->Write();
   sf_v_pcal_debug->Write();
   sf_w_pcal_debug->Write();

   sf_v_ecalIN_a_debug->Write();
   sf_w_ecalIN_a_debug->Write();
   sf_v_ecalOUT_a_debug->Write();
   sf_w_ecalOUT_a_debug->Write();
   sf_v_pcal_a_debug->Write();
   sf_w_pcal_a_debug->Write();

   cd_particles_a->Write();
   cd_particles_b->Write();

   pid_proton_fd_debug->Write();
   pid_deuteron_fd_debug->Write();
   pid_piplus_fd_debug->Write();
   pid_piminus_fd_debug->Write();
   pid_kplus_fd_debug->Write();
   pid_kminus_fd_debug->Write();
   pid_neutrals_fd_debug->Write();

   pid_proton_cd_debug->Write();
   pid_deuteron_cd_debug->Write();
   pid_piplus_cd_debug->Write();
   pid_piminus_cd_debug->Write();
   pid_kplus_cd_debug->Write();
   pid_kminus_cd_debug->Write();
   pid_neutrals_cd_debug->Write();

   pid_cd_debug->Write();
   pid_fd_debug->Write();

   el_vz_debug->Write();
   el_vz_p_debug->Write();

   f_debugOut->Close();
 }

