#include "clas12debug.h"

double tof_difference(const clas12::region_part_ptr &p, double mass)
{
  double c_light = 29.9792; //cm/ns
  double mom  = p->par()->getP();
  double path = p->getPath();
  double beta_expected = mom/sqrt(pow(mom,2) + pow(mass,2));
  double tof_expected = path/(beta_expected*c_light);
  double tof_meas = path/(p->par()->getBeta()*c_light);
  //  std::cout<<tof_meas<<" "<<tof_expected<<std::endl;
  return tof_meas - tof_expected;
}

void clas12debug::fillBeforeEl(const clas12::region_part_ptr &el)
{
  int sector = el->getSector();
  double el_mom = el->getP();
  double el_sf = getSF(el);
  double el_pcal_energy = el->cal(clas12::PCAL)->getEnergy();
  
  double ecin_v = el->cal(clas12::ECIN)->getLv();
  double ecin_w = el->cal(clas12::ECIN)->getLw();
  double ecout_v = el->cal(clas12::ECOUT)->getLv();
  double ecout_w = el->cal(clas12::ECOUT)->getLw();
  double pcal_v = el->cal(clas12::PCAL)->getLv();
  double pcal_w = el->cal(clas12::PCAL)->getLw();
  double chi2_ndf = el->trk(clas12::DC)->getChi2()/el->trk(clas12::DC)->getNDF();

  auto traj_index_1 = el->traj(clas12::DC,6) ->getIndex(); //layer 1                                                                                                                        
  auto traj_index_2 = el->traj(clas12::DC,18)->getIndex(); //layer 2                                                                                                                         
  auto traj_index_3 = el->traj(clas12::DC,36)->getIndex(); //layer 3                                                                                                                         
  
  auto traj_edge_1  = el->traj(clas12::DC,6) ->getFloat("edge",traj_index_1);
  auto traj_edge_2  = el->traj(clas12::DC,18)->getFloat("edge",traj_index_2);
  auto traj_edge_3  = el->traj(clas12::DC,36)->getFloat("edge",traj_index_3);
  
  pcal_energy_b_debug->Fill(el_pcal_energy,el->cal(clas12::ECIN)->getEnergy()+el->cal(clas12::ECOUT)->getEnergy());

  sf_v_ecalIN_debug->Fill(ecin_v,el_sf);
  sf_w_ecalIN_debug->Fill(ecin_w,el_sf);
  sf_v_ecalOUT_debug->Fill(ecout_v,el_sf);
  sf_w_ecalOUT_debug->Fill(ecout_w,el_sf);
  sf_v_pcal_debug->Fill(pcal_v,el_sf);
  sf_w_pcal_debug->Fill(pcal_w,el_sf);
  
  if( sector <= 6 && sector >= 1)
    {
      sf_e_debug_b[sector-1]->Fill(el_pcal_energy,el_sf); //0 indexed vector
      sf_p_debug_b[sector-1]->Fill(el_mom,el_sf); 

      dc_edge_el_chi2_r1[sector-1]->Fill(traj_edge_1,chi2_ndf);
      dc_edge_el_chi2_r2[sector-1]->Fill(traj_edge_2,chi2_ndf);
      dc_edge_el_chi2_r3[sector-1]->Fill(traj_edge_3,chi2_ndf);

      dc_edge_el_r1[sector-1]->Fill(traj_edge_1);
      dc_edge_el_r2[sector-1]->Fill(traj_edge_2);
      dc_edge_el_r3[sector-1]->Fill(traj_edge_3);
    }

  el_vz_b_debug->Fill( el->par()->getVz());  
  fillDCdebug(el,dc_hit_map_b);


}

void clas12debug::fillAfterEl(const clas12::region_part_ptr &el)
{
  int sector = el->getSector();
  double el_mom = el->getP();
  double el_sf = getSF(el);
  double el_pcal_energy = el->cal(clas12::PCAL)->getEnergy();
  
  double ecin_v = el->cal(clas12::ECIN)->getLv();
  double ecin_w = el->cal(clas12::ECIN)->getLw();
  double ecout_v = el->cal(clas12::ECOUT)->getLv();
  double ecout_w = el->cal(clas12::ECOUT)->getLw();
  double pcal_v = el->cal(clas12::PCAL)->getLv();
  double pcal_w = el->cal(clas12::PCAL)->getLw();
  
  pcal_energy_a_debug->Fill(el_pcal_energy,el->cal(clas12::ECIN)->getEnergy()+el->cal(clas12::ECOUT)->getEnergy());
  //DEBUG plots
  if(debug_plots && sector <= 6 && sector >= 1)
    {
      sf_e_debug_a[sector-1]->Fill(el->cal(clas12::PCAL)->getEnergy(),el_sf);
      sf_p_debug_a[sector-1]->Fill(el_mom,el_sf);
    }
  
  el_vz_a_debug->Fill( el->par()->getVz());
  
  sf_v_ecalIN_a_debug->Fill(ecin_v,el_sf);
  sf_w_ecalIN_a_debug->Fill(ecin_w,el_sf);
  sf_v_ecalOUT_a_debug->Fill(ecout_v,el_sf);
  sf_w_ecalOUT_a_debug->Fill(ecout_w,el_sf);
  sf_v_pcal_a_debug->Fill(pcal_v,el_sf);
  sf_w_pcal_a_debug->Fill(pcal_w,el_sf);
  
  
  fillDCdebug(el,dc_hit_map_a); //electron DC hit debug maps
}

void clas12debug::fillBeforePart(const clas12::region_part_ptr &p)
{
  double mass_proton    = 0.938272; //GeV/c2                                                                                                                                                 

  int pid = p->par()->getPid();

  if( pid == 2212)
    fillDCdebug(p,dc_hit_map_b_proton);
  if( pid == 211)
    fillDCdebug(p,dc_hit_map_b_pion);
  
  double par_mom  = p->par()->getP();
  double par_beta = p->par()->getBeta();
  
  bool is_cd = ( p->getRegion()==clas12::CD);
  bool is_fd = ( p->getRegion()==clas12::FD);
  
  int sector = p->getSector();
  double chi2_ndf = p->trk(clas12::DC)->getChi2()/p->trk(clas12::DC)->getNDF();

  auto traj_index_1 = p->traj(clas12::DC,6) ->getIndex(); //layer 1                                                                                                                          
  auto traj_index_2 = p->traj(clas12::DC,18)->getIndex(); //layer 2                                                                                                                          
  auto traj_index_3 = p->traj(clas12::DC,36)->getIndex(); //layer 3                                                                                                                          
  
  auto traj_edge_1  = p->traj(clas12::DC,6) ->getFloat("edge",traj_index_1);
  auto traj_edge_2  = p->traj(clas12::DC,18)->getFloat("edge",traj_index_2);
  auto traj_edge_3  = p->traj(clas12::DC,36)->getFloat("edge",traj_index_3);


  //DEBUG plots
  if(debug_plots && ( p->par()->getCharge() >= 1) && ( pid != 11) )
    {
      if(is_cd)
	{
	  pid_cd_debug->Fill(par_mom,par_beta);
	  cd_particles_b->Fill(p->getPhi()*180/pi,sqrt( pow(p->par()->getPx(),2) + pow(p->par()->getPy(),2)));

	  if(pid == 2212)
	    {
	      pid_proton_chi2_cd_debug->Fill(p->par()->getChi2Pid());
	      pid_proton_tof_cd_b_debug->Fill(par_mom,tof_difference(p,mass_proton));
	    }
	}
      else if(is_fd)
	{
	  if(pid == 2212)
	    {
	      pid_fd_debug->Fill(par_mom,par_beta);
	      pid_proton_chi2_fd_debug->Fill(p->par()->getChi2Pid());
	      pid_proton_tof_fd_b_debug->Fill(par_mom,tof_difference(p,mass_proton));
	    }
	}      

      if( sector <= 6 && sector >= 1)
	{
	  dc_edge_p_chi2_r1[sector-1]->Fill(traj_edge_1,chi2_ndf);
	  dc_edge_p_chi2_r2[sector-1]->Fill(traj_edge_2,chi2_ndf);
	  dc_edge_p_chi2_r3[sector-1]->Fill(traj_edge_3,chi2_ndf);
	  
	  dc_edge_p_r1[sector-1]->Fill(traj_edge_1);
	  dc_edge_p_r2[sector-1]->Fill(traj_edge_2);
	  dc_edge_p_r3[sector-1]->Fill(traj_edge_3);
	}
    }
  
}

void clas12debug::fillAfterPart(const clas12::region_part_ptr &p)
{
  double mass_proton    = 0.938272; //GeV/c2                                                                                                                                             
  int pid = p->par()->getPid();
  double par_mom  = p->par()->getP();
  double par_beta = p->par()->getBeta();
      
  if(p->par()->getCharge() >= 1 && p->par()->getPid() != 11 && p->par()->getPid() == 2212)
    {
      debugByPid(p);
      if(pid == 2212)
	fillDCdebug(p,dc_hit_map_a_proton);
      if(pid == 211)
	fillDCdebug(p,dc_hit_map_a_pion);
      
      bool is_cd = ( p->getRegion()==clas12::CD);
      bool is_fd = ( p->getRegion()==clas12::FD);
      
      if(is_cd)
	{
	  cd_particles_a->Fill(p->getPhi()*180/pi,sqrt( pow(p->par()->getPx(),2) + pow(p->par()->getPy(),2)));
	  pid_proton_tof_cd_a_debug->Fill(par_mom,tof_difference(p,mass_proton));
	}
      else if(is_fd)
	{
	  pid_proton_tof_fd_a_debug->Fill(par_mom,tof_difference(p,mass_proton));
	}
    }
}

double clas12debug::getSF(const clas12::region_part_ptr &p)
{
  if(p->par()->getPid() == 11)
    return (p->cal(clas12::PCAL)->getEnergy() +  p->cal(clas12::ECIN)->getEnergy() +  p->cal(clas12::ECOUT)->getEnergy()) / p->par()->getP();
  else
    return -9999.;
}

void clas12debug::fillDCdebug(const clas12::region_part_ptr &p, std::vector<std::unique_ptr<TH2D> > &h)
{
  h.at(0)->Fill(p->traj(clas12::DC,6)->getX(),p->traj(clas12::DC,6)->getY());
  h.at(1)->Fill(p->traj(clas12::DC,18)->getX(),p->traj(clas12::DC,18)->getY());
  h.at(2)->Fill(p->traj(clas12::DC,36)->getX(),p->traj(clas12::DC,36)->getY());
}


void clas12debug::plotDebug()
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

void clas12debug::debugByPid(const clas12::region_part_ptr &p)
{
  int pid = p->par()->getPid();
  double par_mom  = p->par()->getP();
  double par_beta = p->par()->getBeta();
  
  bool is_cd = (p->getRegion()==clas12::CD);
  bool is_fd = (p->getRegion()==clas12::FD);
  
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


 void clas12debug::InitDebugPlots()
 {

   hists_2D.push_back(pid_cd_debug.get());
   hists_2D.push_back(pid_fd_debug.get());

   for(int i = 1; i <= 6; i++)
     {

       sf_p_debug_b.push_back(std::make_unique<TH2D>(Form("sf_p_debug_b_sector_%d",i),Form("Sampling Fraction Before Cuts Sector_%d;Momentum (GeV/c);Sampling Fraction",i),100,0,6,100,0,.4));
       sf_p_debug_a.push_back(std::make_unique<TH2D>(Form("sf_p_debug_a_sector_%d",i),Form("Sampling Fraction  After Cuts Sector_%d;Momentum (GeV/c);Sampling Fraction",i),100,0,6,100,0,.4));

       sf_e_debug_b.push_back(std::make_unique<TH2D>(Form("sf_e_debug_b_sector_%d",i),Form("Sampling Fraction Before Cuts Sector_%d;Energy (GeV);Sampling Fraction",i),100,0,1.5,100,0,.4));
       sf_e_debug_a.push_back(std::make_unique<TH2D>(Form("sf_e_debug_a_sector_%d",i),Form("Sampling Fraction  After Cuts Sector_%d;Energy (GeV);Sampling Fraction",i),100,0,1.5,100,0,.4));

       dc_edge_el_chi2_r1.push_back(std::make_unique<TH1D>(Form("dc_edge_el_chi2_r1_%d",i),Form("Region 1 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_el_chi2_r2.push_back(std::make_unique<TH1D>(Form("dc_edge_el_chi2_r2_%d",i),Form("Region 2 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_el_chi2_r3.push_back(std::make_unique<TH1D>(Form("dc_edge_el_chi2_r3_%d",i),Form("Region 3 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));

       dc_edge_p_chi2_r1.push_back(std::make_unique<TH1D>(Form("dc_edge_p_chi2_r1_%d",i),Form("Region 1 Protons DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_p_chi2_r2.push_back(std::make_unique<TH1D>(Form("dc_edge_p_chi2_r2_%d",i),Form("Region 2 Protons DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_p_chi2_r3.push_back(std::make_unique<TH1D>(Form("dc_edge_p_chi2_r3_%d",i),Form("Region 3 Protons DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));

       dc_edge_el_r1.push_back(std::make_unique<TH1D>(Form("dc_edge_el_r1_%d",i),Form("Region 1 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_el_r2.push_back(std::make_unique<TH1D>(Form("dc_edge_el_r2_%d",i),Form("Region 2 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_el_r3.push_back(std::make_unique<TH1D>(Form("dc_edge_el_r3_%d",i),Form("Region 3 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));

       dc_edge_p_r1.push_back(std::make_unique<TH1D>(Form("dc_edge_p_r1_%d",i),Form("Region 1 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_p_r2.push_back(std::make_unique<TH1D>(Form("dc_edge_p_r2_%d",i),Form("Region 2 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
       dc_edge_p_r3.push_back(std::make_unique<TH1D>(Form("dc_edge_p_r3_%d",i),Form("Region 3 DC edge %d;Distance to Edge (cm);#chi^{2}/DOF",i),100,0,50));
     }


   //DC hit maps
   for(int i = 1; i <=3 ; i++)
     {
       dc_hit_map_b.push_back(std::make_unique<TH2D>(Form("dc_hitmap_before_%d",i), Form("Region %d Before Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));
       dc_hit_map_a.push_back(std::make_unique<TH2D>(Form("dc_hitmap_after_%d",i), Form("Region %d After Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));

       dc_hit_map_a_proton.push_back(std::make_unique<TH2D>(Form("dc_hitmap_after_proton_%d",i), Form("Region %d After Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));
       dc_hit_map_b_proton.push_back(std::make_unique<TH2D>(Form("dc_hitmap_before_proton_%d",i), Form("Region %d Before Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));
       dc_hit_map_a_pion.push_back(std::make_unique<TH2D>(Form("dc_hitmap_after_pion_%d",i), Form("Region %d After Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));
       dc_hit_map_b_pion.push_back(std::make_unique<TH2D>(Form("dc_hitmap_before_pion_%d",i), Form("Region %d Before Cuts;x-position (cm);y-position (cm)",i),600,-300,300,600,-300,300));

      //       dc_hit_map_b[i] = new TH2D(Form("dc_hitmap_before_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);
       //       dc_hit_map_a[i] = new TH2D(Form("dc_hitmap_after_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       /*
       dc_hit_map_a_proton[i] = new TH2D(Form("dc_hitmap_after_proton_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_b_proton[i] = new TH2D(Form("dc_hitmap_before_proton_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_a_pion[i] = new TH2D(Form("dc_hitmap_after_pion_%d",i), Form("Region %d After Cuts",i),600,-300,300,600,-300,300);
       dc_hit_map_b_pion[i] = new TH2D(Form("dc_hitmap_before_pion_%d",i), Form("Region %d Before Cuts",i),600,-300,300,600,-300,300);
       */
     }


 }


 void clas12debug::WriteDebugPlots()
 {
   TFile f_debugOut(debug_fileName,"RECREATE");

   std::for_each(sf_p_debug_b.begin(),sf_p_debug_b.end(), [](auto &el) {el->Write();} );
   std::for_each(sf_p_debug_a.begin(),sf_p_debug_a.end(), [](auto &el) {el->Write();} );
   std::for_each(sf_e_debug_a.begin(),sf_e_debug_a.end(), [](auto &el) {el->Write();} );
   std::for_each(sf_e_debug_b.begin(),sf_e_debug_b.end(), [](auto &el) {el->Write();} );

   std::for_each(dc_edge_el_chi2_r1.begin(),dc_edge_el_chi2_r1.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_el_chi2_r2.begin(),dc_edge_el_chi2_r2.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_el_chi2_r3.begin(),dc_edge_el_chi2_r3.end(),[](auto &p) {p->Write();} );

   std::for_each(dc_edge_el_r1.begin(),dc_edge_el_r1.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_el_r2.begin(),dc_edge_el_r2.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_el_r3.begin(),dc_edge_el_r3.end(),[](auto &p) {p->Write();} );

   std::for_each(dc_edge_p_chi2_r1.begin(),dc_edge_p_chi2_r1.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_p_chi2_r2.begin(),dc_edge_p_chi2_r2.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_p_chi2_r3.begin(),dc_edge_p_chi2_r3.end(),[](auto &p) {p->Write();} );

   std::for_each(dc_edge_p_r1.begin(),dc_edge_p_r1.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_p_r2.begin(),dc_edge_p_r2.end(),[](auto &p) {p->Write();} );
   std::for_each(dc_edge_p_r3.begin(),dc_edge_p_r3.end(),[](auto &p) {p->Write();} );


   for(int i = 0; i <= 2; i++)
       dc_hit_map_b[i]->Write();

   for(int i = 0; i <= 2; i++)
       dc_hit_map_a[i]->Write();

   for(int i = 0; i <= 2; i++)
       dc_hit_map_b_proton[i]->Write();

   for(int i = 0; i <= 2; i++)
       dc_hit_map_a_proton[i]->Write();

   for(int i = 0; i <= 2; i++)
       dc_hit_map_b_pion[i]->Write();

   for(int i = 0; i <= 2; i++)
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

   el_vz_b_debug->Write();
   el_vz_a_debug->Write();
   el_vz_p_debug->Write();

   pcal_energy_b_debug->Write();
   pcal_energy_a_debug->Write();

   pid_proton_chi2_fd_debug->Write();      
   pid_proton_chi2_cd_debug->Write();      

   pid_proton_tof_fd_b_debug->Write();
   pid_proton_tof_cd_b_debug->Write();
   pid_proton_tof_fd_a_debug->Write();
   pid_proton_tof_cd_a_debug->Write();

   f_debugOut.Close();
 }

