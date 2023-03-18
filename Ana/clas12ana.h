 #ifndef CLAS12ANA_HH
 #define CLAS12ANA_HH

 #include <iostream>
 #include <vector>
 #include <TF1.h>
 #include <math.h>
 #include "clas12reader.h"
 #include "region_particle.h"
 #include <map>

 using namespace std;
 using namespace clas12;

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



 //#############
 //Analysis suite for CLAS12 analysis
 //#############

 class clas12ana : public clas12reader
 {

  public:
   clas12ana() 
     {
       Init();

     };

   void Init();
   void InitCuts();
   void readInputParam(const char* inFile);
   void readEcalPar(const char* filename);
   void printParams();
   void InitDebugPlots();
   void WriteDebugPlots();
   void Clear();
   void Run(const std::unique_ptr<clas12::clas12reader>& c12);
   void plotDebug();


   double getSF(region_part_ptr p);

   //   void pidCuts(std::vector<region_part_ptr> &particles);
   //   void pidCuts(std::vector<std::vector<region_part_ptr>> &particles);

   void setEcalSFCuts(bool flag = true)    {f_ecalSFCuts = flag;}; //option to have several cuts
   void setDCEdgeCuts(bool flag = true)    {f_DCEdgeCuts = flag;};
   void setEcalEdgeCuts(bool flag = true){f_ecalEdgeCuts = flag;};
   void setPidCuts(bool flag = true)     {f_pidCuts = flag;};
   void setVertexCuts(bool flag = true)  {f_vertexCuts = flag;};
   void setVertexCorrCuts(bool flag = true)  {f_corr_vertexCuts = flag;};

   TVector3 getCOM(TLorentzVector l, TLorentzVector r, TLorentzVector q);

   std::vector<region_part_ptr> getByPid(int pid)
     {
       if( pid == 11)
	 return electrons;
       else if( pid == 2212)
	 return protons;
       else if( pid == 2112)
	 return neutrons;
       else if( pid == 45)
	 return deuterons;
       else if( pid == 211)
	 return piplus;
       else if( pid == -211)
	 return piminus;
       else if( pid == 321)
	 return kplus;
       else if( pid == -321)
	 return kminus;
       else if( pid == 0)
	 return neutrals;
       else
	 return otherpart;
     }

   void setByPid(region_part_ptr p)
     {
       int pid = p->par()->getPid();
       if( pid == 11)
	 electrons.push_back(p);
       else if( pid == 2212)
	 protons.push_back(p);
       else if( pid == 2112)
	 neutrons.push_back(p);
       else if( pid == 45)
	  deuterons.push_back(p);
       else if( pid == 211)
	  piplus.push_back(p);
       else if( pid == -211)
	  piminus.push_back(p);
       else if( pid == 321)
	  kplus.push_back(p);
       else if( pid == -321)
	  kminus.push_back(p);
       else if( pid == 0)
	 neutrals.push_back(p);
       else 
	 otherpart.push_back(p);

     }


   double getEventMult(){return event_mult;};

   void debugByPid(region_part_ptr p);

   bool EcalEdgeCuts(region_part_ptr p);
   bool checkEcalCuts(region_part_ptr p);
   bool checkPidCut(region_part_ptr p);
   bool checkVertex(region_part_ptr p);
   bool checkVertexCorrelation(region_part_ptr el,region_part_ptr p);
   bool DCEdgeCuts(region_part_ptr p);

   void setVxcuts(double min, double max){vertex_x_cuts.at(0)=min; vertex_x_cuts.at(1)=max;};
   void setVycuts(double min, double max){vertex_y_cuts.at(0)=min; vertex_y_cuts.at(1)=max;};
   void setVzcuts(double min, double max){vertex_z_cuts.at(0)=min; vertex_z_cuts.at(1)=max;};
   void setVertexCorrCuts(double min, double max){vertex_corr_cuts.at(0) = min; vertex_corr_cuts.at(1) = max; };

   void fillDCdebug(region_part_ptr p,TH2D **h);

   void getLeadRecoilSRC(TLorentzVector beam, TLorentzVector target, TLorentzVector el);
   std::vector<region_part_ptr> getLeadSRC(){return lead_proton;};
   std::vector<region_part_ptr> getRecoilSRC(){return recoil_proton;};

   //  void getByChi2Pid(std::vector<region_part_ptr> &p,double mean, double sigma);
   std::vector<region_part_ptr> getByPid(std::vector<region_part_ptr> particles, int pid);
   //  std::vector<region_part_ptr> getByPidChi2(int pid, double chi2);


  private:

   std::vector<region_part_ptr> electrons;
   std::vector<region_part_ptr> protons;
   std::vector<region_part_ptr> deuterons;
   std::vector<region_part_ptr> neutrals;
   std::vector<region_part_ptr> neutrons;
   std::vector<region_part_ptr> piplus;
   std::vector<region_part_ptr> piminus;
   std::vector<region_part_ptr> kplus;
   std::vector<region_part_ptr> kminus;
   std::vector<region_part_ptr> otherpart;


   //SRC 
   std::vector<region_part_ptr> lead_proton;
   std::vector<region_part_ptr> recoil_proton;

   //prototype function for fitting ECAL electron cuts
   TF1 *ecal_fcn[2][7]; //0 upper 1 lower fiducial
   double ecal_fcn_par[7][6]; //sector, parameter
   int sigma_cut = 4;
   //vector to hold constants of function for each particle type
   //  vector<cutpar> dc_cuts; 
   //  vector<cutpar> ecal_cuts; 

   //pid quality cuts e.g. {2212, {0.,2}} for a proton with chisq2pid mean = 0. and sigma = 2
   //  vector<cutpar> pid_cuts;
   //vector<cutpar> vertex_cuts;

   bool f_ecalSFCuts         = false;
   bool f_ecalEdgeCuts       = false;
   bool f_DCEdgeCuts         = false;
   bool f_pidCuts            = false;
   bool f_vertexCuts         = false;
   bool f_corr_vertexCuts    = false;


   //  map<int,vector<double> > test_cuts;
   //  map<int,int> test_cuts;
   map<int,vector<double> > pid_cuts_cd; // map<pid, {min,max cut}> Central Detector (CD)
   map<int,vector<double> > pid_cuts_fd; // map<pid, {min,max cut}> Forward Detector (FD)

   vector<double> vertex_x_cuts = {-99,99};
   vector<double> vertex_y_cuts = {-99,99};
   vector<double> vertex_z_cuts = {-99,99};
   map<string,vector<double> > vertex_cuts; //map< x,y,z, {min,max}> 
   vector<double> vertex_corr_cuts = {-99,99}; //electron vertex <-> particle vertex correlation cuts

   double ecal_edge_cut = 14;
   double dc_edge_cut   = 5;

   //SRC Cuts
   double q2_cut = 1.5; //Q^2 cut
   double xb_cut = 1.2; //x-borken
   double pmiss_cut = .2; //missing momentum cut
   double recoil_mom_cut = .3; //missing momentum cut
   double mmiss_cut[2] = {.84,1.2}; //missing mass cut
   double pq_cut[2] = {0.6,0.96}; //|p|/|q| cut
   double theta_pq_cut = 25; //degrees angle between pLead & q

   //SRC variables
   double q2_e = 0;
   double xb_e = 0;
   double pmiss_e = 0; //missing momentum
   double mmiss_e = 0;
   double pq_e = 0;
   double theta_pq_e = 0;


   //  map<int,vector<int> > pid_cuts_stats;
   //  map<string,vector<int> > vertex_cuts_stats;

   //constants
   double mass_proton    = 0.938272; //GeV/c2
   double mass_neutron   = 0.939565;
   double mass_pion      = 0.13957;
   double mass_deuterium = 1.8756;

   double event_mult = 0; //charged particle multiplicity 

   //debugging tools
   TString debug_fileName = "./debugOutputFile.root";
   bool debug_plots = true;
   TH1D *ecal_sf[7]; //ECAL sampling fraction
   TH1D *dc[4][7];   //DC hit map



   TH2D *sf_debug_b[7] = {nullptr};
   TH2D *sf_debug_a[7] = {nullptr};

   TH2D *pid_cd_debug = new TH2D("pid_cd_debug","PID Uncut CD",100,0,3,100,0,1.2);
   TH2D *pid_fd_debug = new TH2D("pid_fd_debug","PID Uncut FD",100,0,5,100,0,1.2);

   TH2D *sf_v_ecalIN_debug = new TH2D("sf_v_ecalIN_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_ecalIN_debug = new TH2D("sf_w_ecalIN_debug","",100,0,30,100,0,.4);

   TH2D *sf_v_ecalOUT_debug = new TH2D("sf_v_ecalOUT_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_ecalOUT_debug = new TH2D("sf_w_ecalOUT_debug","",100,0,30,100,0,.4);

   TH2D *sf_v_pcal_debug = new TH2D("sf_v_pcal_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_pcal_debug = new TH2D("sf_w_pcal_debug","",100,0,30,100,0,.4);

   TH2D *sf_v_ecalIN_a_debug = new TH2D("sf_v_ecalIN_a_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_ecalIN_a_debug = new TH2D("sf_w_ecalIN_a_debug","",100,0,30,100,0,.4);

   TH2D *sf_v_ecalOUT_a_debug = new TH2D("sf_v_ecalOUT_a_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_ecalOUT_a_debug = new TH2D("sf_w_ecalOUT_a_debug","",100,0,30,100,0,.4);

   TH2D *sf_v_pcal_a_debug = new TH2D("sf_v_pcal_a_debug","",100,0,30,100,0,.4);
   TH2D *sf_w_pcal_a_debug = new TH2D("sf_w_pcal_a_debug","",100,0,30,100,0,.4);

   TH2D *pid_proton_fd_debug = new TH2D("pid_proton_fd_debug","PID Cut Proton FD",100,0,5,100,0,1.2);
   TH2D *pid_proton_cd_debug = new TH2D("pid_proton_cd_debug","PID Cut Proton CD",100,0,5,100,0,1.2);
   TH2D *pid_piplus_fd_debug = new TH2D("pid_piplus_fd_debug","PID Cut #pi + FD",100,0,5,100,0,1.2);
   TH2D *pid_piplus_cd_debug = new TH2D("pid_piplus_cd_debug","PID Cut #pi + CD",100,0,5,100,0,1.2);
   TH2D *pid_kplus_fd_debug = new TH2D("pid_kplus_fd_debug","PID Cut K+ FD",100,0,5,100,0,1.2);
   TH2D *pid_kplus_cd_debug = new TH2D("pid_kplus_cd_debug","PID Cut K+ CD",100,0,5,100,0,1.2);

   TH2D *pid_piminus_fd_debug = new TH2D("pid_piminus_fd_debug","PID Cut #pi + FD",100,0,5,100,0,1.2);
   TH2D *pid_piminus_cd_debug = new TH2D("pid_piminus_cd_debug","PID Cut #pi + CD",100,0,5,100,0,1.2);
   TH2D *pid_kminus_fd_debug = new TH2D("pid_kminus_fd_debug","PID Cut K+ FD",100,0,5,100,0,1.2);
   TH2D *pid_kminus_cd_debug = new TH2D("pid_kminus_cd_debug","PID Cut K+ CD",100,0,5,100,0,1.2);
   TH2D *pid_neutrals_fd_debug = new TH2D("pid_neutrals_fd_debug","PID Cut neutrals FD",100,0,5,100,0,1.2);
   TH2D *pid_neutrals_cd_debug = new TH2D("pid_neutrals_cd_debug","PID Cut neutrals CD",100,0,5,100,0,1.2);
   TH2D *pid_deuteron_fd_debug = new TH2D("pid_deuteron_fd_debug","PID Cut deuteron FD",100,0,5,100,0,1.2);
   TH2D *pid_deuteron_cd_debug = new TH2D("pid_deuteron_cd_debug","PID Cut deutereon CD",100,0,5,100,0,1.2);

   TH1D *el_vz_debug = new TH1D("el_vz_debug","El vertex ",100,-20,10);
   TH1D *el_vz_p_debug = new TH1D("el_vz_p_debug","El-proton vertex ",100,-10,10);


   TH2D *dc_hit_map_a[4]; //3 regions
   TH2D *dc_hit_map_b[4]; //3 regions

   TH2D *dc_hit_map_a_proton[4]; //3 regions
   TH2D *dc_hit_map_b_proton[4]; //3 regions

   TH2D *dc_hit_map_a_pion[4]; //3 regions
   TH2D *dc_hit_map_b_pion[4]; //3 regions

 };


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

   for (auto el = electrons_det.begin(); el != electrons_det.end();)
	 {

	   int sector = (*el)->getSector();
	   double el_mom = (*el)->getP();
	   double el_sf = getSF(*el);

	   double ecin_v = (*el)->cal(ECIN)->getLv();
	   double ecin_w = (*el)->cal(ECIN)->getLw();
	   double ecout_v = (*el)->cal(ECOUT)->getLv();
	   double ecout_w = (*el)->cal(ECOUT)->getLw();
	   double pcal_v = (*el)->cal(PCAL)->getLv();
	   double pcal_w = (*el)->cal(PCAL)->getLw();

	   //DEBUG plots
	   if(debug_plots)
	     {
	       sf_v_ecalIN_debug->Fill(ecin_v,el_sf);
	       sf_w_ecalIN_debug->Fill(ecin_w,el_sf);
	       sf_v_ecalOUT_debug->Fill(ecout_v,el_sf);
	       sf_w_ecalOUT_debug->Fill(ecout_w,el_sf);
	       sf_v_pcal_debug->Fill(pcal_v,el_sf);
	       sf_w_pcal_debug->Fill(pcal_w,el_sf);

	       if(sector <= 6 && sector >= 1)
		 sf_debug_b[sector]->Fill(el_mom,el_sf); 
	       fillDCdebug(*el,dc_hit_map_b);
	     }


	   //HTCC photo electron cuts
	   if( (*el)->che(HTCC)->getNphe() <= 2)
	       el = electrons_det.erase(el);


	   if(!checkEcalCuts(*el) && f_ecalSFCuts)       //ECAL SF cuts
	     {
	       el = electrons_det.erase(el);
	     }
	   else if(!EcalEdgeCuts(*el) && f_ecalEdgeCuts) //ECAL edge cuts
	     {
	       el = electrons_det.erase(el);
	     }
	   else if(!checkVertex(*el)  && f_vertexCuts)  //Vertex cut
	     {
	       el = electrons_det.erase(el);
	     }
	   else if(!DCEdgeCuts(*el)   && f_DCEdgeCuts)  //DC edge cut
	     {
	       el = electrons_det.erase(el);
	     }
	   else
	     {
	       //DEBUG plots
	       if(debug_plots && sector <= 6 && sector >= 1)
		 sf_debug_a[sector]->Fill(el_mom,el_sf);

	       el_vz_debug->Fill( (*el)->par()->getVz());

	       sf_v_ecalIN_a_debug->Fill(ecin_v,el_sf);
	       sf_w_ecalIN_a_debug->Fill(ecin_w,el_sf);
	       sf_v_ecalOUT_a_debug->Fill(ecout_v,el_sf);
	       sf_w_ecalOUT_a_debug->Fill(ecout_w,el_sf);
	       sf_v_pcal_a_debug->Fill(pcal_v,el_sf);
	       sf_w_pcal_a_debug->Fill(pcal_w,el_sf);

	       ++el; //itterate
	     }

	 }


   if(electrons_det.size() == 1) //good trigger electron
     {

       setByPid(electrons_det[0]); //set good trigger electron

       if(debug_plots)
	 fillDCdebug(electrons_det[0],dc_hit_map_a); //electron DC hit debug maps

       //DON'T FORGET TO ADD ++p ITTERATOR in this loop, it's not added in the for statement for a reason
       for (auto p = particles.begin(); p != particles.end();)
	 {

	   if(debug_plots)
	     {
	       if( (*p)->par()->getPid() == 2212)
		 fillDCdebug(*p,dc_hit_map_b_proton);
	       if( (*p)->par()->getPid() == 211)
		 fillDCdebug(*p,dc_hit_map_b_pion);
	     }


	   //neutrals don't follow cuts below, skip them 
	   if((*p)->par()->getCharge() == 0)
	     {
	       setByPid(*p);
	       ++p; //itterate 
	       continue; 
	     }
	   else
	     {
	       if((*p)->par()->getPid() != 11)
		 event_mult++;
	     }

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
	     }


	   //	   bool pid_cut    = checkPidCut(*p); 
	   //	   bool vertex_cut = checkVertex(*p);
	   //	   bool vertex_corr_cut = checkVertexCorrelation(electrons_det[0],*p); //correlation between good electron and particles vertex
	   //	   bool dc_edge_cut     = DCEdgeCuts(*p); 

	   if(!checkPidCut(*p) && f_pidCuts)	     //PID cuts
	     {
	       p = particles.erase(p);
	     }
	   else if(!checkVertex(*p) && f_vertexCuts) //Vertex cut
	     {
	       p = particles.erase(p);
	     }
	   else if(!DCEdgeCuts(*p) && f_DCEdgeCuts)  //DC edge cut
	     {
	       p = particles.erase(p);
	     }
	   else if(!checkVertexCorrelation(electrons_det[0],*p) && f_corr_vertexCuts) //Vertex correlation cut between electron
	     {
	       p = particles.erase(p);
	     }
	   else 	   //itterate 
	     {
	       setByPid(*p);
	       if(debug_plots)
		 {
		   if((*p)->par()->getCharge() != 0 && (*p)->par()->getPid() != 11)
		     el_vz_p_debug->Fill( (*p)->par()->getVz() - electrons_det[0]->par()->getVz() );

		   debugByPid(*p);
		   if( (*p)->par()->getPid() == 2212)
		     fillDCdebug(*p,dc_hit_map_a_proton);
		   if( (*p)->par()->getPid() == 211)
		     fillDCdebug(*p,dc_hit_map_a_pion);
		 }
	
	       ++p;
	     }

	 }//particle loop
     }//good electron loop


}


void clas12ana::fillDCdebug(region_part_ptr p,TH2D **h)
{
  //  if(p->par()->getPid() == 11)
  //    {
  h[1]->Fill(p->traj(DC,6)->getX(),p->traj(DC,6)->getY());
  h[2]->Fill(p->traj(DC,18)->getX(),p->traj(DC,18)->getY());
  h[3]->Fill(p->traj(DC,36)->getX(),p->traj(DC,36)->getY());
      //    }
}


void clas12ana::plotDebug()
{

  TCanvas *c1 = new TCanvas("c1","c1",1000,2000);
  c1->Divide(2,6);

  for(int i = 1; i <=6; i++)
    {   
      c1->cd(i);
      sf_debug_b[i]->Draw("colz");
    }
  for(int i = 7; i <=12; i++)
    {   
      c1->cd(i);
      sf_debug_a[i]->Draw("colz");
    }


}


void clas12ana::InitCuts()
{
 
   for(int i = 1; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  cout<<"sector "<<i <<" j "<<j<<" par "<<ecal_fcn_par[i][j]<<endl;
	  ecal_fcn[0][i]->SetParameter(j,ecal_fcn_par[i][j]);
	  ecal_fcn[1][i]->SetParameter(j,ecal_fcn_par[i][j]);
	}

      ecal_fcn[1][i]->SetParameter(6,sigma_cut);
      ecal_fcn[1][i]->SetParameter(6,sigma_cut);
    }

}



void clas12ana::Init()
{
      for(int i = 0; i < 7; i++)
	{
	  ecal_fcn[0][i] = new TF1(Form("ecal_fcn_0_%d",i),"[0] + [1]/sqrt(x) + [2]/x - [6]* ([3] + [4]/sqrt(x) + [5]/x)",0,10);
	  ecal_fcn[1][i] = new TF1(Form("ecal_fcn_1_%d",i),"[0] + [1]/sqrt(x) + [2]/x + [6]* ([3] + [4]/sqrt(x) + [5]/x)",0,10);
	}

  for(int i = 0; i < 7; i++)
    {
      for(int j = 0; j < 6; j++)
	{
	  if(j == 3)
	    ecal_fcn_par[i][j] = 9999; 

	  ecal_fcn_par[i][j] = 0;

	}
    }

  if(debug_plots)
    InitDebugPlots();

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


bool clas12ana::checkEcalCuts(region_part_ptr p)
{
  //true if inside cut

  if(p->par()->getPid() == 11)
    {
      double sampling_frac = getSF(p);
      int sector = p->getSector();
      
      //Turn on for functional form 
      //      double sf_max_cut = ecal_fcn[1][sector]->Eval(p->par()->getP() );
      //      double sf_min_cut = ecal_fcn[0][sector]->Eval(p->par()->getP() );
      //      cout<<"sf cut "<<sf_max_cut<<" "<<sf_min_cut<< " "<< sampling_frac <<" mom "<<p->par()->getP()<<endl;
      //      cout<<ecal_fcn[0][sector]->GetParameter(0)<<" "<<ecal_fcn[0][sector]->GetParameter(1)<<" "<<ecal_fcn[0][sector]->GetParameter(2)<<" sector "<<sector<<endl;

      double sf_max_cut = .28;
      double sf_min_cut = .2;
      
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

/*
void clas12ana::pidCuts(std::vector<region_part_ptr> &particles)
{


  for(auto &p : particles)
    {
      cout<<" Part ID " << p->par()->getPid() <<" " << endl;
      cout<<" Part ID " << p->par()->getChi2Pid() << " "<< endl;
    }  


  particles.erase(std::remove_if(particles.begin(), particles.end(), [this](const region_part_ptr& p) {
	auto itter = pid_cuts.find(p->par()->getPid());
	if(itter != pid_cuts.end())
	  return (abs(p->par()->getChi2Pid() - itter->second.at(0)) > itter->second.at(1));
	else
	  return false;
      }), particles.end());


  cout <<"After remove "<<endl;
  for(auto &p : particles)
    {
      cout<<" Part ID "<<p->par()->getPid()<<" "<<endl;
      cout<<" Part ID "<<p->par()->getChi2Pid()<<" "<<endl;
    }

}
*/

void clas12ana::readEcalPar(const char* filename)
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
	      ecal_fcn_par[i][j] = parameter;
	    }
	}


      InitCuts();
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
          else if(parameter == "target_A")
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

  if(abs(com.Dot(vy)) < 20)
    cout<<com.Mag()<<" "<<com.Dot(vx)<<" "<<com.Dot(vy)<<" "<<com.Dot(vz)<<" "<<lead.P()<<" "<<recoil.P()<<endl;

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
       sf_debug_b[i] = new TH2D(Form("sf_debug_b_sector_%d",i),Form("Sampling Fraction Before Cuts Sector_%d",i),100,0,6,100,0,.4);
       sf_debug_a[i] = new TH2D(Form("sf_debug_a_sector_%d",i),Form("Sampling Fraction  After Cuts Sector_%d",i),100,0,6,100,0,.4);

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
       sf_debug_b[i]->Write();
       sf_debug_a[i]->Write();

     }

   for(int i = 1; i <= 6; i++)
     {
       ecal_fcn[0][i]->Write();
       ecal_fcn[1][i]->Write();
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






#endif

