 #ifndef CLAS12DEBUG_HH
 #define CLAS12DEBUG_HH

 #include <iostream>
 #include <vector>
 #include <TF1.h>
 #include <math.h>
 #include "clas12reader.h"
 #include "region_particle.h"
 #include <map>
 #include "TH2D.h"
 #include "TH1D.h"
 #include "TVector3.h"
 #include "TLorentzVector.h"
 #include "TFile.h"
 #include "TCanvas.h"
 #include <sstream>

 #define CLAS12DEBUG_DIR _CLAS12DEBUG_DIR

 //#############
 //Debug clas for plotting before and after 
 //#############

 class clas12debug
 {

  public:
   clas12debug() 
     {
       //       InitDebugPlots();
     };

   ~clas12debug(){};

   void InitDebugPlots();
   void WriteDebugPlots();
   void plotDebug();
   void debugByPid(const clas12::region_part_ptr &p);
   void fillDCdebug(const clas12::region_part_ptr &p,std::vector<std::unique_ptr<TH2D>> &h);
   void fillBeforeEl(const clas12::region_part_ptr &p);
   void fillAfterEl(const clas12::region_part_ptr &p);
   double getSF(const clas12::region_part_ptr &p);
   void fillAfterPart(const clas12::region_part_ptr &p);
   void fillBeforePart(const clas12::region_part_ptr &p);

  private:

   double pi = 3.1415926;

   //debugging tools
   TString debug_fileName = "./debugOutputFile.root";
   bool debug_plots = true;

   std::vector<const TH2D *> hists_2D;

   std::vector<std::unique_ptr<TH2D>> sf_e_debug_b;
   std::vector<std::unique_ptr<TH2D>> sf_e_debug_a;
   std::vector<std::unique_ptr<TH2D>> sf_p_debug_b;
   std::vector<std::unique_ptr<TH2D>> sf_p_debug_a;

   /*
   TH2D *sf_e_debug_b[7] = {nullptr};
   TH2D *sf_e_debug_a[7] = {nullptr};
   TH2D *sf_p_debug_b[7] = {nullptr};
   TH2D *sf_p_debug_a[7] = {nullptr};
   */

   std::unique_ptr<TH2D> pid_cd_debug = std::make_unique<TH2D>("pid_cd_debug","PID Uncut CD",100,0,3,100,0,1.2);
   std::unique_ptr<TH2D> pid_fd_debug = std::make_unique<TH2D>("pid_fd_debug","PID Uncut FD",100,0,5,100,0,1.2);

   std::unique_ptr<TH2D>sf_v_ecalIN_debug = std::make_unique<TH2D>("sf_v_ecalIN_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D>sf_w_ecalIN_debug = std::make_unique<TH2D>("sf_w_ecalIN_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalOUT_debug = std::make_unique<TH2D>("sf_v_ecalOUT_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalOUT_debug = std::make_unique<TH2D>("sf_w_ecalOUT_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_pcal_debug = std::make_unique<TH2D>("sf_v_pcal_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_pcal_debug = std::make_unique<TH2D>("sf_w_pcal_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalIN_a_debug = std::make_unique<TH2D>("sf_v_ecalIN_a_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalIN_a_debug = std::make_unique<TH2D>("sf_w_ecalIN_a_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalOUT_a_debug = std::make_unique<TH2D>("sf_v_ecalOUT_a_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalOUT_a_debug = std::make_unique<TH2D>("sf_w_ecalOUT_a_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_pcal_a_debug = std::make_unique<TH2D>("sf_v_pcal_a_debug","",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_pcal_a_debug = std::make_unique<TH2D>("sf_w_pcal_a_debug","",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> pid_proton_fd_debug = std::make_unique<TH2D>("pid_proton_fd_debug","PID Cut Proton FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_proton_cd_debug = std::make_unique<TH2D>("pid_proton_cd_debug","PID Cut Proton CD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_piplus_fd_debug = std::make_unique<TH2D>("pid_piplus_fd_debug","PID Cut #pi + FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_piplus_cd_debug = std::make_unique<TH2D>("pid_piplus_cd_debug","PID Cut #pi + CD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_kplus_fd_debug  = std::make_unique<TH2D>("pid_kplus_fd_debug","PID Cut K+ FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_kplus_cd_debug  = std::make_unique<TH2D>("pid_kplus_cd_debug","PID Cut K+ CD",100,0,5,100,0,1.2);

   std::unique_ptr<TH2D> pid_piminus_fd_debug  = std::make_unique<TH2D>("pid_piminus_fd_debug","PID Cut #pi + FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_piminus_cd_debug  = std::make_unique<TH2D>("pid_piminus_cd_debug","PID Cut #pi + CD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_kminus_fd_debug   = std::make_unique<TH2D>("pid_kminus_fd_debug","PID Cut K+ FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_kminus_cd_debug   = std::make_unique<TH2D>("pid_kminus_cd_debug","PID Cut K+ CD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_neutrals_fd_debug = std::make_unique<TH2D>("pid_neutrals_fd_debug","PID Cut neutrals FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_neutrals_cd_debug = std::make_unique<TH2D>("pid_neutrals_cd_debug","PID Cut neutrals CD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_deuteron_fd_debug = std::make_unique<TH2D>("pid_deuteron_fd_debug","PID Cut deuteron FD",100,0,5,100,0,1.2);
   std::unique_ptr<TH2D> pid_deuteron_cd_debug = std::make_unique<TH2D>("pid_deuteron_cd_debug","PID Cut deutereon CD",100,0,5,100,0,1.2);

   std::unique_ptr<TH1D> el_vz_debug   = std::make_unique<TH1D>("el_vz_debug","El vertex ",100,-20,10);
   std::unique_ptr<TH1D> el_vz_p_debug = std::make_unique<TH1D>("el_vz_p_debug","El-proton vertex ",100,-10,10);

   std::unique_ptr<TH2D> cd_particles_b = std::make_unique<TH2D>("cd_particles_b","Pt;phi(deg);CD protons before edge cut",100,-180,180,100,0,1.);
   std::unique_ptr<TH2D> cd_particles_a = std::make_unique<TH2D>("cd_particles_a","Pt;phi(deg);CD protons after edge cut ",100,-180,180,100,0,1.);


   std::vector<std::unique_ptr<TH2D>> dc_hit_map_a; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_b; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_a_proton; //3 regions
   std::vector<std::unique_ptr<TH2D>>  dc_hit_map_b_proton; //3 regions

   std::vector<std::unique_ptr<TH2D>> dc_hit_map_a_pion; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_b_pion; //3 regions

   /*
   //   TH2D *dc_hit_map_a[4]; //3 regions
   //   TH2D *dc_hit_map_b[4]; //3 regions
   TH2D *dc_hit_map_a_proton[4]; //3 regions
   TH2D * dc_hit_map_b_proton[4]; //3 regions

   TH2D *dc_hit_map_a_pion[4]; //3 regions
   TH2D *dc_hit_map_b_pion[4]; //3 regions
   */
 };

#endif

