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
   void WriteDebugPlots(TString file);
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

   //regions 1,2,3 for each sector
   std::vector<std::unique_ptr<TH1D>> dc_edge_el_r1;
   std::vector<std::unique_ptr<TH1D>> dc_edge_el_r2;
   std::vector<std::unique_ptr<TH1D>> dc_edge_el_r3;

   std::vector<std::unique_ptr<TH1D>> dc_edge_p_r1; 
   std::vector<std::unique_ptr<TH1D>> dc_edge_p_r2;
   std::vector<std::unique_ptr<TH1D>> dc_edge_p_r3;

   std::vector<std::unique_ptr<TH1D>> dc_edge_el_chi2_r1;
   std::vector<std::unique_ptr<TH1D>> dc_edge_el_chi2_r2;
   std::vector<std::unique_ptr<TH1D>> dc_edge_el_chi2_r3;

   std::vector<std::unique_ptr<TH1D>> dc_edge_p_chi2_r1;
   std::vector<std::unique_ptr<TH1D>> dc_edge_p_chi2_r2;
   std::vector<std::unique_ptr<TH1D>> dc_edge_p_chi2_r3;


   std::unique_ptr<TH2D> pid_cd_debug = std::make_unique<TH2D>("pid_cd_debug","PID Uncut CD;Momentum (GeV/c);#Beta (v/c)",1000,0,3,1000,0,1.2);
   std::unique_ptr<TH2D> pid_fd_debug = std::make_unique<TH2D>("pid_fd_debug","PID Uncut FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);

   std::unique_ptr<TH2D> pcal_energy_b_debug = std::make_unique<TH2D>("pcal_energy_b_debug",";PCAL Energy (GeV);ECAL Inner + Outer (GeV)",1000,0,.6,100,0,.6);
   std::unique_ptr<TH2D> pcal_energy_a_debug = std::make_unique<TH2D>("pcal_energy_a_debug",";PCAL Energy (GeV);ECAL Inner + Outer (GeV)",1000,0,.6,100,0,.6);

   std::unique_ptr<TH2D>sf_v_ecalIN_debug = std::make_unique<TH2D>("sf_v_ecalIN_debug",";ECAL IN V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D>sf_w_ecalIN_debug = std::make_unique<TH2D>("sf_w_ecalIN_debug",";ECAL IN W (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalOUT_debug = std::make_unique<TH2D>("sf_v_ecalOUT_debug",";ECAL OUT V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalOUT_debug = std::make_unique<TH2D>("sf_w_ecalOUT_debug",";ECAL OUT W (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_pcal_debug = std::make_unique<TH2D>("sf_v_pcal_debug",";PCAL V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_pcal_debug = std::make_unique<TH2D>("sf_w_pcal_debug",";PCAL W (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalIN_a_debug = std::make_unique<TH2D>("sf_v_ecalIN_a_debug",";ECAL IN V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalIN_a_debug = std::make_unique<TH2D>("sf_w_ecalIN_a_debug",";ECAL IN W (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_ecalOUT_a_debug = std::make_unique<TH2D>("sf_v_ecalOUT_a_debug",";ECAL OUT V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_ecalOUT_a_debug = std::make_unique<TH2D>("sf_w_ecalOUT_a_debug",";ECAL OUT W (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> sf_v_pcal_a_debug = std::make_unique<TH2D>("sf_v_pcal_a_debug",";PCAL V (cm);Sampling Fraction",100,0,30,100,0,.4);
   std::unique_ptr<TH2D> sf_w_pcal_a_debug = std::make_unique<TH2D>("sf_w_pcal_a_debug",";PCAL V (cm);Sampling Fraction",100,0,30,100,0,.4);

   std::unique_ptr<TH2D> pid_proton_fd_debug = std::make_unique<TH2D>("pid_proton_fd_debug","PID Cut Proton FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_proton_cd_debug = std::make_unique<TH2D>("pid_proton_cd_debug","PID Cut Proton CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_piplus_fd_debug = std::make_unique<TH2D>("pid_piplus_fd_debug","PID Cut #pi + FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_piplus_cd_debug = std::make_unique<TH2D>("pid_piplus_cd_debug","PID Cut #pi + CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_kplus_fd_debug  = std::make_unique<TH2D>("pid_kplus_fd_debug","PID Cut K+ FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_kplus_cd_debug  = std::make_unique<TH2D>("pid_kplus_cd_debug","PID Cut K+ CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);

   std::unique_ptr<TH2D> pid_piminus_fd_debug  = std::make_unique<TH2D>("pid_piminus_fd_debug","PID Cut #pi + FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_piminus_cd_debug  = std::make_unique<TH2D>("pid_piminus_cd_debug","PID Cut #pi + CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_kminus_fd_debug   = std::make_unique<TH2D>("pid_kminus_fd_debug","PID Cut K+ FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_kminus_cd_debug   = std::make_unique<TH2D>("pid_kminus_cd_debug","PID Cut K+ CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_neutrals_fd_debug = std::make_unique<TH2D>("pid_neutrals_fd_debug","PID Cut neutrals FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_neutrals_cd_debug = std::make_unique<TH2D>("pid_neutrals_cd_debug","PID Cut neutrals CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_deuteron_fd_debug = std::make_unique<TH2D>("pid_deuteron_fd_debug","PID Cut deuteron FD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);
   std::unique_ptr<TH2D> pid_deuteron_cd_debug = std::make_unique<TH2D>("pid_deuteron_cd_debug","PID Cut deutereon CD;Momentum (GeV/c);#Beta (v/c)",1000,0,5,1000,0,1.2);

   std::unique_ptr<TH1D> pid_proton_chi2_fd_debug = std::make_unique<TH1D>("pid_proton_chi2_fd_debug","chi2PID Proton FD;chi2Pid;Counts",100,-10,10);
   std::unique_ptr<TH1D> pid_proton_chi2_cd_debug = std::make_unique<TH1D>("pid_proton_chi2_cd_debug","chi2PID Proton CD;chi2Pid;Counts",100,-10,10);

   std::unique_ptr<TH2D> pid_proton_tof_fd_b_debug = std::make_unique<TH2D>("pid_proton_tof_fd_b_debug","chi2PID Proton FD;Momentum (GeV/c);TOF_{Measured}-TOF_{Expected} (ns)",1000,0,3.5,100,-1,1);
   std::unique_ptr<TH2D> pid_proton_tof_cd_b_debug = std::make_unique<TH2D>("pid_proton_tof_cd_b_debug","chi2PID Proton CD;Momentum (GeV/c);TOF_{Measured}-TOF_{Expected} (ns)",1000,0,3.5,100,-1,1);

   std::unique_ptr<TH2D> pid_proton_tof_fd_a_debug = std::make_unique<TH2D>("pid_proton_tof_fd_a_debug","chi2PID Proton FD;Momentum (GeV/c);TOF_{Measured}-TOF_{Expected} (ns)",1000,0,3.5,100,-1,1);
   std::unique_ptr<TH2D> pid_proton_tof_cd_a_debug = std::make_unique<TH2D>("pid_proton_tof_cd_a_debug","chi2PID Proton CD;Momentum (GeV/c);TOF_{Measured}-TOF_{Expected} (ns)",1000,0,3.5,100,-1,1);

   std::unique_ptr<TH1D> el_vz_b_debug = std::make_unique<TH1D>("el_vz_b_debug","El vertex;z-vertex (cm);Counts",100,-20,10);
   std::unique_ptr<TH1D> el_vz_a_debug = std::make_unique<TH1D>("el_vz_a_debug","El vertex;z-vertex (cm);Counts",100,-20,10);

   std::unique_ptr<TH1D> p_vz_cd_debug = std::make_unique<TH1D>("p_vz_cd_debug","Proton vertex;z-vertex (cm);Counts",100,-20,10);
   std::unique_ptr<TH1D> p_vz_fd_debug = std::make_unique<TH1D>("p_vz_fd_debug","Proton vertex;z-vertex (cm);Counts",100,-20,10);

   std::unique_ptr<TH1D> el_vz_p_debug = std::make_unique<TH1D>("el_vz_p_debug","El-proton vertex;z-vertex (cm);Counts ",100,-10,10);

   std::unique_ptr<TH2D> cd_particles_b = std::make_unique<TH2D>("cd_edge_before","CD protons before edge cut;#Phi angle (deg);Transverse Momentum P_t (GeV/c)",100,-180,180,100,0,1.);
   std::unique_ptr<TH2D> cd_particles_a = std::make_unique<TH2D>("cd_edge_after","CD protons after edge cut;#Phi angle (deg);Transverse Momentum P_t (GeV/c)",100,-180,180,100,0,1.);

   std::vector<std::unique_ptr<TH2D>> dc_hit_map_a; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_b; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_a_proton; //3 regions
   std::vector<std::unique_ptr<TH2D>> dc_hit_map_b_proton; //3 regions

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

