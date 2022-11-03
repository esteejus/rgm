#source ../../environment.csh

# naming convention: target_experiment_energy (e.g. d_rgm_2gev)


# HYDROGEN

# Run h(e,e'pi+ n) - RGM 6 GeV
mkdir -p h_epin_rgm_6gev
./neff_h_epin 0 5.98636 h_epin_rgm_6gev/h_rgm_6gev.root h_epin_rgm_6gev/h_rgm_6gev.pdf ../cutfile_np_nosrc.txt ../../../data/h_6gev/h_*_e.hipo

# Run h(e,e'pi+ n) - RGK 7.5 GeV
mkdir -p h_epin_rgk
./neff_h_epin 0 7.546 h_epin_rgk/h_rgk_75gev.root h_epin_rgk/h_rgk_75gev.pdf ../cutfile_np_nosrc.txt ../../../data/rgk/*_e.hipo 


# DEUTERIUM

# Run d(e,e'pn) - RGB 4.5 GeV, p in FD
mkdir -p d_epFDn_rgb
./neff_d_pfdn 0 4.244 d_epFDn_rgb/dFD_rgb_42gev.root d_epFDn_rgb/dFD_rgb_42gev.pdf ../cutfile_d_neff_pfd.txt ../../../data/rgb_ler/d_*_e.hipo

# Run d(e,e'pn) - RGM 2 GeV, p in FD
mkdir -p d_epFDn_rgm_2gev
./neff_d_pfdn 0 2.07052 d_epFDn_rgm_2gev/dFD_rgm_2gev.root d_epFDn_rgm_2gev/dFD_rgm_2gev.pdf ../cutfile_d_neff_pfd.txt ../../../data/d_2gev/e_pFD_*.hipo

# Run d(e,e'pn) - RGM 6 GeV, p in FD, simulation
mkdir -p d_epFDn_sim_6gev
./neff_d_pfdn 1 6.0 d_epFDn_sim_6gev/dFD_rgm_6gev_sim.root d_epFDn_sim_6gev/dFD_rgm_6gev_sim.pdf ../cutfile_d_neff_pfd.txt ../../Simulation/reconhipo/recon_qe_d_6gev_*.hipo

# Run d(e,e'pn) - RGM 2 GeV, p in CD
mkdir -p d_epCDn_rgm_2gev
./neff_d_pcdn 0 2.07052 d_epCDn_rgm_2gev/dCD_rgm_2gev.root d_epCDn_rgm_2gev/dCD_rgm_2gev.pdf ../cutfile_d_neff_pcd.txt ../../../data/d_2gev/e_pCD_*.hipo

# Run d(e,e'pn) - RGM 6 GeV, p in CD, simulation
mkdir -p d_epCDn_sim_6gev
./neff_d_pcddn 1 6.0 d_epCDn_sim_6gev/dCD_rgm_6gev_sim.root d_epCDn_sim_6gev/dCD_rgm_6gev_sim.pdf ../cutfile_d_neff_pcd.txt ../../Simulation/reconhipo/recon_qe_d_6gev_*.hipo
