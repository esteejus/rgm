#source ../../environment.csh

# naming convention: target_experiment_energy (e.g. d_rgm_2gev)


# HYDROGEN

# Run h(e,e'pi+ n) - RGM 6 GeV
./neff_h_epin 0 5.98636 neff_out/h_rgm_6gev.root neff_out/h_rgm_6gev.pdf ../cutfile_np_nosrc.txt ../../../data/h_6gev/h_015*_e.hipo

# Run h(e,e'pi+ n) - RGK 7.5 GeV
./neff_h_epin 7.546 neff_out/h_rgk_75gev.root neff_out/h_rgk_75gev.pdf ../cutfile_np_nosrc.txt ../../../data/rgk/*_e.hipo 


# DEUTERIUM

# Run d(e,e'pn) - RGB 4.5 GeV, p in FD
./neff_d_pfdn 0 4.244 neff_out/dFD_rgb_42gev.root neff_out/dFD_rgb_42gev.pdf ../cutfile_d_neff_pfd.txt ../../../data/rgb_ler/d_*_e.hipo

# Run d(e,e'pn) - RGM 2 GeV, p in FD
./neff_d_pfdn 0 2.07052 neff_out/dFD_rgm_2gev.root neff_out/dFD_rgm_2gev.pdf ../cutfile_d_neff_pfd.txt ../../../data/d_2gev/e_pFD_*.hipo

# Run d(e,e'pn) - RGM 6 GeV, p in FD, simulation
./neff_d_pfdn 0 6.0 neff_out/dFD_rgm_6gev_sim.root neff_out/dFD_rgm_6gev_sim.pdf ../cutfile_d_neff_pfd.txt ../../Simulation/reconhipo/recon_qe_d_6gev_*.hipo

# Run d(e,e'pn) - RGM 2 GeV, p in CD
./neff_d_pcdn 0 2.07052 neff_out/dCD_rgm_2gev.root neff_out/dCD_rgm_2gev.pdf ../cutfile_d_neff_pcd.txt ../../../data/d_2gev/e_pCD_*.hipo

# Run d(e,e'pn) - RGM 6 GeV, p in CD, simulation
./neff_d_pcddn 0 6.0 neff_out/dCD_rgm_6gev_sim.root neff_out/dCD_rgm_6gev_sim.pdf ../cutfile_d_neff_pcd.txt ../../Simulation/reconhipo/recon_qe_d_6gev_*.hipo
