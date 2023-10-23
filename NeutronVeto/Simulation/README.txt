# Training Charged Particle Veto on Simulations

1. Simulate e'p and e'n events using a flat generator

The following bash scripts run generate_neutrons.sh and generate_protons.sh in series to simulate events with one electron and one neutron/proton with momentum from 0-1 GeV/c. The output format is lund.

```
bash ex_isotropic_CD_neutrons.sh
bash ex_isotropic_CD_protons.sh
```

2. Simulate detector response with clas12 geant4 (gemc)

The following script submits slurm jobs to run the e'n and e'p events through gemc. You'll have to change the OUTPUT_FILE variable and run it twice (once for protons and once for neutrons). The output format is hipo.

```
sbatch submit_GEMC_nobkg.sh
```

3. Add RGM background

The following python scripts run add_en_bkg.sh and add_ep_bkg.sh in parallel slurm jobs to add clas12 background (using a 6 GeV beam on Carbon) to each hipo file.

```
python submit_en_bkg.py
python submit_ep_bkg.py
```

4. Run CLAS12 reconstruction

The following script submits slurm jobs to perform clas12 reconstruction on each hipo file with clas12 background. You'll have to change the OUTPUT_FILE variable and run it twice.

```
sbatch submit_recon_bkg.sh
```

5. Analyze the events to get ML features

```
sbatch submit_analyze.sh
```
The output format will be root and txt files.


6. Train the Machine Learning algorithm

```
root -l TrainNeutronVeto_TMVA.C\(\"MLP,BDT\"\)
```
The trained model will be in xml format in the datasets directory. See veto_test.cpp for an example of how to implement the trained model.
