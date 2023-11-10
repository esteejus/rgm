# Creating signal and background samples (simulation)

N_getfeatures.cpp defines the features for simulated events to create signal and background events, depending on which value is selected for the input parameter called charged (0 for neutrons, 1 for protons).

The submit_sim.sh script shows an example of how to use N_getfeatures.cpp for slurm jobs.

# Creating signal and background samples (background)

D_getfeatures.cpp defines the features for signal events from Deuterium using the d(e,e'pn) channel.

submit_6gev_goodn.sh shows an example of how to use D_getfeatures.cpp for slurm jobs.

D_getfeatures_ppipn.sh defines the features for background events from Deuterium using the d(e,e'p pi- p) channel.

submit_6gev_ppipn.sh shows an example of how to use D_getfeatures_ppipn.sh for slurm jobs.

# Training charged particle veto ML algorithm for CND

Train the ML model using the following macro, using whichever ML models you want. (e.g. The following example uses MLP and BDT.) Before running, uncomment the relevant lines of code in the macro to choose which dataset to train on: simulation, D at 2 GeV, or D at 6 GeV. (Note: I'm not really using 2 GeV for training the model to use on nuclear targets. This is a holdover from early development.)

```
root -l TrainNeutronVeto_TMVA.C\(\"MLP,BDT\"\)
```

Depending on which dataset you want to train on, change the following settings:
- uncomment lines that read in data from either simulation, D 2 GeV, or D 6 GeV
- in the TMVA::DataLoader line, change the name to either dataset_sim (for simulation), dataset_2gev_pCD (for training on D 2 GeV), or dataset_6gev_pCD (for training on D 6 GeV) ... this determines which directory the trained model gets saved in
- in the PrepareTrainingAndTestTree line, choose how many signal and background events you want in your training sample (if you choose more than are available, you'll get an error and the macro will fail)

The trained models are xml files stored in their respective directories.

# Visualizing output plots

To combine the root plots from the signal and background files analyzed in parallel slurm jobs, execute the following for each type of event (signal or background, simulation or data).

```
$ROOTSYS/bin/hadd /path/to/target.root /path/to/source/files/file*.root
```

Run the following files to produce an output pdf with lots of histograms using the big root file.

```
root -q print_nsim_plots.c
root -q print_goodn_plots.c
root -q print_ppipn_plots.c
```

The first root script is for simulation (update the name to reflect choice of neutrons/protons). The second is for good neutrons from deuterium. The third is for protons misreconstructed as neutrons from deuterium. Make sure the file name in the script matches the file name from the combined root files.

# Applying the neutron veto

For an example application, use veto_test.cpp as a template.
