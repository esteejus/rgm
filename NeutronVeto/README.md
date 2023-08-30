# Training charged particle veto ML algorithm for CND

Train the ML model using this macro, using whichever ML models you want. (e.g. The following example uses MLP and BDT.) Before running, uncomment the relevant lines of code in the macro to choose which dataset to train on: simulation, D at 2 GeV, or D at 6 GeV.

```
root -l TrainNeutronVeto_TMVA.C\(\"MLP,BDT\"\)
```

The trained models are xml files stored in separate directories depending on the training dataset: dataset_2gev_pCD, dataset_6gev_pCD, and dataset_sim_eN_bknd.

# Applying the neutron veto

For an example application, use veto_test.cpp as a template.
