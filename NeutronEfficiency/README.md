# Measuring Neutron Efficiency

Run the submit scripts for each channel [ h(e,e'p pi-), d(e,e'pFDn), and d(e,e'pCDn) ] to perform the efficiency analysis for each run.

Combine the output from the efficiency analysis into one root file, e.g. using hadd.

```
$ROOTSYS/bin/hadd allevents_2gev_pCDn.root /output_directory/hepin_6gev_0150*.root
$ROOTSYS/bin/hadd allevents_6gev_epin.root /output_directory/dCD_2gev_0150*.root
```

Run the following macros to combine the neutron candidates and detected neutrons from all runs to get a measure of the neutron efficiency as a function of both momentum and polar angle.

```
combine_h_epin_neff.c
combine_d_pcdn_neff.c
```

The efficiency will be saved in a root file and also printed to the screen.

```
neff_6gev_epin.root
neff_2gev_pCDn.root
```

The workflow above illustrates the neutron efficiency calculation for the h(e,e'p pi-) and d(e,e'pCDn) channels. The parallelization for the d(e,e'pFDn) channel is still in development.
