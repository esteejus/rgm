# Measuring Neutron Efficiency

Each channel [ i.e. p(e,e' pi+ n) and d(e,e'pn) ] has a program for calculating the neutron efficiency using that channel [ neff_h_epin.cpp and neff_d_pn.cpp, respectively ].

The hydrogen program has an option to use the event selection that was used to obtain the RGK neutron efficiency measured that is reported in the NIM paper, or to use event selection for RGM.

The deuterium program has two options for dealing with the inelastic background: cutting on xB to get good separation between the signal and background and then cutting on the missing mass, or fitting the signal and background each to Gaussian distributions and then subtracting the background using the fits.

I have created skims for each channel at the following locations.

```
/lustre19/expphy/volatile/clas12/users/erins/neutron-efficiency/skims/h_epin
/lustre19/expphy/volatile/clas12/users/erins/neutron-efficiency/skims/d_pn
/lustre19/expphy/volatile/clas12/users/erins/neutron-efficiency/skims/d_pn_xb1
```

For the h_epin directory, I use RG-M hydrogen 6 GeV data, and the skim requires the predicted momentum of the neutron (q - pion momentum) to fall within 40-140 degrees and 0.2-1.2 GeV/c.

For the d_pn directory, I use RG-M deuterium 2 GeV data, and the skim requires the predicted momentum of the neutron (p - proton momentum) to fall within 40-140 degrees and 0.25-1.25 GeV/c. It also applies a cut requiring xB>0.3.

For the d_pn_xb1 directory, the event selection is the same as for d_pn, except that the xB cut is xB>0.1.
