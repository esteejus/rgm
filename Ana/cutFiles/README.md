# Cut files for data analysis

The .par files handle particle PID cuts (for FD and CD), vertex cuts, vertex correlation cuts:
ana.par      - generic example file
ana_he4.par  - 4He target, Deuterium, Hydrogen standard liquid cell targets
ana_cx4.par  - 4-foil 6GeV carbon target
ana_ca40.par - 40Ca target cell
ana_ca48.par - 48Ca target cell
ana_snx4.par - 120Sn 4-foil target 


The sampling fraction cuts are done as a function of momentum and deposited energy of the electon in the calorimeter. 

These files handle each sector parameterization for the momentum cuts:
```
paramsPI_40Ca_x2.dat
paramsPI_LD2_x2.dat
```

These files handle each sector parameterization for the deposited energy cuts:
```
paramsSF_40Ca_x2.dat
paramsSF_LD2_x2.dat
```
