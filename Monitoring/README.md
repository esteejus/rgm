# Install Monitoring Code

Setup environment 
```
module load cmake
module load sqlite/4.4.1
module load clas12/pro
module switch clas12root/1.7.3
```

```
mkdir build
cd build
cmake /path/to/Monitoring/ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
make
```

# Running Monitoring Code

Running the compiled code will output a message showing the inputs required and the required order.

```
./code <MC =1,Data = 0> <Ebeam(GeV)> <path/to/ouput.root> <path/to/ouput.pdf> <path/to/cutfile.txt> <path/to/input.hipo>
```

# Cut file for skimmer and monitoring

To create a cutfile, structure the cut as follows:

```
cutname: minimum maximum
```

For the PID cuts you can use either:
```
cutname: 2212
cutname: Proton
```

For the choice of scintiallator, use either:
```
cutname: FTOF1
cutname: FTOF2
cutname: FTOF
cutname: CTOF
cutname: TOF
```

Finally, all cuts are grouped together into "Electron Cuts", "Lead Cuts", "Lead SRC Cuts", and "Recoil SRC Cuts". If "Electron Cuts" is turned off, all cuts related to the electron will be turned off. Turn these cuts on and off using either:

```
Electron Cuts: ON
Electron Cuts: OFF
```

To exlude a cut, remove it from your cut file. The follow is a list of the current cutnames. The order of the cuts listed in the file are irrelevant. They are seperated into different cut groups:

```
(e,e') Cuts:
Electron Cuts (If this is off, none of the below will cut)
e_nphe (Number of photo electrons in HTCC for electron)
e_calv (Calorimiter hit positon of electron)
e_calw (Calorimiter hit positon of electron)
e_SF (Calorimiter sampling fraction of electron)
e_mom (Electron momentum)
e_vtze (Electron vertex position in z coordinate)

(e,e'N_{Lead}) Cuts:
Lead Cuts (If this is off, none of the below will cut)
l_pid (Lead nucleon PID)
l_scint (Lead nucleon scintillator)
l_theta (Lead nucleon angle)
l_thetalq (Angle between lead nucleon and q)
l_chipid (Lead nucleon PID chi squared)
l_vtzdiff (Electron minus lead nucleon vertex position in z coordinate)
l_phidiff (Electron minus lead nucleon vertex phi)

(e,e'N_{Lead,SRC}) Cuts:
Lead SRC Cuts (If this is off, none of the below will cut)
lsrc_Q2 (Q2)
lsrc_xB (xB)
lsrc_pmiss (pmiss)
lsrc_mmiss (mmiss)
lsrc_loq (Lead nucleon momentum over |q|)

(e,e'N_{Lead,SRC}N_{Recoil,SRC}) Cuts:
Recoil SRC Cuts (If this is off, none of the below will cut)
rsrc_pid (Recoil nucleon PID)
rsrc_mom (Recoil nucleon momentum)
rsrc_chipid (Recoil nucleon PID chi squared)
```