To create a cutfile, structure the cut as follows:
cutname: minimum maximum

For the PID cuts you can use either:
cutname: 2212
or 
cutname: Proton

For the choice of scintiallator, use:
cutname: FTOF

To exlude a cut, remove it from your cut file


The follow is a list of the current cutnames.
They are seperated into different cut groups:

(e,e') Cuts:
e_nphe (Number of photo electrons in HTCC for electron)
e_calv (Calorimiter hit positon of electron)
e_calw (Calorimiter hit positon of electron)
e_SF (Calorimiter sampling fraction of electron)
e_mom (Electron momentum)

(e,e'N_{Lead}) Cuts:
l_pid (Lead nucleon PID)
l_scint (Lead nucleon scintillator)
l_theta (Lead nucleon angle)
l_thetalq (Angle between lead nucleon and q)
l_chipid (Lead nucleon PID chi squared)

(e,e'N_{Lead,SRC}) Cuts:
lsrc_Q2 (Q2)
lsrc_xB (xB)
lsrc_pmiss (pmiss)
lsrc_mmiss (mmiss)
lsrc_loq (Lead nucleon momentum over |q|)

(e,e'N_{Lead,SRC}N_{Recoil,SRC}) Cuts:
rsrc_pid (Recoil nucleon PID)
rsrc_mom (Recoil nucleon momentum)
