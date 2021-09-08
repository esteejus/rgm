# Setting up environment 
```
module load clas12/pro
```

# Fast Monte Carlo Roads
/work/clas12/users/devita/roads/coatjava-roads/bin/dict-maker

Options:
```
   -t      : torus scale factor (defualt = -1.0)
   -s      : solenoid scale factor (defualt = -1.0)
   -q      : charge (default = -1)
   -n      : number of events 
   -pmin   : minimum momentum (default = 0.3)
   -pmin   : minimum momentum (default = 11.0)
   -thmin  : theta minimum  (default = 5.0)
   -thmax  : theta maximum  (default = 40)
   -phimin : phi minimum  (default = -30)
   -phimax : phi maximum  (default = -30)
   -var    : remove duplicates in dictionary creation, 0=false, 1=true (default = 1)
   -seed   : RNG seed
   -vzmin  : z-vertex min (-5.0)
   -vzmax  : z-vertex max (5.0)

```



# Roads from GEMC
Run GEMC simulations and pass to dict-validate through the -create flag. 

/work/clas12/users/devita/roads/coatjava-roads/bin/dict-validate

```
dict-validate -create <simulation.hipo>
```


# Validate Roads
To validate roads use the -dict flag to point to the dictionary file, -i for the input file for validating the dictionary along with other options desired. 
/work/clas12/users/devita/roads/coatjava-roads/bin/dict-validate

Options :
```
   -charge : select particle charge for new dictionary, 0: no selection (default = 0)
   -create : select filename for new dictionary created from event file (default = )
     -dict : dictionary file name (default = )
    -dupli : remove duplicates in dictionary creation, 0=false, 1=true (default = 1)
        -i : event file for dictionary test (default = )
     -mode : select test mode, 0: DC only, 1: DC-FTOF-pcalU, 2: DC-FTOF-pcalUVW, 3: DC-FTOF-pcalUVW-HTCC (default = 0)
        -n : maximum number of events to process for validation (default = -1)
      -pid : select particle PID for new dictonary, 0: no selection, (default = 0)
   -sector : sector dependent roads, 0=false, 1=true) (default = 0)
    -strip : pcal strip smearing in road finding (default = 0)
-threshold : select roads momentum threshold in GeV (default = 1)
     -wire : dc wire smearing in road finding (default = 0)
```


# Merge Roads

To merge roads while removing duplicate roads (slow)
/work/clas12/users/devita/roads/coatjava-roads/bin/dict-merge

```
dict-merge <input/file1> <input/file2> ... <input/fileN> -o <output/file>
```

To combine road maps without removing duplicates (fast)
```
awk 'FNR==1{print ""}1' *.txt > finalfile.txt
```