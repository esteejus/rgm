# Setup Environment on Farm
```
source ../environment.csh
```

# GCF Code Github
If you want to create GCF simulation files, you can pull the GCF code here:
```
cd /path/to/new/GCF/directory
git clone git@github.com:JacksonPybus/GCF_Generator_Suite.git
```

# Generating LUND files
Target types are listed in the targets.h file. These strings listed below have the correct nominal target positions and can be passed as a parameter in the LUND file scripts below.
```
4-foil - multi-foil targets such as C or Sn
1-foil - single foil targets at 2.5 cm downstream of center (part of Ar cell)
Ar     - Argon target -2.5 cm upstream of center
Ca     - Calcium targets at center of hall
liquid - liquid target cell extends -2.5 to 2.5 cm around center of hall

*note center of Hall is -3cm upstream in GEMC which is accounted for here.
```

Documentation on the LUND file format can be found here:
https://gemc.jlab.org/gemc/html/documentation/generator/lund.html?highlight=lund

Convert GENIE simulations to LUND
Takes root output of GENIE simulations and converts to LUND format. Typically GENIE files are given as a single large file which will be split into nFiles LUND files using the script.
```
root 'GENIE_to_LUND("inputFile","outputFile",nFiles,"target-type",A,Z)'
```

Convert GCF simulations to LUND
Takes root output of GCF simulations and converts to LUND format. Note that you need a root tree with branches for pLead, pRed, pe, pRel, q, and weight.
```
root 'GCF_to_LUND("inputFile","outputFile","target-type",A,Z)'
```

# Submitting Simulations on the Farm
Use ./submit/submit_GEMC.sh for submitting batch jobs on the farm. Keep the number of jobs reasonable and always use sqlite. This script will submit an array of jobs which is controled by tthe parameter ```#SBATCH --array=min-max```

LUND files have the prefix ```../lundfiles/lund_${FILE_PREFIX}_${SLURM_ARRAY_TASK_ID}.txt``` where ```SLURM_ARRAY_TASK_ID``` is the array value that you set above.
Also change
```
NEVENTS - number of events in each file
FILE_PREFIX - lund file prefix
```

Submit to Slurm on the JLab farm using
```
sbatch submit_GEMC.sh
```

Alternatively, I have also included a script to run GCF, GEMC, and reconstruction all in the same submit. Make sure that you already have the GCF code installed and running.
```
submit_GCF_GEMC_RECON.sh
```

Check the status of the simulations using
```
squeue -u yourUserName
```

*General tips for submitting on the farm. Submit an array of 1-1 just to test one file works without crashing. Log outputs can be found in /farm_out/userName/. Submit a reasonable amount of jobs i.e. no more than 1000 and try to keep your output file size < 1Tb 

# GCards and Yaml Files
The yaml files and gcards can be found in the submit folder. The torus field is set to -1 for most files. Make sure that this is consistent with the data files.
