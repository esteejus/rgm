# Run Group M repository for reconstruction (a.k.a. cooking). 
This is a brief document of the reconstruction workflow specific to RGM. 

# Setting up workflow environment on the JLab Farm

```
source /group/clas12/packages/setup.csh #consider putting in ~/.cshrc file
module load workflow
```

# Tags for the workflow input
clas12-workflow.py takes care of setting up the workflow and submitting to the farm. Running clas12-workflow.py will show you the tags that can be input. Several key tags are listed below where some are optional. 

```
--model  ("dec","rec","ana" are decode,reconstruct, and analysis (produce trains) can be combined i.e. "decrecana" in any combination)
--runGroup (rgm)
--tag  (optional tag for labeling workflow)
--threads (32)
--reconSize (2)
--runs (accepts a list of runs i.e. "15359,15360,..." or a .txt file as an input)
--inputs (raw input of data usually from tape(/mss) or /cache i.e. /mss/clas12/rg-m/data/)
--outDir (output directory)
--clara  (location of clas12 reconstruction software (clara))
--reconYaml (reconstruction yaml file)
--trainYaml (optional train yaml file)
```

# Checking the workflow status 

Checking the workflow status is made easy by running 

```
swif-status.py 

```

for all the optional flags check 

```
swif-status.py --help

```

# Resubmitting failed jobs

For jobs that failed in the workflow resubmission is as simple as 

```
swif-status.py --retry

```

You should see an output of each failed job which was resubmitted separated by the reason why each job has failed. 


# Updating the ccdb data base for beam offset parameters


Example beamoffset.txt file:
```
#
  0            0            0            x-value          y-value          0            0     

```
All beam offset values must be set to 0 for the initial cook prior to doing the beam offset calibration which the calibration suite can be found here https://github.com/jeffersonlab/clas12-beamspot 

To update the ccdb database with a corresponding beamoffset.txt file one needs to submit from a user account and not clas12-5. Updating the database for an example #runNumber is done by submitting or by editing the file ccdb/modify_table.sh in this repo.

```
ccdb -c mysql://clas12writer:geom3try@clasdb.jlab.org/clas12 add /geometry/beam/position -r #runNumber -v rgm_fall2021 beamoffset.txt

```

The -r #runNumber flag needs to be specified for a range of run numbers.

```
-r 15000-15000 (for only run 15000)
-r 15000-16000 (for a run range 15000-16000)
-r 15000-      (for run 15000 - to infinity typically used in the initial calibration and is constantly overwritten with progressive calibrations)

```

