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


