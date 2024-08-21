# Run Group M 
A repository for Run Group M. 

# Setting up your environment on the JLab Farm

```
module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module purge
module load sqlite/dev
module load clas12
module switch coatjava/10.0.2
```

alternatively 

```
source environment.csh
```

# Installing 

```
mkdir build
cd build
cmake /path/to/rgm/repo/
make
```