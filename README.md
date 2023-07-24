# Run Group M 
A repository for Run Group M. 


If you want to also clone the latest GCF generator to do simulations, clone using:

```
git clone --recurse-submodules https://github.com/esteejus/rgm.git
```

# Setting up your environment on the JLab Farm

```
module load cmake
module load clas12/pro
module load sqlite/dev
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
cmake /path/to/rgm/repo/ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
make
```