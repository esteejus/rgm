# Setup Environment on Farm
```
module load cmake
module load sqlite/4.4.1
module load clas12/pro
module switch clas12root/1.7.3
```

# Compile GCF generator
```
cd Simulation/GCF_Generator_Suite/
cmake ./
make
```