
# Install Monitoring Code

```
module load cmake
module load clas12/dev
mkdir build
cd build
cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
make
```
# Running Monitoring Code

Running the compiled code will output a message showing the inputs required and the required order .
```
./monitorPID  
```

# Reading many files within a directory
You can read these files by passing the string /path/to/files/file_\*.hipo as the input file. The TChain inside will read all the events with that placeholder. The \ is needed in front of the * to pass the character along to main. 

```
/path/to/files/file_0.hipo
/path/to/files/file_1.hipo
/path/to/files/file_2.hipo
....
/path/to/files/file_N.hipo
```
