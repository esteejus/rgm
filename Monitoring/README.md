
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
Run PrepareDatabases.C first to setup a local copy of RCDB and SQLite as not to overload the lab servers with queries. Add your hipo files of interest into the 

```
chain.Add(...)
```


Running the compiled code will output a message showing the inputs required and the required order .
```
/monitorPID <path/to/ouput.root> <path/to/ouput.pdf>  <path/to/input.hipo>
```

Optional flags:

```
-S         (optional flag to turn on MC simulation mode, DEFUALT is data and will try to read values from RCDB  
-E <value> (for manually entering beam energy)
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
