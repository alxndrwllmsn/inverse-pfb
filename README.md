# Inverse PFB

A new implementation of the inverse PFB.

To compile:
1. go into src directory and type:
```bash
make
```

To Run in series:
1. Set up the data in a directory, the code assumes the input of HTR-VCS for now.
2. This means that there should be a certain number of input files equal to the number of channels. Specify the parameters in a text file, an example can be found in the testing directory.
3. A file called [datadir].info needs to be created containing the list of files to include in the processing (make sure these are in channel order).
4. for series:
```bash
./ipfb [parameterfile]
```
for parallel (using multiprocess, run logistics.py without paramters for more information on the inputs):
```bash
python logistics.py [parameterfile] -t [tile range]
```
for parallel (using mpi4py):
```bash
mpiexec -n [nprocs] python logistics.py [parameterfile] -t [tile range] -m
```

To clean up:
1. go into src directory and type:
```bash
make clean
```
