# Inverse PFB

A new implementation of the inverse PFB.

## Compilation
1. go into src directory and type:
```bash
make
```
## Set up
1. Set up the data in a directory, use the "-v" flag to read in vcs input, omit this for HTR-VCS format.
2. Specify the parameters in a text file, an example can be found in the testing directory (as either fine or coarse parfile).
3. A file called \[datadir\].info needs to be created containing the list of files to include in the processing (make sure these are in channel order).
## Running
### Single processor
```bash
./ipfb [parameterfile]
```
### multiple processors
Using multiprocess, run logistics.py without paramters for more information on the inputs. ( this may not work in current version, use mpi if possible):
```bash
python logistics.py [parameterfile] -t [tile range]
```
Using MPI:
```bash
mpiexec -n [nprocs] python logistics.py [parameterfile] -t [tile range] -m
```
The overall_logistics.py script will facilitate a full vcs -> htr-vcs -> baseband inversion run, this is run via:
```bash
mpirun -n [nprocs] python overall_logistics.py [parameter_file]
```
This script requires 3 different parameter files: 
1. a master parameter file
2. a fine inversion paramater file
3. a coarse inversion parameter file

Examples of these can be found in the testing directory.

### Cleaning up
Go into src directory and type:
```bash
make clean
```
