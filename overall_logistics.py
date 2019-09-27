import argparse
import subprocess as sp
import os
import numpy as np
from logistics import readpars
from rearrangedata import rearrange_as_module


def fine_inversion(datadir, fprefix, parfile, fchanC, nchanC, nchanF, srun, vcs):
    for k in range(fchanC, fchanC + nchanC):
        write_info_f(datadir, fprefix, k, nchanF)
        pars = write_parfile(parfile, k)
        setup_out_dir(pars, k)
        run_logistics(srun, parfile, pars, vcs)
    return pars


def write_info_f(directory, prefix, chan, nchanF):
    owd = os.getcwd()
    os.chdir("{}/../".format(directory))
    fname = directory.split[-1]
    file = open("{}.info".format(fname), "w")
    textout = "{}_{}.dat".format(prefix, chan)      # Check format of input filenames
    for i in range(nchanF):
        file.write("{}\n".format(textout))
    file.close()
    os.chdir(owd)


def write_parfile(parfile, chanC):
    pars = readpars(parfile)
    file = open("tmpparfileF.txt", "w")
    file.write("""{}    {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}
{}  {}""".format('datadir', pars['datadir'],
                 'filterfile', pars['filterfile'],
                 'outputdir', "{}/{}".format(pars['outputdir'],chanC),
                 'filter_length', pars['filter_length'],
                 'filter_chans', pars['filter_chans'],
                 'amplification', pars['amplification'],
                 'nchannels', pars['nchannels'],
                 'firstchan', pars['firstchan'],
                 'ntiles', pars['ntiles'],
                 'tile', 0,
                 'pol', 0))

    file.close()
    return pars


def setup_out_dir(outdir, chanC):
    owd = os.getcwd()
    os.chdir(outdir)
    os.mkdir("{}/{}".format(outdir, chanC))
    os.chdir(owd)


def run_logistics(srun, parfile, pars, vcs):
    nprocs = int(pars['ntiles']*2)
    nnodes = int(np.ceil(nprocs/24))
    logrun = ["mpirun", "-n", "{}".format(nprocs), "-N", "{}".format(nnodes), "-c", "1", "python",
              "logistics.py", parfile, "-t", "0,{}".format(pars['ntiles']), "-m"]
    if srun:
        logrun[0] = "srun"
    if vcs:
        logrun.append("-v")

    out = sp.check_output(logrun, stderr=sp.STDOUT)
    print(out)


def rearrange(fchanC, nchanC, pars, prefix, datadir): # note: nsamples is set for 1 second files
    for i in range(fchanC, fchanC + nchanC):
        directory = "{}/{}".format(pars["outputdir"], nchanC)
        rearrange_as_module(directory, 1280000, pars["ntiles"], "{}/{}_{}.sub".format(datadir, prefix, i))


def coarse_inversion(pars, fchanC, nchanC, cPrefix, srun, parfile, vcs):
    write_info_c(pars["datadir"], fchanC, nchanC, cPrefix)
    run_logistics(srun, parfile, pars, vcs)


def write_info_c(directory, fchanC, nchanC, cPrefix):
    owd = os.getcwd()
    os.chdir("{}/..".format(directory))
    file = open("{}.info".format(directory.split('/')[-1]), "w")
    for i in range(fchanC, fchanC + nchanC):
        file.write("{}_{}.sub\n".format(cPrefix, i))
    file.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("parfile", help="The master parfile for entire run.")

    args = parser.parse_args()

    mpars = readpars(args.parfile)

    pars = fine_inversion(mpars["datadir"], mpars["fine_prefix"], mpars["fine_parfile"], mpars["coarse_first_chan"],
                          mpars["coarse_nchans"], mpars["fine_nchans"], mpars["srun"], 1)

    rearrange(mpars["coarse_first_chan"], pars["coarse_nchans"], pars, mpars["coarse_prefix"], mpars["datadir"])

    coarse_inversion(pars, mpars["coarse_first_chan"], pars["coarse_nchans"], mpars["coarse_prefix"], mpars["srun"],
                     mpars["coarse_parfile"], 0)

