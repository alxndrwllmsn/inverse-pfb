import argparse
import subprocess as sp
import os
import numpy as np
from logistics import readpars
from rearrangedata import rearrange_as_module


def fine_inversion(datadir, fprefix, parfile, fchanC, nchanC, nchanF, srun, vcs):
    for k in range(fchanC, fchanC + nchanC):
        print("Coarse Channel: {}".format(k))
        write_info_f(datadir, fprefix, k, nchanF)
        pars = write_parfile(parfile, k)
        setup_out_dir(pars["outputdir"], k)
        run_logistics(srun, "tmpparfileF.txt", pars, vcs)
    return pars


def write_info_f(directory, prefix, chan, nchanF):
    print("Writing .info file")
    owd = os.getcwd()
    os.chdir("{}/../".format(directory))
    fname = directory.split('/')[-1]
    file = open("{}.info".format(fname), "w")
    textout = "{}_ch{}.dat".format(prefix, chan)      # Check format of input filenames
    for i in range(nchanF):
        file.write("{}\n".format(textout))
    file.close()
    os.chdir(owd)


def write_parfile(parfile, chanC):
    print("writing temporary parameter file F")
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
    print("setting up output directories: {}/{}".format(outdir, chanC))
    owd = os.getcwd()
    os.chdir(outdir)
    os.mkdir("{}".format(chanC))
    os.chdir(owd)


def run_logistics(srun, parfile, pars, vcs):
    print("running ipfb logistics")
    nprocs = int(pars['ntiles'])*2
    nnodes = int(np.ceil(nprocs/24))
    logrun = ["mpirun", "-n", "{}".format(nprocs), "python3", "logistics.py", "{}".format(parfile), "-t",
              "0,{}".format(pars['ntiles']), "-m"]
    if srun:
        logrun[0] = "srun"
        logrun.insert(3, "-N")
        logrun.insert(4, "{}".format(nnodes))
        logrun.insert(5, "-c")
        logrun.insert(6, "1")
    if vcs:
        logrun.append("-v")

    try:
        output = sp.check_output(logrun, stderr=sp.STDOUT).decode()
        print(output)
    except sp.CalledProcessError as e:
        output = e.output.decode()
        print(output)
        raise sp.CalledProcessError


def rearrange(fchanC, nchanC, pars, prefix, datadir): # note: nsamples is set for 1 second files
    owd = os.getcwd()
    for i in range(fchanC, fchanC + nchanC):
        directory = "{}/{}".format(pars["outputdir"], i)
        rearrange_as_module(directory, 1280000, pars["ntiles"], "{}/{}/{}_{}.sub".format(owd, datadir, prefix, i))


def coarse_inversion(pars, fchanC, nchanC, cPrefix, srun, parfile, vcs):
    write_info_c(pars["datadir"], fchanC, nchanC, cPrefix, pars["outputdir"])
    run_logistics(srun, parfile, pars, vcs)


def write_info_c(directory, fchanC, nchanC, cPrefix, outdir):
    print("writing info file")
    owd = os.getcwd()
    os.chdir("{}/..".format(directory))
    file = open("{}.info".format(directory.split('/')[-1]), "w")
    for i in range(fchanC, fchanC + nchanC):
        file.write("{}_{}.sub\n".format(cPrefix, i))
    file.close()
    os.chdir(owd)
    os.chdir(outdir)
    os.mkdir("all")
    os.chdir(owd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("parfile", help="The master parfile for entire run.")

    args = parser.parse_args()

    print("Reading parameter file: {}".format(args.parfile))
    mpars = readpars(args.parfile)

    print("running fine inversion")
    pars = fine_inversion(mpars["datadir"], mpars["fine_prefix"], mpars["fine_parfile"], int(mpars["coarse_first_chan"]),
                          int(mpars["coarse_nchans"]), int(mpars["fine_nchans"]), int(mpars["srun"]), 1)

    print("rearranging data for re - input")
    rearrange(int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), pars, mpars["coarse_prefix"], mpars["datadir"])

    print("running coarse inversion")
    coarse_inversion(pars, int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), mpars["coarse_prefix"], int(mpars["srun"]),
                     mpars["coarse_parfile"], 0)

