import argparse
import subprocess as sp
import os
import numpy as np
from logistics import readpars
from rearrangedata import rearrange_as_module
from time import time


def fine_inversion(datadir, fprefix, pars, fchanC, nchanC, nchanF, srun, vcs):
    for k in range(fchanC, fchanC + nchanC):
        print("Coarse Channel: {}".format(k))
        write_info_f(datadir, fprefix, k, nchanF)
        write_parfile(pars, k)
        setup_out_dir(pars["outputdir"], k)
        setup_out_dir("tmppars", k)
        if k == fchanC +nchanC - 1:
            nowait = False
        else:
            nowait = True
        run_logistics(srun, "tmppars/tmpparfileF_{}.txt".format(k), pars, vcs, "tmppars/{}".format(k), nowait)


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


def write_parfile(pars, chanC):
    print("writing temporary parameter file F")
    file = open("tmppars/tmpparfileF_{}.txt".format(chanC), "w")
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


def setup_out_dir(outdir, chanC):
    print("setting up output directories: {}/{}".format(outdir, chanC))
    owd = os.getcwd()
    os.chdir(outdir)
    newpath = "{}".format(chanC)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
    os.chdir(owd)


def run_logistics(srun, parfile, pars, vcs, parprefix, nowait):
    print("running ipfb logistics")
    nprocs = int(pars['ntiles'])*2
    nnodes = int(np.ceil(nprocs/24))
    logrun = ["mpirun", "-n", "{}".format(nprocs), "python3", "logistics.py", "{}".format(parfile), "{}".format(parprefix), "-t",
              "0,{}".format(pars['ntiles']), "-m"]
    if srun:
        logrun[0] = "srun"
        logrun.insert(1, "--export=all")
        logrun.insert(4, "-N")
        logrun.insert(5, "{}".format(nnodes))
        logrun.insert(6, "-c")
        logrun.insert(7, "1")
    if vcs:
        logrun.append("-v")
    if nowait:
        logrun.append("-n")

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
        rearrange_as_module(directory, 1280000, pars["ntiles"], "{}/{}_{}.sub".format(datadir, prefix, i))


def coarse_inversion(pars, fchanC, nchanC, cPrefix, srun, parfile, vcs):
    write_info_c(pars["datadir"], fchanC, nchanC, cPrefix, pars["outputdir"])
    run_logistics(srun, parfile, pars, vcs, "tmppars", False)


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
    if not os.path.exists("all"):
        os.mkdir("all")
    os.chdir(owd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("parfile", help="The master parfile for entire run.")
    parser.add_argument("-f", "--fine", help="Run the fine inversion (default is to run everything)",
                        action="store_true")
    parser.add_argument("-r", "--rearrange", help="Run data rearrangement (default is to run everything)",
                        action="store_true")
    parser.add_argument("-c", "--coarse", help="Run the coarse inversion (default is to run everything)",
                        action="store_true")

    args = parser.parse_args()
    mcontrol = False

    if args.fine or args.rearrange or args.coarse:
        mcontrol = True

    if not mcontrol:
        args.fine = True
        args.rearrange = True
        args.coarse = True

    print("Reading parameter file: {}".format(args.parfile))
    start = time()
    mpars = readpars(args.parfile)
    pars = readpars(mpars["fine_parfile"])
    print(time()-start)

    print("running fine inversion")
    start = time()
    if args.fine:
        fine_inversion(mpars["datadir"], mpars["fine_prefix"], pars, int(mpars["coarse_first_chan"]),
                       int(mpars["coarse_nchans"]), int(mpars["fine_nchans"]), int(mpars["srun"]), 1)
    print(time()-start)

    print("rearranging data for re - input")
    start = time()
    if args.rearrange:
        rearrange(int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), pars, mpars["coarse_prefix"],
                  mpars["datadir"])
    print(time() - start)

    print("running coarse inversion")
    start = time()
    if args.coarse:
        coarse_inversion(pars, int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), mpars["coarse_prefix"],
                         int(mpars["srun"]), mpars["coarse_parfile"], 0)
    print(time()-start)

