import argparse
import os
import numpy as np
from logistics import readpars, module_parser, run_MPI, runMultiProcess
from rearrangedata import rearrange_as_module
from time import time
import multiprocessing as mp


def fine_inversion(datadir, fprefix, pars, fchanC, nchanC, nchanF, vcs, args):
    """
    Runs the fine pfb inversion
    :param datadir: directory of input data
    :type datadir: str
    :param fprefix: prefix for input vcs files
    :type fprefix: str
    :param pars: a set of parameters read from a file
    :type pars: dict
    :param fchanC: the first coarse channel number
    :type fchanC: int
    :param nchanC: the number of coarse channels
    :type nchanC: int
    :param nchanF: the number of fine channels
    :type nchanF: int
    :param vcs: specifies vcs flag for ./ipfb
    :type vcs: bool
    :param args: the command line arguments
    """
    amp = args.multiprocessing
    if not amp:
        comm = MPI.COMM_WORLD
    for k in range(fchanC, fchanC + nchanC):
        print("Coarse Channel: {}".format(k))
        if not amp:
            if comm.Get_rank() == 0:
                write_info_f(datadir, fprefix, k, nchanF)
                write_parfile(pars, k)
                setup_out_dir(pars["outputdir"], k)
                setup_out_dir("tmppars", k)
        else:
            write_info_f(datadir, fprefix, k, nchanF)
            write_parfile(pars, k)
            setup_out_dir(pars["outputdir"], k)
            setup_out_dir("tmppars", k)
            run_logistics("tmppars/tmpparfileF_{}.txt".format(k), pars, vcs, "tmppars/{}".format(k), False, args)
        if not amp:
            comm.barrier()
            run_logistics("tmppars/tmpparfileF_{}.txt".format(k), pars, vcs, "tmppars/{}".format(k), True, args)
    if not amp:
        comm.barrier()
    else:
        for p in mp.active_children():
            p.join()


def write_info_f(directory, prefix, chan, nchanF):
    """
    Writes the .info file for fine inversion
    :param directory: directory of input data
    :type directory: str
    :param prefix: the prefix for vcs files
    :type prefix: str
    :param chan: the current coarse channel
    :type chan: int
    :param nchanF: the number of fine channels
    :type nchanF: int
    """
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
    """
    Writes a temporary parameter file for the fine inversion
    :param pars: the parameters to transcribe
    :type pars: dict
    :param chanC: the current coarse channel
    :type chanC: int
    """
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
    """
    Sets up the output directories for the fine inversion
    :param outdir: the output directory name
    :type outdir: str
    :param chanC: the current coarse channel
    :type chanC: int
    """
    print("setting up output directories: {}/{}".format(outdir, chanC))
    owd = os.getcwd()
    os.chdir(outdir)
    newpath = "{}".format(chanC)
    if not os.path.exists(newpath):
        os.mkdir(newpath)
    os.chdir(owd)


def run_logistics(parfile, pars, vcs, parprefix, nowait, args):
    """
    Sends the parameters to the inverse pfb code and runs it
    :param parfile: the parameter file name
    :type parfile: str
    :param pars: the parameters from the master parameter file
    :type pars: dict
    :param vcs: sends vcs flag to ./ipfb code
    :type vcs: bool
    :param parprefix: the prefix for the parameter files
    :type parprefix: str
    :param nowait: sends nowait flag to run_MPI or MP
    :type nowait: bool
    :param args: the command line arguments
    """
    print("running ipfb logistics")
    logrun = ["{}".format(parfile), "{}".format(parprefix), "-t", "0,{}".format(pars['ntiles'])]

    if vcs:
        logrun.append("-v")
    if nowait:
        logrun.append("-n")

    logargs = module_parser(logrun)
    if not args.multiprocessing:
        run_MPI(logargs, np.array([0, int(pars['ntiles'])], dtype=np.int))
    else:
        runMultiProcess(logargs, np.array([0, int(pars['ntiles'])], dtype=np.int))


def rearrange(fchanC, nchanC, pars, prefix, datadir, args): # note: nsamples is set for 1 second files
    """
    Run's the rearrangedata.py as a module
    :param fchanC: the first coarse channel
    :type fchanC: int
    :param nchanC: the number of coarse channels
    :type nchanC: int
    :param pars: the parameters read from the master parameter file
    :type pars: dict
    :param prefix: the prefix for the output file
    :type prefix: str
    :param datadir: the directory of the output files
    :type datadir: str
    :param args: the command line arguments
    """
    if not args.multiprocessing:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocs = comm.Get_size()
    for i in range(fchanC, fchanC + nchanC):
        if not args.multiprocessing:
            if rank == (i-fchanC)% nprocs:
                directory = "{}/{}".format(pars["outputdir"], i)
                rearrange_as_module(directory, 1280000, pars["ntiles"], "{}/{}_{}.sub".format(datadir, prefix, i))
        else:
            directory = "{}/{}".format(pars["outputdir"], i)
            rearrange_as_module(directory, 1280000, pars["ntiles"], "{}/{}_{}.sub".format(datadir, prefix, i))
    if not args.multiprocessing:
        comm.barrier()


def coarse_inversion(pars, fchanC, nchanC, cPrefix, parfile, vcs, args):
    """
    Facilitates the coarse channel inversion
    :param pars: the paramters read from the master parameter file
    :type pars: dict
    :param fchanC: the first coarse channel
    :type fchanC: int
    :param nchanC: the number of coarse channels
    :type nchanC: int
    :param cPrefix: the prefix of the coarse channel input files
    :type cPrefix: str
    :param parfile: the parameter file name for coarse inversion
    :type parfile: str
    :param vcs: the vcs flag to send to run_logistics
    :type vcs: bool
    :param args: the command line arguments
    """
    if not args.multiprocessing:
        comm = MPI.COMM_WORLD
        if comm.Get_rank() == 0:
            write_info_c(pars["datadir"], fchanC, nchanC, cPrefix, pars["outputdir"])
        comm.barrier()
    else:
        write_info_c(pars["datadir"], fchanC, nchanC, cPrefix, pars["outputdir"])
    run_logistics(parfile, pars, vcs, "tmppars", False, args)


def write_info_c(directory, fchanC, nchanC, cPrefix, outdir):
    """
    Writes info file for coarse inversion
    :param directory: location of input files
    :type directory: str
    :param fchanC: the first coarse channel
    :type fchanC: int
    :param nchanC: the number of coarse channels
    :type nchanC: int
    :param cPrefix: the prefix of the coarse channel files
    :type cPrefix: str
    :param outdir: the output directory of the coarse inversion
    :type outdir: str
    """
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
    parser.add_argument("-p", "--multiprocessing", help="run using python multiprocessing rather than MPI",
                        action="store_true")

    args = parser.parse_args()
    mcontrol = False

    if args.fine or args.rearrange or args.coarse:
        mcontrol = True

    if not mcontrol:
        args.fine = True
        args.rearrange = True
        args.coarse = True
    if not args.multiprocessing:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nprocs = comm.Get_size()

    print("Reading parameter file: {}".format(args.parfile))
    start = time()
    mpars = readpars(args.parfile)
    pars = readpars(mpars["fine_parfile"])
    print(time()-start)

    print("running fine inversion")
    start = time()
    if args.fine:
        fine_inversion(mpars["datadir"], mpars["fine_prefix"], pars, int(mpars["coarse_first_chan"]),
                       int(mpars["coarse_nchans"]), int(mpars["fine_nchans"]), True, args)
    print(time()-start)

    print("rearranging data for re - input")
    start = time()
    if args.multiprocessing:
        for p in mp.active_children():
            p.join()
    if args.rearrange:
        rearrange(int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), pars, mpars["coarse_prefix"],
                  mpars["datadir"], args)
    print(time() - start)

    print("running coarse inversion")
    start = time()
    if args.coarse:
        coarse_inversion(pars, int(mpars["coarse_first_chan"]), int(mpars["coarse_nchans"]), mpars["coarse_prefix"],
                         mpars["coarse_parfile"], False, args)
    print(time()-start)

