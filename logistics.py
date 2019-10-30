import numpy as np
import subprocess as sp
import argparse


def genParFile(pars, t, p, pref):
    """
    Generates a temporary parameter file to be used by ./ipfb
    :param pars: A set of parameters to transcribe into a new file
    :param t: The tile number
    :type t: int
    :param p: The polarisation number (either 0 or 1)
    :type p: int
    :param pref: The prefix for the directory to store the parameter files
    :type pref: str
    """
    file = open("{}/tmppar{}_{}.txt".format(pref, t, p), "w")
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
                 'outputdir', pars['outputdir'],
                 'filter_length', pars['filter_length'],
                 'filter_chans', pars['filter_chans'],
                 'amplification', pars['amplification'],
                 'nchannels', pars['nchannels'],
                 'firstchan', pars['firstchan'],
                 'ntiles', pars['ntiles'],
                 'tile', t,
                 'pol', p))

    file.close()


def clipCheck(pars, t, p):
    """
    Unused Definition
    :param pars:
    :param t:
    :param p:
    :rtype: bool
    """
    n = np.loadtxt("{}/norms_{}_{}.txt".format(pars['outputdir'], t, p))
    if n.sum() > 0:
        return True
    else:
        return False


def worker(rank, pars, t, p, args):
    """
    This is the main worker that runs the ./ipfb program with the correct arguments
    :param rank: the process rank
    :type rank: int
    :param pars: the parameters for input to genParFile
    :type pars: dict
    :param t: the tile number
    :type t: int
    :param p: the polarisation (0 or 1)
    :type p: int
    :param args: the script input arguments
    """
    repeat = True
    count = 0
    while repeat:
        # create parfile
        genParFile(pars, t, p, args.tmppardir)
        amp = 1
        print("amplification: {}, processor: {}".format(pars['amplification'], rank))
        # run ipfb on parfile
        if args.vcs:
            ipfbcall = ['./ipfb', '{}/tmppar{}_{}.txt'.format(args.tmppardir, t, p), '1']
        else:
            ipfbcall = ['./ipfb', '{}/tmppar{}_{}.txt'.format(args.tmppardir, t, p), '0']
        try:
            sp.check_output(ipfbcall, stderr=sp.STDOUT)
        except sp.CalledProcessError as e:
            print(e.output.decode('utf-8'))
            out = e.output.decode('utf-8').split('\n')[-2].split(' ')
            amp = float(out[1])
            # check for clipping
        if amp != 1:
            # print("Value clipped, repeating with lower amplification.")
            if int(float(pars['amplification'])) == 1:
                repeat = False
                print("amplification is already 1, cannot reduce anymore (there is an issue).")
            else:
                repeat = True
            pars['amplification'] = str(int(float(pars['amplification'])) * 100 // amp)
        else:
            repeat = False
        # remove parfile
        # os.remove('tmppars/tmppar{}_{}.txt'.format(t, p))


def readpars(parfile):
    """
    Reads a set of parameter files and inputs to a dictionary
    :param parfile: the name of the parameter file
    :type parfile: str
    :returns: the parameters and names as keys
    :rtype: dict
    """
    file = open(parfile, 'r')
    pardict = {}
    for line in file:
        words = line.split()
        if len(words) != 0:
            pardict[words[0]] = words[1]
    file.close()
    return pardict


def runMultiProcess(args, trange):
    """
    This runs the code using multiprocess, currently untested.
    :param args: the script input args
    :param trange: a range of tiles [min_t,max_t]
    :type trange: list, int

    """
    import multiprocessing as mp

    # read parameters from parfile
    pars = readpars(args.parfile)
    # loop over each stream
    jobs = []
    for t in range(trange[0], trange[1]):
        for p in range(2):
            proc = mp.Process(target=worker, args=(t*2+p, pars, t, p, args,))
            jobs.append(proc)
            proc.start()
    if not args.nowait:
        for proc in jobs:
            proc.join()


def run_MPI(args, trange):
    """
    Run's MPI using the worker defined above
    :param args: script input arguments
    :param trange: a range of tiles [min_t, max_t]
    :type trange: list, int
    """
    from mpi4py import MPI
    pars = readpars(args.parfile)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    nstreams = 2 * (trange[1] - trange[0])
    if nstreams != nprocs:
        raise ValueError("Please set the number of processors to be equal to the number of voltage streams"
                         " (i.e. 2 times the number of tiles).")
    else:
        t, p = np.mgrid[trange[0]:trange[1]:1, 0:2:1]
        t = t.reshape(nstreams)
        p = p.reshape(nstreams)
        print("tile: {}, pol: {}".format(t[rank], p[rank]))
        worker(rank, pars, t[rank], p[rank], args)
        # print("processor {} complete\n".format(rank))
    if not args.nowait:
        comm.barrier()


def module_parser(argarr):
    """
    Parses a set of parameters for if the script is run as a module
    :param argarr: an array of input arguments to parse
    :type argarr: list, str
    :return: a list of arguments as an object
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("parfile", help="The parameter file from which to read.")
    parser.add_argument("tmppardir", help="The location for storing temp parameter files.")
    parser.add_argument("-n", "--nowait", help="Do not wait for process to finish", action="store_true")
    parser.add_argument("-t", "--tiles", help="The range of tiles to process (as present within the input file,"
                                              " eg. '-t 0,12')", default=None)
    parser.add_argument("-m", "--mpi", help="Use mpi rather than multiprocess", action="store_true")
    parser.add_argument("-v", "--vcs", help="The input is vcs format", action="store_true")
    args = parser.parse_args(argarr)
    return args


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("parfile", help="The parameter file from which to read.")
    parser.add_argument("tmppardir", help="The location for storing temp parameter files.")
    parser.add_argument("-t", "--tiles", help="The range of tiles to process (as present within the input file,"
                                              " eg. '-t 0,12')", default=None)
    parser.add_argument("-m", "--mpi", help="Use mpi rather than multiprocess", action="store_true")
    parser.add_argument("-v", "--vcs", help="The input is vcs format", action="store_true")
    parser.add_argument("-n", "--nowait", help="Do not wait for process to finish", action="store_true")
    args = parser.parse_args()

    if args.tiles is None:
        raise ValueError("Please specify a range of tiles (eg. 0,12). Note that there are 24 cores on a magnus node, "
                         "therefore this range should have a multiple of 12 tiles (each tile has 2 polarisations)")
    else:
        trange = np.array(args.tiles.split(','), dtype=np.int)

    if args.mpi:
        print("Running with MPI")
        run_MPI(args, trange)
    else:
        print("Running with Multiprocess")
        runMultiProcess(args, trange)
