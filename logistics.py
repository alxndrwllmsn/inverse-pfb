import numpy as np
import os
import subprocess as sp
import argparse


def genParFile(pars, t, p):
    file = open("tmppars/tmppar{}_{}.txt".format(t, p), "w")
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
    n = np.loadtxt("{}/norms_{}_{}.txt".format(pars['outputdir'], t, p))
    if n.sum() > 0:
        return True
    else:
        return False


def worker(rank, pars, t, p):
    repeat = True
    count = 0
    while repeat:
        # create parfile
        genParFile(pars, t, p)
        amp = 1
        print("amplification: {}, processor: {}".format(pars['amplification'], rank))
        # run ipfb on parfile
        if args.vcs:
            ipfbcall = ['./ipfb', 'tmppars/tmppar{}_{}.txt'.format(t, p), '1']
        else:
            ipfbcall = ['./ipfb', 'tmppars/tmppar{}_{}.txt'.format(t, p), '0']
        try:
            sp.check_output(ipfbcall, stderr=sp.STDOUT)
        except sp.CalledProcessError as e:
            out = e.output.decode('utf-8').split('\n')[-2].split(' ')
            print(e.output.decode('utf-8'))
            amp = float(out[1])
            # check for clipping
        if amp != 1:
            # print("Value clipped, repeating with lower amplification.")
            if int(pars['amplification']) == 1:
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
    file = open(parfile, 'r')
    pardict = {}
    for line in file:
        words = line.split()
        if len(words) != 0:
            pardict[words[0]] = words[1]
    file.close()
    return pardict


def runMultiProcess(args, trange):
    import multiprocessing as mp

    # read parameters from parfile
    pars = readpars(args.parfile)
    # loop over each stream
    jobs = []
    for t in range(trange[0], trange[1]):
        for p in range(2):
            proc = mp.Process(target=worker, args=(pars, t, p,))
            jobs.append(proc)
            proc.start()


def run_MPI(args, trange):
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
        worker(rank, pars, t[rank], p[rank])
        # print("processor {} complete\n".format(rank))

    comm.barrier()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("parfile", help="The parameter file from which to read.")
    parser.add_argument("-t", "--tiles", help="The range of tiles to process (as present within the input file,"
                                              " eg. '-t 0,12')", default=None)
    parser.add_argument("-m", "--mpi", help="Use mpi rather than multiprocess", action="store_true")
    parser.add_argument("-v", "--vcs", help="The input is vcs format", action="store_true")
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
