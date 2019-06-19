import numpy as np
import os
import subprocess as sp
import argparse


def genParFile(pars,t,p):
    file = open("tmppars/tmppar{}_{}.txt".format(t, p),"w")
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


def worker(pars, t, p):
    repeat = True
    count = 0
    while repeat:
        # create parfile
        genParFile(pars, t, p)
        # run ipfb on parfile
        out = sp.check_output(['./ipfb', 'tmppars/tmppar{}_{}.txt'.format(t, p)], stderr=sp.STDOUT)
        print(out.decode('utf-8'))
        # check for clipping
        if clipCheck(pars, t, p):
            print("Value clipped, repeating with lower amplification.")
            if pars['amplification'] == 1:
                repeat = False
                print("amplification is already 1, cannot reduce anymore (there is an issue).")
            else:
                repeat = True
            pars['amplification'] = str(int(pars['amplification'])//2)
        # remove parfile
        # os.remove('tmppars/tmppar{}_{}.txt'.format(t, p))



def readpars(parfile):
    file = open(parfile, 'r')
    pardict = {}
    for line in file:
        words = line.split()
        pardict[words[0]] = words[1]
    file.close()
    return pardict


def runMultiProcess(args,trange):
    import multiprocessing as mp

    # read parameters from parfile
    pars = readpars(args.parfile)
    # loop over each stream
    jobs = []
    for t in range(trange[0],trange[1]):
        for p in range(2):
            proc = mp.Process(target=worker, args=(pars, t, p, ))
            jobs.append(proc)
            proc.start()


def run_MPI(args,trange):
    from mpi4py import MPI
    pars = readpars(args.parfile)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    nstreams = 2*(trange[1] - trange[0])
    if nstreams != nprocs:
        raise ValueError("Please set the number of processors to be equal to the number of voltage streams"
                         " (i.e. 2 times the number of tiles).")
    else:
        t, p = np.mgrid[trange[0]:trange[1]:1, 0:2:1]
        t = t.reshape(nstreams)
        p = p.reshape(nstreams)
        worker(pars, t[rank], p[rank])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("parfile", help="The parameter file from which to read.")
    parser.add_argument("-t", "--tiles", help="The range of tiles to process (as present within the input file,"
                                              " eg. '-t 0,12')",default=None)
    parser.add_argument("-m","--mpi",help="Use mpi rather than multiprocess",action="store_true")
    args = parser.parse_args()
    if args.tiles == None:
        raise ValueError("Please specify a range of tiles (eg. 0,12). Note that there are 24 cores on a magnus node, "
                         "therefore this range should have a multiple of 12 tiles (each tile has 2 polarisations)")
    else:
        trange = np.array(args.tiles.split(','), dtype=np.int)

    if args.mpi:
        print("Running with MPI")
        run_MPI(args,trange)
    else:
        print("Running with Multiprocess")
        runMultiProcess(args,trange)
