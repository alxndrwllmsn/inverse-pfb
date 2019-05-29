import multiprocessing as mp
import numpy as np
import os
import subprocess as sp
import argparse


def genParFile(pars,t,p):
    file = open("tmppar{}_{}.txt".format(t, p),"w")
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
                 'nsamples', pars['nsamples'],
                 'nchannels', pars['nchannels'],
                 'firstchan', pars['firstchan'],
                 'ntiles', pars['ntiles'],
                 'tile', t,
                 'pol', p))

    file.close()


def worker(pars, t, p):
    #create parfile
    genParFile(pars,t,p)
    #run ipfb on parfile
    sp.call(['./ipfb', 'tmppar{}_{}.txt'.format(t, p)])
    #remove parfile
    os.remove('tmppar{}_{}.txt'.format(t, p))



def readpars(parfile):
    file = open(parfile, 'r')
    pardict = {}
    for line in file:
        words = line.split()
        pardict[words[0]] = words[1]
    file.close()
    return pardict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("parfile", help="The parameter file from which to read.")
    parser.add_argument("-t", "--tiles", help="The range of tiles to process (as present within the input file,"
                                              " eg. '-t 0,12')",default=None)
    args = parser.parse_args()
    if args.tiles == None:
        raise ValueError("Please specify a range of tiles (eg. 0,12). Note that there are 24 cores on a magnus node, "
                         "therefore this range should have a multiple of 12 tiles (each tile has 2 polarisations)")
    else:
        trange = np.array(args.tiles.split(','), dtype=np.int)

    # read parameters from parfile
    pars = readpars(args.parfile)
    # loop over each stream
    jobs = []
    for t in range(trange[0],trange[1]):
        for p in range(2):
            proc = mp.Process(target=worker, args=(pars, t, p, ))
            jobs.append(proc)
            proc.start()
