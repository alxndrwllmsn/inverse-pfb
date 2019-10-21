import numpy as np
import argparse as ap
import os


def rearrange_as_module(directory, nsamples, ntiles, output_file):
    """
    Runs the script to be called from another
    :param directory: the directory of the input files to be reformatted
    :type directory: str
    :param nsamples: the number of samples in the files to be reformatted
    :type nsamples: int
    :param ntiles: the number of tiles in the files
    :type ntiles: int
    :param output_file: the output file name after reformatting
    :type output_file: str
    """
    parser = ap.ArgumentParser()
    parser.add_argument("directory", help="The directory of the out files or the name of the vcs input files")
    parser.add_argument("nsamples", help="The number of samples in the input file", type=np.int)
    parser.add_argument("ntiles", help="The number of tiles to process", type=np.int)
    parser.add_argument("output_file", help="The name of the output file")
    parser.add_argument("-v", "--vcs", help="use this flag to change the mode to read in vcs data (set directory to be "
                                            "the filename of the vcs file", action="store_true")

    args = parser.parse_args([directory, "{}".format(nsamples), "{}".format(ntiles), output_file])
    if not args.vcs:
        out_to_htr_vcs(args)
    else:
        vcs_to_htr(args)


def out_to_htr_vcs(args):
    """
    Rearranges the files within the directory specified to a single file rearranged for input into ./ipfb
    :param args: the list of arguments specified above
    """
    data = np.zeros((args.nsamples//51200, args.ntiles, 2, 51200, 2), dtype=np.int8)
    owd = os.getcwd()
    print(args.directory)
    os.chdir(args.directory)

    for tile in range(args.ntiles):
        for pol in range(2):
            data[:, tile, pol, :, :] = np.fromfile("out_{}_{}.dat".format(tile, pol), dtype=np.int8).reshape(args.nsamples//51200, 51200, 2)

    header = np.zeros(4096+102400*2*args.ntiles, np.int8)
    print(owd, args.output_file)
    file = open("{}/{}".format(owd,args.output_file), "w")
    header.tofile(file)
    data.tofile(file)
    file.close()
    os.chdir(owd)


def vcs_to_htr(args):
    """
    Converts the vcs file into one that can be read into ./ipfb more efficiently, unfortunately it is too slow.
    :param args: the script input arguments
    """
    data = np.fromfile(args.directory, dtype=np.int8)
    filesize = data.shape[0]
    print(filesize)
    nchans = ((filesize//args.nsamples)//args.ntiles)//2
    print(nchans)
    data = data.reshape(args.nsamples, nchans, args.ntiles, 2)
    data = np.transpose(data, (1, 2, 3, 0))
    data.tofile(args.output_file)

if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("directory", help="The directory of the out files or the name of the vcs input files")
    parser.add_argument("nsamples", help="The number of samples in the input file", type=np.int)
    parser.add_argument("ntiles", help="The number of tiles to process", type=np.int)
    parser.add_argument("output_file", help="The name of the output file")
    parser.add_argument("-v", "--vcs", help="use this flag to change the mode to read in vcs data (set directory to be "
                                            "the filename of the vcs file", action="store_true")

    args = parser.parse_args()
    if not args.vcs:
        out_to_htr_vcs(args)
    else:
        vcs_to_htr(args)
