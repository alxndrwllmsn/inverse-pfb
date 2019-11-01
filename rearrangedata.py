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
    parser.add_argument("directory")
    parser.add_argument("nsamples", type=np.int)
    parser.add_argument("ntiles", type=np.int)
    parser.add_argument("output_file")

    args = parser.parse_args([directory, "{}".format(nsamples), "{}".format(ntiles), output_file])
    main(args)


def main(args):
    """
    Rearranges the files within the directory specified to a single file rearranged for input into ./ipfb
    :param args: the list of arguments specified above
    :type args: object
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
    if (owd in args.output_file) or (args.output_file[0] == "/"):
        fstring = args.output_file
    else:
        fstring = "{}/{}".format(owd, args.output_file)
    file = open(fstring, "w")
    header.tofile(file)
    data.tofile(file)
    file.close()
    os.chdir(owd)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("nsamples", type=np.int)
    parser.add_argument("ntiles", type=np.int)
    parser.add_argument("output_file")

    args = parser.parse_args()
    main(args)
