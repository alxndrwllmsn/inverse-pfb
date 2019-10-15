import numpy as np
import argparse as ap
import os


def rearrange_as_module(directory, nsamples, ntiles, output_file):
    parser = ap.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("nsamples", type=np.int)
    parser.add_argument("ntiles", type=np.int)
    parser.add_argument("output_file")

    args = parser.parse_args([directory, "{}".format(nsamples), "{}".format(ntiles), output_file])
    main(args)


def main(args):
    data = np.zeros((args.nsamples//51200, args.ntiles, 2, 51200, 2), dtype=np.int8)
    owd = os.getcwd()
    os.chdir(args.directory)

    for tile in range(args.ntiles):
        for pol in range(2):
            data[:, tile, pol, :, :] = np.fromfile("out_{}_{}.dat".format(tile, pol), dtype=np.int8).reshape(args.nsamples//51200, 51200, 2)

    header = np.zeros(4096+102400*2*args.ntiles, np.int8)
    file = open(args.output_file, "w")
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