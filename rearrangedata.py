import numpy as np
import argparse as ap
import os

if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument("directory")
    parser.add_argument("nsamples")
    parser.add_argument("ntiles")

    args = parser.parse_args()

    data = np.zeros((args.nsamples//51200, args.ntiles, 2, 51200, 2), dtype=np.int8)

    os.chdir(args.directory)

    for tile in range(args.ntiles):
        for pol in range(2):
            data[:, tile, pol, :, :] = np.fromfile("out_{}_{}.dat", dtype=np.int8).reshape(args.nsamples//51200, 51200, 2)

    np.tofile(data)
