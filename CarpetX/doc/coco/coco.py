"""COCO - CarpetX OpenPMD COmpanion.

Usage:
  coco.py list <openPMD-file> [--data-dir=<dir>]
  coco.py (-h | --help)
  coco.py --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  --data-dir=<dir>  The directory where data files are located [default: .].
"""

from docopt import docopt
import openpmd_api as opmd

import os


def list_file(args):
    file = args["<openPMD-file>"]

    print("Trying to open", file, "as a openPMD file")

    series = opmd.Series(file, opmd.Access.read_only)

    print("File loaded\n")

    num_iter = len(series.iterations)

    print("The file contains", num_iter, "iterations:")
    for i in series.iterations:
        print(" ", i)
    print("")

    for i in series.iterations:
        iteration = series.iterations[i]
        meshes = iteration.meshes
        particles = iteration.particles

        num_meshes = len(meshes)
        num_particles = len(particles)

        if num_meshes != 0:
            print("Iteration", i, "has", num_meshes,
                  "meshes and", num_particles, "particles")

            for m in meshes:
                print(" ", m)
            print("")

        if num_particles != 0:
            for p in particles:
                print(" ", p)
            print("")


def main():
    args = docopt(__doc__, version="COCO 1.0.0")

    if (args["list"]):
        list_file(args)


if __name__ == '__main__':
    main()
