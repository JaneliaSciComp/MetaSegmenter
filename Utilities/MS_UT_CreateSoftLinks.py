#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# Purpose: in the output dor, create soft links to all the files from the input dir
#

from PIL import Image
import h5py
import numpy
import sys, os, re

if __name__ == "__main__":

    if len(sys.argv) == 3:
        input_dir  = sys.argv[1]
        output_dir = sys.argv[2]
    else:
        sys.exit("usage: MS_UT_CreateSoftLinks.py <input_dir> <output_dir> \n")

    for file in os.listdir(input_dir):
        print "...linking file ", file
        infile_path  = os.path.join(input_dir, file)
        outfile_path = os.path.join(output_dir, file)
        os.system("ln -s " + infile_path + " " + outfile_path)
