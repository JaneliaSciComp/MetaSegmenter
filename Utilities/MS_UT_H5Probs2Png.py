#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os

if __name__ == "__main__":

    if len(sys.argv) in [2, 3]:
        h5_file_name =     sys.argv[1]
        png_file_name = h5_file_name.split(".")[0] + ".png"
        if len(sys.argv) == 3:
            png_file_name = sys.argv[2]
    else:
        sys.exit("usage: MS_UT_H5ToPng.py h5file [ output_png_file ] \n")

    f  = h5py.File(h5_file_name, 'r')   
    data = f["probabilities"]
    img  = Image.fromarray(numpy.uint8(data/numpy.max(data)*255.))
    img.save(png_file_name)               
