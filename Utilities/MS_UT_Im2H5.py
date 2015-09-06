#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os, copy

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"

    # parse input
    if len(sys.argv) == 3:      
        image_file   = sys.argv[1]
        h5_file_name = sys.argv[2]
    else:
        sys.exit("\nusage: MS_UT_Im2H5.py <image_file> <hdf5_file> \n")

    data_stack = numpy.array(Image.open(image_file))
    f  = h5py.File(h5_file_name, 'w')
    f.create_dataset('stack', data_stack.shape, data = numpy.transpose(data_stack))


      
