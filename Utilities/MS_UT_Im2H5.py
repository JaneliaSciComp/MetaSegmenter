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

    data       = numpy.array(Image.open(image_file))
    if 1:
        data_stack = numpy.zeros([1, data.shape[0], data.shape[1]], dtype = numpy.uint32)
        data_stack[0,:,:] = data
    else:
        data_stack = numpy.zeros([   data.shape[0], data.shape[1]], dtype = numpy.uint32)
        data_stack[  :,:] = data
    f  = h5py.File(h5_file_name, 'w')
    f.create_dataset('stack', data_stack.shape, data = data_stack, dtype = numpy.uint32)                       


      
