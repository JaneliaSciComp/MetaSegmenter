#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import re
import matplotlib
from PIL import Image
import h5py
import numpy
import sys, os, copy
from scipy import misc

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"

    # parse input
    if len(sys.argv) in [2, 3]:
        input_png_dir  = sys.argv[1]
        h5_file_name   = "groundtruth.h5"
        if len(sys.argv) == 3:
            h5_file_name = sys.argv[2]
    else:
        sys.exit("\nusage: MS_UT_CreateH5Groundtruth.py  <groundtruth png dir> [ output h5_file ] \n")

    print "len(sys.argv) =", len(sys.argv), " h5_file_name=", h5_file_name
    input_files = []
    for file in os.listdir(input_png_dir):
        path1 = os.path.join(input_png_dir, file)
        input_files.append(path1)

    input_files1 = sorted(input_files)
    f  = h5py.File(h5_file_name, 'w')

    for i in range(0, len(input_files1)): 
        data = copy.copy(numpy.asarray(Image.open(input_files1[i], 'r'))) 
        if i == 0:
            my_shape = [len(input_files1), data.shape[0], data.shape[1]]
            data_stack = numpy.zeros(my_shape, dtype = numpy.int32)
        data_stack[i, :,:] = data
    f.create_dataset('stack', my_shape, data = data_stack)


      
