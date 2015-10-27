#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

from PIL import Image
import h5py
import numpy
import sys, os, copy
from numpy import squeeze

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"

    # parse input
    if len(sys.argv) == 3:      
        NF_probs_file   = sys.argv[1]
        GALA_probs_file = sys.argv[2]
    else:
        sys.exit("\nusage: MS_UT_ProbsNF2GALA.py <NF_probs_file> < GALA_probs_file> \n")

    fin   = h5py.File(NF_probs_file, 'r')
    data_stack = numpy.array(fin['/volume']['predictions'])
    data_stack1 = data_stack.transpose((2,1,0,3))
    print "data_stack1.shape=", data_stack1.shape
    data_stack2 = numpy.zeros(data_stack1.shape, dtype = data_stack.dtype)
    data_stack2[:,:,:,0] = data_stack1[:,:,:,0]
    data_stack2[:,:,:,1] = data_stack1[:,:,:,1]
    data_stack2[:,:,:,2] = data_stack1[:,:,:,2]
    data_stack2[:,:,:,3] = data_stack1[:,:,:,3]
    print "data_stack.dtype=", data_stack.dtype, " data_stack2.dtype=", data_stack2.dtype
    fout  = h5py.File(GALA_probs_file, 'w')
    fout.create_dataset('stack', data_stack2.shape, data = data_stack2)

      
