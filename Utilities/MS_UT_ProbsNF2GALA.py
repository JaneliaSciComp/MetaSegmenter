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
        NF_probs_file   = sys.argv[1]
        GALA_probs_file = sys.argv[2]
    else:
        sys.exit("\nusage: MS_UT_ProbsNF2GALA.py <NF_probs_file> < GALA_probs_file> \n")

    fin   = h5py.File(NF_probs_file, 'r')
    data_stack = numpy.array(fin['/volume']['predictions'])
    fout  = h5py.File(GALA_probs_file, 'w')
    fout.create_dataset('stack', data_stack.shape, data = data_stack)


      
