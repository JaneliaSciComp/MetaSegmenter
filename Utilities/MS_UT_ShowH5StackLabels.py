#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import Image
import h5py
import numpy
import sys, os
from scipy import misc

if __name__ == "__main__":

    # parse input
#   print "len(sys.argv)=", len(sys.argv), "sys.argv=", sys.argv
    if len(sys.argv) == 3:
        h5_file_name =     sys.argv[1]
        layer_id     = int(sys.argv[2])
    else:
        sys.exit("usage: MS_UT_ShowH5StackLabels.py h5file layer_id \n")

    f = h5py.File(h5_file_name, 'r')
    key = f.keys()[0]
    data = numpy.transpose(f[key])
    print "type=", type(data)
    print "type=", type(numpy.array(data))
    print "data.shape=", numpy.array(data).shape
#   data1 = data[0,:,:,layer_id]
#   data2 = data[1,:,:,layer_id]
#   data3 = data[2,:,:,layer_id]
    print "data=", data
    print "max=", numpy.max(data), " min=", numpy.min(data)
#   data1_new = data1;
#   data1 == data2
#   data1 
#:red
#   data1_new(data1 == data2) = "."
