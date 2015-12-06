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

# Purpose: count and report the # of pixels corresponding to each label in h5 stack
    if len(sys.argv) == 2:
        h5_file_name =     sys.argv[1]
    else:
        sys.exit("usage: MS_UT_CountH5StackLabels.py h5file\n")

    f = h5py.File(h5_file_name, 'r')
    key = f.keys()[0]
    data = numpy.transpose(f[key])
    print "type=", data.dtype
    print "data.shape=", numpy.array(data).shape

    data_shape = numpy.array(data).shape
    num_labels_z = []
    if len(data_shape) > 2:
        for i in range (0, data_shape[2]):
            labels_z = numpy.unique(data[:,:,i])
            print "layer=", i, " num_zero labels=", (data[:,:,i] == 0).sum()
            num_labels_z.append(len(labels_z))
        print "\nnum_labels_z=", num_labels_z

    num_labels = {}
    all_labels = list(numpy.unique(data))
    print "\nall_labels=", all_labels
    largest_label = ""
    largest_label_count = 0
    for label in all_labels:
       num_labels[str(label)] = (data == label).sum()
       num_layers = 0
       for i in range(data_shape[2]):
           if (data[:,:,i] == label).sum() > 0:
               num_layers += 1
       if int(label) > 0 and  largest_label_count < num_labels[str(label)]:
           largest_label_count = num_labels[str(label)]
           largest_label       = str(label)
       print "label=", label, " size=", num_labels[str(label)], " num_layers=", num_layers, \
             " largest_label=", largest_label, " largest_size=", largest_label_count
       if len(data_shape) > 2:
           print  " frac_size=", float(largest_label_count)/float(data_shape[0]*data_shape[1]*data_shape[2])
    print "\nlargest_3D_label=", largest_label 
