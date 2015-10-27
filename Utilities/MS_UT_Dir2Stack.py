#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
# Purpose: given a directory of image or HDF5 files, create an HDF5 stack
#

from PIL import Image
import h5py
import numpy
import sys, os, re

if __name__ == "__main__":

    transpose = 1
    
    inverse_probs = 1

    if len(sys.argv) == 3 and os.path.isdir(sys.argv[1]):
        indir_name = sys.argv[1]
        outfile    = sys.argv[2]
        match_string = "m4"
        infiles = []
        print "num input files=", len(os.listdir(indir_name))
        for file in os.listdir(indir_name):
            if re.search(match_string, file):
                infile_path = os.path.join(indir_name, file)
                infiles.append(infile_path)
        infiles = sorted(infiles)
        print "infiles=", infiles
    else:
        sys.exit("usage: MS_UT_Dir2Stack.py <input_dir> <output_h5_file> \n")

    print "len(infiles)=", len(infiles)
    if len(infiles) > 0:
        f  = h5py.File(outfile,    'w')
        f1 = h5py.File(infiles[0], 'r')          
        print "\nfile=", infiles[0], " keys=", str(f1.keys())                       
        key   = f1.keys()[0]
        data1 = numpy.squeeze(f1[key])
        data  = numpy.zeros([len(infiles), data1.shape[0], data1.shape[1]], dtype = data1.dtype)
        dtype = data1.dtype
        print "data.shape=", data.shape
        for i in range(0, len(infiles)):
            if re.search(".h5", infiles[i]):
                print "...Processing ", i+1, "/", len(infiles),  infiles[i]
                f1 = h5py.File(infiles[i], 'r')                               
                key = f1.keys()[0]
                data1 = numpy.array(f1[key])
                data[i,:,:] = data1
        f.create_dataset('main', data.shape, data = data, dtype = dtype)
