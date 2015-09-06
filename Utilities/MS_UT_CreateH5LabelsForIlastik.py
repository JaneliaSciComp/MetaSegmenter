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
    if len(sys.argv) in [3, 4]:
        input_dir_memb = sys.argv[1]
        h5_file_name   = sys.argv[2]
        input_dir_mito = ''
        if len(sys.argv) == 4:
            input_dir_mito = sys.argv[3]
    else:
        sys.exit("\nusage: MS_UT_CreateH5LabelsForIlastik.py <input_dir_memb_labels> <input_dir_mito_labels> { h5_file ] \n")

    print "len(sys.argv) =", len(sys.argv), " h5_file_name=", h5_file_name
    input_files_memb = []
    for file in os.listdir(input_dir_memb):
        path1 = os.path.join(input_dir_memb, file)
        input_files_memb.append(path1)
    input_files1 = sorted(input_files_memb)

    if len(input_dir_mito) > 0:
        input_files_mito = []
        for file in os.listdir(input_dir_mito):
            path1 = os.path.join(input_dir_mito, file)
            input_files_mito.append(path1)
        input_files2 = sorted(input_files_mito)

    f  = h5py.File(h5_file_name, 'w')

    for i in range(0, len(input_files1)): 
        data_memb = copy.copy(numpy.asarray(Image.open(input_files1[i], 'r'))) # 0 = memb, 255 = cells
        data_memb[data_memb == 255] = 1
        if len(input_dir_mito) > 0:
            data_mito =           numpy.asarray(Image.open(input_files2[i], 'r'))  # 255 = mito, 0 = the rest
            data_memb[data_mito == 255] = 2
        if i == 0:
            my_shape = [len(input_files1), data_memb.shape[0], data_memb.shape[1]]
            data_stack = numpy.zeros(my_shape, dtype = numpy.uint8)
        data_stack[i, :,:] = data_memb

#           my_shape = [data_memb.shape[0], data_memb.shape[1], len(input_files1)]
#           data_stack = numpy.zeros(my_shape, dtype = numpy.uint8)
#       data_stack[:,:, i] = data_memb
    print "my_shape=", my_shape
    f.create_dataset('main', my_shape, data = numpy.transpose(data_stack))
        


      
