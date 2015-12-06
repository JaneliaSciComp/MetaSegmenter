#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: given a segmentation stack,
# produce a directory of RGB files

import os, sys, re, h5py
import tifffile as tiff
import numpy
from skimage.morphology import label
from scipy.ndimage import measurements
from skimage import filters            
from skimage.color import label2rgb
from PIL import Image

def read_input(stack_file_name):
    if re.search(".h5", stack_file_name):
        f = h5py.File(stack_file_name, 'r')
        key = f.keys()[0]
        data = numpy.transpose(f[key])
    elif re.search(".tif", stack_file_name):
        data = numpy.transpose(tiff.imread(stack_file_name))
    else:
        sys.exit("Unrecognized stack type")
    return data

# -------------------------------------------------------------------------------

if __name__ == "__main__":

    print "len(sys.argv)=", len(sys.argv)
    if len(sys.argv) in [3,4]:          
        stack_file_name = sys.argv[1]
        out_dir_name    = sys.argv[2]
        if not os.path.isdir(out_dir_name):
            os.mkdir(out_dir_name) 
        first_file_name = 'crop.00004225'
        if len(sys.argv) == 4:
            first_file_name = sys.argv[3]  
    else:
        sys.exit("\nusage: MS_UT_Stack2Dir.py seg_stack_file output_dir_name [ first_file_name] \n")

    data = read_input(stack_file_name)
    print "data_shape =", data.shape, " num_layers=", data.shape[2]

    layer_id = 1
    num_layers = data.shape[2]
    while layer_id <= num_layers:
        data1  = numpy.transpose(numpy.float64(numpy.squeeze(data[:,:,layer_id-1])))
        labels = numpy.unique(data1)
        # If label 0 is missing, detect edges of components and set their label to 0
        print "min_label=", min(labels), " data1.dtype=", data1.dtype
        edges = filters.sobel(data1) # detect e
        data1[edges > 0] = -1

        output_file_name = "0000" + str(layer_id) + "_RGB.png"
        if len(str(layer_id)) == 2:
            output_file_name = "000" + str(layer_id) + "_RGB.png"
        if len(str(layer_id)) == 3:
            output_file_name = "00"  + str(layer_id) + "_RGB.png"
        if len(str(layer_id)) == 4:
            output_file_name = "0"   + str(layer_id) + "_RGB.png"
        out_path = os.path.join(out_dir_name, output_file_name)
        outdata1 = label2rgb(data1)
        outdata3 = numpy.zeros([outdata1.shape[0], outdata1.shape[1], 3], 'uint8')
        for i in [0,1,2]:
            outdata3[:,:,i] = numpy.uint8(outdata1[:,:,i]/numpy.max(outdata1[:,:,i])*255.)
        img = Image.fromarray(outdata3)
        img.save(out_path)
        print "layer_id=", layer_id, " out_path=", out_path
        layer_id = layer_id + 1


