#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, re, h5py
import tifffile as tiff
import numpy as np
from skimage.color import label2rgb
from skimage import filters
from PIL import Image

def read_input(stack_file_name, layer_id):
    if re.search(".h5", stack_file_name):
        f = h5py.File(stack_file_name, 'r')
        key = f.keys()[0]
        data = f[key]
    elif re.search(".tif", stack_file_name):
        data = tiff.imread(stack_file_name)
    else:
        sys.exit("Unrecognized stack type")
    print "data.shape=", data.shape, " layer_id=", layer_id
    return data[layer_id,:,:]

# -------------------------------------------------------------------------------

if __name__ == "__main__":

    print "len(sys.argv)=", len(sys.argv)
    if len(sys.argv) in [3,4]:      
        stack_file_name = sys.argv[1]
        RGB_name        = sys.argv[2]
        layer_id        = 0
        if len(sys.argv) == 4:
            layer_id = int(sys.argv[3])
    else:
        sys.exit("\nusage: MS_UT_SegStack2RGB.py input_stack_file output_RGB_file\n")

    data = np.squeeze(read_input(stack_file_name, layer_id))
    print "data.shape=", data.shape
    max_value = data.max()
    min_value = data.min()
    print "data.max=", max_value, " data.min=", min_value 

    labels = np.unique(data)
    # If label 0 is missing, detect edges of components and set their label to 0
    if not 0 in labels:
        edges = filters.sobel(data) # detect e
        data[edges > 0] = 0

    if round(max_value) == max_value and round(min_value) == min_value and max_value - min_value > 1:
        num_values = 0
        for i in range(min_value, max_value + 1):
            if (data == i).sum() > 0:
                num_values += 1
        print "num_values=", num_values

    outdata = label2rgb(data)
    print "\noutdata.shape=", outdata.shape
    max_value = outdata.max()
    min_value = outdata.min()
    print "outdata.max=", max_value, " outdata.min=", min_value

    outdata2 = np.zeros([outdata.shape[0], outdata.shape[1], 3], 'uint8')
    for i in [0,1,2]:
        outdata2[:,:,i] = np.uint8(outdata[:,:,i]/np.max(outdata[i,:,:])*255.)
    img = Image.fromarray(outdata2)
    img.save(RGB_name)




