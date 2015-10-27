#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, re, h5py
import tifffile as tiff
import numpy
from skimage.morphology import label
from scipy.ndimage import measurements,morphology
import matplotlib.pyplot as plt
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
    if len(sys.argv) == 3:          
        stack_file_name = sys.argv[1]
        dir_name        = sys.argv[2]
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)   
    else:
        sys.exit("\nusage: MS_UT_Stack2Dir.py input_stack_file output_dir_name\n")

    show_image = 0

    data = read_input(stack_file_name)
    print "data_shape =", data.shape, " num_layers=", data.shape[2]

    layer_id = 1
    num_layers = data.shape[2]
    while layer_id <= num_layers:
        print "\nlayer_id=", layer_id
        im = data[:,:,layer_id-1]
        labels1, nbr_objects = measurements.label(im)
        print "Number of objects:", nbr_objects, " labels1.shape=", labels1.shape, \
              " labels1.dtype=", labels1.dtype, \
              " max_label=", numpy.max(labels1), " min_label=", numpy.min(labels1)
        

        if len(stack_file_name2) > 0:
            im2 = data2[:,:,layer_id-1]
            labels2, nbr_objects2 = measurements.label(im2)
            print "    Number of objects2:", nbr_objects2, " labels2.shape=", \
                  labels2.shape, " labels2.dtype=", labels2.dtype, \
                  " max_label2=", numpy.max(labels2), " min_label2=", numpy.min(labels2)
            for i in range(1, nbr_objects2+1):
                labels1[labels2 == i] = nbr_objects + i
        
        labels[:,:,layer_id-1] = labels1
        print "max label=", numpy.max(labels1), " min label=", numpy.min(labels1)

        if show_image:
            img = Image.fromarray(labels1) #.convert('RGB')
            plt.imshow(img)
            plt.show()

            try:
                id = input("Enter layer_id (=0 to exit) or press enter to continue to next layer")
                print "id=", id
                layer_id = int(id)
                if layer_id == 0:
                    sys.exit("Exiting ...")
            except SyntaxError:
                layer_id = layer_id + 1
                pass
        else:
            layer_id = layer_id + 1

    output_stack_file_name = stack_file_name.split('.')[0] + ".h5"
    print "output_stack_file_name=", output_stack_file_name
    f  = h5py.File(output_stack_file_name, 'w')
    output_data = numpy.transpose(labels)
    f.create_dataset('main', output_data.shape, data=output_data)               



