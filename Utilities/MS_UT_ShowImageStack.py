#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, re, h5py
import tifffile as tiff
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, ion, show
ion()

if __name__ == "__main__":

    if len(sys.argv) in [2, 3, 4]:
        stack_file_name =     sys.argv[1]
        layer_id = 1
        image_id = 1
        if len(sys.argv) > 2:
            layer_id     = int(sys.argv[2])
        if len(sys.argv) > 3:
            image_id = int(sys.argv[3]) 
    else:
        sys.exit("\nusage: MS_UT_ShowH5Stack.py h5file [ layer_id [ image_id ]]\n")
if re.search(".h5", stack_file_name):
    f = h5py.File(stack_file_name, 'r')
    key = f.keys()[0]
    if "segmentations" in f.keys():
        key = "segmentations"
    data0 = f[key]
    print "data0.shape=", data0.shape, " data0.dtype=", data0.dtype
    data =                 f[key]
    print "data.shape=", data.shape
    # Layer id is supposed to be the last dimension:
    if len(data.shape) > 2 and data.shape[0] < min(data.shape[1], data.shape[2]):
        data = numpy.transpose(f[key])
elif re.search(".tif", stack_file_name):
    data = numpy.transpose(tiff.imread(stack_file_name))
else:
    sys.exit("Unrecognized stack type")

data_shape = data.shape
print "data_shape =", data.shape
if len(data.shape) > 2:
    num_layers = data.shape[2]
    print " layer_id=", layer_id, " num_layers=", num_layers
else:
    num_layers = 1

while layer_id <= num_layers:
    if len(data.shape) == 2:
        image_data = numpy.uint8(data)
    elif len(data.shape) == 3:
        image_data = numpy.uint8(data[:,:,layer_id-1])
        #image_data = data[layer_id-1,:,:]
        if len(sys.argv) == 4:
            print "\nWarning: ignoring redundant imadge_id\n"
    else:
        image_data = data[:,:,layer_id-1,:]
#           sys.exit("\nLayer contains > 1 image. Please, specify image_id\n")

    print "Layer=", layer_id, " max_label=", numpy.max(image_data), \
          " image_data.shape=", image_data.shape, " image_data.dtype=", image_data.dtype

    plot_3 = 0
    if len(image_data.shape) == 2:
#       print "max=", numpy.amax(image_data), " min=", numpy.amin(image_data)
        image_data[image_data > 127] = 0
        plt.imshow(image_data, cmap = cm.Greys_r) # show as grayscale image
        plt.show()
    elif len(image_data.shape) == 3:
        if plot_3:
            fig = plt.figure(figsize=(12, 9))
            num_images = int(image_data.shape[2])
            for i in range(0, num_images):
               a=fig.add_subplot(1,num_images,i+1)
               image_data1 = image_data[:,:,i]
#              imgplot = plt.imshow(image_data1)                     # RGB
               imgplot = plt.imshow(image_data1, cmap = cm.Greys_r)  # Grayscale
               a.set_title('Image' + str(i+1))
               plt.colorbar(ticks=[0.1,0.3,0.5,0.7], orientation ='horizontal')
            show()
        else:
            image_data1 = image_data[:,:,1]
#           print "max=", numpy.amax(image_data1), " min=", numpy.amin(image_data1)
            plt.imshow(image_data1) 
            plt.show()
    else:
        sys.exit("\nProcessing laters of shape " + str(image_data.shape) + " is not supported\n")
    try:
        id = input("Enter layer_id (=0 to exit) or press enter to continue to next layer")
        print "id=", id
        layer_id = int(id)
        if layer_id == 0:
            sys.exit("Exiting ...")
    except SyntaxError:
        layer_id = layer_id + 1
        pass



