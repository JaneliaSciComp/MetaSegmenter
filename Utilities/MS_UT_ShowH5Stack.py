#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, h5py
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, ion, show
ion()

if __name__ == "__main__":

    if len(sys.argv) in [2, 3, 4]:
        h5_file_name =     sys.argv[1]
        layer_id = 1
        image_id = 1
        if len(sys.argv) > 2:
            layer_id     = int(sys.argv[2])
        if len(sys.argv) > 3:
            image_id = int(sys.argv[3]) 
    else:
        sys.exit("\nusage: MS_UT_ShowH5Stack.py h5file [ layer_id [ image_id ]]\n")

file    = h5py.File(h5_file_name, 'r')
dataset = numpy.transpose(file['/main'])
data_shape = dataset.shape
num_layers = dataset.shape[2]
print "data_shape =", dataset.shape, " layer_id=", layer_id, " num_layers=", num_layers

while layer_id <= num_layers:
    if len(dataset.shape) == 3:
        image_data = dataset[:,:,layer_id-1]
        #image_data = dataset[layer_id-1,:,:]
        if len(sys.argv) == 4:
            print "\nWarning: ignoring redundant imadge_id\n"
    else:
        image_data = dataset[:,:,layer_id-1,:]
#           sys.exit("\nLayer contains > 1 image. Please, specify image_id\n")

    print "Layer=", layer_id

    plot_3 = 0
    if len(image_data.shape) == 2:
#       print "max=", numpy.amax(image_data), " min=", numpy.amin(image_data)
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



