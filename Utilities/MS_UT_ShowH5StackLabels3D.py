#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import os, sys, h5py
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D

if __name__ == "__main__":

    if len(sys.argv) in [3,4]:
        h5_file_name =     sys.argv[1]
        label_id     = int(sys.argv[2])
        layer_id = 0
        if len(sys.argv) == 4:
            layer_id = int(sys.argv[3])
    else:
        sys.exit("\nusage: MS_UT_ShowH5StackLabels3D.py label_id\n")

file    = h5py.File(h5_file_name, 'r')
dataset = numpy.transpose(file['/main'])
num_layers = dataset.shape[2]

if len(sys.argv) == 3: # 3D scatter plot
    X = []
    Y = []
    Z = []
    layer_id = 1
    while layer_id <= num_layers:
        image_data = dataset[:,:,layer_id-1]
        z = layer_id
        xy = numpy.column_stack(numpy.where(image_data == label_id))
        for item in xy:
            X.append(item[0])
            Y.append(item[1])
            Z.append(z)
        print "Collected data from layer", layer_id
        layer_id = layer_id + 1

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(X, Y, Z, s=1, c='b')
    my_title = "Label " + str(label_id)
    plt.show()
    plt.title(my_title)
else:
    while layer_id <= num_layers:
        image_data = dataset[:,:,layer_id-1]
        image_data1 = numpy.zeros(image_data.shape, dtype="int")
        image_data1[image_data == label_id] = 1
        if len(image_data1.shape) == 2:
            print "image_data1.shape == 2"
            plt.imshow(image_data1, cmap = cm.Greys_r) # show as grayscale image
            plt.show()
        elif len(image_data1.shape) == 3:
            print "image_data1.shape == 3"
            image_data1 = image_data[:,:,1]
            plt.imshow(image_data1)
            plt.show()
        else:
            sys.exit("\nProcessing layers of shape " + str(image_data1.shape) + " is not supported\n")
        try:
            print "Input:"
            id = input("Enter layer_id (=0 to exit) or press enter to continue to next layer")
            print "id=", id
            layer_id = int(id)
            if layer_id == 0:
                sys.exit("Exiting ...")
        except SyntaxError:
            layer_id = layer_id + 1
            pass
