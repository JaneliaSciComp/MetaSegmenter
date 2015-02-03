#! /usr/local/python-2.7.8/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import sys, h5py
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.pyplot import figure, show

if __name__ == "__main__":

    if len(sys.argv) in [3, 4]:
        h5_file_name =     sys.argv[1]
        layer_id     = int(sys.argv[2])
        if len(sys.argv) == 4:
            image_id = int(sys.argv[3]) 
    else:
        sys.exit("\nusage: MS_UT_ShowH5Stack.py h5file layer_id [ image_id ]\n")

file    = h5py.File(h5_file_name, 'r')
dataset = numpy.transpose(file['/main'])
data_shape = dataset.shape
if len(dataset.shape) == 3:
    image_data = dataset[:,:,layer_id-1]
    #image_data = dataset[layer_id-1,:,:]
    if len(sys.argv) == 4:
        print "\nWarning: ignoring redundant imadge_id\n"
else:
    if len(sys.argv) == 4:
        image_data = dataset[:,:,layer_id-1,image_id-1]
    else:
        image_data = dataset[:,:,layer_id-1,:]
#       sys.exit("\nLayer contains > 1 image. Please, specify image_id\n")

print "Stack shape=", dataset.shape
print "Layer shape=", image_data.shape

if len(image_data.shape) == 2:
    plt.imshow(image_data, cmap = cm.Greys_r) # show as grayscale image
    plt.show()
elif len(image_data.shape) == 3:
    fig = plt.figure(figsize=(12, 9))
    num_images = int(image_data.shape[2])
    for i in range(0, num_images):
       a=fig.add_subplot(1,num_images,i+1)
       image_data1 = image_data[:,:,i]
#      imgplot = plt.imshow(image_data1)                     # RGB
       imgplot = plt.imshow(image_data1, cmap = cm.Greys_r)  # Grayscale
       a.set_title('Image' + str(i+1))
       plt.colorbar(ticks=[0.1,0.3,0.5,0.7], orientation ='horizontal')
    show()
else:
    sys.exit("\nProcessing laters of shape " + str(image_data.shape) + " is not supported\n")
exit()


