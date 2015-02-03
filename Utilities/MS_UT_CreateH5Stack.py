#! /usr/local/python-2.7.8/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import Image
import h5py
import numpy
import sys, os, optparse
from scipy import misc

# ----------------------------------------------------------------------

def command_line_parser(parser):
    parser.add_option("-a", "--alpha", dest="alpha", help="smoothing coefficient",  metavar="alpha", default="0.01")
    parser.add_option("-i", "--image_dir",dest="image_dir",help="optional when type=labels')", metavar="imd", default="")
    parser.add_option("-o", "--output_name", dest="output_name", help="name of the output HDF5 file", metavar="out", default="output.h5")
    parser.add_option("-t", "--type", dest="output_type", help="output type (='data','labels' or 'mask')", metavar="ot", default="")
    parser.add_option("-v", "--verbose",action="store_true",dest="verbose",help="increase verbosity of output", default=False)
    return parser

# ----------------------------------------------------------------------

# Set a label at pixel i to 1 if pixels i and i+1 are "connected",
# i.e. are both contained in an object
# (this idea applies separately to x-, y- and z-direction)

def update_labels_connect(old_labels_stack):
    print "...Updating labels by connecting"
    stack_shape = old_labels_stack.shape
    new_labels_stack = numpy.zeros(stack_shape, dtype="int", order='F')
    print "   x-dimension ..."
    for j in range(0, stack_shape[1]):
        for k in range(0, stack_shape[2]):
            for i in range(0, stack_shape[0]-1):
                if old_labels_stack[i,j,k,0] * old_labels_stack[i+1,j,k,0] == 1:
                    new_labels_stack[i,j,k,0] = 1
    print "   y-dimension ..."
    for i in range(0, stack_shape[0]):
        for k in range(0, stack_shape[2]):
            for j in range(0, stack_shape[1]-1):
                if old_labels_stack[i,j,k,1] * old_labels_stack[i,j+1,k,1] == 1:
                    new_labels_stack[i,j,k,1] = 1
    print "   z-dimension ..."
    for i in range(0, stack_shape[0]):
        for j in range(0, stack_shape[1]):
            for k in range(0, stack_shape[2]-1):
                if old_labels_stack[i,j,k,2] * old_labels_stack[i,j,k+1,2] == 1:
                    new_labels_stack[i,j,k,2] = 1
    return new_labels_stack

# ----------------------------------------------------------------------

def update_labels_smooth(labels_stack, image_stack, mean_signal, alpha):
    print "...Updating labels by smoothing"
    num_iter = 3
    stack_shape = labels_stack.shape
    print "old_labels_stack.shape=", stack_shape, " image_stack.shape=", image_stack.shape    
    for l in range(3):
        if l == 0:
            # Update in x-direction
            for j in range(stack_shape[1]):                # y-direction
                for k in range(stack_shape[2]) :           # z-direction
                    for iter in range(num_iter):
                        for i in range(stack_shape[0]-1):  # x-direction
                            lab0 = labels_stack[i  ,j,k,l]
                            lab1 = labels_stack[i+1,j,k,l]
                            if lab0 + lab1 == 1:
                                s0 = image_stack[i  ,j,k]
                                s1 = image_stack[i+1,j,k]
                                der = abs(s1 - s0)
#                               if der > 0:
#                                  print "der=" der " sig*alpha=" mean_signal*alpha ' update=' (der < mean_signal*alpha)
#                               end
                                if der < mean_signal*alpha:
                                    labels_stack[i+1,j,k,l] = 1;
                                    labels_stack[i  ,j,k,l] = 1;
        elif l == 1:
            # Update in y-direction
            for i in range(stack_shape[0]):                # x-direction
                for k in range(stack_shape[2]):            # z-direction
                    for iter in range(num_iter):
                        for j in range(stack_shape[1]-1):  # x-direction
                            lab0 = labels_stack[i,j  ,k,l]
                            lab1 = labels_stack[i,j+1,k,l]
                            if lab0 + lab1 == 1:
                                s0 = image_stack[i,j  ,k]
                                s1 = image_stack[i,j+1,k]
                                der = abs(s1 - s0)
#                               if der > 0:
#                                  print "der=" der " sig*alpha=" mean_signal*alpha ' update=' (der < mean_signal*alpha)
#                               end
                                if der < mean_signal*alpha:
                                    labels_stack[i,j  ,k,l] = 1;
                                    labels_stack[i,j+1,k,l] = 1;
        else: # l == 2
            # Update in z-direction
            for i in range(stack_shape[0]):                # x-direction
                for j in range(stack_shape[1]):            # y-direction
                    for iter in range(num_iter):
                        for k in range(stack_shape[2]-1):  # z-direction
                            lab0 = labels_stack[i,j,k+1,l]
                            lab1 = labels_stack[i,j,k+1,l]
                            if lab0 + lab1 == 1:
                                s0 = image_stack[i,j,k  ]
                                s1 = image_stack[i,j,k+1]
                                der = abs(s1 - s0)
#                               if der > 0:
#                                  print "der=" der " sig*alpha=" mean_signal*alpha ' update=' (der < mean_signal*alpha)
#                               end
                                if der < mean_signal*alpha:
                                    labels_stack[i,j,k  ,l] = 1;
                                    labels_stack[i,j,k+1,l] = 1;
    return labels_stack
 
# ----------------------------------------------------------------------

# In the default stack, pixel values = 1 for data type=mask and 0 otherwise
def create_default_stack(image_shape, num_images, output_type):
    print "...Creating a default stack for output_type=", output_type
    stack_shape = (image_shape[0],image_shape[1],num_images)
    if len(image_shape) > 2 or not output_type == "data":
        stack_shape = (image_shape[0],image_shape[1],num_images,3)
    print "image_shape=", image_shape, " stack_shape=", stack_shape
    dtype = "float"
    if not output_type == "data":
        dtype = "int" 
    if output_type == "mask":
        nd_stack = numpy.ones(stack_shape, dtype=dtype)  # 1 for mask
    else:
        nd_stack = numpy.zeros(stack_shape, dtype=dtype)
    print "Default stack shape=", nd_stack.shape
    return nd_stack

# ----------------------------------------------------------------------

def scale_image(orig_image):
    std  = numpy.std( orig_image)
    scaled_image = orig_image/ std
    return scaled_image

# ----------------------------------------------------------------------

def normalize_image(orig_image):
    mean = numpy.mean(orig_image)
    std  = numpy.std( orig_image)
    norm_image = (orig_image - mean)/ std
    return norm_image

# ----------------------------------------------------------------------

def populate_stack(nd_stack, input_dir, input_files, output_type, imdir):
    print "...Populating stack for output_type=", output_type
    i = 0
    for ifile in input_files:
        i = i + 1
        ifile_path = input_dir + "/" + ifile
        im = misc.imread(ifile_path)
#       if output_type == "data" and len(imdir) == 0:
#           im = normalize_image(im)
#       elif output_type == "data" and len(imdir) > 0:
#           im = scale_image(im)
        print "Updating stack with ", ifile, " (image shape=", im.shape, ") ..."
        nd_stack = update_nd_stack(nd_stack, im, i-1, output_type)
    return nd_stack

# ----------------------------------------------------------------------

# For data stack, set pixel value t0 intensity of input image
# For labels stack, set pixel value initially to 1 if the pixel is within object
# (this value will be updated later)
# For mask stack, leave pixels unchanged (they have been already set all to 1)
def update_nd_stack(nd_stack, im, k, output_type):
    stack_shape = nd_stack.shape
#   print "stack_shape=", stack_shape, " im.shape=", im.shape, " k=", k
    if output_type == "data":
        for j in range(0, stack_shape[1]):
#           print "j=", j, " im[:,j].shape=", im[:,j].shape, " nd_stack[:,j,k].shape=", nd_stack[:,j,k].shape
            nd_stack[:,j,k] = im[:,j]
    else:
        if len(im.shape) == 3 and output_type == "labels":
            for i in range(0, stack_shape[0]):
                for j in range(0, stack_shape[1]):
                    if im[i,j,0] > 0 or im[i,j,1] > 0 or im[i,j,2] > 0:
                        nd_stack[i,j,k,:] = 1
    return nd_stack

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = command_line_parser(parser)
    (options, args) = parser.parse_args()

    # parse input
#   print "len(sys.argv)=", len(sys.argv), "sys.argv=", sys.argv
    if len(args) == 1:
        input_dir    = args[0]
    else:
        sys.exit("\nusage: MS_UT_CreateH5Stack.py input_dir -t type [ other options ] \n")

    if len(options.output_type) == 0:
        sys.exit("\nPlease, specify the input type (data, labels or mask) with option -t\n")
    elif not options.output_type in ["data", "labels", "mask"]:
        sys.exit("\nThe allowed input types are 'data', 'labels' and 'mask'\n")

    h5_file_name = options.output_name
    input_files = os.listdir(input_dir)
    list_stack = []
    first_file_path = input_dir + "/" + input_files[0]
    image_shape = misc.imread(first_file_path).shape
    nd_stack = create_default_stack(image_shape, len(input_files), options.output_type)
    nd_stack = populate_stack(nd_stack, input_dir, input_files, options.output_type,"")
        
#   print "nd_stack.shape=", nd_stack.shape
    if options.output_type == "labels" and len(image_shape) > 2: # input folder contains labels, not images
        if len(options.image_dir) == 0:
            nd_stack = update_labels_connect(nd_stack)            
        else:
            print "options.input_dir=", options.image_dir
            input_files = os.listdir(options.image_dir)
            first_file_path = options.image_dir + "/" + input_files[0]
            image_shape =  misc.imread(first_file_path).shape
            image_stack = create_default_stack(image_shape, len(input_files), "data")
            image_files = os.listdir(options.image_dir)
            image_stack = populate_stack(image_stack, options.image_dir, image_files, "data", "imdir")
            mean_signal = numpy.mean(image_stack[nd_stack[:,:,:,0] > 0])
            print "mean_signal=", mean_signal
            print "alpha=", options.alpha
            nd_stack    = update_labels_smooth(nd_stack, image_stack, mean_signal, float(options.alpha))

    # Transpose the array to  comply with Matlab
    nd_stack = numpy.transpose(nd_stack)
    print "Final nd_stack.shape=", nd_stack.shape
    # NOTE: for reasons which I don't understand, the Numpy's shape of nd_stack
    #       turns out to be inverted relative to the Matlab's size of h5 stack 
    f  = h5py.File(h5_file_name, 'w')
    f.create_dataset('main', nd_stack.shape, data=nd_stack)    



      
