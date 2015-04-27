#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

import re
import matplotlib
from PIL import Image
import h5py
import numpy
import sys, os, optparse
from scipy import misc

import MS_LIB_Options

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
def create_default_stack(image_data, num_images, output_type):
    print "...Creating a default stack for output_type=", output_type
    stack_shape = (image_data.shape[0],image_data.shape[1],num_images)
    if len(image_data.shape) > 2 or not output_type == "data":
        stack_shape = (image_data.shape[0],image_data.shape[1],num_images,3)
    print "image_data.shape=", image_data.shape, " stack_shape=", stack_shape
    dtype = "float"
    if not output_type == "data":
        dtype = "int" 
    if output_type == "mask":
        nd_stack = numpy.ones(stack_shape, dtype=dtype)  # 1 for mask
    else:
        nd_stack = numpy.zeros(stack_shape, dtype=image_data.dtype)
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

def populate_stack(nd_stack, input_dir, input_files, output_type, imdir, options):
    print "...Populating stack for output_type=", output_type
    i = 0
    for ifile in input_files:
        if i < int(options.zmin) or i >= int(options.zmax):
            continue
        ifile_path = os.path.join(input_dir, ifile)
        print "In populate_stack: i=", i, " ifile_path=", ifile_path
        try:
            print "Trying Image.open ..."
            im = Image.open(ifile_path)
        except:
            try:  
                print "Trying misc.imread ..."
                im = misc.imread(ifile_path)
            except:
                print "Trying matplotlib.image.imread ..."
                matplotlib.image.imread(ifile_path)
#       if output_type == "data" and len(imdir) == 0:
#           im = normalize_image(im)
#       elif output_type == "data" and len(imdir) > 0:
#           im = scale_image(im)
        i = i + 1
        print "Updating stack with ", ifile, " (i=", i, ")..."
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
        im_data = numpy.asarray(im)
        print "nd_stack.shape=", nd_stack.shape, " im_data.shape=", im_data.shape 
        if len(im_data.shape) == 2:
            # BW or Seg
            print "type(nd_stack)=", type(nd_stack[0,0,0]), " type(im_data)=",  type(im_data[0,0])
            im_data = numpy.array(im.getdata(), numpy.uint8).reshape(im.size[0], im.size[1])
            print "type(nd_stack)=", type(nd_stack[0,0,0]), " type(im_data)=",  type(im_data[0,0])
            nd_stack[:,:,k] = im_data
        else:
            # RGB image
            print "type(nd_stack)=", type(nd_stack[0,0,0,0]), " type(im_data)=",  type(im_data[0,0,0])
            im_data = numpy.array(im.getdata(), numpy.uint8).reshape(im.size[0], im.size[1], 3)
            nd_stack[:,:,k,:] = im_data
        print "...done"
    else:
        if len(im.shape) == 3 and output_type == "labels":
            for i in range(0, stack_shape[0]):
                for j in range(0, stack_shape[1]):
                    if im[i,j,0] > 0 or im[i,j,1] > 0 or im[i,j,2] > 0:
                        nd_stack[i,j,k,:] = 1
    return nd_stack

# ----------------------------------------------------------------------

def function_z_cmp(my_str):
    return int(my_str.split("_z")[1].split("_")[0])

# ----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.CreateH5Stack_command_line_parser(parser)
    (options, args) = parser.parse_args()

    # parse input
    if len(args) == 1:
        input_dir = args[0]
    else:
        sys.exit("\nusage: MS_UT_CreateH5Stack.py input_dir -t type [ other options ] \n")

    if len(options.output_type) == 0:
        sys.exit("\nPlease, specify the input type (data, labels or mask) with option -t\n")
    elif not options.output_type in ["data", "labels", "mask"]:
        sys.exit("\nThe allowed input types are 'data', 'labels' and 'mask'\n")

    h5_file_name = options.output_name
    print "input_dir=", input_dir, " match_string=", options.match_string, " unmatch_string=", options.unmatch_string
    input_files = []
    # Allow using match_string as a wild card
    match_string_list = options.match_string.split('*')
    for file in os.listdir(input_dir):
        match = True
        for m in match_string_list:
            if not re.search(m, file):
                match = False
        if match and \
           (len(options.unmatch_string) == 0 or not re.search(options.unmatch_string, file)):
            input_files.append(file)
    if re.search("_z", input_files[0]):
        sorted_input_files = sorted(input_files, key=function_z_cmp)
    else:
        sorted_input_files = sorted(input_files)
    imin = max(0, int(options.zmin))
    imax = min(len(input_files), int(options.zmax))
    input_files = sorted_input_files[imin:imax]
    print "sorted input_files=", input_files

    list_stack = []
    print "Opening file ", os.path.join(input_dir, input_files[0]), " ... "
    myfile = os.path.join(input_dir, input_files[0])
    try:
        # Using PIL
        image_data = numpy.asarray(Image.open(myfile))
    except:
        try:
            # Using scipy
            image_data = numpy.asarray(misc.imread(myfile))
        except:
            try:
                # Using Matplotlib
                image_data = numpy.asarray(matplotlib.image.imread(myfile, format = 'png'))
            except:
                print "Unable to read image file", myfile
                sys.exit(2)
    num_images = min(len(input_files), int(options.zmax)-int(options.zmin))
    nd_stack = create_default_stack(image_data, num_images, options.output_type)
    nd_stack = populate_stack(nd_stack, input_dir, input_files, options.output_type,"", options)
    if options.chunked:
        chunk_shape = (1, image_shape[0], image_shape[1])
        
#   print "nd_stack.shape=", nd_stack.shape
    if options.output_type == "labels" and len(image_shape) > 2: # input folder contains labels, not images
        if len(options.image_dir) == 0:
            nd_stack = update_labels_connect(nd_stack)            
        else:
            print "options.input_dir=", options.image_dir
            files = os.listdir(options.image_dir)
            if re.search("_z", files[0]):
                input_file_paths = sorted(files, key=function_z_cmp)
            else:
                input_file_paths = sorted(files)
            first_file_path  = input_file_paths[0]
            image_shape = misc.imread(first_file_path).shape
            num_images  = min(len(input_files), int(options.zmax)-int(options.zmin))
            image_stack = create_default_stack(image_shape, num_images, "data")
            image_files = os.listdir(options.image_dir)
            image_stack = populate_stack(image_stack, options.image_dir, image_files, "data", "imdir", options)
            mean_signal = numpy.mean(image_stack[nd_stack[:,:,:,0] > 0])
            print "mean_signal=", mean_signal
            print "alpha=", options.alpha
            nd_stack    = update_labels_smooth(nd_stack, image_stack, mean_signal, float(options.alpha))

    # Transpose the array to  comply with Matlab
    nd_stack = numpy.transpose(nd_stack)
    print "Final nd_stack.shape=", nd_stack.shape, " h5_file_name=", h5_file_name 
    # NOTE: for reasons which I don't understand, the Numpy's shape of nd_stack
    #       turns out to be inverted relative to the Matlab's size of h5 stack 
    f  = h5py.File(h5_file_name, 'w')
    if not options.chunked:
        f.create_dataset('main', nd_stack.shape, data=nd_stack)    
    else:
        if options.verbose:
            print "Chunk shape=", chunk_shape
        f.create_dataset('chunked', nd_stack.shape, data=nd_stack, chunks=chunk_shape)              



      
