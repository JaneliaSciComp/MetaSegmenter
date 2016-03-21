#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os
import sys, re, optparse
import tifffile as tiff
import numpy
import h5py
from scipy import misc
import httplib
from pydvid import voxels, general
import mimetypes
import Image

import MS_LIB_Dict
import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------------

def get_input_type(input_data, options):
    input_type = ""
    if   os.path.isfile(input_data) or \
         os.path.isfile(os.path.join(ms_data,input_data)) or \
         os.path.isfile(os.path.join(ms_temp,input_data)):
        input_type = "file"
        # Check if a text file
        if mimetypes.guess_type(input_data)[0] == 'text/plain':
            input_type = "fileoffiles"
    elif os.path.isdir(input_data) or \
         os.path.isdir(os.path.join(ms_data,input_data)) or \
         os.path.isdir(os.path.join(ms_temp,input_data)):
        input_type = "directory"
#   elif is_DVID_subvolume(input_data, options):
#       input_type = "DVID"
    return input_type

# -----------------------------------------------------------------------------

def get_data_dimensions(input_data, input_type, options):
    input_dim = [0,0,0]
    if input_type == "file":
        input_data_path = os.path.join(ms_data,input_data)
        if re.search(".h5", input_data):
            try:
                f  = h5py.File(input_data_path, 'r')
                key = f.keys()[0]
                image_data = numpy.transpose(f[key])
                input_dim = image_data.shape[0:3]
            except:
                print "Input file",input_data,"is not an HDF5 stack\n"
                sys.exit(2)
        elif re.search(".tif", input_data):
            try:
                image_data = numpy.transpose(tiff.imread(input_data_path))
                input_dim = image_data.shape
                print "In get_data_dimensions: input_dim=", input_dim
            except:
                print "Input file",input_data_path,"is not a tiff stack\n"
                sys.exit(2)
        else:
            print "Unsupported input file\n"
            sys.exit(2)
    elif input_type == "directory":
        files = sorted(os.listdir(os.path.join(ms_data,input_data)))
        if options.verbose:
            print "len(files)=", len(files)
        num_files = 0
        image_size = []
        for i in range(0, len(files)):
            ifile = files[i]
            myfile = os.path.join(os.path.join(ms_data,input_data), ifile)
            try:
                image_shape = Image.open(myfile).size
                input_dim[0:2] = [image_shape[1],image_shape[0]]
                num_files = num_files + 1
            except:
                print "Cannot open ", myfile
                try:
                    if re.search(".png",  myfile) or \
                       re.search(".tif",  myfile) or \
                       re.search(".jpeg", myfile):
                        num_files = num_files + 1
                        if num_files == 1:
#                           print "Reading first file", myfile," ..."
#                           image_shape = misc.imread(myfile).shape
                            try:
                                # Using PIL
                                image_shape = Image.open(myfile).size
                            except:
                                try:
                                    # Using scipy
                                    image_shape = misc.imread(myfile).shape
                                except:
                                    try:
                                        # Using Matplotlib
                                        image_shape = matplotlib.image.imread(myfile, format = 'png').shape
                                    except:
                                        print "Unable to read image file", myfile
                                        sys.exit(2)
                            input_dim[0:2] = [image_shape[1],image_shape[0]] # transpose
                except:
                    continue
        input_dim[2] = min(num_files, int(options.zmax)-int(options.zmin))
        print "num_files=", num_files, " int(options.zmax)-int(options.zmin=", int(options.zmax)-int(options.zmin), " input_dim[2]=", input_dim[2]

    elif input_type == "DVID":
        ranges = input_data[1:len(input_data)-1].split(",")
        for i in range(0,3):
            coords = ranges[i].split(":")
            input_dim[i] = int(coords[1]) - int(coords[0])
    return input_dim

# -----------------------------------------------------------------------

def read_input(input_data,input_type,ymin,ymax,xmin,xmax,zmin,zmax,options):
    # Extract image data
    if input_type == "file":
        input_data_path = os.path.join(ms_data, input_data)
        if re.search(".h5", input_data): 
            f  = h5py.File(input_data_path, 'r')
            key = f.keys()[0]
            data = numpy.transpose(f[key])

            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
        elif re.search(".tif", input_data):
            if options.verbose:
                print "tiff input_data=", input_data_path
            data = numpy.transpose(tiff.imread(input_data_path))
            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
            print "Input data shape=", image_data.shape
        else:
            sys.exit("\nUnsupported type of input file. Exit.")
    elif input_type == "directory":
        # Create a list of files
        files = []
        i = 0
        dir_path = os.path.join(ms_data,input_data)
        files_list = sorted(os.listdir(dir_path))
        files = []
        for i in range(0, len(files_list)):
            file = files_list[i]
            file_path = os.path.join(dir_path, file)
            if i in range(int(zmin), int(zmax)):
                files.append(os.path.join(input_data, file_path))
            i = i+1
        if re.search(".tif", files[0]):
            im1 = tiff.imread(files[0])   
        else:
            im1 = misc.imread(files[0])
        data_shape = im1.shape
        image_data = numpy.zeros((ymax-ymin,xmax-xmin), dtype="float")
        # Read files
        i = 0
        for i in range(0, len(files)):
            file = files[i]
            if re.search(".tif", file):
                im = tiff.imread(file)
            else:
                im = misc.imread(file)
            if options.verbose:
                print "file=", file, " data_shape=", data_shape, " im.shape=", im.shape
            if len(im.shape) == 2:
                image_data[:,:] = im[ymin:ymax, xmin:xmax]
            else:
                image_data[:,:] = im[ymin:ymax, xmin:xmax,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmin:zmax])
        image_data = image_data.copy(order='C')
    return image_data

# -------------------------------------------------------------------------------

def read_input1(seg_stack_h5_path, layer_id):
    if re.search(".h5", seg_stack_h5_path):
        f    = h5py.File(seg_stack_h5_path, 'r')
        key  = f.keys()[0]
        data = f[key]
    elif re.search(".tif", seg_stack_h5_path):
        data = tiff.imread(seg_stack_h5_path)
    else:
        sys.exit("Unrecognized stack type")
    print "data.shape=", data.shape, " layer_id=", layer_id
    if layer_id >= 0:
        return data[layer_id,:,:]
    return data

