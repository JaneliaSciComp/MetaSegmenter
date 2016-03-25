#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#

# Purpose: given segmentation stack, alighn a couple of layers
#          using ANTS                   

# ------------------------- imports -------------------------

import os, glob, gc
import shutil, commands
import sys, re, optparse
import numpy
import h5py
from scipy import misc
import imghdr
import matplotlib
from PIL import Image
import matlab.engine
import MS_LIB_Options

ms_home = os.environ['MS_HOME'] # where source code is located
ms_data = os.environ['MS_DATA'] # dir for initial input data and final output 
ms_temp = os.environ['MS_TEMP'] # dir for intermediate data files
nii_home= os.environ['NIFTI_HOME'] 

eng = matlab.engine.start_matlab()
eng.addpath(nii_home)

# -----------------------------------------------------------------------

def is_int(my_sting):     
    try: 
        int(my_sting)
        return True
    except ValueError:
        return False

# -----------------------------------------------------------------------------

def is_DVID_subvolume(input_data, options):
    is_subvolume = True
    if not input_data[0] == "[" or not input_data[len(input_data)-1] == "]":
        is_subvolume = False
        print "Input does not start with '[' or does not end with ']'"
        return is_subvolume
    ranges = input_data[1:len(input_data)-1].split(",")
    if not len(ranges) == 3:
        is_subvolume = False
        print "Input does not contain 3 ranges"
        return is_subvolume
    for range in ranges:
        coords = range.split(":")
        if not len(coords) == 2:
            is_subvolume = False
            print "Input coordinates", coords, " are missing ':'"
            return is_subvolume
        if not is_int(coords[0]) or \
           not is_int(coords[1]) or \
           int(coords[0]) > int(coords[1]):
            is_subvolume = False
            print "Input coordinates", coords, " are in bad order"
            return is_subvolume
    return is_subvolume

# -----------------------------------------------------------------------------

def get_input_type(input_data, options):
    input_type = ""
    if   os.path.isfile(input_data) or \
         os.path.isfile(os.path.join(ms_data,input_data)) or \
         os.path.isfile(os.path.join(ms_temp,input_data)):
        input_type = "file"
    elif os.path.isdir(input_data) or \
         os.path.isdir(os.path.join(ms_data,input_data)) or \
         os.path.isdir(os.path.join(ms_temp,input_data)):
        input_type = "directory"
    elif is_DVID_subvolume(input_data, options):
        input_type = "DVID"
    return input_type

# ----------------------------------------------------------------------

def get_data_dimensions(input_data, input_type, options):
    input_dim = [0,0,0]
    if input_type == "file":
        try:
            f  = h5py.File(os.path.join(ms_data,input_data), 'r')
            key = f.keys()[0]
            data = numpy.transpose(f[key])
            input_dim = data.shape[0:3]
        except:
            print "Input file",input_data,"is not an HDF5 stack\n"
            sys.exit(2)
    elif input_type == "directory":
        files = os.listdir(os.path.join(ms_data,input_data))
        num_files = 0
        image_size = []
        for ifile in files:
            myfile = os.path.join(os.path.join(ms_data,input_data), ifile)
            try:
                image_shape = Image.open(myfile).size
            except:
                myfile = os.path.join(os.path.join(ms_temp,input_data), ifile)
            try:
                suffix3 = myfile[(len(myfile)-3):len(myfile)]
                suffix4 = myfile[(len(myfile)-4):len(myfile)]
#               print "    suffix3=", suffix3, " suffix4=", suffix4
#               print "suffix3 in ['png','tif']:", suffix3 in ['png','tif']
                if  suffix3 in ['png', 'tif' ] or \
                    suffix4 in ['tiff','jpeg'] or \
                    imghdr.what(myfile) in ['png','tiff','jpeg']:
                    num_files = num_files + 1
                    if num_files == 1:
#                       print "Reading first file", myfile," ..."
#                       image_shape = misc.imread(myfile).shape
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
                    input_dim[0:2] = image_shape[0:2]
            except:
                continue
        input_dim[2] = min(num_files, int(options.zmax)-int(options.zmin))
    elif input_type == "DVID":
        ranges = input_data[1:len(input_data)-1].split(",")
        for i in range(0,3):
            coords = ranges[i].split(":")
            input_dim[i] = int(coords[1]) - int(coords[0]) 
    return input_dim

# ----------------------------------------------------------------------

def read_input_data(input_data, input_type, ymin, ymax, xmin, xmax, \
                    zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
#       key = f.keys()[0]
#       image_data1 = numpy.transpose(f[key][zmin]) # read one chunk
#       image_data2 = numpy.transpose(f[key][zmax-1])
        key = f.keys()[0]
        dataset = numpy.transpose(f[key])
        data_shape = dataset.shape
        num_layers = dataset.shape[2]
        image_data1= dataset[:,:,zmin]
        image_data2= dataset[:,:,zmin+1]
        if options.verbose:
            print "num_layers=", num_layers
            print "numpy.max(image_data1)=", numpy.max(image_data1)
            print "numpy.max(image_data2)=", numpy.max(image_data2)
    elif input_type == "directory":
        files = []
        i = 0
        dir_path = os.path.join(ms_data,input_data)
        files_list = os.listdir(dir_path)
        files_list.sort()
        files = []
        for file in files_list:
            file_path = os.path.join(dir_path, file)
            if i in range(int(zmin), int(zmax)):
                files.append(os.path.join(input_data, file_path))
            i = i+1
        im1 =  misc.imread(files[0])
        data_shape = misc.imread(files[0]).shape
        image_data1 = numpy.zeros(data_shape, dtype="float")
        image_data2 = numpy.zeros(data_shape, dtype="float")
        i = 0
        for file in files:
            im = misc.imread(file)
            image_data1[:,:,i] = im[ymin:ymax, xmin:xmax,1]
            i = i + 1
    elif input_type == "DVID":
        connection = httplib.HTTPConnection("emdata2.int.janelia.org:80",\
                     timeout=5.0 )
        dvid_volume = voxels.VoxelsAccessor(connection, options.uuid, \
                                            options.dataset)
        image_data1 = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmin])
        image_data2 = numpy.asarray(dvid_volume[:, ymin:ymax, xmin:xmax, zmax-1])
        # Put data in C-contiguous order (last-index varies the fastest).
        image_data1 = image_data1.copy(order='C')
        image_data2 = image_data2.copy(order='C')
    image_fragment1 = image_data1[ymin:ymax, xmin:xmax]
    image_fragment2 = image_data2[ymin:ymax, xmin:xmax]
    return (image_fragment1, image_fragment2)

# ----------------------------------------------------------------------

# Create and submit all cluster jobs
def process_inputs(input_data, input_type, input_dim, options):
 
    output_path = options.output_name 
    if len(options.output_name) == 0:
        output_path = "BW_align_12.png"
    
    if input_type == 'DVID':
        input_label = "ms3_DVID"
    elif input_type == "directory":
        input_label = "ms3_" + input_data[4:]
    else:
        input_label = "ms3_" + input_data.split('.')[0][4:]

    ymin = 0
    ymax = int(input_dim[0])
    xmin = 0
    xmax = int(input_dim[1])
    zmin = int(options.zmin)
    zmax = int(options.zmax)

    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, z=",\
              [y, ymin, ymax, x, xmin, xmax, z]

    # Extract image data
    labels1, labels2 = \
        read_input_data(input_data, input_type, \
        ymin, ymax, xmin, xmax, zmin, zmax, options)

    if options.verbose:
        print "image_fragment1.shape=", image_fragment1.shape, \
             " image_fragment2.shape=", image_fragment2.shape
        print "\n...Generating overlap matrix..."

    bw1 = numpy.zeros(labels1.shape, 'uint8')
    bw1[labels1 == 0] = 255
#   BW1 = matlab.double(numpy.ndarray.tolist(bw1))
#   eng.make_and_save_nii(BW1, 'BW1.nii')
    bw1_path = "BW1.png"
    img = Image.fromarray(bw1)
    img.save(bw1_path)
    nii1_path = eng.bw2nii(bw1_path)

    bw2 = numpy.zeros(labels2.shape, 'uint8')
    bw2[labels2 == 0] = 255
#   BW2 = matlab.double(numpy.ndarray.tolist(bw2))
#   eng.make_and_save_nii(BW2, 'BW2.nii')
    bw2_path = "BW2.png"
    img = Image.fromarray(bw2)
    img.save(bw2_path)
    nii2_path = eng.bw2nii(bw2_path)

    command_prefix = "ANTS 2 -m CC[" + nii2_path + ", " + nii1_path + ", 1, 8]"

    command1 = command_prefix + "-o BW_align_12 -i 0 --number-of-affine-iterations " \
                              + " 10000x10000x10000"
    print "command1=", command1
    os.system(command1)
    
    initaffine="BW_align_12Affine.txt"
    diffeom=   "BW_align_12_diffeom"

    command2 = command_prefix + " -t SyN[0.25] -r Gauss[3,0] -o   " + diffeom \
                              + " -i 100x100x100 --initial-affine " + initaffine
    print "command2=", command2
    os.system(command2)

    affmat=diffeom + "Affine.txt"
    warpfd=diffeom + "Warp.nii.gz"

    command3 = "WarpImageMultiTransform 2 " + nii1_path + " BW_align_12.nii  -R " + nii2_path \
             + " " + warpfd + " " + affmat + " --use-BSpline"
    print "command3=", command3
    os.system(command3)

#   nii = eng.load_nii('BW_align_12.nii')
#   eng.imwrite(nii.img, 'BW_align_12.png', 'png')
    bw_path = eng.nii2bw('BW_align_12.nii')


# -----------------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_data [options (-h to list)] \nwhere\n\
    input_data=\n\
    1) path to image stack file in HDF5 format, or \n\
    2) path to folder containing image files, or \n\
    3) coordinates of DVID subvolume: [xmin:xmax,ymin:ymax,zmin:zmax]"   

#   Steps of processing:
#   1) MS_F3D_GenerateMatrices.py
#   2) MS_F3D_MergeMatrices.py
#   3) MS_F3D_TraverseFusionTrees.py
#   4) MS_F3D_RelabelSegmentedData.py
#   5) produce ouput (in HDF5 or other format)

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.Align3D_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if len(args) == 1:
        input_data = args[0]
        if not os.path.exists(os.path.join(ms_data,input_data)):
            sys.exit("\nInput data is not found")
        input_type = get_input_type(input_data, options)

        if input_type in ["file", "directory", "DVID"]:
            if options.verbose and int(options.node) == 0:
                print "Input type=", input_type
        else:
            print "\nIncorrectly specified input data", input_data, "\n"
            parser.print_usage()
            sys.exit(2)         
        input_dim = get_data_dimensions(input_data, input_type, options)
        input_dim1 = [input_dim[0], input_dim[1], input_dim[2]-1]
        if options.verbose:
            print "Input data dimensions: ", input_dim, " num_nodes=", num_nodes
            print "\nProcessing input data ..."
        process_inputs(input_data, input_type, input_dim, options)
    else:
        parser.print_usage()
        sys.exit(2)
