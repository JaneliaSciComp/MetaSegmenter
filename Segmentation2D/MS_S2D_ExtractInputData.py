#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os
import sys, re, optparse
import tifffile as tiff
import numpy
from scipy import misc
import httplib
from pydvid import voxels, general

import MS_Dict
import MS_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# -----------------------------------------------------------------------

def extract_image_data(input_data, input_type, ymin, ymax,\
                       xmin, xmax, zmin, zmax, options):
    # Extract image data
    if input_type == "file":
        if re.search(".h5", input_data): 
            f  = h5py.File(input_data, 'r')
            key = f.keys()[0]
            data = numpy.transpose(f[key])
            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
        elif re.search(".tif", input_data):
            if options.verbose:
                print "tiff input_data=", input_data
            data = numpy.transpose(tiff.imread(input_data))
            data_shape = data.shape
            if len(data.shape) == 3:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax]
            else:
                image_data = data[ymin:ymax, xmin:xmax, zmin:zmax, 1]
            print "Input data shape=", image_data.shape
        else:
            sys.exit("\nUnsupported type of input file. Exit.")
    elif input_type == "directory":
        files = []
        i = 0
        dir_path = os.path.join(ms_data,input_data)
        files_list = sorted(os.listdir(dir_path))
        files = []
        for file in files_list:
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
        i = 0
        for file in files:
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

# -----------------------------------------------------------------------

def process_inputs(input_data, input_type, input_label, dict_node_xyz, options):

    node = int(options.node)
    if node > 0:
        y, ymin, ymax, x, xmin, xmax, z = dict_node_xyz[node]
        zmin = z
        zmax = z + 1
    else:
        ymin = 0
        ymax = int(options.ydim)
        xmin = 0
        xmax = int(options.xdim)
        zmin = max(0,               int(options.zmin))
        zmax = min(int(options.zdim,int(options.zmax)))
    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, z=",\
              [y, ymin, ymax, x, xmin, xmax, z]
        print "    ysize=", ymax - ymin

    # Extract image data
    image_data = extract_image_data(input_data, input_type, ymin, ymax,\
                                    xmin, xmax, zmin, zmax, options)
    # Extract output path
    if len(options.output_path) == 0:
        if   int(options.nx) > 1 and int(options.ny) > 1:
            output_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                       "_x" + str(x+1) + "_z" + str(z+1) + ".png")
        elif int(options.nx) == 1 and int(options.ny) > 1:
            output_path = os.path.join(ms_temp, input_label + "_y" + str(y+1) +\
                                                         "_z" + str(z+1) + ".png")
        elif int(options.nx) > 1 and int(options.ny) == 1:
            output_path = os.path.join(ms_temp, input_label + \
                                       "_x" + str(x+1) + "_z" + str(z+1) + ".png")
        else:
            output_path = os.path.join(ms_temp, input_label + \
                                                         "_z" + str(z+1) + ".png")
    else:
        output_path = options.output_path

    # Ourput data
    if options.verbose:
        print "Image_data shape=", image_data.shape
        print "Writing image data to file", output_path
    misc.imsave(output_path, image_data)   

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [ options ]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_Options.ExtractInputData_command_line_parser(parser)
    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        input_data, input_type = args
        input_dim = [int(options.ydim), int(options.xdim), int(options.zdim)]
        input_label = "ms2_DVID"
        if input_type in ["file", "directory"]:
            input_label = "ms2_" + input_data.split('.')[0]
        num_nodes, dict_node_xyz = \
            MS_Dict.map_node_to_xyz(input_dim, input_label, ".png", options)
        if options.verbose:
            print "num_nodes=", num_nodes
            print "input_dim=", input_dim
        if len(options.output_path) == 0 and len(options.node) == 0:
            print "\nPlease, specify the output name (with option -o)"
            sys.exit(2) 
        if num_nodes > 0:
            process_inputs(input_data, input_type, input_label, dict_node_xyz, options)
    else:
        print usage
        sys.exit(2)

