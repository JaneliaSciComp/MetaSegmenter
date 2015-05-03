#! /usr/local/python-2.7.6/bin/python
#
# Copyright (C) 2015 by Howard Hughes Medical Institute.
#
import os
import sys, re, optparse
import numpy
from scipy import misc
import httplib
import h5py
from pydvid import voxels, general

import MS_LIB_Dict
import MS_LIB_Options

ms_home = os.environ['MS_HOME']
ms_data = os.environ['MS_DATA']
ms_temp = os.environ['MS_TEMP']

# ----------------------------------------------------------------------

def generate_overlap_matrix(image_fragment1, image_fragment2, options):
    fragment_matrix1 = numpy.matrix(image_fragment1)      
    fragment_matrix2 = numpy.matrix(image_fragment2)
    max_label1 = fragment_matrix1.max()
    max_label2 = fragment_matrix2.max()
    if options.verbose:
        print "max_label1=", max_label1, " max_label2=", max_label2    
    labels1 = []
    labels2 = []
    for i in range(1, int(max_label1)+1):
        count = (fragment_matrix1 == i).sum()
        if count > 0:
            labels1.append(i)
    if options.verbose:
        print "len(labels1)=", len(labels1)
    for i in range(1, int(max_label2)+1):
        count = (fragment_matrix2 == i).sum()
        if count > 0:
            labels2.append(i)
    overlap_matrix = numpy.zeros((len(labels1)+1, len(labels2)+1), dtype=int)
    if options.verbose:
        print "len(labels2)=", len(labels2)
    overlap_matrix[0 , 0] = 0
    overlap_matrix[1:, 0] = labels1
    overlap_matrix[0 ,1:] = labels2 
    for i in range(0, len(labels1)):
        for j in range(0, len(labels2)):
            area1 = (fragment_matrix1 == int(labels1[i])).sum()
            area2 = (fragment_matrix2 == int(labels2[j])).sum()
            overlap_matrix[i+1, j+1] = \
                ((fragment_matrix1 == int(labels1[i]))  & \
                 (fragment_matrix2 == int(labels2[j]))).sum()
            if  overlap_matrix[i+1, j+1] < float(options.overlap_fraction) * min(area1, area2):
                overlap_matrix[i+1, j+1] = 0.
            if  overlap_matrix[i+1, j+1] < float(options.overlap_area):
                overlap_matrix[i+1, j+1] = 0.
            if options.verbose and options.debug and overlap_matrix[i+1, j+1] > 0:
                print "i=", i, "/", len(labels1), " j=", j, "/", len(labels2), \
                      " overlap=", overlap_matrix[i+1, j+1]
    return overlap_matrix

# -----------------------------------------------------------------------

def read_input_data(input_data, input_type, ymin, ymax, xmin, xmax, \
                    zmin, zmax, options):
    if input_type == "file":
        input_path = os.path.join(ms_data, input_data)
        f  = h5py.File(input_path, 'r')
#       key = f.keys()[0]
#       image_data1 = numpy.transpose(f[key][zmin]) # read one chunk
#       image_data2 = numpy.transpose(f[key][zmax-1])
        dataset = numpy.transpose(f['/main'])
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

# -----------------------------------------------------------------------

def output_overlap_matrix(overlap_matrix, input_label, x, y, z, options):
    # Extract output path
    if len(options.output_path) == 0:
        output_path = os.path.join(ms_temp, input_label + "_overlap" +\
                                   "_y" + str(y+1) +\
                                   "_x" + str(x+1) + "_z" + str(z+1) + ".txt")
    else:
        output_path = options.output_path

    # Ourput data
    if options.verbose:
        print "Storing the overlap matrix in file ", output_path
    numpy.savetxt(output_path, overlap_matrix, fmt='%10u', delimiter='\t')

# -----------------------------------------------------------------------

def convert_to_sparse_format(omatrix):
    sparse_omatrix = [ omatrix.shape[0]-1,  omatrix.shape[1]-1, 0 ]
    for i in range(1,omatrix.shape[0]):
        for j in range(1, omatrix.shape[1]):
            if omatrix[i][j] == 1:
                sparse_omatrix.append([omatrix[i,0], omatrix[0,j], omatrix[i,j]])
    return numpy.array(sparse_omatrix)

# -----------------------------------------------------------------------

def process_inputs(input_data, input_type, dict_node_xyz, options):
    if input_type == 'DVID':
        input_label = "ms3_DVID" 
    elif input_type == "directory":
        input_label = "ms3_" + input_data[4:] 
    else:
        input_label = "ms3_" + input_data.split('.')[0][4:]

    node = int(options.node)
    if node > 0:
        y, ymin, ymax, x, xmin, xmax, z = dict_node_xyz[node]
        zmin = z
        zmax = z+1
    else:
        ymin = 0
        ymax = int(options.ydim)
        xmin = 0
        xmax = int(options.xdim)
        zmin = int(options.zmin)
        zmax = int(options.zmax)

    if options.verbose:
        print "    node=", node, "y, ymin, ymax, x, xmin, xmax, z=",\
              [y, ymin, ymax, x, xmin, xmax, z]

    # Extract image data
    image_fragment1, image_fragment2 = \
        read_input_data(input_data, input_type, \
        ymin, ymax, xmin, xmax, zmin, zmax, options)

    if options.verbose:
        print "image_fragment1.shape=", image_fragment1.shape, \
             " image_fragment2.shape=", image_fragment2.shape
        print "\n...Generating overlap matrix..."
    overlap_matrix = generate_overlap_matrix(image_fragment1, image_fragment2, \
                                             options) 
#   print "overlap_matrix=\n", overlap_matrix
    sparse_overlap_matrix = convert_to_sparse_format(overlap_matrix)
    output_overlap_matrix(sparse_overlap_matrix, input_label, \
                          x, y, zmin, options)

# -----------------------------------------------------------------------

def convert_to_sparse_format(omatrix):
    sparse_matrix = [[ omatrix.shape[0]-1, omatrix.shape[1]-1, 0]]
    for i in range(1,omatrix.shape[0]):
        for j in range(1, omatrix.shape[1]):
            if omatrix[i][j] > 0: 
                sparse_matrix.append([omatrix[i,0],omatrix[0,j],omatrix[i,j]])
    print "sparse_matrix=\n", sparse_matrix
    return numpy.array(sparse_matrix)

# -----------------------------------------------------------------------

if __name__ == "__main__":

    usage = "\nUsage: %prog input_data input_type [ options ]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = MS_LIB_Options.GenerateMatrices_command_line_parser(parser)
    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        input_data, input_type = args
        input_dim = [int(options.ydim), int(options.xdim), int(options.zdim)]

        if input_type == 'DVID':
            input_label = "ms3_DVID"
        elif input_type == "directory":
            input_label = "ms3_" + os.path.basename(input_data)[4:]
        else:
            input_label = "ms3_" + os.path.basename(input_data).split('.')[0][4:]

        input_dim1 = [input_dim[0], input_dim[1], input_dim[2]-1];
        num_nodes, dict_node_xyz = \
            MS_LIB_Dict.map_node_to_xyz(input_dim1, input_label+"_fusion", ".txt", options)
        if options.verbose:
#           print "dict_node_xyz=", dict_node_xyz
            print "In MS_F3D_GenerateOverlapMatrix.py: num_nodes=", num_nodes
        if len(options.output_path) == 0 and len(options.node) == 0:
            print "\nPlease, specify the output name (with option -o)"
            sys.exit(2) 
        process_inputs(input_data, input_type, dict_node_xyz, options)
    else:
        print usage
        sys.exit(2)

